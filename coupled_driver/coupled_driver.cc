/*
 =============================================================================
 coupled_driver.cc

 External driver coupling Cytosim (filaments) and Mem3DG (membrane).

 Neither simulator drives the other: this program owns the loop and alternates
 between them (Lie-Trotter / Strang operator splitting).  Each iteration:

   1. gather the filament segments from Cytosim
   2. hand them to Mem3DG, which computes the contacts and keeps the membrane
      reaction internally      (mem3dg_system_filament_repulsion)
   3. step Cytosim; the filament-side forces are injected inside its implicit
      solve via Simul::externalForceCallback -> Meca::addForce
   4. step Mem3DG, which applies the reaction it stored in step 2

 WHY A C ABI
 -----------
 Cytosim is C++11, Mem3DG is C++20 and exposes Eigen/geometry-central in its
 headers.  This file therefore talks to Mem3DG only through <mem3dg/capi/mem3dg_c.h>,
 which is plain C over an opaque handle.  Nothing here needs C++20.

 STERIC MODEL
 ------------
 The contact geometry lives in Mem3DG's edge-face kernel, reached through the C
 ABI: for every (filament segment, membrane triangle) pair it finds the true
 minimum distance -- including the piercing case -- and decomposes a harmonic
 repulsion with DDG barycentric weights, so Newton's third law holds exactly.

 An earlier version computed contacts here instead, comparing each filament
 POINT to each membrane VERTEX.  That is blind to a filament crossing the
 middle of a triangle away from its corners, which is exactly how filaments
 tunnel through a membrane.

 Cost is still O(nSegments * nFaces) per call, with a centroid quick-reject
 inside Mem3DG.  A BVH or spatial grid is the next optimisation if filament
 counts grow.
 =============================================================================
*/

#include <mem3dg/capi/mem3dg_c.h>

#include "simul.h"
#include "parser.h"
#include "meca.h"
#include <map>
#include "fiber.h"
#include "fiber_set.h"
#include "single.h"
#include "single_set.h"
#include "couple.h"
#include "couple_set.h"
#include "glossary.h"
#include "exceptions.h"
#include "spaces/space_mem3dg.h"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>

//------------------------------------------------------------------------------

/// one filament-side force to be added inside Cytosim's solve
struct FiberForce {
    Mecable const* mecable;   ///< the fiber owning the point
    unsigned       point;     ///< LOCAL point index within that fiber
    Vector         force;
};

/// parameters of the toy steric interaction
struct StericParams {
    double cutoff   = 0.05;   ///< interaction range; in confine mode, the shell
    double stiff    = 1.0;    ///< harmonic stiffness
    /*
     0 = symmetric repulsion (mem3dg_system_filament_repulsion)
         pushes a segment away from the surface on whichever side it is, within
         `cutoff`.  Works on open OR closed meshes, but does NOT confine: a
         segment that gets outside is pushed further out and never returns.

     1 = one-sided wall (mem3dg_system_filament_confinement)
         keeps filaments INSIDE and pushes the membrane OUTWARD.  `cutoff` is
         read as the shell clearance.  Requires a CLOSED mesh.
     */
    int    mode     = 0;
};

//------------------------------------------------------------------------------

/// which fiber point a segment endpoint belongs to
struct SegmentRef {
    Mecable const* mecable;   ///< the fiber owning this segment
    unsigned       k;         ///< LOCAL index of p0; p1 is point k+1
};

/// scratch buffers reused across steps, so the loop does not reallocate
struct ContactWork {
    std::vector<double>     segP0, segP1;   ///< nSeg*3, flat xyz (repulsion mode)
    std::vector<double>     segF0, segF1;   ///< nSeg*3, forces returned
    std::vector<SegmentRef> segRef;
    // confinement mode works per POINT, so it needs its own arrays
    std::vector<double>     ptP, ptF;       ///< nPts*3
    std::vector<SegmentRef> ptRef;          ///< .k is the point's local index
};


/**
 Compute filament-membrane steric forces for the current configuration.

 The contact geometry is evaluated INSIDE Mem3DG, via
 mem3dg_system_filament_repulsion().  That kernel measures the true minimum
 distance between each filament segment and each membrane triangle -- including
 the case where a segment pierces a triangle, which is precisely how filaments
 tunnel through a membrane.  It then decomposes the force with DDG barycentric
 weights, so Newton's third law holds exactly.

 This replaces an earlier vertex-to-vertex loop written here in the driver,
 which compared each filament POINT against each membrane VERTEX.  That test
 sees nothing when a filament crosses the middle of a triangle without coming
 near any of its three corners, so filaments could pass straight through; it
 also treated the membrane as a cloud of points rather than a surface.

 Membrane side: accumulated into Mem3DG's external force buffer by the call
 itself (add_to_membrane=1), after we clear it, so it holds exactly this
 contribution and is applied on the next mem3dg_system_step().

 Filament side: returned per segment ENDPOINT and folded back onto fiber points
 here.  A fiber of N points has N-1 segments, and point k is shared by segment
 k-1 (as its p1) and segment k (as its p0), so both contributions are summed.

 Returns the number of contacting (segment, face) pairs.
 */
static long contactForces(Simul const&             simul,
                          mem3dg_system*           mem,
                          const StericParams&      p,
                          ContactWork&             w,
                          std::vector<FiberForce>& fiberForces)
{
    // ---- gather every filament segment as a flat array ----------------------
    w.segP0.clear(); w.segP1.clear(); w.segRef.clear();

    for ( Fiber const* F = simul.fibers.first(); F; F = F->next() )
    {
        const unsigned n = F->nbPoints();
        for ( unsigned k = 0; k + 1 < n; ++k )
        {
            Vector a = F->posPoint(k), b = F->posPoint(k+1);
            w.segP0.push_back(a.XX); w.segP0.push_back(a.YY); w.segP0.push_back(a.ZZ);
            w.segP1.push_back(b.XX); w.segP1.push_back(b.YY); w.segP1.push_back(b.ZZ);
            SegmentRef r; r.mecable = F; r.k = k;
            w.segRef.push_back(r);
        }
    }

    fiberForces.clear();
    const size_t nSeg = w.segRef.size();
    if ( nSeg == 0 )
    {
        mem3dg_system_clear_external_forces(mem);
        return 0;
    }

    w.segF0.assign(nSeg*3, 0.0);
    w.segF1.assign(nSeg*3, 0.0);

    // ---- the actual contact computation, inside Mem3DG ----------------------
    mem3dg_system_clear_external_forces(mem);

    if ( p.mode == 1 )
    {
        /*
         Confinement works per POINT, not per segment.  Passing each filament
         vertex exactly once means every one carries its own full wall force.
         The earlier segment form had to halve each endpoint (consecutive
         segments share interior points), which left each fiber's two TERMINAL
         points at half strength -- and those are precisely the points that lead
         growth into the membrane, so they were the ones that leaked out.
         */
        w.ptP.clear(); w.ptRef.clear();
        for ( Fiber const* F = simul.fibers.first(); F; F = F->next() )
        {
            const unsigned n = F->nbPoints();
            for ( unsigned k = 0; k < n; ++k )
            {
                Vector a = F->posPoint(k);
                w.ptP.push_back(a.XX); w.ptP.push_back(a.YY); w.ptP.push_back(a.ZZ);
                SegmentRef r; r.mecable = F; r.k = k;
                w.ptRef.push_back(r);
            }
        }
        const size_t nPts = w.ptRef.size();
        if ( nPts == 0 )
            return 0;
        w.ptF.assign(nPts*3, 0.0);

        long n = mem3dg_system_filament_confinement(
            mem, &w.ptP[0], nPts, /*shell=*/p.cutoff, p.stiff,
            &w.ptF[0], /*add_to_membrane=*/1);

        if ( n < 0 )
        {
            std::fprintf(stderr,
                "error: filament_confinement failed (needs a CLOSED mesh)\n");
            return 0;
        }

        for ( size_t i = 0; i < nPts; ++i )
        {
            const Vector f(w.ptF[3*i+0], w.ptF[3*i+1], w.ptF[3*i+2]);
            if ( f.normSqr() > 0 )
            {
                FiberForce ff;
                ff.mecable = w.ptRef[i].mecable;
                ff.point   = w.ptRef[i].k;
                ff.force   = f;
                fiberForces.push_back(ff);
            }
        }
        return n;
    }

    long nContacts = mem3dg_system_filament_repulsion(
        mem, &w.segP0[0], &w.segP1[0], nSeg,
        p.cutoff, p.stiff,
        &w.segF0[0], &w.segF1[0], /*add_to_membrane=*/1);

    if ( nContacts < 0 )
    {
        std::fprintf(stderr, "error: filament_repulsion failed\n");
        return 0;
    }

    return nContacts;
}

//------------------------------------------------------------------------------

/**
 Write one membrane frame as OBJ into `dir` with a zero-padded index, matching
 the sequence the Blender importer expects:  <dir>/membrane_0000.obj, ...

 NOTE the extension is .obj, not .ply.  mem3dg_system_save_ply() is a misnomer:
 it forwards to geometry-central's writer, which only implements OBJ output --
 asking for .ply detects the type, opens the file, then throws, leaving a
 0-byte file.  OBJ imports into Blender natively, so this is the right target
 anyway.  mem3dg_system_save_ply() returns 0 on success, non-zero on failure.

 The Cytosim side is exported separately, offline, by the `export_fibers` tool
 reading objects.cmo -- see COORDINATES.md.  Both use the SAME frame index so
 Blender overlays them: membrane frame f and fiber frame f are the same instant.
 */
static void saveMembraneFrame(mem3dg_system* mem, const std::string& dir, int frame)
{
    char path[512];
    std::snprintf(path, sizeof(path), "%s/membrane_%04d.obj", dir.c_str(), frame);
    if ( mem3dg_system_save_ply(mem, path) != 0 )
        std::fprintf(stderr, "warning: could not write %s\n", path);
}


/**
 Write one filament frame as polylines, in the SAME format and directory as the
 membrane frames:  <dir>/fibers_0000.txt, ...

 Why here rather than via the `export_fibers` tool: this driver parses the
 config with Run=0, Write=0, so Cytosim never writes objects.cmo -- there is no
 trajectory for export_fibers to read.  Exporting from inside the coupled loop
 also makes synchronisation exact: fibers_000K.txt and membrane_000K.obj are
 written at the same instant of the same run, by construction.

 (export_fibers remains the right tool for a standalone Cytosim trajectory; the
 file format written here is identical, so the Blender importer handles both.)
 */
static void saveFiberFrame(Simul const& simul, const std::string& dir, int frame)
{
    char path[512];
    std::snprintf(path, sizeof(path), "%s/fibers_%04d.txt", dir.c_str(), frame);
    std::FILE* f = std::fopen(path, "w");
    if ( ! f )
    {
        std::fprintf(stderr, "warning: could not write %s\n", path);
        return;
    }

    std::fprintf(f, "# cytosim fibers  frame %d  time %g\n", frame, simul.time());

    for ( Fiber const* F = simul.fibers.first(); F; F = F->next() )
    {
        // steric_radius gives a natural on-screen thickness; fall back if unset
        double radius = ( F->prop && F->prop->steric_radius > 0 )
                      ? F->prop->steric_radius : 0.005;

        const unsigned n = F->nbPoints();
        std::fprintf(f, "fiber %lu %g %u\n",
                     (unsigned long)F->identity(), radius, n);

        for ( unsigned k = 0; k < n; ++k )
        {
            Vector p = F->posPoint(k);
#if ( DIM >= 3 )
            std::fprintf(f, "%g %g %g\n", p.XX, p.YY, p.ZZ);
#elif ( DIM == 2 )
            std::fprintf(f, "%g %g %g\n", p.XX, p.YY, 0.0);
#else
            std::fprintf(f, "%g %g %g\n", p.XX, 0.0, 0.0);
#endif
        }
    }
    std::fclose(f);
}

/**
 Write one singles/couples frame as point clouds, alongside the membrane and
 fiber frames:  <dir>/points_0000.txt, ...

 Points are grouped into CLASSES -- one block per (kind, property name, binding
 state) -- rather than written one per line with attributes.  Two reasons:

  - Blender renders a class as a single vertex-instanced mesh.  One object per
    class costs one draw call; one object per single would cost hundreds per
    frame and make scrubbing unusable.
  - binding state is what you actually want to see (does the couple bridge two
    filaments, or is it just diffusing?), and having it in the block header lets
    the importer give each state its own material without per-vertex attributes.

 The state comes for free from Cytosim's own list structure: SingleSet keeps
 free/attached in separate lists, CoupleSet keeps all four combinations, so no
 per-object queries are needed.

 A couple's position() is its centre: the diffusing position while free, and the
 midpoint of the two hands once it bridges.
 */
static void savePointFrame(Simul const& simul, const std::string& dir, int frame)
{
    char path[512];
    std::snprintf(path, sizeof(path), "%s/points_%04d.txt", dir.c_str(), frame);
    std::FILE* f = std::fopen(path, "w");
    if ( ! f )
    {
        std::fprintf(stderr, "warning: could not write %s\n", path);
        return;
    }

    std::fprintf(f, "# cytosim points  frame %d  time %g\n", frame, simul.time());

    // group positions by "<kind> <name> <state>"; std::map keeps the block
    // order stable from frame to frame, so Blender materials stay put
    std::map<std::string, std::vector<Vector> > classes;

    /*
     A single has TWO meaningful positions once it is attached, and exporting
     only one hides the thing you usually want to see:

       posFoot() -- the anchor.  For a `fixed' single (Picket) this is the point
                    it was created at, e.g. on the membrane.
       posHand() -- where its hand grips the filament.

     While free the two coincide, so one point is enough and it is emitted as
     `free'.  Once attached they separate by the link stretch, and BOTH are
     emitted -- as `anchor' and `hand' -- so the viewer shows the anchor on the
     membrane and the grip on the filament as two distinct clouds.  With
     single:stiffness = 0 the separation is unbounded, which is exactly the
     failure this export makes visible.
    */
    for ( Single const* S = simul.singles.firstF(); S; S = S->next() )
        classes[std::string("single ") + S->prop->name() + " free"]
            .push_back(S->position());

    for ( Single const* S = simul.singles.firstA(); S; S = S->next() )
    {
        const std::string base = std::string("single ") + S->prop->name();
        classes[base + " anchor"].push_back(S->posFoot());
        classes[base + " hand"].push_back(S->posHand());
    }

    // free-free, one hand bound (either way), and bridging both filaments
    Couple const* heads[4] = { simul.couples.firstFF(), simul.couples.firstAF(),
                               simul.couples.firstFA(), simul.couples.firstAA() };
    const char* states[4] = { "free", "bound", "bound", "bridge" };
    for ( int c = 0; c < 4; ++c )
        for ( Couple const* C = heads[c]; C; C = C->next() )
            classes[std::string("couple ") + C->prop->name() + " " + states[c]]
                .push_back(C->position());

    for ( std::map<std::string, std::vector<Vector> >::const_iterator
          it = classes.begin(); it != classes.end(); ++it )
    {
        std::fprintf(f, "class %s %lu\n", it->first.c_str(),
                     (unsigned long)it->second.size());
        for ( size_t k = 0; k < it->second.size(); ++k )
        {
            Vector p = it->second[k];
#if ( DIM >= 3 )
            std::fprintf(f, "%g %g %g\n", p.XX, p.YY, p.ZZ);
#elif ( DIM == 2 )
            std::fprintf(f, "%g %g %g\n", p.XX, p.YY, 0.0);
#else
            std::fprintf(f, "%g %g %g\n", p.XX, 0.0, 0.0);
#endif
        }
    }
    std::fclose(f);
}

//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    Glossary arg;

    // --- options -------------------------------------------------------------
    std::string config = "config.cym";
    std::string mesh   = "membrane.ply";
    double dt          = 1e-4;
    double nSteps      = 100;
    int    strang      = 1;    ///< 1 = Strang (error ~dt^2), 0 = Lie-Trotter (~dt)
    int    saveEvery   = 10;   ///< export a membrane frame every N steps (0 = off)
    std::string frameDir = "frames";  ///< output directory for the frame sequence
    double membraneKb  = 8.22e-5;     ///< membrane bending modulus
    double Ksg         = 0.0;         ///< stretching modulus (0 = no tension)
    double areaFrac    = 0.98;        ///< preferred area as a fraction of the
                                      ///  initial area; <1 pre-tenses the
                                      ///  membrane, which is what makes it
                                      ///  resist local deformation at all
    double Kv          = 0.0;         ///< osmotic modulus (0 = off)
    double Kse         = 0.0;         ///< edge-spring regularization; the LOCAL
                                      ///  constraint that stops the mesh from
                                      ///  concentrating stretch into a needle
    int    remesh      = 0;           ///< 1 = adaptive remeshing (resolves tubes
                                      ///  properly instead of resisting them)
    int    remeshEvery = 50;          ///< remesh every N steps
    double maxEdge     = 0.0;         ///< split edges longer than this (0 = auto)
    StericParams sp;

    if ( argc > 1 ) config    = argv[1];
    if ( argc > 2 ) mesh      = argv[2];
    if ( argc > 3 ) nSteps    = std::atof(argv[3]);
    if ( argc > 4 ) sp.cutoff = std::atof(argv[4]);
    if ( argc > 5 ) sp.stiff  = std::atof(argv[5]);
    if ( argc > 6 ) strang    = std::atoi(argv[6]);
    if ( argc > 7 ) saveEvery = std::atoi(argv[7]);
    if ( argc > 8 ) frameDir  = argv[8];
    if ( argc > 9 ) sp.mode   = std::atoi(argv[9]);
    if ( argc > 10) membraneKb= std::atof(argv[10]);
    if ( argc > 11) Ksg       = std::atof(argv[11]);
    if ( argc > 12) Kv        = std::atof(argv[12]);
    if ( argc > 13) areaFrac  = std::atof(argv[13]);
    if ( argc > 14) Kse       = std::atof(argv[14]);
    if ( argc > 15) remesh    = std::atoi(argv[15]);
    if ( argc > 16) maxEdge   = std::atof(argv[16]);
    if ( argc > 17) remeshEvery = std::atoi(argv[17]);

    // --- set up Mem3DG -------------------------------------------------------
    mem3dg_system* mem = mem3dg_system_create_from_file(
        mesh.c_str(), membraneKb, membraneKb, /*H0c*/0.0, dt);
    if ( ! mem )
    {
        std::fprintf(stderr, "Error: could not load membrane mesh `%s'\n", mesh.c_str());
        return 1;
    }

    /*
     Material properties.  Without these the membrane has ONLY bending, which is
     a purely local resistance: a filament dents one spot and the rest of the
     surface never learns about it.

     `areaFrac` < 1 sets the preferred area BELOW the initial area, so the
     membrane starts under tension.  That matters more than the modulus itself:
     tension = Ksg*(A-At)/At is a global scalar, and a local dent changes total
     area by ~0.01%, so with At == A0 the tension stays ~0 and nothing pulls
     back.  Pre-tensed, the surface is taut everywhere and responds at once.
     */
    if ( Ksg > 0 )
        mem3dg_system_set_tension(mem, Ksg, areaFrac * mem3dg_system_area(mem));
    if ( Kv > 0 )
        mem3dg_system_set_osmotic(mem, Kv, -1.0);   // defend the initial volume
    /*
     Edge springs.  Tension only sees TOTAL area, so the mesh can dump all its
     stretch into a few edges for free -- which is precisely how a filament
     turns into a needle instead of a bulge.  Kse penalises each edge on its
     own, forcing the deformation to be shared with the neighbours.
     */
    if ( Kse > 0 )
        mem3dg_system_set_regularization(mem, Kse, 0.0, 0.0);

    /*
     Adaptive remeshing.  Where Kse RESISTS the stretch that a growing tube
     causes, remeshing RESOLVES it: edges longer than maxEdge are split, so the
     tube gains vertices as it extends.  The two are alternatives -- Kse
     measures each edge against a reference mesh that remeshing invalidates.
     */
    if ( remesh )
    {
        if ( Kse > 0 )
            std::fprintf(stderr, "warning: remeshing with Kse>0; edge springs "
                                 "reference a mesh that remeshing changes\n");
        // default: split anything much longer than the initial mean edge
        double mx = maxEdge;
        if ( mx <= 0 )
        {
            // mean edge length ~ sqrt(2*A/(sqrt(3)*nF)) for a triangulation
            double A = mem3dg_system_area(mem);
            double nF = (double)mem3dg_system_num_faces(mem);
            mx = 1.5 * std::sqrt(2.0*A/(1.7320508*nF));
        }
        mem3dg_system_set_remeshing(mem, 1, remeshEvery, 0.2*mx, mx);
        std::printf("remeshing: every %d steps, maxEdge=%.4f\n", remeshEvery, mx);
    }

    std::printf("membrane:  Kb=%g  Ksg=%g (At=%.3f*A0)  Kv=%g  Kse=%g\n",
                membraneKb, Ksg, areaFrac, Kv, Kse);

    const size_t nVert = mem3dg_system_num_vertices(mem);
    std::printf("membrane: %zu vertices, %zu faces\n",
                nVert, mem3dg_system_num_faces(mem));
    std::printf("splitting: %s\n",
                strang ? "Strang (error ~dt^2)" : "Lie-Trotter (error ~dt)");
    std::printf("bending:   Kb = %g\n", membraneKb);
    std::printf("steric:    %s\n", sp.mode==1
                ? "one-sided wall (confine inside, push membrane out)"
                : "symmetric repulsion (does not confine)");

    std::vector<double> vertPos(nVert * 3);
    std::vector<double> membraneForces(nVert * 3, 0.0);
    std::vector<FiberForce> fiberForces;
    ContactWork work;   // reusable scratch for the segment arrays

    // where the PLY sequence lands; created if absent (POSIX mkdir, ignore EEXIST)
    if ( saveEvery > 0 )
    {
        std::string cmd = "mkdir -p '" + frameDir + "'";
        if ( std::system(cmd.c_str()) != 0 )
            std::fprintf(stderr, "warning: could not create %s\n", frameDir.c_str());
    }
    int frameIndex = 0;

    // --- set up Cytosim ------------------------------------------------------
    /*
     Hand the membrane to any `shape = mem3dg' Space BEFORE the config is
     parsed: `new couple ... { position = surface }' is resolved during
     parsing, so the Space must already be able to answer project()/inside().

     SpaceMem3DG borrows this handle; the driver keeps ownership.  That way the
     membrane the filaments push on and the one objects are confined by are the
     same object.
     */
    SpaceMem3DG::setHandle(mem);

    Simul simul;
    try {
        simul.initCytosim();
        // Run=0, Write=0: build the objects declared in the config, but do NOT
        // execute its `run` commands -- this program owns the time loop.
        Parser(&simul, 1, 1, 1, 0, 0).readConfig(config);
        simul.prepare();

        /*
         Adopt Cytosim's configured time_step for the membrane too.

         There are two clocks here -- Cytosim's `simul:time_step' and the dt we
         pass to mem3dg_system_step().  If they differ the two halves of the
         split advance by different amounts and the coupling is silently
         inconsistent.  Taking Cytosim's value makes `time_step' in the .cym the
         single knob, so the two can never drift apart.
         */
        if ( simul.time_step() > 0 )
            dt = simul.time_step();
        std::printf("dt:        %g (adopted from simul:time_step)\n", dt);
    }
    catch ( Exception& e ) {
        std::fprintf(stderr, "Cytosim setup failed: %s\n", e.brief().c_str());
        mem3dg_system_destroy(mem);
        return 2;
    }

    /*
     Injected inside Simul::solve_meca(), after Cytosim has registered all of
     its own interactions and before Meca::solve().  Adding here means the
     membrane reaction forces take part in the same implicit solve as
     everything else, rather than being applied as a separate explicit kick.
     */
    simul.externalForceCallback = [&fiberForces](Meca& meca) {
        for ( size_t i = 0; i < fiberForces.size(); ++i )
            meca.addForce(fiberForces[i].mecable,
                          fiberForces[i].point,
                          fiberForces[i].force);
    };

    /*
     --- coupled loop ---------------------------------------------------------

     Two splitting schemes are available.

     LIE-TROTTER (strang = 0), error ~ dt:
         contacts -> Cytosim(dt) -> Mem3DG(dt)
     The two sub-systems simply alternate.  Each one advances a full step using
     the other's position from the START of the step, so the information it acts
     on is one step stale, and that staleness accumulates linearly.

     STRANG (strang = 1), error ~ dt^2:
         contacts -> Mem3DG(dt/2) -> Cytosim(dt) -> recompute contacts -> Mem3DG(dt/2)
     The membrane is advanced half a step, Cytosim takes its full step against
     that half-advanced (i.e. mid-interval) membrane, and the membrane then
     completes its second half.  The scheme is symmetric under time reversal,
     which cancels the leading error term and leaves O(dt^2).

     The extra cost is one additional contact evaluation per step, NOT a second
     Cytosim solve -- Cytosim is still stepped exactly once.  The contacts must
     genuinely be recomputed between the half-steps: reusing the first set
     destroys the symmetry and silently degrades the scheme back to first order.
     */
    for ( int s = 0; s < (int)nSteps; ++s )
    {
        /*
         1+2+3. Contacts on the frozen configuration.

         One call now covers what used to be three steps: it reads the current
         membrane geometry itself, computes the contacts, writes the membrane
         reaction straight into Mem3DG's external force buffer, and hands back
         the filament-side forces.  No vertex positions cross the ABI.
         */
        long nContacts = contactForces(simul, mem, sp, work, fiberForces);

        // 3b. Strang: advance the membrane by half a step FIRST, so that
        //     Cytosim below sees the membrane at the middle of the interval
        if ( strang )
            mem3dg_system_step(mem, 0.5*dt);

        // 4. filament side (fiberForces consumed by the callback)
        try {
            /*
             Both calls are required, in this order.  Cytosim's own run loop
             does exactly this (Interface::step_simul):

                 (sim_->*FUNC)();   // e.g. solve_meca()
                 sim_->steps();

             solve_meca() is what assembles the interactions, fires
             externalForceCallback, solves, and applies the result -- i.e. it is
             the mechanics.  steps() only performs the stochastic/object updates
             (binding, polymerisation, object bookkeeping).

             Calling steps() alone -- as this driver originally did -- leaves the
             mechanics unsolved: the callback never fires and the filaments never
             move under any force.
             */
            simul.solve_meca();
            simul.steps();
        }
        catch ( Exception& e ) {
            std::fprintf(stderr, "Cytosim step failed: %s\n", e.brief().c_str());
            break;
        }

        // 5. membrane
        if ( strang )
        {
            /*
             Second half-step.  Recompute the contacts against the filaments'
             NEW positions first -- reusing the forces from step 2 would break
             the time-reversal symmetry that buys the extra order.
             */
            nContacts = contactForces(simul, mem, sp, work, fiberForces);
            mem3dg_system_step(mem, 0.5*dt);
        }
        else
            mem3dg_system_step(mem, dt);

        // export both sides at the same instant, under the same frame index
        if ( saveEvery > 0 && s % saveEvery == 0 )
        {
            saveMembraneFrame(mem, frameDir, frameIndex);
            saveFiberFrame(simul, frameDir, frameIndex);
            savePointFrame(simul, frameDir, frameIndex);
            ++frameIndex;
        }

        if ( s % 10 == 0 )
        {
            // net magnitudes, to confirm force actually crosses the boundary
            double fFib = 0, fMem = 0;
            for ( size_t i = 0; i < fiberForces.size(); ++i )
                fFib += fiberForces[i].force.norm();

            /*
             Re-query the vertex count: remeshing changes it, and this buffer
             was sized at startup.  Reading into a stale buffer would overflow.
             */
            const size_t nv = mem3dg_system_num_vertices(mem);
            if ( membraneForces.size() != nv*3 )
                membraneForces.assign(nv*3, 0.0);
            mem3dg_system_get_external_forces(mem, &membraneForces[0]);
            for ( size_t i = 0; i + 2 < membraneForces.size(); i += 3 )
                fMem += std::sqrt(membraneForces[i]*membraneForces[i]
                                + membraneForces[i+1]*membraneForces[i+1]
                                + membraneForces[i+2]*membraneForces[i+2]);

            std::printf("step %4d  contacts %5ld  |F_fib| %.4g  |F_mem| %.4g  t_mem %.6g\n",
                        s, nContacts, fFib, fMem, mem3dg_system_time(mem));
        }
    }

    // final frame, so the last state is always captured regardless of cadence
    if ( saveEvery > 0 )
    {
        saveMembraneFrame(mem, frameDir, frameIndex);
        saveFiberFrame(simul, frameDir, frameIndex);
        savePointFrame(simul, frameDir, frameIndex);
        ++frameIndex;
    }
    /*
     No separate "final" file: the block above already captures the last state
     as the last frame of the sequence.  (An earlier version wrote
     membrane_final.ply here, which always produced a 0-byte file -- see the
     note on save_ply in saveMembraneFrame.)
     */
    mem3dg_system_destroy(mem);
    std::printf("done: wrote %d membrane frame(s) to %s/\n",
                saveEvery > 0 ? frameIndex : 0, frameDir.c_str());
    return 0;
}
