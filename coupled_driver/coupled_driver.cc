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

/**


   Cytosim  : (s, um, pN)   
              energy = pN.um = 1e-18 J

   Mem3DG   : lengths are um (the meshes are um-scale), and energies are
              1e-15 J.  That second half is not guessed: pymem3dg/util.py's
              femtoJToKBT() documents the energy unit as femtojoules, and it is
              confirmed by the canonical Kb = 8.22e-5 used throughout the Mem3DG
              tests --

                  8.22e-5 * 1e-15 J = 8.22e-20 J
                  kT (298 K)        = 1.38e-23 * 298 = 4.11e-21 J
                  ratio             = 20.0  -> exactly 20 kT, as intended.

              So Mem3DG's force unit is 1e-15 J / 1e-6 m = 1e-9 N = nN.

 */
static const double PN_PER_M3F = 1000.0;   ///< pN per Mem3DG force unit (nN)


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
         Confinement works per POINT.
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
            mem, &w.ptP[0], nPts, /*shell=*/p.cutoff, p.stiff / PN_PER_M3F,
            &w.ptF[0], /*add_to_membrane=*/1);

        if ( n < 0 )
        {
            std::fprintf(stderr,
                "error: filament_confinement failed (needs a CLOSED mesh)\n");
            return 0;
        }

        for ( size_t i = 0; i < nPts; ++i )
        {
            const Vector f(w.ptF[3*i+0] * PN_PER_M3F,
                           w.ptF[3*i+1] * PN_PER_M3F,
                           w.ptF[3*i+2] * PN_PER_M3F);
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

    /*
     Same pN/um -> Mem3DG conversion as the confinement branch above.
     */
    long nContacts = mem3dg_system_filament_repulsion(
        mem, &w.segP0[0], &w.segP1[0], nSeg,
        p.cutoff, p.stiff / PN_PER_M3F,
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
/**
 Total magnitude of the membrane-side external force, in pN.

 Must be called while the external force buffer still holds the current
 contact reaction -- i.e. after contactForces() and BEFORE the membrane step,
 which zeroes the buffer whenever it remeshes.

 `scratch` is resized as needed: remeshing changes the vertex count, and reading
 into a buffer sized at startup would overflow.
 */
static double sumExternalForces(mem3dg_system* mem, std::vector<double>& scratch)
{
    const size_t nv = mem3dg_system_num_vertices(mem);
    if ( scratch.size() != nv*3 )
        scratch.assign(nv*3, 0.0);
    mem3dg_system_get_external_forces(mem, &scratch[0]);

    double s = 0;
    for ( size_t i = 0; i + 2 < scratch.size(); i += 3 )
        s += std::sqrt(scratch[i]*scratch[i] + scratch[i+1]*scratch[i+1]
                     + scratch[i+2]*scratch[i+2]);
    return s * PN_PER_M3F;      // Mem3DG force units -> pN
}


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
    int    memEvery    = 1;           ///< advance the membrane once every N
                                      ///  Cytosim steps, by N*dt.  1 = every
                                      ///  step (the original behaviour)
    int    anchorFollow= 0;           ///< 1 = re-project anchored single feet
                                      ///  onto the membrane every step, so they
                                      ///  track it instead of staying pinned in
                                      ///  the lab frame
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
    if ( argc > 18) anchorFollow= std::atoi(argv[18]);
    if ( argc > 19) memEvery    = std::atoi(argv[19]);

    if ( memEvery < 1 )
        memEvery = 1;

    // --- set up Mem3DG -------------------------------------------------------
    mem3dg_system* mem = mem3dg_system_create_from_file(
        mesh.c_str(), membraneKb, membraneKb, /*H0c*/0.0, dt);
    if ( ! mem )
    {
        std::fprintf(stderr, "Error: could not load membrane mesh `%s'\n", mesh.c_str());
        return 1;
    }

    if ( Ksg > 0 )
        mem3dg_system_set_tension(mem, Ksg, areaFrac * mem3dg_system_area(mem));
    if ( Kv > 0 )
        mem3dg_system_set_osmotic(mem, Kv, -1.0);   // defend the initial volume
    if ( Kse > 0 )
        mem3dg_system_set_regularization(mem, Kse, 0.0, 0.0);

    if ( remesh )
    {
        if ( Kse > 0 )
            std::fprintf(stderr, "warning: remeshing with Kse>0; edge springs "
                                 "reference a mesh that remeshing changes\n");
        // default: split anything much longer than the initial mean edge
        double mx = maxEdge;
        if ( mx <= 0 )
        {

            double A = mem3dg_system_area(mem);
            double nF = (double)mem3dg_system_num_faces(mem);
            mx = 1.5 * std::sqrt(4.0*A/(1.7320508*nF));
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
    std::printf("bending:   Kb = %g  (= %.3g kT at 298 K)\n",
                membraneKb, membraneKb * 1e-15 / (1.38e-23 * 298.0));
    std::printf("units:     Cytosim (s,um,pN); Mem3DG (s,um,%g pN); "
                "steric stiffness read as pN/um\n", PN_PER_M3F);
    std::printf("steric:    %s\n", sp.mode==1
                ? "one-sided wall (confine inside, push membrane out)"
                : "symmetric repulsion (does not confine)");

    std::vector<double> vertPos(nVert * 3);
    std::vector<double> membraneForces(nVert * 3, 0.0);
    std::vector<FiberForce> fiberForces;
    ContactWork work;   // reusable scratch for the segment arrays

    if ( saveEvery > 0 )
    {
        std::string cmd = "mkdir -p '" + frameDir + "'";
        if ( std::system(cmd.c_str()) != 0 )
            std::fprintf(stderr, "warning: could not create %s\n", frameDir.c_str());
    }
    int frameIndex = 0;

    // --- set up Cytosim ------------------------------------------------------
    SpaceMem3DG::setHandle(mem);

    Simul simul;
    try {
        simul.initCytosim();
        Parser(&simul, 1, 1, 1, 0, 0).readConfig(config);
        simul.prepare();

        if ( simul.time_step() > 0 )
            dt = simul.time_step();
        std::printf("dt:        %g (adopted from simul:time_step)\n", dt);

        if ( memEvery > 1 )
        {
            std::printf("membrane:  sub-stepped every %d Cytosim steps "
                        "(membrane dt = %g); contacts still every step\n",
                        memEvery, dt*memEvery);
            if ( dt*memEvery > 5e-5 )
                std::fprintf(stderr, "warning: membrane dt %g exceeds 5e-5, the "
                             "largest value this project has run stably; "
                             "check for membrane instability\n", dt*memEvery);
        }
    }
    catch ( Exception& e ) {
        std::fprintf(stderr, "Cytosim setup failed: %s\n", e.brief().c_str());
        mem3dg_system_destroy(mem);
        return 2;
    }

    simul.externalForceCallback = [&fiberForces](Meca& meca) {
        for ( size_t i = 0; i < fiberForces.size(); ++i )
            meca.addForce(fiberForces[i].mecable,
                          fiberForces[i].point,
                          fiberForces[i].force);
    };

    for ( int s = 0; s < (int)nSteps; ++s )
    {
        
        long nContacts = contactForces(simul, mem, sp, work, fiberForces);

        const bool stepMembrane = ( (s % memEvery) == 0 );
        const double memDt = dt * memEvery;

        double fMem = 0;
        if ( s % 10 == 0 )
            fMem = sumExternalForces(mem, membraneForces);

        // 3b. Strang: advance the membrane by half a step FIRST, so that
        //     Cytosim below sees the membrane at the middle of the interval
        if ( strang && stepMembrane )
            mem3dg_system_step(mem, 0.5*memDt);

        // 4. filament side (fiberForces consumed by the callback)
        try {
            simul.solve_meca();
            simul.steps();
        }
        catch ( Exception& e ) {
            std::fprintf(stderr, "Cytosim step failed: %s\n", e.brief().c_str());
            break;
        }

        // 5. membrane
        if ( stepMembrane )
        {
            if ( strang )
            {
                nContacts = contactForces(simul, mem, sp, work, fiberForces);
                if ( s % 10 == 0 )
                    fMem = sumExternalForces(mem, membraneForces);
                mem3dg_system_step(mem, 0.5*memDt);
            }
            else
                mem3dg_system_step(mem, memDt);
        }

        if ( anchorFollow )
        {
            for ( int pass = 0; pass < 2; ++pass )
            {
                Single * sgl = pass ? simul.singles.firstA() : simul.singles.firstF();
                for ( ; sgl; sgl = sgl->next() )
                {
                    if ( ! sgl->hasLink() )
                        continue;
                    Vector p = sgl->position();
                    double prj[3];
                    if ( 0 == mem3dg_system_project(mem, p.XX, p.YY, p.ZZ, prj, nullptr) )
                        sgl->setPosition(Vector(prj[0], prj[1], prj[2]));
                }
            }
        }

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
            double fFib = 0;
            for ( size_t i = 0; i < fiberForces.size(); ++i )
                fFib += fiberForces[i].force.norm();

            std::printf("step %4d  contacts %5ld  |F_fib| %.4g pN  |F_mem| %.4g pN  t_mem %.6g\n",
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

    mem3dg_system_destroy(mem);
    std::printf("done: wrote %d membrane frame(s) to %s/\n",
                saveEvery > 0 ? frameIndex : 0, frameDir.c_str());
    return 0;
}
