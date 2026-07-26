/* ===========================================================================
 * mem3dg_c.h  —  Plain-C ABI facade for Mem3DG.
 *
 * PURPOSE
 * -------
 * Lets a C or C++ program (e.g. Cytosim, compiled as C++11) drive a Mem3DG
 * membrane WITHOUT including any geometry-central / Eigen headers and WITHOUT
 * compiling as C++20.  Everything Mem3DG-specific stays sealed inside
 * libmem3dg.dylib / libmem3dg.so.
 *
 * Only opaque handles and plain-old-data cross this boundary, so:
 *   - the consumer's C++ standard is irrelevant (this is a C linkage layer);
 *   - no template / STL types are shared across the ABI (no ODR hazard);
 *   - the consumer needs exactly ONE thing on its link line: the Mem3DG shared
 *     library.  geometry-central is statically baked into that library.
 *
 * MEMORY / OWNERSHIP
 * ------------------
 * The library owns the System.  The caller owns any buffer it passes in.
 * "get" functions copy into caller-provided buffers sized per the nVertices /
 * nFaces queries.  Nothing returned by this API needs to be freed except the
 * handle itself, via mem3dg_system_destroy().
 *
 * THREADING
 * ---------
 * A handle is not thread-safe; serialize calls on a given handle.
 * ==========================================================================*/
#ifndef MEM3DG_C_H
#define MEM3DG_C_H

#include <stddef.h>

#ifdef _WIN32
#define MEM3DG_C_API __declspec(dllexport)
#else
#define MEM3DG_C_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handle to a Mem3DG membrane (System + integrator + force buffer). */
typedef struct mem3dg_system mem3dg_system;

/* ---------------------------------------------------------------------------
 * Lifecycle
 * ------------------------------------------------------------------------- */

/* Create a membrane from a surface mesh file (.ply/.obj) with a harmonic
 * bending model.  Returns NULL on failure (bad path, non-manifold mesh, ...).
 *   Kb  : bending modulus
 *   Kbc : (spatially-uniform) constrained bending modulus
 *   H0c : spontaneous curvature
 *   characteristicTimeStep : fixed dt used by mem3dg_system_step()
 */
MEM3DG_C_API mem3dg_system *
mem3dg_system_create_from_file(const char *mesh_path, double Kb, double Kbc,
                               double H0c, double characteristicTimeStep);

/* Create a membrane from raw arrays instead of a file.
 *   vertices : nV * 3 doubles, row-major (x0,y0,z0, x1,y1,z1, ...)
 *   faces    : nF * 3 size_t,  row-major vertex indices (triangles only)
 */
MEM3DG_C_API mem3dg_system *mem3dg_system_create_from_arrays(
    const double *vertices, size_t nV, const size_t *faces, size_t nF,
    double Kb, double Kbc, double H0c, double characteristicTimeStep);

/* Destroy a handle previously returned by a create function.  Safe on NULL. */
MEM3DG_C_API void mem3dg_system_destroy(mem3dg_system *h);

/* ---------------------------------------------------------------------------
 * Geometry queries  (for the caller to build its own steric contacts)
 * ------------------------------------------------------------------------- */

MEM3DG_C_API size_t mem3dg_system_num_vertices(const mem3dg_system *h);
MEM3DG_C_API size_t mem3dg_system_num_faces(const mem3dg_system *h);

/* Copy current vertex positions into out_xyz (must hold 3*nVertices doubles),
 * row-major (x0,y0,z0, x1,y1,z1, ...). */
MEM3DG_C_API void mem3dg_system_get_vertex_positions(const mem3dg_system *h,
                                                     double *out_xyz);

/* Copy triangle connectivity into out_ijk (must hold 3*nFaces size_t),
 * row-major vertex indices.  Connectivity is constant unless the mesh mutates. */
MEM3DG_C_API void mem3dg_system_get_faces(const mem3dg_system *h,
                                          size_t *out_ijk);

/* Overwrite the vertex positions (3*nVertices doubles, row-major).  Intended
 * for diagnostics: impose a known deformation, then read the forces back with
 * mem3dg_system_get_mechanical_forces() to see how the membrane responds. */
MEM3DG_C_API void
mem3dg_system_set_vertex_positions(mem3dg_system *h, const double *xyz);

/* Recompute the membrane's own forcing (bending + tension + osmotic + ...) for
 * the CURRENT geometry, then copy the per-vertex total into out_xyz
 * (3*nVertices doubles).  Does NOT step, and does NOT include the external
 * force buffer -- this is the membrane's internal response alone. */
MEM3DG_C_API void
mem3dg_system_get_mechanical_forces(mem3dg_system *h, double *out_xyz);

/* ---------------------------------------------------------------------------
 * Membrane material properties
 *
 * A membrane created by the functions above has ONLY bending stiffness.  That
 * is enough to smoke-test the coupling but it is not a membrane: bending is a
 * purely LOCAL resistance to curvature, so a filament pressing on one spot
 * deforms that spot and nothing else -- neighbouring faces have no way to feel
 * it and the surface cannot relax.
 *
 * Surface tension and osmotic pressure are what couple distant parts of the
 * surface: tension makes any local stretch cost area globally, and osmotic
 * pressure gives a closed vesicle a volume it defends.  Enable at least one for
 * physically meaningful relaxation.
 * ------------------------------------------------------------------------- */

/* Preferred-area surface tension:
 *
 *     tension(A) = Ksg * (A - At) / At
 *     energy (A) = Ksg * (A - At)^2 / (2*At)
 *
 * Ksg              : stretching modulus; 0 DISABLES surface tension
 * preferred_area   : At.  Pass <= 0 to use the membrane's CURRENT area, i.e.
 *                    it defends the shape it was loaded with.
 */
MEM3DG_C_API void mem3dg_system_set_tension(mem3dg_system *h, double Ksg,
                                            double preferred_area);

/* Preferred-volume osmotic pressure:
 *
 *     pressure(V) = Kv * (Vt - V) / Vt      (negative, i.e. inward, when V > Vt)
 *     energy  (V) = Kv * (V - Vt)^2 / (2*Vt)
 *
 * Kv               : osmotic modulus; 0 DISABLES osmotic pressure
 * preferred_volume : Vt.  Pass <= 0 to use the membrane's CURRENT volume.
 *
 * Only meaningful for a closed surface.
 */
MEM3DG_C_API void mem3dg_system_set_osmotic(mem3dg_system *h, double Kv,
                                            double preferred_volume);

/* Mesh regularization: springs that hold each EDGE, FACE AREA and TRIANGLE
 * SHAPE near its value in the reference (initial) mesh.
 *
 * WHY THIS MATTERS -- and why surface tension is not a substitute.
 * Tension is a GLOBAL constraint: it penalises the total area only.  The mesh
 * can therefore concentrate an arbitrary amount of stretch into a handful of
 * edges at almost no cost, as long as the total is preserved.  In practice a
 * filament pushing on one face drags those few vertices out into a NEEDLE with
 * edges stretched several-fold, while the mean edge length is unchanged and
 * tension registers nothing wrong.
 *
 * Kse penalises each edge individually, so stretch has to be shared with the
 * neighbours -- which is what turns a needle into a smooth bulge.
 *
 *   Kse : edge spring    (per-edge length, relative to the reference mesh)
 *   Ksl : local area     (per-face area)
 *   Kst : triangle ratio (per-triangle shape / aspect)
 *
 * Pass 0 for any of them to leave that term off.
 */
MEM3DG_C_API void mem3dg_system_set_regularization(mem3dg_system *h, double Kse,
                                                   double Ksl, double Kst);

/* Adaptive remeshing.
 *
 * WHY: a filament pushing out a tube stretches the triangles it drags along,
 * because a fixed mesh has no material to add.  Edge springs (Kse) suppress
 * that by resisting the stretch, but they do so by resisting the TUBE -- they
 * trade an artefact for a stiffer membrane.  Remeshing resolves it properly:
 * long edges are split, so the tube gains vertices as it grows and the
 * triangles stay well shaped.
 *
 * enable   : 0 turns all mesh mutation off (the default)
 * period   : mutate every `period` calls to mem3dg_system_step(); 1 = every step
 * min_edge : edges shorter than this are collapsed
 * max_edge : edges longer than this are split -- this is the knob that actually
 *            resolves a growing tube, so set it near the mesh's initial edge
 *            length
 *
 * IMPORTANT: remeshing CHANGES THE VERTEX AND FACE COUNT.  After any step that
 * mutates, re-query mem3dg_system_num_vertices() / num_faces() and resize any
 * buffers before calling the get/set functions again.  The internal external
 * force buffer is resized automatically, but it is also CLEARED, because the
 * old per-vertex forces no longer correspond to anything.
 *
 * NOTE: combining remeshing with edge springs (mem3dg_system_set_regularization)
 * is not recommended -- Kse measures each edge against a REFERENCE mesh, and
 * remeshing invalidates that reference for every new edge.  Use one or the
 * other.
 */
MEM3DG_C_API void mem3dg_system_set_remeshing(mem3dg_system *h, int enable,
                                              int period, double min_edge,
                                              double max_edge);

/* Current total surface area and enclosed volume (diagnostics). */
MEM3DG_C_API double mem3dg_system_area(const mem3dg_system *h);
MEM3DG_C_API double mem3dg_system_volume(const mem3dg_system *h);

/* ---------------------------------------------------------------------------
 * Coupling: inject external forces, then advance
 * ------------------------------------------------------------------------- */

/* Overwrite the per-vertex external force buffer (e.g. filament->membrane
 * steric forces the caller computed).  force_xyz holds 3*nVertices doubles,
 * row-major.  The buffer persists until overwritten; it is applied on every
 * subsequent step().  Pass all-zeros (or call mem3dg_system_clear_external_forces)
 * to stop injecting. */
MEM3DG_C_API void
mem3dg_system_set_external_forces(mem3dg_system *h, const double *force_xyz);

MEM3DG_C_API void mem3dg_system_clear_external_forces(mem3dg_system *h);

/* Add to the per-vertex external force buffer instead of overwriting it, so
 * several force sources can be combined before stepping. */
MEM3DG_C_API void
mem3dg_system_add_external_forces(mem3dg_system *h, const double *force_xyz);

/* Read the per-vertex external force buffer back into out_xyz (3*nVertices
 * doubles).  Useful for diagnostics when the buffer was filled internally, e.g.
 * by mem3dg_system_filament_repulsion(). */
MEM3DG_C_API void
mem3dg_system_get_external_forces(const mem3dg_system *h, double *out_xyz);

/* ---------------------------------------------------------------------------
 * Filament-membrane steric repulsion
 * ------------------------------------------------------------------------- */

/* Compute excluded-volume repulsion between a set of filament segments and the
 * membrane surface, using the edge-face kernel inside Mem3DG.
 *
 * WHY THIS EXISTS
 * ---------------
 * The kernel is C++20 and speaks Eigen, so a C++11 host (Cytosim) cannot call
 * it directly.  This entry point takes and returns plain doubles, keeping the
 * ABI clean while letting the caller use the real contact geometry instead of
 * reimplementing it.
 *
 * WHAT IT COMPUTES
 * ----------------
 * For every (segment, membrane triangle) pair it finds the TRUE minimum
 * distance -- including the case where the segment pierces the triangle, which
 * a vertex-to-vertex test misses entirely and which is how filaments tunnel
 * through a membrane.  Pairs closer than `cutoff` get a harmonic repulsion
 *
 *     |F| = stiffness * (cutoff - distance)
 *
 * decomposed with DDG barycentric weights onto the two segment endpoints and
 * the three face vertices.  Newton's third law holds pair-by-pair: the total
 * force on the filament side is exactly minus the total on the membrane side.
 *
 * PARAMETERS
 * ----------
 * seg_p0, seg_p1  : nSeg*3 doubles each, row-major xyz; segment k runs from
 *                   seg_p0[3k..3k+2] to seg_p1[3k..3k+2]
 * nSeg            : number of segments
 * cutoff          : interaction range (same length units as the mesh)
 * stiffness       : harmonic stiffness
 * out_f0, out_f1  : nSeg*3 doubles each, receive the force on each segment's
 *                   p0 and p1 endpoints.  OVERWRITTEN (not accumulated).
 *                   Either may be NULL to skip that output.
 * add_to_membrane : if nonzero, the equal-and-opposite membrane reaction is
 *                   ACCUMULATED into the external force buffer, ready for the
 *                   next mem3dg_system_step().  Call
 *                   mem3dg_system_clear_external_forces() first if the buffer
 *                   should hold only this contribution.
 *
 * Returns the number of contacting (segment, face) pairs, or -1 on error.
 */
MEM3DG_C_API long mem3dg_system_filament_repulsion(
    mem3dg_system *h, const double *seg_p0, const double *seg_p1, size_t nSeg,
    double cutoff, double stiffness, double *out_f0, double *out_f1,
    int add_to_membrane);

/* Is the point (x,y,z) inside the membrane?
 *
 * Ray casting: shoot a ray and count surface crossings -- odd means inside.
 * This is exact for any closed, oriented manifold, including strongly
 * non-convex shapes, which is why it is used instead of testing the sign of
 * the nearest face's normal.  That local test is only reliable when the
 * closest point lies in a face INTERIOR; on an edge or vertex, or when the
 * surface folds back on itself, it can report the wrong side.
 *
 * Returns 1 inside, 0 outside, -1 on error (null handle or open mesh).
 * Cost is O(nFaces) per query.
 */
MEM3DG_C_API int mem3dg_system_inside(const mem3dg_system *h, double x, double y,
                                      double z);

/* Closest point on the membrane surface to (x,y,z).
 *
 * out_point  : the closest surface point            (3 doubles, may be NULL)
 * out_normal : OUTWARD unit normal there            (3 doubles, may be NULL)
 *
 * The normal is the angle-weighted pseudonormal, so
 * dot(query - out_point, out_normal) is a correct SIGNED distance even when the
 * closest point lands on an edge or a vertex.
 *
 * Accelerated by a KD-tree over the vertices, so this is not O(nFaces).
 * Returns 0 on success, -1 on error.
 */
MEM3DG_C_API int mem3dg_system_project(const mem3dg_system *h, double x,
                                       double y, double z, double *out_point,
                                       double *out_normal);

/* Axis-aligned bounding box of the membrane (3 doubles each). */
MEM3DG_C_API void mem3dg_system_bounds(const mem3dg_system *h, double *out_inf,
                                       double *out_sup);

/* Confine filaments INSIDE a closed membrane, and push the membrane outward.
 *
 * Takes filament POINTS, not segments -- one entry per filament vertex.
 *
 * HOW THIS DIFFERS FROM mem3dg_system_filament_repulsion()
 * --------------------------------------------------------
 * That function is SYMMETRIC: it pushes a segment away from the surface on
 * whichever side the segment already lies, and only within `cutoff`.  It
 * therefore cannot confine -- a point that drifts outside is pushed further
 * out, and beyond `cutoff` it feels nothing and never returns.
 *
 * This function is a ONE-SIDED WALL.  For each point it finds the closest
 * point on the surface, determines the side by RAY CASTING (see
 * mem3dg_system_inside -- the sign of the nearest face's normal is unreliable
 * when the closest point lands on an edge or vertex, or where the surface
 * folds), and forms the signed distance `s` (positive outside).  For
 * s > -shell it applies
 *
 *     force on filament =  -stiffness * (s + shell) * n_out   (inward)
 *     force on membrane =  +stiffness * (s + shell) * n_out   (outward)
 *
 * so the wall sits at s = -shell.  Deeper inside there is no force.  The
 * magnitude grows with penetration and is NOT capped by a cutoff, so a point
 * that gets outside is always pulled back.  The membrane reaction is spread
 * over the closest face's vertices with DDG barycentric weights, so Newton's
 * third law holds point-by-point.
 *
 * WHY POINTS RATHER THAN SEGMENTS
 * -------------------------------
 * An earlier segment-based version evaluated the wall once per segment, at the
 * segment's closest point to the surface.  For a segment straddling the
 * membrane that point lies essentially ON the surface, so the reported
 * penetration was only `shell` and the restoring force saturated at
 * stiffness*shell however far the outer endpoint had escaped (measured 17x too
 * weak).  Weighting each endpoint by 1/2, to compensate for consecutive
 * segments sharing interior points, then left each fiber's two TERMINAL points
 * at half force -- and those are exactly the points that lead growth into the
 * membrane.  Passing points directly removes both problems.
 *
 * PARAMETERS
 * ----------
 * points          : nPts*3 doubles, row-major xyz
 * nPts            : number of points
 * shell           : clearance held between filament and surface (>= 0); the
 *                   wall is this far inside the membrane
 * stiffness       : force per unit penetration
 * out_f           : nPts*3 doubles, force on each point.  OVERWRITTEN.  May be
 *                   NULL.
 * add_to_membrane : nonzero -> accumulate the outward reaction into the
 *                   external force buffer
 *
 * REQUIRES A CLOSED MESH; returns -1 if the membrane has boundary edges.
 * Returns the number of points that received a force, or -1 on error.
 */
MEM3DG_C_API long mem3dg_system_filament_confinement(
    mem3dg_system *h, const double *points, size_t nPts, double shell,
    double stiffness, double *out_f, int add_to_membrane);

/* Advance the membrane by one dt:
 *   recompute conservative + nonconservative forcing,
 *   add the external force buffer,
 *   integrate one forward-Euler step.
 * dt <= 0 uses the characteristicTimeStep set at creation. */
MEM3DG_C_API void mem3dg_system_step(mem3dg_system *h, double dt);

/* Current simulation time of the membrane. */
MEM3DG_C_API double mem3dg_system_time(const mem3dg_system *h);

/* ---------------------------------------------------------------------------
 * Output (optional)
 * ------------------------------------------------------------------------- */

/* Write the current mesh to a .ply file. Returns 0 on success, nonzero on error. */
MEM3DG_C_API int mem3dg_system_save_ply(mem3dg_system *h, const char *out_path);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* MEM3DG_C_H */
