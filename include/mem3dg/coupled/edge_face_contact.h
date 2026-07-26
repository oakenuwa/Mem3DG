#pragma once

// =============================================================================
// edge_face_contact.h
//
// Minimum-distance contact between a filament EDGE (line segment between two
// adjacent Cytosim fiber points) and a membrane FACE (triangle of the Mem3DG
// surface mesh).
//
// WHY edge-face instead of vertex-segment?
// ─────────────────────────────────────────
// The vertex-based StericKernel (steric_kernel.h) computes, for each membrane
// vertex, the closest point on the nearest filament segment.  This works well
// when the filament is thin relative to the mesh resolution, but misses the
// case where a filament segment crosses a large membrane triangle without
// passing near any of the triangle's three vertices.
//
// The edge-face kernel fixes this by evaluating the MINIMUM distance between
// the filament segment (as a geometric edge) and the whole membrane triangle
// (as a geometric face).  That minimum is attained by one of five primitive
// sub-cases:
//
//   ① Endpoint p0 of segment  vs  membrane triangle   (vertex-face)
//   ② Endpoint p1 of segment  vs  membrane triangle   (vertex-face)
//   ③ Filament segment        vs  face edge v0–v1     (edge-edge)
//   ④ Filament segment        vs  face edge v1–v2     (edge-edge)
//   ⑤ Filament segment        vs  face edge v2–v0     (edge-edge)
//
// DDG FORCE DECOMPOSITION
// ────────────────────────
// Given the closest-point pair (pSeg on filament, pFace on membrane face):
//
//   Filament side — split by parametric coord s on the segment:
//     force_on_p0 = (1−s) · force_on_segment
//     force_on_p1 =    s  · force_on_segment
//
//   Membrane side — split by barycentric coords {w0,w1,w2} of pFace:
//     force_on_face_vertex_k = wFace[k] · force_on_face
//
// Newton's 3rd law:  force_on_face = −force_on_segment  (enforced by kernel).
// See steric_kernel_edge_face.h for the kernel stub where you fill in the
// actual force magnitude and direction.
// =============================================================================

#include "mem3dg/macros.h"

#include <Eigen/Dense>
#include <array>
#include <limits>

namespace coupled {

// ---------------------------------------------------------------------------
// Contact descriptor for one (filament_segment, membrane_face) pair.
//
// pSeg   : closest point on the filament segment in world space.
// pFace  : closest point on the membrane face in world space.
// s      : parametric coord on segment  (s=0 → p0, s=1 → p1).
// wFace  : barycentric coords of pFace in face's vertex order {v0, v1, v2}
//           — follows the halfedge order: he, he.next(), he.next().next().
// distance : ‖pSeg − pFace‖  (set to ∞ if pair is invalid / beyond cutoff).
// ---------------------------------------------------------------------------
struct EdgeFaceContact {
    Eigen::Vector3d      pSeg;
    Eigen::Vector3d      pFace;
    double               s        = 0.0;
    std::array<double,3> wFace    = {};
    double               distance = std::numeric_limits<double>::infinity();
};

// ---------------------------------------------------------------------------
// Computes the minimum-distance EdgeFaceContact between filament segment
// [p0, p1] and membrane triangle {fv[0], fv[1], fv[2]}.
//
// Evaluates all 5 primitive sub-cases (2 vertex-face + 3 edge-edge) and
// returns the one with the smallest distance.
//
// The returned s and wFace are barycentric coords ready for DDG decomposition.
// ---------------------------------------------------------------------------
// DLL_PUBLIC is required: this project compiles with
// CMAKE_CXX_VISIBILITY_PRESET=hidden, so without it the symbol is emitted as
// local and cannot be linked against from outside libmem3dg.
DLL_PUBLIC EdgeFaceContact
minimumEdgeFaceContact(const Eigen::Vector3d&               p0,
                       const Eigen::Vector3d&               p1,
                       const std::array<Eigen::Vector3d,3>& fv);

} // namespace coupled
