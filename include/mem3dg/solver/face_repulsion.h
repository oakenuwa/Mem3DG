#pragma once

#include "mem3dg/macros.h"

#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/utilities/vector3.h>

#include <array>
#include <functional>
#include <limits>

// =============================================================================
// Face-to-face and edge-to-face repulsion via discrete differential geometry.
//
// MOTIVATION
// ──────────
// Vertex-vertex repulsion (existing computeSelfAvoidanceForce) misses tunneling:
// a thin triangle can pierce another face even when all six vertices are far apart.
// Face-level primitives close this gap.
//
// APPROACH
// ────────
// For each non-adjacent face pair (Fi, Fj), the minimum Euclidean distance
// between the two triangles is attained by exactly one of two primitive cases:
//
//   ① Vertex–Face (VF):  a vertex of Fi (or Fj) is closest to the interior
//                         or boundary of the opposing face.
//   ② Edge–Edge   (EE):  an edge of Fi and an edge of Fj have closest-point
//                         pair that lies in the interior of both segments.
//
// A face–face pair has 6 VF sub-cases and 9 EE sub-cases; the minimum over all
// 15 is the triangle–triangle distance.
//
// FORCE DECOMPOSITION (DDG / virtual work)
// ─────────────────────────────────────────
// Given potential E(d) with derivative E'(d), and the closest-point pair
// (pA on face A, pB on face B) at distance d = ‖pA − pB‖:
//
//   n̂  =  (pA − pB) / d          (unit normal, from B toward A)
//
//   Force on vertex k of face A  =  −E'(d) · wA[k] · n̂
//   Force on vertex k of face B  =  +E'(d) · wB[k] · n̂
//
//   where wA[k], wB[k] are the barycentric coordinates of pA and pB in their
//   respective triangles.  Newton's 3rd law is satisfied automatically since
//   ΣwA[k] = ΣwB[k] = 1.
//
//   For repulsion E'(d) < 0, so forces on A point in +n̂ (away from B) ✓
//   For repulsion E'(d) < 0, forces on B point in −n̂ (away from A)   ✓
//
// See face_repulsion.cpp for the derivation of each barycentric weight.
// =============================================================================

namespace mem3dg {
namespace solver {

namespace gc  = geometrycentral;
namespace gcs = geometrycentral::surface;

// ---------------------------------------------------------------------------
// Contact descriptor returned by minimumFacePairContact().
//
// pA / pB    : closest world-space points on face A and face B.
// wA / wB    : barycentric coordinates within {v0, v1, v2} of face A / face B.
//              Layout follows the halfedge order: he, he.next(), he.next().next().
// distance   : ‖pA − pB‖.  Set to ∞ if the pair should be skipped.
// ---------------------------------------------------------------------------
struct FacePairContact {
    gc::Vector3          pA;
    gc::Vector3          pB;
    double               distance = std::numeric_limits<double>::infinity();
    std::array<double,3> wA{};   // barycentric coords of pA in face A
    std::array<double,3> wB{};   // barycentric coords of pB in face B
};

// ---------------------------------------------------------------------------
// Closest-point helpers (implemented in face_repulsion.cpp).
// These are standard computational geometry algorithms; their only purpose
// here is to produce barycentric coordinates suitable for DDG decomposition.
// ---------------------------------------------------------------------------

// Closest point on triangle {a,b,c} to query point p.
// Returns (closest_world_point, {λa, λb, λc}) using Voronoi-region projection.
// The barycentric coordinates are always non-negative and sum to 1.
DLL_PUBLIC std::pair<gc::Vector3, std::array<double,3>>
closestPointOnTriangle(gc::Vector3 a, gc::Vector3 b, gc::Vector3 c,
                       gc::Vector3 p);

// Closest points between segment [e0,e1] and segment [g0,g1].
// Returns (s*, t*) the parametric coordinates on each segment (both ∈ [0,1]).
// The closest points are  pe = (1−s)·e0 + s·e1  and  pg = (1−t)·g0 + t·g1.
DLL_PUBLIC std::pair<double,double>
closestSegmentSegment(gc::Vector3 e0, gc::Vector3 e1,
                      gc::Vector3 g0, gc::Vector3 g1);

// Minimum-distance contact between two triangles A = {a0,a1,a2} and
// B = {b0,b1,b2}, evaluated over all 6 VF + 9 EE sub-cases.
// Returns the FacePairContact with the smallest distance.
DLL_PUBLIC FacePairContact
minimumFacePairContact(const std::array<gc::Vector3,3>& A,
                       const std::array<gc::Vector3,3>& B);

} // namespace solver
} // namespace mem3dg
