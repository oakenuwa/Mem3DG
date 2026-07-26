// =============================================================================
// face_repulsion.cpp
//
// Implements face-to-face and edge-to-face excluded-volume repulsion for
// Mem3DG membranes, with DDG-consistent force decomposition to vertices.
//
// Call System::computeFaceFaceRepulsionForce() from computeConservativeForcing()
// after the existing self-avoidance block (or instead of it, if replacing).
// =============================================================================

#include "mem3dg/solver/face_repulsion.h"
#include "mem3dg/solver/system.h"

#include <geometrycentral/utilities/vector3.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace mem3dg {
namespace solver {

namespace gc  = geometrycentral;
namespace gcs = geometrycentral::surface;

// ═════════════════════════════════════════════════════════════════════════════
// §1  GEOMETRY HELPERS
// ═════════════════════════════════════════════════════════════════════════════

// -----------------------------------------------------------------------------
// §1.1  closestPointOnTriangle
//
// Given triangle {a, b, c} and point p, returns the closest point q on the
// (filled) triangle and its barycentric coordinates (λa, λb, λc) such that
//   q = λa·a + λb·b + λc·c,   λa+λb+λc = 1,   λi ≥ 0.
//
// Algorithm: Voronoi-region projection (Ericson, "Real-Time Collision Detection",
// §5.1.5).  Each of the 7 Voronoi regions (3 vertex, 3 edge, 1 interior)
// corresponds to a different analytic formula for the barycentric coords.
//
// WHY barycentric coordinates matter:
//   The DDG force decomposition requires ∂(closest_point)/∂vertex_i = λi · I,
//   which holds exactly when the closest point is expressed as a convex
//   combination of the triangle's vertices.  Any motion of vertex i shifts the
//   closest point by λi times that motion (to first order), so the chain rule
//   gives force_on_vertex_i = λi · (force_at_closest_point).
// -----------------------------------------------------------------------------
std::pair<gc::Vector3, std::array<double,3>>
closestPointOnTriangle(gc::Vector3 a, gc::Vector3 b, gc::Vector3 c,
                       gc::Vector3 p)
{
    gc::Vector3 ab = b - a;
    gc::Vector3 ac = c - a;
    gc::Vector3 ap = p - a;

    double d1 = gc::dot(ab, ap);
    double d2 = gc::dot(ac, ap);

    // Voronoi region of vertex a
    if (d1 <= 0.0 && d2 <= 0.0)
        return { a, {1.0, 0.0, 0.0} };

    gc::Vector3 bp = p - b;
    double d3 = gc::dot(ab, bp);
    double d4 = gc::dot(ac, bp);

    // Voronoi region of vertex b
    if (d3 >= 0.0 && d4 <= d3)
        return { b, {0.0, 1.0, 0.0} };

    gc::Vector3 cp = p - c;
    double d5 = gc::dot(ab, cp);
    double d6 = gc::dot(ac, cp);

    // Voronoi region of vertex c
    if (d6 >= 0.0 && d5 <= d6)
        return { c, {0.0, 0.0, 1.0} };

    // Voronoi region of edge ab
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        double v = d1 / (d1 - d3);
        gc::Vector3 q = a + v * ab;
        return { q, {1.0 - v, v, 0.0} };
    }

    // Voronoi region of edge ac
    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        double w = d2 / (d2 - d6);
        gc::Vector3 q = a + w * ac;
        return { q, {1.0 - w, 0.0, w} };
    }

    // Voronoi region of edge bc
    double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        gc::Vector3 q = b + w * (c - b);
        return { q, {0.0, 1.0 - w, w} };
    }

    // Interior of triangle
    double denom = 1.0 / (va + vb + vc);
    double v     = vb * denom;
    double w     = vc * denom;
    gc::Vector3 q = a + v * ab + w * ac;
    return { q, {1.0 - v - w, v, w} };
}


// -----------------------------------------------------------------------------
// §1.2  closestSegmentSegment
//
// Finds the parameters (s*, t*) ∈ [0,1]² that minimise
//   ‖ pe(s) − pg(t) ‖²    where   pe(s) = (1−s)·e0 + s·e1
//                                  pg(t) = (1−t)·g0 + t·g1
//
// Uses the closed-form analytic solution (Ericson §5.1.9), with clamping and
// re-projection at the boundary of [0,1]².
//
// WHY the clamping matters for force decomposition:
//   If s or t hits the boundary, the closest-point pair transitions into a
//   vertex-edge (or vertex-vertex) sub-case.  The barycentric weights below
//   handle this correctly: at s=0 the full force goes to e0 (weight 1-s=1).
// -----------------------------------------------------------------------------
std::pair<double,double>
closestSegmentSegment(gc::Vector3 e0, gc::Vector3 e1,
                      gc::Vector3 g0, gc::Vector3 g1)
{
    constexpr double eps = 1e-12;
    gc::Vector3 d1 = e1 - e0;
    gc::Vector3 d2 = g1 - g0;
    gc::Vector3 r  = e0 - g0;
    double a = gc::dot(d1, d1);   // ‖d1‖²
    double e = gc::dot(d2, d2);   // ‖d2‖²
    double f = gc::dot(d2, r);

    double s, t;

    if (a < eps && e < eps) { return {0.0, 0.0}; }   // both degenerate

    if (a < eps) {           // edge e is a point
        s = 0.0;
        t = std::clamp(f / e, 0.0, 1.0);
    } else {
        double c = gc::dot(d1, r);
        if (e < eps) {       // edge g is a point
            t = 0.0;
            s = std::clamp(-c / a, 0.0, 1.0);
        } else {
            double b     = gc::dot(d1, d2);
            double denom = a * e - b * b;   // > 0 unless edges are parallel

            if (std::abs(denom) > eps) {
                s = std::clamp((b * f - c * e) / denom, 0.0, 1.0);
            } else {
                s = 0.0;   // parallel: arbitrary, handle below
            }

            t = (b * s + f) / e;

            // Re-project if t is outside [0,1]
            if (t < 0.0) {
                t = 0.0;
                s = std::clamp(-c / a, 0.0, 1.0);
            } else if (t > 1.0) {
                t = 1.0;
                s = std::clamp((b - c) / a, 0.0, 1.0);
            }
        }
    }
    return {s, t};
}


// ═════════════════════════════════════════════════════════════════════════════
// §2  TRIANGLE–TRIANGLE MINIMUM CONTACT
// ═════════════════════════════════════════════════════════════════════════════

// -----------------------------------------------------------------------------
// §2.1  minimumFacePairContact
//
// Evaluates all 15 primitive sub-cases for triangles A = {a0,a1,a2} and
// B = {b0,b1,b2} and returns the one with the smallest distance.
//
// Sub-cases:
//   VF (3+3 = 6): each vertex of A vs triangle B, each vertex of B vs triangle A.
//   EE (3×3 = 9): each directed edge of A vs each directed edge of B.
//
// Barycentric weights for the winning sub-case are stored in contact.wA / wB.
// These weights are used by the caller to distribute the contact force to the
// six triangle vertices according to the DDG virtual-work principle.
// -----------------------------------------------------------------------------
FacePairContact
minimumFacePairContact(const std::array<gc::Vector3,3>& A,
                       const std::array<gc::Vector3,3>& B)
{
    FacePairContact best;    // distance = ∞ by default

    auto update = [&](gc::Vector3 pA, std::array<double,3> wA,
                      gc::Vector3 pB, std::array<double,3> wB)
    {
        double d = gc::norm(pA - pB);
        if (d < best.distance) {
            best.pA = pA;  best.wA = wA;
            best.pB = pB;  best.wB = wB;
            best.distance = d;
        }
    };

    // ── Vertex-of-A  vs  face B ──────────────────────────────────────────────
    for (int i = 0; i < 3; ++i) {
        auto [pB, wB] = closestPointOnTriangle(B[0], B[1], B[2], A[i]);
        // Barycentric weight of the query vertex on face A: it is a vertex, so
        // weight = 1 at position i and 0 elsewhere.
        std::array<double,3> wA{};
        wA[i] = 1.0;
        update(A[i], wA, pB, wB);
    }

    // ── Vertex-of-B  vs  face A ──────────────────────────────────────────────
    for (int j = 0; j < 3; ++j) {
        auto [pA, wA] = closestPointOnTriangle(A[0], A[1], A[2], B[j]);
        std::array<double,3> wB{};
        wB[j] = 1.0;
        update(pA, wA, B[j], wB);
    }

    // ── Edge-of-A  vs  Edge-of-B ─────────────────────────────────────────────
    // Edges of triangle A: (0,1), (1,2), (2,0)
    // Edges of triangle B: (0,1), (1,2), (2,0)
    constexpr int edgeA[3][2] = {{0,1},{1,2},{2,0}};
    constexpr int edgeB[3][2] = {{0,1},{1,2},{2,0}};

    for (auto [ia, ja] : edgeA) {
        for (auto [ib, jb] : edgeB) {
            auto [s, t] = closestSegmentSegment(A[ia], A[ja], B[ib], B[jb]);

            gc::Vector3 pA = (1.0 - s) * A[ia] + s * A[ja];
            gc::Vector3 pB = (1.0 - t) * B[ib] + t * B[jb];

            // Barycentric weights on face A: only vertices ia and ja participate.
            // The edge from A[ia] to A[ja] has parametric point at s:
            //   pA = (1-s)·A[ia] + s·A[ja]  →  wA[ia] = 1-s, wA[ja] = s.
            std::array<double,3> wA{}, wB{};
            wA[ia] = 1.0 - s;   wA[ja] = s;
            wB[ib] = 1.0 - t;   wB[jb] = t;

            update(pA, wA, pB, wB);
        }
    }

    return best;
}


// ═════════════════════════════════════════════════════════════════════════════
// §3  SYSTEM METHOD:  computeFaceFaceRepulsionForce
// ═════════════════════════════════════════════════════════════════════════════

// -----------------------------------------------------------------------------
// §3.1  USER-DEFINED REPULSION POTENTIAL
//
// Replace this function with your specific steric potential.
// It receives the signed gap  gap = d - d0  (negative when overlapping).
//
// Requirements:
//   • dEdd(gap) < 0  when gap < 0  (repulsion pushes surfaces apart).
//   • dEdd(gap) = 0  when gap ≥ 0  (no force beyond cutoff).
//
// Examples:
//   Harmonic:      E(gap) = 0.5·μ·gap²  →  dEdd = μ·gap   (gap < 0, so < 0 ✓)
//   Inverse-power: E(gap) = −μ/gap      →  dEdd = μ/gap²  (but μ must be neg)
//   IPC barrier:   E = −κ·gap²·ln(d/d0) →  dEdd = ... (see §3.2 below)
// -----------------------------------------------------------------------------
static double repulsionDeriv(double gap, double mu)
{
    // Harmonic repulsion: E(gap) = 0.5·μ·gap², dE/dd = μ·gap
    // gap = d − d0 < 0 in contact, so dEdd < 0 → force pushes apart ✓
    return mu * gap;   // USER: replace with your potential derivative here
}


// -----------------------------------------------------------------------------
// §3.2  (Optional) IPC-style smooth barrier derivative for reference.
//
// The Incremental Potential Contact barrier (Li et al. 2020) gives a smooth,
// C² potential that goes to ∞ as d→0:
//
//   E(d) = −κ · (d − d̂)² · ln(d / d̂)     for 0 < d < d̂
//   dE/dd = −κ · [ 2(d−d̂)·ln(d/d̂) + (d−d̂)²/d ]
//
// Use this instead of the harmonic kernel if you need guaranteed non-penetration
// (harmonic repulsion can be overcome by large forces in stiff problems).
// Uncomment and substitute into repulsionDeriv() above.
//
// static double ipcBarrierDeriv(double d, double d_hat, double kappa)
// {
//     if (d >= d_hat) return 0.0;
//     double diff = d - d_hat;
//     return -kappa * (2.0 * diff * std::log(d / d_hat) + diff * diff / d);
// }
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// §3.3  System::computeFaceFaceRepulsionForce
//
// Main loop.  For every non-adjacent face pair (Fi, Fj):
//   1. Quick-reject by face-centre distance (optional; marked below).
//   2. Compute minimumFacePairContact to find (pA, wA, pB, wB, d).
//   3. If d < cutoff, evaluate the repulsion force via repulsionDeriv.
//   4. Distribute force to the six involved vertices using the DDG formula:
//
//       n̂ = (pA − pB) / d      [unit normal, B → A]
//       force_on_vert_k_of_A  +=  −dE/dd · wA[k] · n̂
//       force_on_vert_k_of_B  +=  +dE/dd · wB[k] · n̂
//
// Adjacency exclusion:
//   Faces that share an edge or a vertex ("1-ring") are skipped.  The threshold
//   is controlled by parameters.selfAvoidance.n (number of vertex rings).
//   Currently we use a simple shared-vertex check (n=1 ring equivalent).
//
// Complexity: O(F²) brute-force.  For large meshes (N > ~2000 vertices,
// F ≈ 2N faces) replace the inner loop with a BVH or spatial hash grid—
// the interface (FacePairContact + wA/wB) is unchanged.
// -----------------------------------------------------------------------------
void System::computeFaceFaceRepulsionForce()
{
    const double d0  = parameters.selfAvoidance.d;   // contact distance / cutoff
    const double mu  = parameters.selfAvoidance.mu;  // stiffness coefficient
    const std::size_t n_excl = parameters.selfAvoidance.n;  // topological exclusion rings

    if (mu == 0.0 || d0 <= 0.0) return;

    auto& vpg  = *geometry.vpg;
    auto& mesh = *geometry.mesh;

    // Zero the force accumulator (re-use the existing vector field).
    forces.selfAvoidanceForceVec.fill({0.0, 0.0, 0.0});

    // Cache vertex positions as an array for efficient lookup.
    const std::size_t nV = mesh.nVertices();

    // ── Adjacency mask ───────────────────────────────────────────────────────
    // For each face, collect the vertex index set (i,j,k).
    // Two faces are excluded if they share any vertex within n_excl rings.
    // Here we implement the simplest case (n_excl = 1 ring: shared vertex).
    // For n_excl > 1, use meshMutator.markVertices() as in computeSelfAvoidanceForce.
    //
    // Note: for the common case n_excl = 1 we just check vertex-set intersection,
    // which is 9 comparisons per face pair — negligible vs. the geometry query.

    auto faceVertexIndices = [&](gcs::Face f) -> std::array<std::size_t,3>
    {
        gcs::Halfedge he = f.halfedge();
        return { he.vertex().getIndex(),
                 he.next().vertex().getIndex(),
                 he.next().next().vertex().getIndex() };
    };

    auto shareVertex = [](const std::array<std::size_t,3>& a,
                          const std::array<std::size_t,3>& b) -> bool
    {
        for (std::size_t i : a)
            for (std::size_t j : b)
                if (i == j) return true;
        return false;
    };

    // ── Main loop ────────────────────────────────────────────────────────────
    // Index-based iteration matches the pattern used throughout force.cpp
    // (e.g. computeSelfAvoidanceForce uses mesh->vertex(i) in the same style).
    const std::size_t nF = mesh.nFaces();

    for (std::size_t i = 0; i < nF; ++i) {
        gcs::Face faceA{mesh.face(i)};
        auto idxA = faceVertexIndices(faceA);

        // Vertex positions of face A.
        gcs::Halfedge heA = faceA.halfedge();
        std::array<gc::Vector3,3> posA = {
            vpg.inputVertexPositions[heA.vertex()],
            vpg.inputVertexPositions[heA.next().vertex()],
            vpg.inputVertexPositions[heA.next().next().vertex()]
        };

        // Face-centre for quick-reject.
        gc::Vector3 centreA = (posA[0] + posA[1] + posA[2]) / 3.0;

        for (std::size_t j = i + 1; j < nF; ++j) {
            gcs::Face faceB{mesh.face(j)};
            auto idxB = faceVertexIndices(faceB);

            // Skip topologically adjacent faces.
            if (shareVertex(idxA, idxB)) continue;

            gcs::Halfedge heB = faceB.halfedge();
            std::array<gc::Vector3,3> posB = {
                vpg.inputVertexPositions[heB.vertex()],
                vpg.inputVertexPositions[heB.next().vertex()],
                vpg.inputVertexPositions[heB.next().next().vertex()]
            };

            // ── Quick-reject: centre-to-centre distance ───────────────────
            // (Optional acceleration: replace with BVH for large meshes.)
            gc::Vector3 centreB = (posB[0] + posB[1] + posB[2]) / 3.0;
            double centreDist = gc::norm(centreA - centreB);
            // Max radius of a triangle from its centre ≤ max edge length.
            // Use a generous bound: if centres are > 2·d0 apart, faces cannot
            // be closer than d0. Tune the multiplier for your mesh resolution.
            if (centreDist > 3.0 * d0) continue;

            // ── Minimum triangle–triangle contact ─────────────────────────
            FacePairContact contact = minimumFacePairContact(posA, posB);
            const double d = contact.distance;
            if (d >= d0 || d < 1e-10) continue;

            // ── Repulsion force magnitude ─────────────────────────────────
            // gap = d − d0  (< 0 in contact)
            // dE/dd = repulsionDeriv(gap, mu)  (< 0 for repulsion)
            // |F| = −dE/dd  (> 0)
            const double gap  = d - d0;
            const double dEdd = repulsionDeriv(gap, mu);

            // ── Unit normal from B toward A ───────────────────────────────
            gc::Vector3 nhat = (contact.pA - contact.pB) / d;

            // ── DDG force decomposition ───────────────────────────────────
            // Force on vertex k of face A:  f_k = −dE/dd · wA[k] · n̂
            // Force on vertex k of face B:  f_k = +dE/dd · wB[k] · n̂
            //
            // Both signs together satisfy Newton's 3rd:
            //   Σ_k f_k(A) + Σ_k f_k(B)
            //   = −dE/dd·n̂·Σ wA[k]  +  dE/dd·n̂·Σ wB[k]
            //   = −dE/dd·n̂ + dE/dd·n̂  = 0  ✓   (since Σwi = 1)

            gcs::Halfedge heA_v = faceA.halfedge();
            for (int k = 0; k < 3; ++k) {
                forces.selfAvoidanceForceVec[heA_v.vertex()] += -dEdd * contact.wA[k] * nhat;
                heA_v = heA_v.next();
            }

            gcs::Halfedge heB_v = faceB.halfedge();
            for (int k = 0; k < 3; ++k) {
                forces.selfAvoidanceForceVec[heB_v.vertex()] += +dEdd * contact.wB[k] * nhat;
                heB_v = heB_v.next();
            }
        }
    }

    // Project onto vertex normals for the scalar force component used by the
    // integrator's error norm (matches the convention in computeSelfAvoidanceForce).
    forces.selfAvoidanceForce = forces.ontoNormal(forces.selfAvoidanceForceVec);
}


} // namespace solver
} // namespace mem3dg
