// =============================================================================
// multi_membrane.cpp
//
// Implementation of inter-membrane excluded-volume repulsion and the
// MultiMembraneSystem driver.  See multi_membrane.h for the design rationale.
// =============================================================================

#include "mem3dg/solver/multi_membrane.h"

#include <geometrycentral/utilities/vector3.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace mem3dg {
namespace solver {

namespace gc = geometrycentral;
namespace gcs = geometrycentral::surface;

// ─────────────────────────────────────────────────────────────────────────────
// §1  Local helpers
// ─────────────────────────────────────────────────────────────────────────────

namespace {

// Harmonic repulsion energy derivative dE/dd.
//   E(d)   = (mu/2) (d - d0)^2   for d < d0,  else 0
//   dE/dd  = mu (d - d0)         (< 0 when overlapping → repulsive)
// Kept inline and local so the potential is easy to swap for a different law
// (e.g. inverse-power) without touching the decomposition logic below.
inline double repulsionDeriv(double gap, double mu) { return mu * gap; }

// Extract the three vertex positions of a face in halfedge order
// {he.vertex(), he.next().vertex(), he.next().next().vertex()}, matching the
// wA/wB barycentric layout produced by minimumFacePairContact().
inline std::array<gc::Vector3, 3>
facePositions(gcs::Face f, const gcs::VertexPositionGeometry &vpg) {
  gcs::Halfedge he = f.halfedge();
  return {vpg.inputVertexPositions[he.vertex()],
          vpg.inputVertexPositions[he.next().vertex()],
          vpg.inputVertexPositions[he.next().next().vertex()]};
}

// Circumradius of a triangle about its own centroid: the smallest r such that a
// sphere of radius r at the centroid contains the triangle.  Needed to make the
// face-pair reject CONSERVATIVE (see the note in computeInterMembraneRepulsion-
// Force) rather than a heuristic multiple of d0.
inline double centroidRadius(const std::array<gc::Vector3, 3> &p,
                             const gc::Vector3 &centre) {
  return std::max({gc::norm(p[0] - centre), gc::norm(p[1] - centre),
                   gc::norm(p[2] - centre)});
}

struct BoundingSphere {
  gc::Vector3 centre{0, 0, 0};
  double radius = 0;
};

/// Bounding sphere of a whole membrane: area-weighted centre and the farthest
/// vertex from it.  Used for the organelle-level broad phase.
BoundingSphere membraneBounds(const System &S) {
  const gcs::VertexPositionGeometry &vpg = *S.geometry.vpg;
  gcs::ManifoldSurfaceMesh &mesh = *S.geometry.mesh;
  gc::Vector3 acc{0, 0, 0};
  double wsum = 0;
  for (gcs::Vertex v : mesh.vertices()) {
    const double w = vpg.vertexDualAreas[v];
    acc += w * vpg.inputVertexPositions[v];
    wsum += w;
  }
  BoundingSphere b;
  b.centre = (wsum > 0) ? acc / wsum : acc;
  b.radius = 0;
  for (gcs::Vertex v : mesh.vertices())
    b.radius =
        std::max(b.radius, gc::norm(vpg.inputVertexPositions[v] - b.centre));
  return b;
}

} // namespace

// ─────────────────────────────────────────────────────────────────────────────
// §2  Pairwise membrane-membrane repulsion
// ─────────────────────────────────────────────────────────────────────────────

void computeInterMembraneRepulsionForce(System &A, System &B,
                                        const InterMembraneParams &params) {
  computeInterMembraneRepulsionForce(A, B, params, nullptr);
}

void computeInterMembraneRepulsionForce(System &A, System &B,
                                        const InterMembraneParams &params,
                                        const FaceBVH *bvhB) {
  const double d0 = params.d;
  const double mu = params.mu;
  if (mu == 0 || d0 <= 0)
    return;

  gcs::ManifoldSurfaceMesh &meshA = *A.geometry.mesh;
  gcs::ManifoldSurfaceMesh &meshB = *B.geometry.mesh;
  gcs::VertexPositionGeometry &vpgA = *A.geometry.vpg;
  gcs::VertexPositionGeometry &vpgB = *B.geometry.vpg;

  /*
   Face-pair reject.

   The bound used here is CONSERVATIVE: a pair is skipped only when
       ‖cA − cB‖ − rA − rB > d0
   where rA, rB are the triangles' circumradii about their own centroids.  Since
   every point of a triangle lies within r of its centroid, that expression is a
   true lower bound on the triangle-triangle distance, so nothing within d0 is
   ever discarded.

   This replaces an earlier test, `‖cA − cB‖ > cutoffFactor * d0`, which was a
   heuristic and not a bound: it ignored triangle size entirely.  On the
   1280-face icosphere at R=0.8 the circumradius is ≈0.070, so rA+rB ≈ 0.139 and
   a correct reject needs ‖cA−cB‖ ≤ d0 + 0.139; with d0=0.1 the old
   cutoffFactor=3 gave 0.3, which happened to pass with only 1.25x margin.  At
   the d0=0.0035 used by the Cytosim coupling it gave 0.0105 — far inside
   rA+rB — and silently threw away real contacts.

   cutoffFactor is still honoured as an ADDITIONAL user-tunable slack so
   existing configs keep working, but it can no longer cause a missed contact.
  */
  const double slack = std::max(0.0, params.cutoffFactor) * d0;

  const std::size_t nFA = meshA.nFaces();
  const std::size_t nFB = meshB.nFaces();

  // Candidate list for one face of A.  Declared outside the loop so its buffer
  // is reused across all nFA queries instead of reallocating each time.
  std::vector<std::size_t> cand;

  for (std::size_t i = 0; i < nFA; ++i) {
    gcs::Face faceA{meshA.face(i)};
    std::array<gc::Vector3, 3> posA = facePositions(faceA, vpgA);
    gc::Vector3 centreA = (posA[0] + posA[1] + posA[2]) / 3.0;
    const double rA = centroidRadius(posA, centreA);

    /*
     Narrow-phase candidates.

     With a BVH: box(faceA) padded by (d0 + slack), queried against B's tree.
     That is conservative -- if two triangles are within t of each other then
     their boxes padded by t must overlap on every axis -- so no contact is
     lost.  Pruned pairs are farther apart than d0, where the harmonic force is
     identically zero, which is why the forces come out bit-identical to the
     exhaustive scan rather than merely similar.

     Without one: every face of B, as before.
    */
    cand.clear();
    if (bvhB) {
      Aabb q;
      q.reset();
      for (int k = 0; k < 3; ++k)
        q.expand(posA[k]);
      q.pad(d0 + slack);
      bvhB->query(q, cand);
    } else {
      cand.reserve(nFB);
      for (std::size_t j = 0; j < nFB; ++j)
        cand.push_back(j);
    }

    for (std::size_t j : cand) {
      gcs::Face faceB{meshB.face(j)};
      std::array<gc::Vector3, 3> posB = facePositions(faceB, vpgB);

      // Two membranes never share vertices, so — unlike the single-mesh case —
      // there is no adjacency skip to perform here.
      gc::Vector3 centreB = (posB[0] + posB[1] + posB[2]) / 3.0;
      const double rB = centroidRadius(posB, centreB);
      if (gc::norm(centreA - centreB) - rA - rB > d0 + slack)
        continue;

      // Minimum triangle-triangle contact (6 vertex-face + 9 edge-edge cases).
      FacePairContact contact = minimumFacePairContact(posA, posB);
      const double d = contact.distance;
      if (d >= d0 || d < 1e-10)
        continue;

      const double gap = d - d0;
      const double dEdd = repulsionDeriv(gap, mu);

      // Unit normal from B toward A.
      gc::Vector3 nhat = (contact.pA - contact.pB) / d;

      // DDG decomposition (identical convention to
      // System::computeFaceFaceRepulsionForce):
      //   face-A vertex k :  f_k = −dE/dd · wA[k] · n̂   (pushes A along +n̂)
      //   face-B vertex k :  f_k = +dE/dd · wB[k] · n̂   (pushes B along −n̂)
      // Σ over both faces cancels because Σ wA = Σ wB = 1  → momentum conserved.
      gcs::Halfedge heA = faceA.halfedge();
      for (int k = 0; k < 3; ++k) {
        gc::Vector3 f = -dEdd * contact.wA[k] * nhat;
        A.forces.mechanicalForceVec[heA.vertex()] += f;
        A.forces.selfAvoidanceForceVec[heA.vertex()] += f;
        heA = heA.next();
      }

      gcs::Halfedge heB = faceB.halfedge();
      for (int k = 0; k < 3; ++k) {
        gc::Vector3 f = +dEdd * contact.wB[k] * nhat;
        B.forces.mechanicalForceVec[heB.vertex()] += f;
        B.forces.selfAvoidanceForceVec[heB.vertex()] += f;
        heB = heB.next();
      }
    }
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// §3  MultiMembraneSystem driver
// ─────────────────────────────────────────────────────────────────────────────

MultiMembraneSystem::MultiMembraneSystem(std::vector<Member> members,
                                         InterMembraneParams params)
    : members_(std::move(members)), params_(params) {}

void MultiMembraneSystem::applyPairwiseRepulsion() {
  const std::size_t n = members_.size();

  /*
   BROAD PHASE (organelle level).

   computeInterMembraneRepulsionForce is O(nF_a * nF_b); wrapping it in the
   C(N,2) pair loop makes the whole step O(N^2 * nF^2).  At N=5 with 1280 faces
   that is ~16M face-pair tests EVERY step, to find a few hundred real contacts.

   Most of those pairs are nowhere near each other.  One bounding sphere per
   membrane rules them out for the cost of a single distance: two closed
   surfaces can only come within d0 if
       ‖c_a − c_b‖ − r_a − r_b < d0
   which is a true lower bound on their surface-surface distance, so this
   discards no contact.  Cost drops to O(N^2) cheap tests plus the narrow phase
   on only the pairs that can actually touch.

   The bounds are computed ONCE per step here, not per pair, so this is O(N)
   mesh passes regardless of how many pairs survive.
  */
  std::vector<BoundingSphere> bounds;
  bounds.reserve(n);
  for (std::size_t k = 0; k < n; ++k)
    bounds.push_back(membraneBounds(*members_[k].system));

  /*
   NARROW-PHASE ACCELERATION (face level).

   One BVH per membrane, built once and refitted every step.  Refit is O(nF)
   with no allocation; a rebuild is only triggered when the face count changes,
   which is how remeshing is detected -- the tree stores face INDICES, so
   refitting a stale tree would reference faces that no longer exist.
  */
  bvh_.resize(n);
  for (std::size_t k = 0; k < n; ++k) {
    System &S = *members_[k].system;
    if (!bvh_[k])
      bvh_[k] = std::make_unique<FaceBVH>();
    if (bvh_[k]->needsRebuild(*S.geometry.mesh))
      bvh_[k]->build(*S.geometry.mesh, *S.geometry.vpg);
    else
      bvh_[k]->refit(*S.geometry.vpg);
  }

  const double reach = params_.d * (1.0 + std::max(0.0, params_.cutoffFactor));

  for (std::size_t a = 0; a < n; ++a) {
    for (std::size_t b = a + 1; b < n; ++b) {
      const double gap = gc::norm(bounds[b].centre - bounds[a].centre) -
                         bounds[a].radius - bounds[b].radius;
      if (gap > reach)
        continue;
      computeInterMembraneRepulsionForce(*members_[a].system,
                                         *members_[b].system, params_,
                                         bvh_[b].get());
    }
  }
}

void MultiMembraneSystem::step() {
  // 1. Recompute each membrane's own (conservative + nonconservative) forcing
  //    on the current configuration.  This zeroes and refills mechanicalForceVec
  //    and selfAvoidanceForceVec, so it MUST precede the repulsion injection.
  for (Member &m : members_) {
    m.system->computeConservativeForcing();
    m.system->addNonconservativeForcing(m.integrator->characteristicTimeStep);
  }

  // 2. Add cross-membrane repulsion evaluated on the frozen snapshot.
  applyPairwiseRepulsion();

  // 3. Refresh the scalar (onto-normal) mechanical force so the integrator's
  //    error norm reflects the injected repulsion, then advance each membrane.
  for (Member &m : members_) {
    m.system->forces.mechanicalForce =
        m.system->forces.ontoNormal(m.system->forces.mechanicalForceVec);
    m.integrator->march();
  }
}

void MultiMembraneSystem::run(double totalTime, double saveInterval) {
  if (members_.empty())
    return;

  // Drive off the first membrane's clock; all share the same characteristic dt
  // in the typical setup.  Each membrane saves its own trajectory.
  double nextSave = members_.front().system->time;

  while (members_.front().system->time < totalTime) {
    if (members_.front().system->time >= nextSave) {
      for (Member &m : members_)
        m.integrator->saveData(m.integrator->ifOutputTrajFile,
                               m.integrator->ifOutputMeshFile, false);
      nextSave += saveInterval;
    }
    step();
  }

  // Final save.
  for (Member &m : members_)
    m.integrator->saveData(m.integrator->ifOutputTrajFile,
                           m.integrator->ifOutputMeshFile, false);
}

} // namespace solver
} // namespace mem3dg
