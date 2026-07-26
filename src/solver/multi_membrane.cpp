// =============================================================================
// multi_membrane.cpp
//
// Implementation of inter-membrane excluded-volume repulsion and the
// MultiMembraneSystem driver.  See multi_membrane.h for the design rationale.
// =============================================================================

#include "mem3dg/solver/multi_membrane.h"

#include <geometrycentral/utilities/vector3.h>

#include <array>
#include <cstddef>

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

} // namespace

// ─────────────────────────────────────────────────────────────────────────────
// §2  Pairwise membrane-membrane repulsion
// ─────────────────────────────────────────────────────────────────────────────

void computeInterMembraneRepulsionForce(System &A, System &B,
                                        const InterMembraneParams &params) {
  const double d0 = params.d;
  const double mu = params.mu;
  if (mu == 0 || d0 <= 0)
    return;

  gcs::ManifoldSurfaceMesh &meshA = *A.geometry.mesh;
  gcs::ManifoldSurfaceMesh &meshB = *B.geometry.mesh;
  gcs::VertexPositionGeometry &vpgA = *A.geometry.vpg;
  gcs::VertexPositionGeometry &vpgB = *B.geometry.vpg;

  const double rejectRadius = params.cutoffFactor * d0;

  const std::size_t nFA = meshA.nFaces();
  const std::size_t nFB = meshB.nFaces();

  for (std::size_t i = 0; i < nFA; ++i) {
    gcs::Face faceA{meshA.face(i)};
    std::array<gc::Vector3, 3> posA = facePositions(faceA, vpgA);
    gc::Vector3 centreA = (posA[0] + posA[1] + posA[2]) / 3.0;

    for (std::size_t j = 0; j < nFB; ++j) {
      gcs::Face faceB{meshB.face(j)};
      std::array<gc::Vector3, 3> posB = facePositions(faceB, vpgB);

      // Quick-reject on centroid distance.  Two membranes never share
      // vertices, so — unlike the single-mesh case — there is no adjacency
      // skip to perform here.
      gc::Vector3 centreB = (posB[0] + posB[1] + posB[2]) / 3.0;
      if (gc::norm(centreA - centreB) > rejectRadius)
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
  for (std::size_t a = 0; a < n; ++a) {
    for (std::size_t b = a + 1; b < n; ++b) {
      computeInterMembraneRepulsionForce(*members_[a].system,
                                         *members_[b].system, params_);
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
