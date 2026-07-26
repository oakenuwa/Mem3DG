#pragma once

// =============================================================================
// multi_membrane.h
//
// Model several independent 3-D membranes that interact through excluded-volume
// (face-face) repulsion.
//
// Each membrane is a full Mem3DG `System` with its own mesh, geometry, physical
// parameters, and time integrator.  The membranes are otherwise independent —
// they only "feel" each other through a short-range repulsive force applied
// wherever a triangle of one membrane approaches a triangle of another.
//
// WHY A SEPARATE MODULE?
// ──────────────────────
// System::computeSelfAvoidanceForce() and System::computeFaceFaceRepulsionForce()
// handle a membrane repelling ITSELF (one mesh).  Between two DISTINCT meshes the
// bookkeeping is actually simpler — no two faces ever share a vertex, so the
// topological-adjacency skip is unnecessary — but the force must be written into
// two different System force buffers.  That cross-System step is what this module
// adds.
//
// INTEGRATION STRATEGY (operator splitting)
// ──────────────────────────────────────────
// Euler::march() reads System::forces.mechanicalForceVec directly; it does NOT
// recompute forces (that happens in Euler::status(), which the internal
// integrate() loop calls).  So to inject a cross-membrane force we drive the
// integrators manually, mirroring the CoupledSimulator pattern:
//
//   for each membrane m:
//       m.computeConservativeForcing();        // resets + fills mechanicalForceVec
//       m.addNonconservativeForcing(dt);
//   computeInterMembraneRepulsion(all pairs);  // += into mechanicalForceVec
//   for each membrane m:
//       m.integrator.march();                  // consumes mechanicalForceVec
//
// The repulsion is evaluated with FROZEN positions (all membranes see the same
// snapshot), then every membrane advances one step — a symmetric Lie–Trotter
// split that keeps Newton's 3rd law intact within each synchronized step.
// =============================================================================

#include "mem3dg/macros.h"
#include "mem3dg/solver/face_repulsion.h"
#include "mem3dg/solver/integrator/forward_euler.h"
#include "mem3dg/solver/system.h"

#include <vector>

namespace mem3dg {
namespace solver {

// ---------------------------------------------------------------------------
// Parameters controlling the pairwise membrane-membrane repulsion.
// ---------------------------------------------------------------------------
struct InterMembraneParams {
  /// Contact distance d0.  Faces closer than this repel; faces farther apart
  /// feel nothing.  Typically ~ the sum of the two membranes' half-thicknesses
  /// (or a fraction of the mean edge length).
  double d = 0;

  /// Penalty coefficient μ (stiffness of the harmonic repulsion).
  double mu = 0;

  /// Face-centre quick-reject multiplier.  A face pair is skipped outright when
  /// their centroids are farther apart than cutoffFactor * d.  Larger = safer
  /// for coarse meshes with big triangles, slower for fine meshes.
  double cutoffFactor = 3.0;
};

// ---------------------------------------------------------------------------
// Add mutual excluded-volume repulsion between two DISTINCT membranes A and B.
//
// For every non-trivially-distant (face_of_A, face_of_B) pair, computes the
// minimum triangle-triangle distance and, when it is below params.d, applies a
// harmonic repulsion decomposed to the six involved vertices via DDG barycentric
// weights (see face_repulsion.h for the decomposition).
//
// Forces are ACCUMULATED into:
//   A.forces.mechanicalForceVec   (and A.forces.selfAvoidanceForceVec for output)
//   B.forces.mechanicalForceVec   (and B.forces.selfAvoidanceForceVec for output)
//
// Call AFTER each System's computeConservativeForcing() (which zeroes the
// buffers) and BEFORE the integrator's march().  Momentum is conserved: the
// total force added to A equals minus the total added to B.
//
// PRECONDITION: A and B must be different System objects (different meshes).
// For a membrane repelling itself, use System::computeFaceFaceRepulsionForce().
// ---------------------------------------------------------------------------
DLL_PUBLIC void computeInterMembraneRepulsionForce(
    System &A, System &B, const InterMembraneParams &params);

// ---------------------------------------------------------------------------
// Driver that advances N membranes in lock-step with mutual face repulsion.
//
// The MultiMembraneSystem does NOT own the System or Euler objects — the caller
// constructs and owns them, and is responsible for keeping them alive for the
// lifetime of the driver.  This mirrors how CoupledSimulator references an
// externally-owned System and Simul.
//
// SYNCHRONIZED-STEP REQUIREMENT
// ──────────────────────────────
// The frozen-snapshot split is only consistent if all membranes advance by the
// SAME dt each step, so that they share a common clock.  Configure every
// member's Euler integrator with a FIXED time step:
//
//     integrator.isBacktrack    = false;   // no per-membrane line search
//     integrator.ifAdaptiveStep = false;   // no per-membrane mesh-size scaling
//
// Leaving backtracking or adaptive stepping on lets each membrane choose its own
// timeStep inside march(), which desynchronizes their clocks and makes the
// repulsion snapshot stale.  Use identical characteristicTimeStep across members.
// ---------------------------------------------------------------------------
class DLL_PUBLIC MultiMembraneSystem {
public:
  // One membrane: a System paired with the integrator that advances it.
  struct Member {
    System *system;
    integrator::Euler *integrator;
  };

  MultiMembraneSystem(std::vector<Member> members, InterMembraneParams params);

  // Advance every membrane by one integrator step, applying pairwise repulsion
  // over all C(N,2) membrane pairs on the frozen pre-step configuration.
  void step();

  // Convenience loop: repeatedly step() until every membrane has reached
  // totalTime, saving each membrane's trajectory every saveInterval.
  void run(double totalTime, double saveInterval);

  // Access the members (e.g., to read positions or change parameters).
  std::vector<Member> &members() { return members_; }
  const std::vector<Member> &members() const { return members_; }

private:
  // Compute repulsion over all unordered membrane pairs and accumulate into
  // each System's mechanicalForceVec.  Assumes each System's conservative
  // forcing has already been computed for this step.
  void applyPairwiseRepulsion();

  std::vector<Member> members_;
  InterMembraneParams params_;
};

} // namespace solver
} // namespace mem3dg
