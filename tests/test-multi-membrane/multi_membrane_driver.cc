/*
 =============================================================================
 multi_membrane_driver.cc

 Reusable driver + physics harness for Mem3DG's multi-membrane framework
 (mem3dg/solver/multi_membrane.h): N independent membranes ("organelles") that
 interact only through short-range face-face excluded-volume repulsion.

 This links libmem3dg's C++20 API DIRECTLY -- unlike coupled_driver.cc, which
 goes through the plain-C ABI because it has to coexist with Cytosim's C++11.
 There is no Cytosim here, so there is no reason to pay that indirection.

 SCENARIOS  (argv[1])
 --------------------
   approach : two equal organelles driven toward each other by a constant
              external body force.  Tests that the repulsion can HOLD a load --
              the forcing is what keeps pushing, the repulsion is what stops
              them merging.

   inflate  : two adjacent organelles, one given a preferred volume larger than
              its current one via an osmotic term.  It swells and displaces its
              neighbour.  Tests that the repulsion transmits a quasi-static
              push rather than only resisting an impulse.

 WHAT IS CHECKED  (every scenario, every step)
 ---------------------------------------------
   1. Newton's 3rd law.  The repulsion is measured in isolation by snapshotting
      mechanicalForceVec, calling computeInterMembraneRepulsionForce(), and
      differencing.  sum(dF_A) + sum(dF_B) must be 0 to round-off.  This is the
      strongest single check on the DDG barycentric decomposition: it holds only
      if the per-face weights sum to 1 on BOTH faces.

   2. Non-interpenetration.  The true minimum surface-surface distance is
      recomputed by brute force (all face pairs, exact triangle-triangle) rather
      than trusting the solver's own accounting.  It must stay > 0.

   3. Separation vs contact distance.  Once in contact, the gap should sit near
      params.d and not collapse -- a harmonic wall, not a soft suggestion.

   4. Volume/area bookkeeping per membrane, so `inflate` can be shown to
      actually inflate.

 OUTPUT
 ------
   <outDir>/membrane<K>_%04d.obj    one sequence per organelle, same frame index
   <outDir>/metrics.csv             the per-step physics record

 The membrane<K>_ prefix is what blender/import_trajectory.py groups on, so a
 multi-membrane run drops into the existing viewer unchanged.
 =============================================================================
*/

#include "mem3dg/solver/face_bvh.h"
#include "mem3dg/solver/geometry.h"
#include "mem3dg/solver/integrator/forward_euler.h"
#include "mem3dg/solver/multi_membrane.h"
#include "mem3dg/solver/system.h"

#include <geometrycentral/surface/surface_mesh.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace ms = mem3dg::solver;
namespace gc = geometrycentral;
namespace gcs = geometrycentral::surface;

using mem3dg::EigenVectorX3dr;
using mem3dg::EigenVectorX3sr;

// ---------------------------------------------------------------------------
// One organelle: geometry + system + integrator, kept alive together.
// ---------------------------------------------------------------------------
struct Organelle {
  std::unique_ptr<ms::Geometry> geometry;
  std::unique_ptr<ms::System> system;
  std::unique_ptr<ms::integrator::Euler> integrator;
};

/**
 Build one organelle from a base mesh, rigidly translated and uniformly scaled.

 The mesh is read once by the caller and passed in as (V,F) so that every
 organelle shares the same topology without re-parsing the file, and so the
 offset can be applied before Geometry construction -- moving vertices after
 System::initialize() would leave cached geometric quantities stale.
 */
static std::unique_ptr<Organelle>
makeOrganelle(const EigenVectorX3dr &V0, const EigenVectorX3sr &F0,
              gc::Vector3 offset, double scale, double Kb, double dt) {
  auto org = std::make_unique<Organelle>();

  EigenVectorX3dr V = V0;
  for (Eigen::Index i = 0; i < V.rows(); ++i) {
    V(i, 0) = V0(i, 0) * scale + offset.x;
    V(i, 1) = V0(i, 1) * scale + offset.y;
    V(i, 2) = V0(i, 2) * scale + offset.z;
  }
  EigenVectorX3sr F = F0;

  org->geometry = std::make_unique<ms::Geometry>(F, V);

  ms::Parameters p;
  p.variation.isShapeVariation = true;
  p.variation.isProteinVariation = false;
  p.proteinMobility = 0;
  p.temperature = 0;
  p.bending.Kb = Kb;
  p.bending.Kbc = Kb;
  p.bending.H0c = 0;
  p.selfAvoidance.mu = 0; // self-avoidance is a different feature; keep it off
  p.dpd.gamma = 0;

  const std::size_t nV = org->geometry->mesh->nVertices();
  Eigen::VectorXd protein = Eigen::VectorXd::Constant(nV, 1.0);
  EigenVectorX3dr velocity = EigenVectorX3dr::Zero(nV, 3);

  org->system = std::make_unique<ms::System>(*org->geometry, protein, velocity,
                                             p, 0.0);
  org->system->initialize();

  org->integrator = std::make_unique<ms::integrator::Euler>(
      *org->system, dt, std::numeric_limits<double>::max(), 1.0, 1e-4, ".");
  // REQUIRED by MultiMembraneSystem: per-membrane backtracking or adaptive
  // stepping would let each organelle pick its own dt inside march(), which
  // desynchronizes their clocks and makes the frozen repulsion snapshot stale.
  org->integrator->isBacktrack = false;
  org->integrator->ifAdaptiveStep = false;
  org->integrator->ifOutputTrajFile = false;
  org->integrator->ifOutputMeshFile = false;
  org->integrator->ifPrintToConsole = false;
  return org;
}

// ---------------------------------------------------------------------------
// Exact triangle-triangle minimum distance, by brute force over all face pairs.
// Deliberately independent of the solver's own contact search: if both agree,
// the solver is not missing contacts.
// ---------------------------------------------------------------------------
static double pointSegDist2(const gc::Vector3 &p, const gc::Vector3 &a,
                            const gc::Vector3 &b) {
  gc::Vector3 ab = b - a;
  double t = dot(p - a, ab);
  double d = dot(ab, ab);
  t = (d > 0) ? std::min(1.0, std::max(0.0, t / d)) : 0.0;
  gc::Vector3 q = a + t * ab;
  return dot(p - q, p - q);
}

static double pointTriDist2(const gc::Vector3 &p, const gc::Vector3 &a,
                            const gc::Vector3 &b, const gc::Vector3 &c) {
  gc::Vector3 ab = b - a, ac = c - a, n = cross(ab, ac);
  double nn = dot(n, n);
  if (nn > 1e-30) {
    // barycentric of the in-plane projection
    gc::Vector3 ap = p - a;
    double u = dot(cross(ab, ap), n) / nn; // weight of c
    double v = dot(cross(ap, ac), n) / nn; // weight of b
    if (u >= 0 && v >= 0 && u + v <= 1) {
      double s = dot(ap, n) / std::sqrt(nn);
      return s * s;
    }
  }
  return std::min({pointSegDist2(p, a, b), pointSegDist2(p, b, c),
                   pointSegDist2(p, c, a)});
}

static double segSegDist2(const gc::Vector3 &p1, const gc::Vector3 &q1,
                          const gc::Vector3 &p2, const gc::Vector3 &q2) {
  gc::Vector3 d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
  double a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
  double s, t;
  const double EPS = 1e-30;
  if (a <= EPS && e <= EPS)
    return dot(r, r);
  if (a <= EPS) {
    s = 0;
    t = std::min(1.0, std::max(0.0, f / e));
  } else {
    double c = dot(d1, r);
    if (e <= EPS) {
      t = 0;
      s = std::min(1.0, std::max(0.0, -c / a));
    } else {
      double b = dot(d1, d2), denom = a * e - b * b;
      s = (denom > EPS) ? std::min(1.0, std::max(0.0, (b * f - c * e) / denom))
                        : 0.0;
      t = (b * s + f) / e;
      if (t < 0) {
        t = 0;
        s = std::min(1.0, std::max(0.0, -c / a));
      } else if (t > 1) {
        t = 1;
        s = std::min(1.0, std::max(0.0, (b - c) / a));
      }
    }
  }
  gc::Vector3 c1 = p1 + s * d1, c2 = p2 + t * d2;
  return dot(c1 - c2, c1 - c2);
}

/**
 Does segment pq pass through triangle abc?  (Moller-Trumbore.)

 REQUIRED for correctness, not an optimization.  Vertex-triangle and edge-edge
 distances do NOT detect two triangles that cross: an edge of A can pass clean
 through the interior of B with every vertex far from the other triangle and no
 edge-edge pair touching, so the minimum over those cases is POSITIVE while the
 true distance is 0.  Without this, an interpenetrating pair of surfaces reports
 a small positive separation and the run is scored as passing.

 This is the same piercing sub-case that minimumEdgeFaceContact() handles in the
 solver's own kernel; the verifier has to handle it too or the two disagree.
 */
static bool segmentPiercesTriangle(const gc::Vector3 &p, const gc::Vector3 &q,
                                   const gc::Vector3 &a, const gc::Vector3 &b,
                                   const gc::Vector3 &c) {
  const gc::Vector3 dir = q - p;
  const gc::Vector3 e1 = b - a, e2 = c - a;
  const gc::Vector3 h = cross(dir, e2);
  const double det = dot(e1, h);
  if (std::fabs(det) < 1e-30)
    return false; // segment parallel to the triangle plane
  const double inv = 1.0 / det;
  const gc::Vector3 s = p - a;
  const double u = inv * dot(s, h);
  if (u < 0.0 || u > 1.0)
    return false;
  const gc::Vector3 sc = cross(s, e1);
  const double v = inv * dot(dir, sc);
  if (v < 0.0 || u + v > 1.0)
    return false;
  const double t = inv * dot(e2, sc);
  return (t >= 0.0 && t <= 1.0); // within the segment, not just the line
}

static double triTriDist(const gc::Vector3 *A, const gc::Vector3 *B) {
  // crossing triangles are at distance 0; check that FIRST, since the
  // feature-distance minima below cannot see it
  for (int i = 0; i < 3; ++i) {
    if (segmentPiercesTriangle(A[i], A[(i + 1) % 3], B[0], B[1], B[2]))
      return 0.0;
    if (segmentPiercesTriangle(B[i], B[(i + 1) % 3], A[0], A[1], A[2]))
      return 0.0;
  }
  double best = std::numeric_limits<double>::max();
  for (int i = 0; i < 3; ++i) {
    best = std::min(best, pointTriDist2(A[i], B[0], B[1], B[2]));
    best = std::min(best, pointTriDist2(B[i], A[0], A[1], A[2]));
  }
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      best = std::min(best, segSegDist2(A[i], A[(i + 1) % 3], B[j],
                                        B[(j + 1) % 3]));
  return std::sqrt(best);
}

/// Brute-force minimum surface-surface distance between two organelles.
static double minSurfaceDistance(const ms::System &A, const ms::System &B) {
  auto &mA = *A.geometry.mesh;
  auto &vA = *A.geometry.vpg;
  auto &mB = *B.geometry.mesh;
  auto &vB = *B.geometry.vpg;
  double best = std::numeric_limits<double>::max();
  for (gcs::Face fa : mA.faces()) {
    gc::Vector3 TA[3];
    int k = 0;
    for (gcs::Vertex v : fa.adjacentVertices())
      TA[k++] = vA.inputVertexPositions[v];
    // centroid quick-reject against the running best
    gc::Vector3 ca = (TA[0] + TA[1] + TA[2]) / 3.0;
    for (gcs::Face fb : mB.faces()) {
      gc::Vector3 TB[3];
      int j = 0;
      for (gcs::Vertex v : fb.adjacentVertices())
        TB[j++] = vB.inputVertexPositions[v];
      gc::Vector3 cb = (TB[0] + TB[1] + TB[2]) / 3.0;
      // triangles cannot be closer than centroid distance minus both radii
      double rA = std::max({norm(TA[0] - ca), norm(TA[1] - ca), norm(TA[2] - ca)});
      double rB = std::max({norm(TB[0] - cb), norm(TB[1] - cb), norm(TB[2] - cb)});
      if (norm(ca - cb) - rA - rB > best)
        continue;
      best = std::min(best, triTriDist(TA, TB));
    }
  }
  return best;
}

/// Total force currently sitting in a System's mechanical force buffer.
/// Differencing this across computeInterMembraneRepulsionForce() isolates the
/// repulsion's own contribution from bending, tension and osmotic terms.
static gc::Vector3 netForce(const ms::System &S) {
  gc::Vector3 sum{0, 0, 0};
  for (gcs::Vertex v : S.geometry.mesh->vertices())
    sum += S.forces.mechanicalForceVec[v];
  return sum;
}

// ---------------------------------------------------------------------------
static gc::Vector3 centroid(const ms::System &S) {
  auto &m = *S.geometry.mesh;
  auto &v = *S.geometry.vpg;
  gc::Vector3 c{0, 0, 0};
  for (gcs::Vertex x : m.vertices())
    c += v.inputVertexPositions[x];
  return c / (double)m.nVertices();
}

static void saveObj(const ms::System &S, const std::string &dir, int member,
                    int frame) {
  char path[512];
  std::snprintf(path, sizeof(path), "%s/membrane%d_%04d.obj", dir.c_str(),
                member, frame);
  std::FILE *f = std::fopen(path, "w");
  if (!f)
    return;
  auto &m = *S.geometry.mesh;
  auto &vpg = *S.geometry.vpg;
  for (gcs::Vertex v : m.vertices()) {
    gc::Vector3 p = vpg.inputVertexPositions[v];
    std::fprintf(f, "v %g %g %g\n", p.x, p.y, p.z);
  }
  for (gcs::Face fa : m.faces()) {
    std::fprintf(f, "f");
    for (gcs::Vertex v : fa.adjacentVertices())
      std::fprintf(f, " %llu", (unsigned long long)(v.getIndex() + 1));
    std::fprintf(f, "\n");
  }
  std::fclose(f);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// `cluster' scenario: N >= 3 organelles jammed together.
//
// WHY THIS EXISTS SEPARATELY FROM approach/inflate
// ------------------------------------------------
// Those two hand-roll the step loop and call computeInterMembraneRepulsionForce
// on the single pair directly, so MultiMembraneSystem::step() -- and with it
// applyPairwiseRepulsion()'s C(N,2) loop -- never runs.  With only two
// membranes that loop has exactly one iteration, so it cannot show the failure
// mode that matters at N >= 3: a membrane must ACCUMULATE force from N-1
// partners into one buffer, and any pair that is skipped, double-counted, or
// written to the wrong System is invisible in a 2-body test.
//
// So this scenario deliberately calls multi.step() and drives the push through
// parameters.external.form -- the only hook that fires INSIDE step(), via
// addNonconservativeForcing() -> prescribeExternalForce().  Injecting the force
// after the fact would mean bypassing the code under test again.
//
// GEOMETRY: N organelles equally spaced on a ring in the xy-plane, each pushed
// toward the origin.  They converge and jam.  The arrangement is symmetric on
// purpose: the exact repulsion force must then be the same magnitude on every
// member and sum to zero over the ring, which is a check no single pair can
// provide.
// ---------------------------------------------------------------------------
static int runCluster(const EigenVectorX3dr &V0, const EigenVectorX3sr &F0,
                      double R, int nOrg, double Kb, double dt, int nSteps,
                      int saveEvery, double contactD, double mu, double push,
                      const std::string &outDir) {
  if (nOrg < 3) {
    std::fprintf(stderr, "cluster needs at least 3 organelles\n");
    return 2;
  }

  // Ring radius: place them so the ADJACENT surface gap starts just clear of
  // contact.  Adjacent centre-to-centre distance on a ring of radius rho is
  // 2*rho*sin(pi/N), so solve that for the gap we want.
  const double startClear = 1.5 * contactD;
  const double chord = 2.0 * R + startClear;
  const double rho = chord / (2.0 * std::sin(M_PI / nOrg));

  std::vector<std::unique_ptr<Organelle>> orgs;
  for (int k = 0; k < nOrg; ++k) {
    const double th = 2.0 * M_PI * k / nOrg;
    gc::Vector3 c{rho * std::cos(th), rho * std::sin(th), 0.0};
    orgs.push_back(makeOrganelle(V0, F0, c, 1.0, Kb, dt));
    if (!orgs.back())
      return 2;
    // uniform body force toward the origin; direction frozen at the initial
    // centre so the field is rigid rather than re-converging every step
    const gc::Vector3 dir = -c / norm(c);
    orgs.back()->system->parameters.external.form =
        [dir, push](EigenVectorX3dr pos, Eigen::VectorXd, double,
                    Eigen::VectorXd) -> EigenVectorX3dr {
      EigenVectorX3dr f(pos.rows(), 3);
      for (Eigen::Index i = 0; i < pos.rows(); ++i) {
        f(i, 0) = push * dir.x;
        f(i, 1) = push * dir.y;
        f(i, 2) = push * dir.z;
      }
      return f;
    };
  }

  ms::InterMembraneParams params;
  params.d = contactD;
  params.mu = mu;
  params.cutoffFactor = 3.0;

  std::vector<ms::MultiMembraneSystem::Member> members;
  for (auto &o : orgs)
    members.push_back({o->system.get(), o->integrator.get()});
  ms::MultiMembraneSystem multi(members, params);

  const int nPairs = nOrg * (nOrg - 1) / 2;
  std::printf("scenario:  cluster\n");
  std::printf("organelles: %d on a ring rho=%.4f, %zu faces each (R=%.4f)\n",
              nOrg, rho, orgs[0]->system->geometry.mesh->nFaces(), R);
  std::printf("contact:   d=%.4g  mu=%.4g   dt=%.4g  steps=%d  pairs=%d\n",
              params.d, params.mu, dt, nSteps, nPairs);
  std::printf("stepping:  MultiMembraneSystem::step()  [N-body path]\n");

  char mpath[512];
  std::snprintf(mpath, sizeof(mpath), "%s/metrics.csv", outDir.c_str());
  std::FILE *csv = std::fopen(mpath, "w");
  if (csv)
    std::fprintf(csv, "step,time,min_sep,contact_pairs,momentum_residual,"
                      "max_member_force,force_spread,ring_radius\n");

  double worstMomentum = 0.0, minSep = 1e30, worstSpread = 0.0;
  int frame = 0;

  for (int s = 0; s <= nSteps; ++s) {
    // ---- the code under test -------------------------------------------
    multi.step();

    // CHECK 1 -- global momentum conservation over ALL pairs.
    //
    // step() has just filled each member's selfAvoidanceForceVec with its
    // inter-membrane contribution only (selfAvoidance.mu is 0, and
    // computeConservativeForcing() zeroed the buffer at the top of step()).
    // Summing over every member therefore totals the whole N-body repulsion,
    // which must vanish.  This is the check that a 2-body test cannot make:
    // it fails if any pair writes into the wrong System or is counted twice.
    gc::Vector3 total{0, 0, 0};
    double maxMember = 0.0, minMember = 1e30;
    for (auto &o : orgs) {
      gc::Vector3 fk{0, 0, 0};
      for (gcs::Vertex v : o->system->geometry.mesh->vertices())
        fk += o->system->forces.selfAvoidanceForceVec[v];
      total += fk;
      const double n = norm(fk);
      maxMember = std::max(maxMember, n);
      minMember = std::min(minMember, n);
    }
    const double residual = norm(total);
    const double relResidual = (maxMember > 1e-14) ? residual / maxMember : 0.0;
    worstMomentum = std::max(worstMomentum, relResidual);

    // CHECK 2 -- symmetry, as a DIAGNOSTIC, not an invariant.
    //
    // Equally spaced members under a radial push should carry similar repulsion
    // magnitudes.  But they will not match exactly: the organelles are
    // TRANSLATED copies of one mesh, not rotated ones, so an icosphere's
    // triangles do not line up with the ring's N-fold axis and each contact
    // patch samples a different set of faces.  A few percent is discretization,
    // not a bug -- so this is reported, not asserted.
    //
    // Gate on minMember, not maxMember: at the step where the first pair
    // engages, one member is loaded and another is still at exactly zero, which
    // is a spread of 1.0 that says nothing about symmetry.
    const double spread =
        (maxMember > 1e-14) ? (maxMember - minMember) / maxMember : 0.0;

    if (s % saveEvery == 0) {
      // CHECK 3 -- non-interpenetration over every pair, brute force.
      // Bounding-sphere reject at the ORGANELLE level.  minSurfaceDistance is
      // O(nF^2) per pair (1.6M triangle pairs here), and C(N,2) grows fast --
      // at N=5 the unfiltered scan does not finish in reasonable time.  Most
      // pairs on a ring are nowhere near each other, and a centroid/radius
      // bound rules them out for the cost of one distance.
      std::vector<gc::Vector3> cen(nOrg);
      std::vector<double> rad(nOrg, 0.0);
      for (int i = 0; i < nOrg; ++i) {
        cen[i] = centroid(*orgs[i]->system);
        for (gcs::Vertex v :
             orgs[i]->system->geometry.mesh->vertices())
          rad[i] = std::max(
              rad[i], norm(orgs[i]->system->geometry.vpg->inputVertexPositions[v]
                           - cen[i]));
      }
      double sep = 1e30;
      int contactPairs = 0;
      for (int i = 0; i < nOrg; ++i)
        for (int j = i + 1; j < nOrg; ++j) {
          // lower bound on the true surface-surface distance for this pair
          const double lb = norm(cen[j] - cen[i]) - rad[i] - rad[j];
          // cannot be in contact, and cannot improve the running minimum
          if (lb >= params.d && lb > sep)
            continue;
          const double dij = minSurfaceDistance(*orgs[i]->system,
                                                *orgs[j]->system);
          sep = std::min(sep, dij);
          if (dij < params.d)
            ++contactPairs;
        }
      minSep = std::min(minSep, sep);

      // Aggregate the symmetry spread ONLY while the contact set is the
      // adjacent-only ring (contactPairs == nOrg), which is the configuration
      // that actually has N-fold symmetry.  Once diagonals start touching
      // (N >= 4) the contact set itself is asymmetric, so a large spread there
      // is correct physics, not a defect -- aggregating it would turn a right
      // answer into a scary-looking number.
      if (contactPairs == nOrg && minMember > 1e-9)
        worstSpread = std::max(worstSpread, spread);

      // ring radius: mean distance of the member centroids from the origin,
      // the scalar that says how far the cluster has collapsed
      double ringR = 0.0;
      for (auto &o : orgs)
        ringR += norm(centroid(*o->system));
      ringR /= nOrg;

      if (csv)
        std::fprintf(csv, "%d,%g,%g,%d,%g,%g,%g,%g\n", s, s * dt, sep,
                     contactPairs, relResidual, maxMember, spread, ringR);
      std::printf("step %5d  sep %10.4e  contacts %d/%d  ring %7.4f  "
                  "|F|max %9.3e  spread %8.2e  momentum %8.2e\n",
                  s, sep, contactPairs, nPairs, ringR, maxMember, spread,
                  relResidual);
      for (int k = 0; k < nOrg; ++k)
        saveObj(*orgs[k]->system, outDir, k, frame);
      ++frame;
    }
  }
  if (csv)
    std::fclose(csv);

  std::printf("\n--- physics summary ---\n");
  std::printf("N-body momentum  : worst |sum_k F_k| / max_k |F_k| = %.3e\n",
              worstMomentum);
  std::printf("ring symmetry    : worst spread in |F_k| = %.3e "
              "(mesh discretization; not an exact invariant)\n", worstSpread);
  std::printf("min separation   : %.6e  (contact distance d = %.4f)\n", minSep,
              params.d);
  const bool interpenetrated = !(minSep > 1e-12);
  std::printf("interpenetration : %s\n", interpenetrated ? "DETECTED" : "NONE");
  std::printf("frames written   : %d x %d organelles -> %s\n", frame, nOrg,
              outDir.c_str());
  return (!interpenetrated && worstMomentum < 1e-8) ? 0 : 1;
}


// ---------------------------------------------------------------------------
// `bvhcheck' scenario: does the BVH narrow phase return the SAME force as the
// exhaustive nF_A * nF_B scan?
//
// Comparing two full trajectories cannot answer this.  Contact dynamics amplify
// round-off, so any difference at all -- including a different summation ORDER
// into the same vertex -- shows up as a visibly different |F| after a few
// thousand steps.  That tells you the runs diverged, not whether the physics
// differs.
//
// So this evaluates BOTH paths from the SAME configuration, one step at a time,
// and reports the per-vertex difference directly.
// ---------------------------------------------------------------------------
static int runBvhCheck(const EigenVectorX3dr &V0, const EigenVectorX3sr &F0,
                       double R, int nOrg, double Kb, double dt, int nSteps,
                       double contactD, double mu, double push) {
  const double startClear = 1.5 * contactD;
  const double chord = 2.0 * R + startClear;
  const double rho = chord / (2.0 * std::sin(M_PI / nOrg));

  std::vector<std::unique_ptr<Organelle>> orgs;
  for (int k = 0; k < nOrg; ++k) {
    const double th = 2.0 * M_PI * k / nOrg;
    gc::Vector3 c{rho * std::cos(th), rho * std::sin(th), 0.0};
    orgs.push_back(makeOrganelle(V0, F0, c, 1.0, Kb, dt));
    const gc::Vector3 dir = -c / norm(c);
    orgs.back()->system->parameters.external.form =
        [dir, push](EigenVectorX3dr pos, Eigen::VectorXd, double,
                    Eigen::VectorXd) -> EigenVectorX3dr {
      EigenVectorX3dr f(pos.rows(), 3);
      for (Eigen::Index i = 0; i < pos.rows(); ++i) {
        f(i, 0) = push * dir.x; f(i, 1) = push * dir.y; f(i, 2) = push * dir.z;
      }
      return f;
    };
  }

  ms::InterMembraneParams params;
  params.d = contactD; params.mu = mu; params.cutoffFactor = 3.0;

  std::vector<ms::MultiMembraneSystem::Member> members;
  for (auto &o : orgs) members.push_back({o->system.get(), o->integrator.get()});
  ms::MultiMembraneSystem multi(members, params);

  std::printf("bvhcheck: %d organelles, d=%.4g mu=%.4g, %d steps\n",
              nOrg, contactD, mu, nSteps);

  std::vector<ms::FaceBVH> trees(nOrg);
  double worstAbs = 0.0, worstRel = 0.0, peakForce = 0.0;
  long comparisons = 0, nonzeroSteps = 0;

  for (int s = 0; s <= nSteps; ++s) {
    for (int i = 0; i < nOrg; ++i)
      for (int j = 0; j < nOrg; ++j) {
        if (i == j) continue;
        ms::System &A = *orgs[i]->system;
        ms::System &B = *orgs[j]->system;

        // --- exhaustive ---
        A.forces.selfAvoidanceForceVec.fill({0, 0, 0});
        B.forces.selfAvoidanceForceVec.fill({0, 0, 0});
        ms::computeInterMembraneRepulsionForce(A, B, params);
        std::vector<gc::Vector3> refA, refB;
        for (gcs::Vertex v : A.geometry.mesh->vertices())
          refA.push_back(A.forces.selfAvoidanceForceVec[v]);
        for (gcs::Vertex v : B.geometry.mesh->vertices())
          refB.push_back(B.forces.selfAvoidanceForceVec[v]);

        // --- BVH ---
        if (trees[j].needsRebuild(*B.geometry.mesh))
          trees[j].build(*B.geometry.mesh, *B.geometry.vpg);
        else
          trees[j].refit(*B.geometry.vpg);
        A.forces.selfAvoidanceForceVec.fill({0, 0, 0});
        B.forces.selfAvoidanceForceVec.fill({0, 0, 0});
        ms::computeInterMembraneRepulsionForce(A, B, params, &trees[j]);

        std::size_t idx = 0;
        double localMax = 0.0;
        for (gcs::Vertex v : A.geometry.mesh->vertices()) {
          const gc::Vector3 d = A.forces.selfAvoidanceForceVec[v] - refA[idx];
          localMax = std::max(localMax, norm(d));
          peakForce = std::max(peakForce, norm(refA[idx]));
          ++idx;
        }
        idx = 0;
        for (gcs::Vertex v : B.geometry.mesh->vertices()) {
          const gc::Vector3 d = B.forces.selfAvoidanceForceVec[v] - refB[idx];
          localMax = std::max(localMax, norm(d));
          peakForce = std::max(peakForce, norm(refB[idx]));
          ++idx;
        }
        if (localMax > 0) ++nonzeroSteps;
        worstAbs = std::max(worstAbs, localMax);
        ++comparisons;
      }
    multi.step();
  }

  worstRel = (peakForce > 0) ? worstAbs / peakForce : 0.0;
  std::printf("\n--- BVH vs exhaustive ---\n");
  std::printf("pair evaluations compared : %ld\n", comparisons);
  std::printf("evaluations that differed : %ld\n", nonzeroSteps);
  std::printf("worst per-vertex |dF|     : %.6e\n", worstAbs);
  std::printf("peak per-vertex |F|       : %.6e\n", peakForce);
  std::printf("worst relative difference : %.3e\n", worstRel);
  std::printf("verdict: %s\n",
              (worstRel < 1e-12) ? "IDENTICAL to round-off"
                                 : "DIFFERENT -- BVH is missing contacts");
  return (worstRel < 1e-12) ? 0 : 1;
}

int main(int argc, char **argv) {
  if (argc < 4) {
    std::fprintf(stderr,
                 "usage: %s <approach|inflate|cluster> <mesh.obj> <outDir> "
                 "[nSteps] [dt] [d] [mu] [saveEvery] [push] [grow] [Kv] [Ksg] "
                 "[nOrganelles]\n",
                 argv[0]);
    return 2;
  }
  const std::string scenario = argv[1];
  const std::string meshPath = argv[2];
  const std::string outDir = argv[3];
  const int nSteps = (argc > 4) ? std::atoi(argv[4]) : 2000;
  const double dt = (argc > 5) ? std::atof(argv[5]) : 1e-4;
  const double contactD = (argc > 6) ? std::atof(argv[6]) : 0.10;
  const double mu = (argc > 7) ? std::atof(argv[7]) : 500.0;
  const int saveEvery = (argc > 8) ? std::atoi(argv[8]) : 100;
  // `approach': per-vertex body force driving the two together.  Exposed
  // because the right value depends on dt, mesh resolution and how far apart
  // the organelles start -- a hardcoded default cannot be right for all three.
  const double pushArg = (argc > 9) ? std::atof(argv[9]) : 0.0;
  // `inflate': target volume for A as a multiple of its initial volume
  const double growArg = (argc > 10) ? std::atof(argv[10]) : 0.0;
  // `inflate': osmotic and stretching moduli.  These two COMPETE -- the
  // osmotic term drives volume up, surface tension resists the area increase
  // that requires.  Laplace balance (p ~ 2*sigma/R) means nothing inflates
  // unless Kv is a few times Ksg, so they have to be tuned together, not
  // separately.
  const double KvArg  = (argc > 11) ? std::atof(argv[11]) : 0.0;
  const double KsgArg = (argc > 12) ? std::atof(argv[12]) : 0.0;
  // `cluster': how many organelles on the ring
  const int nOrgArg = (argc > 13) ? std::atoi(argv[13]) : 3;
  const double Kb = 8.22e-5;

  // ---- read the base mesh once -------------------------------------------
  EigenVectorX3dr V0;
  EigenVectorX3sr F0;
  {
    ms::Geometry base(meshPath);
    auto &m = *base.mesh;
    auto &vpg = *base.vpg;
    V0.resize(m.nVertices(), 3);
    for (gcs::Vertex v : m.vertices()) {
      gc::Vector3 p = vpg.inputVertexPositions[v];
      V0(v.getIndex(), 0) = p.x;
      V0(v.getIndex(), 1) = p.y;
      V0(v.getIndex(), 2) = p.z;
    }
    F0.resize(m.nFaces(), 3);
    for (gcs::Face fa : m.faces()) {
      int k = 0;
      for (gcs::Vertex v : fa.adjacentVertices())
        F0(fa.getIndex(), k++) = v.getIndex();
    }
  }

  // radius of the base mesh, used to place the two organelles just clear of
  // each other so the run starts out of contact and comes into it
  double R = 0;
  for (Eigen::Index i = 0; i < V0.rows(); ++i)
    R = std::max(R, std::sqrt(V0(i, 0) * V0(i, 0) + V0(i, 1) * V0(i, 1) +
                              V0(i, 2) * V0(i, 2)));

  // `cluster' has its own N-body loop (it exercises MultiMembraneSystem::step()
  // rather than driving the pair force by hand), so it branches out here before
  // the 2-organelle setup below.
  if (scenario == "bvhcheck")
    return runBvhCheck(V0, F0, R, nOrgArg, Kb, dt, nSteps, contactD, mu,
                       (pushArg > 0) ? pushArg : 1.0);

  if (scenario == "cluster")
    return runCluster(V0, F0, R, nOrgArg, Kb, dt, nSteps, saveEvery, contactD,
                      mu, (pushArg > 0) ? pushArg : 1.0, outDir);

  // ---- scenario setup -----------------------------------------------------
  std::vector<std::unique_ptr<Organelle>> orgs;
  double pushForce = 0.0;   // `approach' only
  double gap0 = 0.6 * R;    // initial clear space between the two surfaces

  if (scenario == "approach") {
    // equal organelles, symmetric about the origin, driven together
    double c = R + 0.5 * gap0;
    orgs.push_back(makeOrganelle(V0, F0, {-c, 0, 0}, 1.0, Kb, dt));
    orgs.push_back(makeOrganelle(V0, F0, {+c, 0, 0}, 1.0, Kb, dt));
    pushForce = (pushArg > 0) ? pushArg : 2.0e-3;
  } else if (scenario == "inflate") {
    // B is the smaller neighbour; A will swell into it.
    //
    // Unlike `approach', this scenario starts the two ALREADY CLOSE -- just
    // outside contact -- because the thing under test is the swelling, not the
    // travel.  Starting far apart would spend the whole run drifting with the
    // repulsion switched off and never exercise the growing contact patch.
    // The offset is derived from the two radii and the contact distance so it
    // stays just-clear for any mesh or any d.
    const double scaleB = 0.8;
    const double startClear = 1.5 * contactD;
    double c = 0.5 * (R + scaleB * R + startClear);
    orgs.push_back(makeOrganelle(V0, F0, {-c, 0, 0}, 1.0, Kb, dt));
    orgs.push_back(makeOrganelle(V0, F0, {+c, 0, 0}, scaleB, Kb, dt));

    // Drive A toward 2.2x its current volume.  Sign convention follows
    // mem3dg_system_set_osmotic: pressure = -Kv*(V-Vt)/Vt, and the volume
    // variation vector points outward, so V < Vt gives positive pressure and
    // inflates.
    const double Kv = (KvArg > 0) ? KvArg : 1.0e-2;
    const double Vt =
        orgs[0]->system->geometry.volume * ((growArg > 0) ? growArg : 2.2);
    orgs[0]->system->parameters.osmotic.form = [Kv, Vt](double volume) {
      const double dV = volume - Vt;
      return std::make_tuple(-Kv * dV / Vt, Kv * dV * dV / (2.0 * Vt));
    };
    // Surface tension on both, so inflation stretches rather than crumples,
    // and B holds its shape while being shoved.
    const double Ksg = (KsgArg > 0) ? KsgArg : 1.0e-2;
    for (auto &o : orgs) {
      const double At = o->system->geometry.surfaceArea;
      o->system->parameters.tension.form = [Ksg, At](double area) {
        const double dA = area - At;
        return std::make_tuple(Ksg * dA / At, Ksg * dA * dA / (2.0 * At));
      };
    }
  } else {
    std::fprintf(stderr, "unknown scenario `%s'\n", scenario.c_str());
    return 2;
  }

  // ---- multi-membrane driver ---------------------------------------------
  ms::InterMembraneParams params;
  params.d = contactD;
  params.mu = mu;
  params.cutoffFactor = 3.0;

  std::vector<ms::MultiMembraneSystem::Member> members;
  for (auto &o : orgs)
    members.push_back({o->system.get(), o->integrator.get()});
  ms::MultiMembraneSystem multi(members, params);

  std::printf("scenario:  %s\n", scenario.c_str());
  std::printf("organelles: %zu, %zu faces each (base R=%.4f)\n", orgs.size(),
              orgs[0]->system->geometry.mesh->nFaces(), R);
  std::printf("contact:   d=%.4g  mu=%.4g   dt=%.4g  steps=%d\n", params.d,
              params.mu, dt, nSteps);

  char mpath[512];
  std::snprintf(mpath, sizeof(mpath), "%s/metrics.csv", outDir.c_str());
  std::FILE *csv = std::fopen(mpath, "w");
  if (csv)
    std::fprintf(csv, "step,time,sep,cx0,cx1,center_dist,vol0,vol1,area0,"
                      "area1,newton_residual,fmag\n");

  double worstNewton = 0.0, worstRelNewton = 0.0, minSep = 1e30;
  int frame = 0;

  for (int s = 0; s <= nSteps; ++s) {
    // ---------------------------------------------------------------------
    // PHYSICS CHECK 1 -- Newton's 3rd law, measured on the repulsion alone.
    //
    // computeConservativeForcing() zeroes and fills mechanicalForceVec; we
    // snapshot it, add ONLY the inter-membrane repulsion, and difference.  The
    // residual is the net force the repulsion injected into the pair, which
    // must vanish.
    // ---------------------------------------------------------------------
    for (auto &o : orgs) {
      o->system->computeConservativeForcing();
      o->system->addNonconservativeForcing(dt);
    }
    if (scenario == "approach") {
      // constant body force driving the two together along x
      for (std::size_t k = 0; k < orgs.size(); ++k) {
        double sgn = (k == 0) ? +1.0 : -1.0;
        auto &S = *orgs[k]->system;
        for (gcs::Vertex v : S.geometry.mesh->vertices())
          S.forces.mechanicalForceVec[v].x += sgn * pushForce;
      }
    }

    // snapshot BEFORE the repulsion so it can be measured in isolation
    gc::Vector3 sumBefore0 = netForce(*orgs[0]->system);
    gc::Vector3 sumBefore1 = netForce(*orgs[1]->system);

    ms::computeInterMembraneRepulsionForce(*orgs[0]->system, *orgs[1]->system,
                                           params);

    gc::Vector3 d0 = netForce(*orgs[0]->system) - sumBefore0;
    gc::Vector3 d1 = netForce(*orgs[1]->system) - sumBefore1;
    const double residual = norm(d0 + d1);
    const double fmag = 0.5 * (norm(d0) + norm(d1));
    worstNewton = std::max(worstNewton, residual);
    if (fmag > 1e-12)
      worstRelNewton = std::max(worstRelNewton, residual / fmag);

    // ---- record + save ---------------------------------------------------
    if (s % saveEvery == 0) {
      double sep = minSurfaceDistance(*orgs[0]->system, *orgs[1]->system);
      minSep = std::min(minSep, sep);
      gc::Vector3 c0 = centroid(*orgs[0]->system);
      gc::Vector3 c1 = centroid(*orgs[1]->system);
      double v0 = orgs[0]->system->geometry.volume;
      double v1 = orgs[1]->system->geometry.volume;
      double a0 = orgs[0]->system->geometry.surfaceArea;
      double a1 = orgs[1]->system->geometry.surfaceArea;
      if (csv)
        std::fprintf(csv, "%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n", s, s * dt,
                     sep, c0.x, c1.x, norm(c1 - c0), v0, v1, a0, a1, residual,
                     fmag);
      std::printf("step %5d  sep %10.4e  Cdist %7.4f  V %7.4f/%7.4f  "
                  "|Frep| %9.3e  newton %8.2e\n",
                  s, sep, norm(c1 - c0), v0, v1, fmag, residual);
      for (std::size_t k = 0; k < orgs.size(); ++k)
        saveObj(*orgs[k]->system, outDir, (int)k, frame);
      ++frame;
    }

    // advance: integrators consume the mechanicalForceVec we just assembled
    for (auto &o : orgs)
      o->integrator->march();
  }

  if (csv)
    std::fclose(csv);

  std::printf("\n--- physics summary ---\n");
  std::printf("Newton 3rd law   : worst |sum F_A + sum F_B| = %.3e "
              "(relative %.3e)\n",
              worstNewton, worstRelNewton);
  // scientific notation, not %f: a separation of 1e-9 is interpenetration for
  // all practical purposes but prints as "0.000000" in fixed format, which is
  // how a failing run can look like a passing one
  std::printf("min separation   : %.6e  (contact distance d = %.4f)\n", minSep,
              params.d);
  const bool interpenetrated = !(minSep > 1e-12);
  std::printf("interpenetration : %s\n",
              interpenetrated ? "DETECTED" : "NONE");
  std::printf("frames written   : %d x %zu organelles -> %s\n", frame,
              orgs.size(), outDir.c_str());
  return (!interpenetrated && worstRelNewton < 1e-8) ? 0 : 1;
}
