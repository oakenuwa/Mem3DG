// ===========================================================================
// mem3dg_c.cpp  —  Implementation of the plain-C ABI facade (mem3dg_c.h).
//
// This translation unit is the ONLY place where the C-linkage boundary meets
// the full C++20 / geometry-central / Eigen world.  It is compiled into
// libmem3dg.dylib (see src/capi/CMakeLists.txt).  Consumers include only
// mem3dg_c.h and never see anything below.
// ===========================================================================

#include "mem3dg/capi/mem3dg_c.h"

#include "mem3dg/solver/geometry.h"
#include "mem3dg/solver/integrator/forward_euler.h"
#include "mem3dg/solver/system.h"

#include "mem3dg/coupled/edge_face_contact.h"

#include <geometrycentral/utilities/knn.h>

#include <geometrycentral/utilities/vector3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include <memory>
#include <tuple>
#include <vector>

namespace {

namespace gc = geometrycentral;
namespace gcs = geometrycentral::surface;
namespace ms = mem3dg::solver;

using EigenVectorX1d = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using EigenVectorX3dr =
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using EigenVectorX3sr =
    Eigen::Matrix<std::size_t, Eigen::Dynamic, 3, Eigen::RowMajor>;

// Everything one membrane needs, kept alive behind the opaque handle.
struct Impl {
  std::unique_ptr<ms::Geometry> geometry;
  std::unique_ptr<ms::System> system;
  std::unique_ptr<ms::integrator::Euler> integrator;
  // Per-vertex external force (3*nV), applied every step until overwritten.
  std::vector<double> externalForce;
  double dt = 1e-3;
  // adaptive remeshing
  bool remesh = false;
  int  remeshPeriod = 1;
  long stepCount = 0;

  /*
   Spatial index over the vertex positions.  Both the closest-point query and
   the inside/outside test used to scan every face, which is O(nFaces) per
   query -- and remeshing pushes nFaces up as the simulation runs.  A KD-tree
   over the vertices lets us look at only the faces incident to the few nearest
   vertices instead.

   Invalidated by anything that moves or adds vertices, so it is rebuilt lazily:
   `knnStale` is set by step(), set_vertex_positions() and remeshing.
   */
  std::unique_ptr<gc::NearestNeighborFinder> knn;
  std::vector<gc::Vector3> knnPoints;
  bool knnStale = true;
  // per-query "already tested" marks for faces, stamped by a generation counter
  // so they never need clearing between queries
  std::vector<unsigned> faceStamp;
  unsigned stampGen = 0;

  std::size_t nV() const { return system->geometry.mesh->nVertices(); }
};

// Build a Parameters with a minimal harmonic bending model and everything else
// inert, so a bare membrane just relaxes / responds to injected forces.
ms::Parameters makeParams(double Kb, double Kbc, double H0c) {
  ms::Parameters p;
  p.variation.isShapeVariation = true;
  p.variation.isProteinVariation = false;
  p.proteinMobility = 0;
  p.temperature = 0;
  p.bending.Kb = Kb;
  p.bending.Kbc = Kbc;
  p.bending.H0c = H0c;
  p.selfAvoidance.mu = 0;
  p.dpd.gamma = 0;
  return p;
}

// Finish wiring an Impl once its Geometry exists: build System + integrator.
Impl *finishSetup(std::unique_ptr<ms::Geometry> geo, ms::Parameters params,
                  double dt) {
  auto impl = std::make_unique<Impl>();
  impl->geometry = std::move(geo);
  impl->dt = (dt > 0 ? dt : 1e-3);

  const std::size_t nV = impl->geometry->mesh->nVertices();
  EigenVectorX1d protein = Eigen::VectorXd::Constant(nV, 1.0);
  EigenVectorX3dr velocity = EigenVectorX3dr::Zero(nV, 3);

  impl->system = std::make_unique<ms::System>(*impl->geometry, protein,
                                              velocity, params, 0.0);
  impl->system->initialize();

  // Fixed-step forward Euler; we drive march() manually from step().
  impl->integrator = std::make_unique<ms::integrator::Euler>(
      *impl->system, impl->dt, /*totalTime=*/std::numeric_limits<double>::max(),
      /*savePeriod=*/1.0, /*tolerance=*/1e-4, /*outputDir=*/".");
  impl->integrator->isBacktrack = false;
  impl->integrator->ifAdaptiveStep = false;
  impl->integrator->ifOutputTrajFile = false;
  impl->integrator->ifOutputMeshFile = false;
  impl->integrator->ifPrintToConsole = false;

  impl->externalForce.assign(3 * nV, 0.0);
  return impl.release();
}

} // namespace

// ===========================================================================
// C ABI
// ===========================================================================

extern "C" {

mem3dg_system *mem3dg_system_create_from_file(const char *mesh_path, double Kb,
                                              double Kbc, double H0c,
                                              double dt) {
  try {
    auto geo = std::make_unique<ms::Geometry>(std::string(mesh_path));
    return reinterpret_cast<mem3dg_system *>(
        finishSetup(std::move(geo), makeParams(Kb, Kbc, H0c), dt));
  } catch (...) {
    return nullptr;
  }
}

mem3dg_system *mem3dg_system_create_from_arrays(const double *vertices,
                                                size_t nV, const size_t *faces,
                                                size_t nF, double Kb, double Kbc,
                                                double H0c, double dt) {
  try {
    EigenVectorX3dr V(nV, 3);
    for (size_t i = 0; i < nV; ++i) {
      V(i, 0) = vertices[3 * i + 0];
      V(i, 1) = vertices[3 * i + 1];
      V(i, 2) = vertices[3 * i + 2];
    }
    EigenVectorX3sr F(nF, 3);
    for (size_t i = 0; i < nF; ++i) {
      F(i, 0) = faces[3 * i + 0];
      F(i, 1) = faces[3 * i + 1];
      F(i, 2) = faces[3 * i + 2];
    }
    auto geo = std::make_unique<ms::Geometry>(F, V);
    return reinterpret_cast<mem3dg_system *>(
        finishSetup(std::move(geo), makeParams(Kb, Kbc, H0c), dt));
  } catch (...) {
    return nullptr;
  }
}

void mem3dg_system_destroy(mem3dg_system *h) {
  delete reinterpret_cast<Impl *>(h);
}

size_t mem3dg_system_num_vertices(const mem3dg_system *h) {
  return reinterpret_cast<const Impl *>(h)->system->geometry.mesh->nVertices();
}

size_t mem3dg_system_num_faces(const mem3dg_system *h) {
  return reinterpret_cast<const Impl *>(h)->system->geometry.mesh->nFaces();
}

void mem3dg_system_get_vertex_positions(const mem3dg_system *h,
                                        double *out_xyz) {
  auto *impl = reinterpret_cast<const Impl *>(h);
  gcs::ManifoldSurfaceMesh &mesh = *impl->system->geometry.mesh;
  gcs::VertexPositionGeometry &vpg = *impl->system->geometry.vpg;
  const std::size_t nV = mesh.nVertices();
  for (std::size_t i = 0; i < nV; ++i) {
    gc::Vector3 p = vpg.inputVertexPositions[mesh.vertex(i)];
    out_xyz[3 * i + 0] = p.x;
    out_xyz[3 * i + 1] = p.y;
    out_xyz[3 * i + 2] = p.z;
  }
}

void mem3dg_system_get_faces(const mem3dg_system *h, size_t *out_ijk) {
  auto *impl = reinterpret_cast<const Impl *>(h);
  gcs::ManifoldSurfaceMesh &mesh = *impl->system->geometry.mesh;
  const std::size_t nF = mesh.nFaces();
  for (std::size_t i = 0; i < nF; ++i) {
    gcs::Halfedge he = mesh.face(i).halfedge();
    out_ijk[3 * i + 0] = he.vertex().getIndex();
    out_ijk[3 * i + 1] = he.next().vertex().getIndex();
    out_ijk[3 * i + 2] = he.next().next().vertex().getIndex();
  }
}

void mem3dg_system_set_vertex_positions(mem3dg_system *h, const double *xyz) {
  auto *impl = reinterpret_cast<Impl *>(h);
  gcs::ManifoldSurfaceMesh &mesh = *impl->system->geometry.mesh;
  gcs::VertexPositionGeometry &vpg = *impl->system->geometry.vpg;
  const std::size_t nV = mesh.nVertices();
  for (std::size_t i = 0; i < nV; ++i)
    vpg.inputVertexPositions[mesh.vertex(i)] =
        gc::Vector3{xyz[3 * i + 0], xyz[3 * i + 1], xyz[3 * i + 2]};
  // refresh cached areas/volumes/normals so the next force call sees the change
  impl->system->updateConfigurations();
  impl->knnStale = true;
}

void mem3dg_system_get_mechanical_forces(mem3dg_system *h, double *out_xyz) {
  auto *impl = reinterpret_cast<Impl *>(h);
  impl->system->computeConservativeForcing();
  impl->system->addNonconservativeForcing(impl->dt);
  gcs::ManifoldSurfaceMesh &mesh = *impl->system->geometry.mesh;
  auto &mv = impl->system->forces.mechanicalForceVec;
  const std::size_t nV = mesh.nVertices();
  for (std::size_t i = 0; i < nV; ++i) {
    gc::Vector3 f = mv[mesh.vertex(i)];
    out_xyz[3 * i + 0] = f.x;
    out_xyz[3 * i + 1] = f.y;
    out_xyz[3 * i + 2] = f.z;
  }
}

void mem3dg_system_set_tension(mem3dg_system *h, double Ksg,
                               double preferred_area) {
  auto *impl = reinterpret_cast<Impl *>(h);
  if (Ksg == 0.0) {
    impl->system->parameters.tension.form = NULL;
    return;
  }
  // default to the area the membrane currently has, so it defends its own shape
  const double At =
      (preferred_area > 0.0) ? preferred_area : impl->system->geometry.surfaceArea;
  impl->system->parameters.tension.form = [Ksg, At](double area) {
    const double dA = area - At;
    const double tension = Ksg * dA / At;
    return std::make_tuple(tension, tension * dA / 2.0);
  };
}

void mem3dg_system_set_osmotic(mem3dg_system *h, double Kv,
                               double preferred_volume) {
  auto *impl = reinterpret_cast<Impl *>(h);
  if (Kv == 0.0) {
    impl->system->parameters.osmotic.form = NULL;
    return;
  }
  const double Vt =
      (preferred_volume > 0.0) ? preferred_volume : impl->system->geometry.volume;
  impl->system->parameters.osmotic.form = [Kv, Vt](double volume) {
    /*
     Sign: osmoticForceVec += pressure * volumeVariationVector, and the volume
     variation vector points OUTWARD, so a positive pressure inflates.  When the
     vesicle is over-inflated (V > Vt) this expression is negative and therefore
     pulls back in, as it should.
     */
    const double dV = volume - Vt;
    const double pressure = -Kv * dV / Vt;
    return std::make_tuple(pressure, Kv * dV * dV / (2.0 * Vt));
  };
}

void mem3dg_system_set_regularization(mem3dg_system *h, double Kse, double Ksl,
                                      double Kst) {
  auto *impl = reinterpret_cast<Impl *>(h);
  impl->system->parameters.spring.Kse = Kse;
  impl->system->parameters.spring.Ksl = Ksl;
  impl->system->parameters.spring.Kst = Kst;
}

void mem3dg_system_set_remeshing(mem3dg_system *h, int enable, int period,
                                 double min_edge, double max_edge) {
  auto *impl = reinterpret_cast<Impl *>(h);
  auto &mut = impl->system->meshProcessor.meshMutator;

  impl->remesh = (enable != 0);
  impl->remeshPeriod = (period > 0 ? period : 1);

  if (!impl->remesh) {
    // clear every toggle; summarizeStatus() then derives isMeshMutate = false
    mut.flipNonDelaunay = false;
    mut.flipNonDelaunayRequireFlat = false;
    mut.splitCurved = mut.splitLarge = mut.splitLong = false;
    mut.splitSharp = mut.splitSkinnyDelaunay = false;
    mut.collapseSkinny = mut.collapseSmall = mut.collapseFlat = false;
    mut.isShiftVertex = false;
    mut.isSmoothenMesh = false;
    impl->system->meshProcessor.summarizeStatus();
    return;
  }

  /*
   splitLarge is what actually resolves a growing tube.  Note that splitLong,
   despite the name, does NOT use maximumEdgeLength: its test is RELATIVE,

       edgeLength(e) > edgeLength(next) + edgeLength(twin.next)

   so an edge must exceed the sum of two neighbours before it splits -- roughly
   a doubling relative to the local mesh, which a smoothly stretching tube never
   reaches.  maximumEdgeLength only guards against COLLAPSING long edges.

   splitLarge is the absolute criterion:

       faceArea(he) + faceArea(twin) > 4 * maximumFaceArea

   Two equilateral triangles of side L sum to 2*(sqrt(3)/4)L^2 = 0.866 L^2, so
   setting maximumFaceArea = 0.2165 * max_edge^2 puts the threshold exactly at
   edges of length max_edge.
   */
  mut.splitLarge = true;
  mut.splitCurved = true;          // also refine where the surface bends sharply
  mut.splitSkinnyDelaunay = true;
  mut.collapseSkinny = true;
  mut.collapseSmall = true;
  mut.flipNonDelaunay = true;
  mut.isShiftVertex = true;

  if (max_edge > 0) {
    mut.maximumEdgeLength = max_edge;
    mut.maximumFaceArea = 0.2165 * max_edge * max_edge;
  }
  if (min_edge > 0) {
    mut.minimumEdgeLength = min_edge;
    mut.minimumFaceArea = 0.2165 * min_edge * min_edge;
  }

  impl->system->meshProcessor.summarizeStatus();
}

double mem3dg_system_area(const mem3dg_system *h) {
  return reinterpret_cast<const Impl *>(h)->system->geometry.surfaceArea;
}

double mem3dg_system_volume(const mem3dg_system *h) {
  return reinterpret_cast<const Impl *>(h)->system->geometry.volume;
}

void mem3dg_system_set_external_forces(mem3dg_system *h,
                                       const double *force_xyz) {
  auto *impl = reinterpret_cast<Impl *>(h);
  const std::size_t n = 3 * impl->nV();
  impl->externalForce.assign(force_xyz, force_xyz + n);
}

void mem3dg_system_clear_external_forces(mem3dg_system *h) {
  auto *impl = reinterpret_cast<Impl *>(h);
  std::fill(impl->externalForce.begin(), impl->externalForce.end(), 0.0);
}

void mem3dg_system_add_external_forces(mem3dg_system *h,
                                       const double *force_xyz) {
  auto *impl = reinterpret_cast<Impl *>(h);
  const std::size_t n = 3 * impl->nV();
  for (std::size_t i = 0; i < n; ++i)
    impl->externalForce[i] += force_xyz[i];
}

void mem3dg_system_get_external_forces(const mem3dg_system *h,
                                       double *out_xyz) {
  auto *impl = reinterpret_cast<const Impl *>(h);
  const std::size_t n = 3 * impl->system->geometry.mesh->nVertices();
  std::copy(impl->externalForce.begin(), impl->externalForce.begin() + n,
            out_xyz);
}

namespace {

/**
 Count how many times a ray from `q` along `dir` crosses the surface.
 Moller-Trumbore per triangle. 
 */
int rayCrossings(gcs::ManifoldSurfaceMesh &mesh,
                 gcs::VertexPositionGeometry &vpg, const Eigen::Vector3d &q,
                 const Eigen::Vector3d &dir, bool &grazing) {
  const double EPS = 1e-12;
  const double EDGE_TOL = 1e-7; // barycentric margin considered "on an edge"
  int hits = 0;
  grazing = false;

  for (std::size_t fi = 0; fi < mesh.nFaces(); ++fi) {
    gcs::Halfedge he = mesh.face(fi).halfedge();
    const gc::Vector3 A = vpg.inputVertexPositions[he.vertex()];
    const gc::Vector3 B = vpg.inputVertexPositions[he.next().vertex()];
    const gc::Vector3 C = vpg.inputVertexPositions[he.next().next().vertex()];
    const Eigen::Vector3d v0(A.x, A.y, A.z), v1(B.x, B.y, B.z),
        v2(C.x, C.y, C.z);

    const Eigen::Vector3d e1 = v1 - v0, e2 = v2 - v0;
    const Eigen::Vector3d pv = dir.cross(e2);
    const double det = e1.dot(pv);
    if (std::abs(det) < EPS)
      continue; // ray parallel to the triangle

    const double inv = 1.0 / det;
    const Eigen::Vector3d tv = q - v0;
    const double u = inv * tv.dot(pv);
    if (u < -EDGE_TOL || u > 1.0 + EDGE_TOL)
      continue;

    const Eigen::Vector3d qv = tv.cross(e1);
    const double v = inv * dir.dot(qv);
    if (v < -EDGE_TOL || u + v > 1.0 + EDGE_TOL)
      continue;

    const double t = inv * e2.dot(qv);
    if (t <= EPS)
      continue; // behind the origin

    // a hit within EDGE_TOL of the boundary may be counted 0 or 2 times
    if (u < EDGE_TOL || v < EDGE_TOL || u + v > 1.0 - EDGE_TOL)
      grazing = true;

    ++hits;
  }
  return hits;
}

/// robust inside test: retry with other directions if a cast grazes an edge
bool pointInside(gcs::ManifoldSurfaceMesh &mesh,
                 gcs::VertexPositionGeometry &vpg, const Eigen::Vector3d &q) {
  
  static const Eigen::Vector3d dirs[3] = {
      Eigen::Vector3d(0.5773502692, 0.5772156649, 0.5774996669).normalized(),
      Eigen::Vector3d(0.2617993878, -0.9510565163, 0.1650670970).normalized(),
      Eigen::Vector3d(-0.8090169944, 0.3141592654, 0.4972893300).normalized()};

  for (int k = 0; k < 3; ++k) {
    bool grazing = false;
    const int n = rayCrossings(mesh, vpg, q, dirs[k], grazing);
    if (!grazing)
      return (n % 2) == 1;
  }
  // every direction grazed: fall back to the last count rather than failing
  bool g = false;
  return (rayCrossings(mesh, vpg, q, dirs[0], g) % 2) == 1;
}

} // namespace


namespace {

/// (re)build the vertex KD-tree if the geometry changed since the last query
void refreshKnn(Impl *impl) {
  if (!impl->knnStale && impl->knn)
    return;
  gcs::ManifoldSurfaceMesh &mesh = *impl->system->geometry.mesh;
  gcs::VertexPositionGeometry &vpg = *impl->system->geometry.vpg;
  const std::size_t nV = mesh.nVertices();
  impl->knnPoints.resize(nV);
  for (std::size_t i = 0; i < nV; ++i)
    impl->knnPoints[i] = vpg.inputVertexPositions[mesh.vertex(i)];
  impl->knn.reset(new gc::NearestNeighborFinder(impl->knnPoints));
  impl->knnStale = false;
}

/// result of a closest-point query on the surface
struct SurfaceHit {
  Eigen::Vector3d point;     ///< closest point on the surface
  Eigen::Vector3d normal;    ///< OUTWARD unit normal there (pseudonormal)
  double          distance;  ///< unsigned distance to it
  gcs::Face       face;      ///< face carrying the closest point
  std::array<double, 3> w;   ///< barycentric weights of `point` on that face
  bool            valid = false;
};

/**
 Closest point on the surface to `q`, accelerated by the vertex KD-tree.

 Instead of scanning every face, take the K nearest VERTICES and test only the
 faces incident to them.  For a well-shaped mesh the closest face touches one of
 the closest few vertices, so K=8 is ample; the brute-force path is kept in the
 tests.

 The returned normal is the ANGLE-WEIGHTED PSEUDONORMAL, not the raw face
 normal.
 */
SurfaceHit closestOnSurface(Impl *impl, const Eigen::Vector3d &q,
                            std::size_t K = 8) {
  refreshKnn(impl);
  gcs::ManifoldSurfaceMesh &mesh = *impl->system->geometry.mesh;
  gcs::VertexPositionGeometry &vpg = *impl->system->geometry.vpg;

  SurfaceHit best;
  best.distance = std::numeric_limits<double>::infinity();

  const std::size_t nV = mesh.nVertices();
  std::vector<size_t> near =
      impl->knn->kNearest(gc::Vector3{q.x(), q.y(), q.z()},
                          std::min(K, nV ? nV - 1 : 0));
  // kNearest excludes nothing here, but be safe and always test the very
  // nearest vertex's own faces
  if (near.empty() && nV)
    near.push_back(0);

  // one generation per query; a face is tested at most once even though it is
  // reachable from each of its three vertices
  if (impl->faceStamp.size() < mesh.nFaces())
    impl->faceStamp.assign(mesh.nFaces(), 0);
  const unsigned gen = ++impl->stampGen;

  auto testFace = [&](gcs::Face f) {
    const size_t fi = f.getIndex();
    if (fi < impl->faceStamp.size()) {
      if (impl->faceStamp[fi] == gen)
        return;
      impl->faceStamp[fi] = gen;
    }
    gcs::Halfedge he = f.halfedge();
    const gc::Vector3 A = vpg.inputVertexPositions[he.vertex()];
    const gc::Vector3 B = vpg.inputVertexPositions[he.next().vertex()];
    const gc::Vector3 C = vpg.inputVertexPositions[he.next().next().vertex()];
    const std::array<Eigen::Vector3d, 3> fv = {Eigen::Vector3d(A.x, A.y, A.z),
                                               Eigen::Vector3d(B.x, B.y, B.z),
                                               Eigen::Vector3d(C.x, C.y, C.z)};
    const coupled::EdgeFaceContact ct =
        coupled::minimumEdgeFaceContact(q, q, fv);
    if (ct.distance < best.distance) {
      best.distance = ct.distance;
      best.point = ct.pFace;
      best.face = f;
      best.w = ct.wFace;
      best.valid = true;
    }
  };

  for (size_t vi : near)
    for (gcs::Face f : mesh.vertex(vi).adjacentFaces())
      testFace(f);

  /*
   Robustness pass, run ONLY when it can matter.  The closest face need not
   touch one of the k closest VERTICES -- a long thin triangle can be nearer
   than any of its corners -- and on a remeshed membrane this missed the true
   closest face on ~0.1% of queries.

   A face can only beat the current best if one of its vertices lies within
   (bestDistance + that face's longest edge) of the query.  Rather than track
   edge lengths, expand over the 1-ring of the vertices already visited, but
   only when the nearest vertex is further away than the best surface distance
   -- i.e. when the closest point is plausibly on a face whose corners are all
   comparatively distant.  On well-shaped meshes this almost never triggers.
   */
  if (best.valid && !near.empty()) {
    const gc::Vector3 nv = impl->knnPoints[near.front()];
    const double dNearestVert =
        (Eigen::Vector3d(nv.x, nv.y, nv.z) - q).norm();
    if (dNearestVert > best.distance * 1.000001) {
      for (size_t vi : near)
        for (gcs::Vertex w : mesh.vertex(vi).adjacentVertices())
          for (gcs::Face f : w.adjacentFaces())
            testFace(f);
    }
  }

  if (!best.valid)
    return best;

  /*
   Angle-weighted pseudonormal at the closest point.  The barycentric weights
   tell us where it landed: interior of the face (all three weights positive),
   on an edge (one weight ~0), or at a vertex (two weights ~0).  Averaging the
   incident face normals with the appropriate weights covers all three cases in
   one expression: contributions come from the faces around whichever vertices
   carry weight, scaled by their corner angle.
   */
  const double TOL = 1e-9;
  Eigen::Vector3d n = Eigen::Vector3d::Zero();
  gcs::Halfedge hv = best.face.halfedge();
  int nCarrying = 0;
  for (int k = 0; k < 3; ++k) {
    if (best.w[k] > TOL)
      ++nCarrying;
    hv = hv.next();
  }

  if (nCarrying == 3) {
    // interior of the face: the face normal is unambiguous
    const gc::Vector3 fn = vpg.faceNormals[best.face];
    n = Eigen::Vector3d(fn.x, fn.y, fn.z);
  } else {
    // on an edge or at a vertex: average the incident faces, weighted by the
    // corner angle at the vertices that carry barycentric weight
    hv = best.face.halfedge();
    for (int k = 0; k < 3; ++k) {
      if (best.w[k] > TOL) {
        gcs::Vertex v = hv.vertex();
        for (gcs::Corner c : v.adjacentCorners()) {
          const gc::Vector3 fn = vpg.faceNormals[c.face()];
          const double ang = vpg.cornerAngle(c);
          n += ang * Eigen::Vector3d(fn.x, fn.y, fn.z);
        }
      }
      hv = hv.next();
    }
  }

  const double nn = n.norm();
  if (nn > 1e-30)
    best.normal = n / nn;
  else {
    const gc::Vector3 fn = vpg.faceNormals[best.face];
    best.normal = Eigen::Vector3d(fn.x, fn.y, fn.z).normalized();
  }
  return best;
}

} // namespace

int mem3dg_system_inside(const mem3dg_system *h, double x, double y, double z) {
  if (!h)
    return -1;
  auto *impl = const_cast<Impl *>(reinterpret_cast<const Impl *>(h));
  if (impl->system->geometry.mesh->hasBoundary())
    return -1; // "inside" is undefined for an open surface

  const Eigen::Vector3d q(x, y, z);
  const SurfaceHit hit = closestOnSurface(impl, q);
  if (!hit.valid)
    return -1;

  /*
   Sign of the projection onto the ANGLE-WEIGHTED PSEUDONORMAL at the closest
   point.  This is exact for a closed manifold and costs nothing beyond the
   closest-point query we already did -- unlike ray casting, which needed a
   separate O(nFaces) sweep (and retries when a ray grazed an edge).
   */
  const double sd = (q - hit.point).dot(hit.normal);
  return (sd < 0.0) ? 1 : 0;
}

int mem3dg_system_project(const mem3dg_system *h, double x, double y, double z,
                          double *out_point, double *out_normal) {
  if (!h)
    return -1;
  auto *impl = const_cast<Impl *>(reinterpret_cast<const Impl *>(h));
  const Eigen::Vector3d q(x, y, z);
  const SurfaceHit hit = closestOnSurface(impl, q);
  if (!hit.valid)
    return -1;
  if (out_point)
    for (int k = 0; k < 3; ++k)
      out_point[k] = hit.point[k];
  if (out_normal)
    for (int k = 0; k < 3; ++k)
      out_normal[k] = hit.normal[k];
  return 0;
}

void mem3dg_system_bounds(const mem3dg_system *h, double *out_inf,
                          double *out_sup) {
  auto *impl = reinterpret_cast<const Impl *>(h);
  gcs::ManifoldSurfaceMesh &mesh = *impl->system->geometry.mesh;
  gcs::VertexPositionGeometry &vpg = *impl->system->geometry.vpg;
  const std::size_t nV = mesh.nVertices();
  if (nV == 0) {
    for (int k = 0; k < 3; ++k) { out_inf[k] = -1; out_sup[k] = 1; }
    return;
  }
  gc::Vector3 lo = vpg.inputVertexPositions[mesh.vertex(0)], hi = lo;
  for (std::size_t i = 1; i < nV; ++i) {
    gc::Vector3 p = vpg.inputVertexPositions[mesh.vertex(i)];
    lo = gc::Vector3{std::min(lo.x,p.x), std::min(lo.y,p.y), std::min(lo.z,p.z)};
    hi = gc::Vector3{std::max(hi.x,p.x), std::max(hi.y,p.y), std::max(hi.z,p.z)};
  }
  out_inf[0]=lo.x; out_inf[1]=lo.y; out_inf[2]=lo.z;
  out_sup[0]=hi.x; out_sup[1]=hi.y; out_sup[2]=hi.z;
}

long mem3dg_system_filament_confinement(mem3dg_system *h, const double *points,
                                        size_t nPts, double shell,
                                        double stiffness, double *out_f,
                                        int add_to_membrane) {
  if (!h || !points || shell < 0.0)
    return -1;

  auto *impl = reinterpret_cast<Impl *>(h);
  gcs::ManifoldSurfaceMesh &mesh = *impl->system->geometry.mesh;
  gcs::VertexPositionGeometry &vpg = *impl->system->geometry.vpg;

  // "Inside" only means something on a closed surface.
  if (mesh.hasBoundary())
    return -1;

  const std::size_t nF = mesh.nFaces();

  if (out_f)
    std::fill(out_f, out_f + 3 * nPts, 0.0);

  long nForced = 0;

  for (std::size_t i = 0; i < nPts; ++i) {
    const Eigen::Vector3d q(points[3 * i + 0], points[3 * i + 1],
                            points[3 * i + 2]);

    // ---- closest point on the surface -------------------------------------
    double bestAbs = std::numeric_limits<double>::infinity();
    coupled::EdgeFaceContact bestCt;
    gcs::Face bestFace;
    bool found = false;

    for (std::size_t fi = 0; fi < nF; ++fi) {
      gcs::Face face = mesh.face(fi);
      gcs::Halfedge he = face.halfedge();
      const gc::Vector3 A = vpg.inputVertexPositions[he.vertex()];
      const gc::Vector3 B = vpg.inputVertexPositions[he.next().vertex()];
      const gc::Vector3 C = vpg.inputVertexPositions[he.next().next().vertex()];
      const std::array<Eigen::Vector3d, 3> fv = {
          Eigen::Vector3d(A.x, A.y, A.z), Eigen::Vector3d(B.x, B.y, B.z),
          Eigen::Vector3d(C.x, C.y, C.z)};

      // prune: this face cannot beat the best found so far
      const Eigen::Vector3d ctr = (fv[0] + fv[1] + fv[2]) / 3.0;
      const double faceRad =
          std::max({(fv[0] - ctr).norm(), (fv[1] - ctr).norm(),
                    (fv[2] - ctr).norm()});
      if ((ctr - q).norm() - faceRad > bestAbs)
        continue;

      // degenerate "segment" == the point itself
      const coupled::EdgeFaceContact ct =
          coupled::minimumEdgeFaceContact(q, q, fv);
      if (ct.distance < bestAbs) {
        bestAbs = ct.distance;
        bestCt = ct;
        bestFace = face;
        found = true;
      }
    }
    if (!found)
      continue;

    // ---- which side, by ray casting ---------------------------------------
    const bool isOutside = !pointInside(mesh, vpg, q);

    /*
     Outward axis taken from the closest-point pair rather than a face normal:
     (q - pFace) points from the surface toward the query point, so flipping it
     by the side test gives a consistently outward direction.
     */
    Eigen::Vector3d sep = q - bestCt.pFace;
    double sepLen = sep.norm();
    Eigen::Vector3d nOut;
    if (sepLen > 1e-12) {
      nOut = (isOutside ? 1.0 : -1.0) * (sep / sepLen);
    } else {
      // exactly on the surface: the pair degenerates, use the face normal
      const gc::Vector3 gn = vpg.faceNormals[bestFace];
      nOut = Eigen::Vector3d(gn.x, gn.y, gn.z);
      const double nn = nOut.norm();
      if (nn < 1e-30)
        continue;
      nOut /= nn;
      sepLen = 0.0;
    }

    const double sd = (isOutside ? 1.0 : -1.0) * sepLen; // + outside
    const double pen = sd + shell;                       // wall at sd = -shell
    if (pen <= 0.0)
      continue;                                          // deep inside

    ++nForced;

    const Eigen::Vector3d Fin = (-stiffness * pen) * nOut; // inward

    if (out_f)
      for (int k = 0; k < 3; ++k)
        out_f[3 * i + k] += Fin[k];

    // membrane: equal and opposite, spread by barycentric weights
    if (add_to_membrane) {
      gcs::Halfedge hv = bestFace.halfedge();
      for (int k = 0; k < 3; ++k) {
        const std::size_t vi = hv.vertex().getIndex();
        for (int j = 0; j < 3; ++j)
          impl->externalForce[3 * vi + j] -= bestCt.wFace[k] * Fin[j];
        hv = hv.next();
      }
    }
  }

  return nForced;
}

long mem3dg_system_filament_repulsion(mem3dg_system *h, const double *seg_p0,
                                      const double *seg_p1, size_t nSeg,
                                      double cutoff, double stiffness,
                                      double *out_f0, double *out_f1,
                                      int add_to_membrane) {
  if (!h || !seg_p0 || !seg_p1 || cutoff <= 0.0)
    return -1;

  auto *impl = reinterpret_cast<Impl *>(h);
  gcs::ManifoldSurfaceMesh &mesh = *impl->system->geometry.mesh;
  gcs::VertexPositionGeometry &vpg = *impl->system->geometry.vpg;

  const std::size_t nF = mesh.nFaces();

  // outputs are overwritten
  if (out_f0)
    std::fill(out_f0, out_f0 + 3 * nSeg, 0.0);
  if (out_f1)
    std::fill(out_f1, out_f1 + 3 * nSeg, 0.0);

  long nContacts = 0;

  for (std::size_t s = 0; s < nSeg; ++s) {
    const Eigen::Vector3d p0(seg_p0[3 * s + 0], seg_p0[3 * s + 1],
                             seg_p0[3 * s + 2]);
    const Eigen::Vector3d p1(seg_p1[3 * s + 0], seg_p1[3 * s + 1],
                             seg_p1[3 * s + 2]);

    // cheap bound: nothing on this segment can reach a face whose centroid is
    // farther than (half the segment length + cutoff + face radius) away.
    const Eigen::Vector3d segMid = 0.5 * (p0 + p1);
    const double segReach = 0.5 * (p1 - p0).norm() + cutoff;

    for (std::size_t fi = 0; fi < nF; ++fi) {
      gcs::Face face = mesh.face(fi);
      gcs::Halfedge he = face.halfedge();
      const gc::Vector3 a = vpg.inputVertexPositions[he.vertex()];
      const gc::Vector3 b = vpg.inputVertexPositions[he.next().vertex()];
      const gc::Vector3 c = vpg.inputVertexPositions[he.next().next().vertex()];

      const std::array<Eigen::Vector3d, 3> fv = {Eigen::Vector3d(a.x, a.y, a.z),
                                                 Eigen::Vector3d(b.x, b.y, b.z),
                                                 Eigen::Vector3d(c.x, c.y, c.z)};

      // quick reject on centroid distance; the face's own extent is bounded by
      // its longest edge, so add that before comparing
      const Eigen::Vector3d ctr = (fv[0] + fv[1] + fv[2]) / 3.0;
      const double faceRad =
          std::max({(fv[0] - ctr).norm(), (fv[1] - ctr).norm(),
                    (fv[2] - ctr).norm()});
      if ((ctr - segMid).norm() > segReach + faceRad)
        continue;

      const coupled::EdgeFaceContact ct =
          coupled::minimumEdgeFaceContact(p0, p1, fv);

      const double d = ct.distance;
      if (d >= cutoff)
        continue;

      ++nContacts;

      /*
       Direction: from the membrane contact point towards the segment, so the
       filament is pushed away from the surface.  When the segment pierces the
       face the two contact points coincide (d == 0) and that direction is
       undefined; fall back to the face normal, which is the only meaningful
       "away" there.
       */
      Eigen::Vector3d nhat;
      if (d > 1e-12) {
        nhat = (ct.pSeg - ct.pFace) / d;
      } else {
        const Eigen::Vector3d fn = (fv[1] - fv[0]).cross(fv[2] - fv[0]);
        const double fnN = fn.norm();
        if (fnN < 1e-30)
          continue; // degenerate face
        nhat = fn / fnN;
      }

      // harmonic: full strength at contact, zero at the cutoff
      const Eigen::Vector3d F = (stiffness * (cutoff - d)) * nhat;

      // filament side: split along the segment by the parametric coordinate
      const double w0 = 1.0 - ct.s, w1 = ct.s;
      if (out_f0)
        for (int k = 0; k < 3; ++k)
          out_f0[3 * s + k] += w0 * F[k];
      if (out_f1)
        for (int k = 0; k < 3; ++k)
          out_f1[3 * s + k] += w1 * F[k];

      // membrane side: equal and opposite, split by barycentric weights.
      // Since w0+w1 == 1 and wFace sums to 1, the totals cancel exactly.
      if (add_to_membrane) {
        gcs::Halfedge hv = face.halfedge();
        for (int k = 0; k < 3; ++k) {
          const std::size_t vi = hv.vertex().getIndex();
          for (int j = 0; j < 3; ++j)
            impl->externalForce[3 * vi + j] -= ct.wFace[k] * F[j];
          hv = hv.next();
        }
      }
    }
  }

  return nContacts;
}

void mem3dg_system_step(mem3dg_system *h, double dt) {
  auto *impl = reinterpret_cast<Impl *>(h);

  /*
   Remesh BEFORE forcing, so the step that follows runs on a mesh that is
   already consistent.  mutateMesh() can change the vertex count, which
   invalidates the external force buffer: it is sized per vertex and its entries
   refer to vertices that may no longer exist.  Resize and clear it -- the
   caller re-supplies forces every step anyway.
   */
  if (impl->remesh && (impl->stepCount % impl->remeshPeriod) == 0 &&
      impl->stepCount > 0) {
    impl->system->mutateMesh();
    impl->system->updateConfigurations();
    const std::size_t n = 3 * impl->nV();
    impl->externalForce.assign(n, 0.0);
  }
  ++impl->stepCount;
  // every step moves vertices, so the index built for the previous geometry
  // no longer reflects the mesh
  impl->knnStale = true;
  const double step = (dt > 0 ? dt : impl->dt);
  impl->integrator->characteristicTimeStep = step;

  // 1. System's own forcing (zeroes + refills mechanicalForceVec).
  impl->system->computeConservativeForcing();
  impl->system->addNonconservativeForcing(step);

  // 2. Add the injected external force buffer.
  gcs::ManifoldSurfaceMesh &mesh = *impl->system->geometry.mesh;
  auto &mechVec = impl->system->forces.mechanicalForceVec;
  auto &extVec = impl->system->forces.externalForceVec;
  const std::size_t nV = mesh.nVertices();
  for (std::size_t i = 0; i < nV; ++i) {
    gc::Vector3 f{impl->externalForce[3 * i + 0], impl->externalForce[3 * i + 1],
                  impl->externalForce[3 * i + 2]};
    gcs::Vertex v = mesh.vertex(i);
    mechVec[v] += f;
    extVec[v] += f; // reflected in output / diagnostics
  }

  // 3. One forward-Euler step (consumes mechanicalForceVec).
  impl->integrator->march();
}

double mem3dg_system_time(const mem3dg_system *h) {
  return reinterpret_cast<const Impl *>(h)->system->time;
}

int mem3dg_system_save_ply(mem3dg_system *h, const char *out_path) {
  try {
    reinterpret_cast<Impl *>(h)->system->geometry.saveGeometry(
        std::string(out_path));
    return 0;
  } catch (...) {
    return 1;
  }
}

} // extern "C"
