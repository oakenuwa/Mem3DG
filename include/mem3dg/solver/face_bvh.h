#pragma once

// =============================================================================
// face_bvh.h
//
// Bounding-volume hierarchy over the faces of one membrane, used to accelerate
// membrane-membrane contact search.
//
// WHY A BVH RATHER THAN A UNIFORM GRID
// ────────────────────────────────────
// For a fixed-radius query a uniform grid is normally the better structure --
// O(1) bucket lookup, trivial to build.  It degrades, though, when triangle
// sizes vary a lot, because the cell size has to suit the LARGEST triangle while
// the occupancy is set by the smallest.  Mem3DG's mesh mutator (splitLarge /
// collapseSmall / flipNonDelaunay) deliberately produces exactly that spread:
// after remeshing a protrusion, face areas across one mesh differ by an order of
// magnitude.  A BVH adapts to whatever distribution it is given.
//
// REFIT, NOT REBUILD
// ──────────────────
// Vertex positions change every step but the TOPOLOGY does not (except when the
// mutator runs).  So the tree is built once and thereafter refitted: leaf boxes
// are recomputed from current positions and interior boxes recomputed bottom-up,
// which is O(nFaces) with no allocation and no re-sorting.  A full rebuild is
// O(nFaces log nFaces) and only needed when the face count changes -- which is
// how remeshing is detected here (see needsRebuild()).
//
// The tree stores FACE INDICES, so a refit after remeshing would silently
// reference faces that no longer exist; always check needsRebuild() first.
// =============================================================================

#include "mem3dg/macros.h"

#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <cstddef>
#include <vector>

namespace mem3dg {
namespace solver {

namespace gc = geometrycentral;
namespace gcs = geometrycentral::surface;

/// Axis-aligned bounding box.
struct DLL_PUBLIC Aabb {
  gc::Vector3 lo{0, 0, 0};
  gc::Vector3 hi{0, 0, 0};

  void reset();
  void expand(const gc::Vector3 &p);
  void expand(const Aabb &o);
  /// grow the box outward by `r` on every side
  void pad(double r);
  bool overlaps(const Aabb &o) const;
};

class DLL_PUBLIC FaceBVH {
public:
  /// Build the tree topology and fit it to the current positions.
  void build(gcs::ManifoldSurfaceMesh &mesh,
             const gcs::VertexPositionGeometry &vpg);

  /// Recompute all boxes from current positions, keeping the topology.
  /// Undefined if the mesh has been remeshed since build(); check
  /// needsRebuild() first.
  void refit(const gcs::VertexPositionGeometry &vpg);

  /// True when the tree cannot be trusted for this mesh any more -- i.e. it was
  /// never built, or the face count changed (remeshing).
  bool needsRebuild(const gcs::ManifoldSurfaceMesh &mesh) const;

  /// Append every face whose box overlaps `box` to `out`.  `out` is NOT cleared,
  /// so a caller can accumulate; it is the caller's job to clear between queries.
  void query(const Aabb &box, std::vector<std::size_t> &out) const;

  std::size_t faceCount() const { return faceIndex_.size(); }
  bool empty() const { return nodes_.empty(); }

private:
  struct Node {
    Aabb box;
    /// child indices into nodes_; -1 on both marks a leaf.
    ///
    /// BOTH are stored explicitly.  The right child is NOT at left+1: the left
    /// call emits its whole subtree first, so left+1 is the left child's own
    /// left child.  Assuming adjacency here silently searches the wrong subtree
    /// and loses contacts.
    int left = -1;
    int right = -1;
    /// leaf payload: range [first, first+count) into faceIndex_
    std::size_t first = 0;
    std::size_t count = 0;
  };

  /// recursive build over faceIndex_[first, first+count)
  int buildRange(std::size_t first, std::size_t count);
  void refitNode(int node, const gcs::VertexPositionGeometry &vpg);

  std::vector<Node> nodes_;
  /// face indices, permuted so every leaf owns a contiguous range
  std::vector<std::size_t> faceIndex_;
  /// cached triangle corner indices, so refit does not walk halfedges
  std::vector<std::size_t> corner_;   // 3 per face, indexed by FACE index
  /// centroids, scratch used during build only
  std::vector<gc::Vector3> centroid_;

  static constexpr std::size_t kLeafSize = 4;
};

} // namespace solver
} // namespace mem3dg
