// =============================================================================
// face_bvh.cpp  --  see face_bvh.h for the design rationale.
// =============================================================================

#include "mem3dg/solver/face_bvh.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace mem3dg {
namespace solver {

// ─────────────────────────────────────────────────────────────────────────────
// Aabb
// ─────────────────────────────────────────────────────────────────────────────

void Aabb::reset() {
  const double inf = std::numeric_limits<double>::infinity();
  lo = gc::Vector3{inf, inf, inf};
  hi = gc::Vector3{-inf, -inf, -inf};
}

void Aabb::expand(const gc::Vector3 &p) {
  lo.x = std::min(lo.x, p.x);
  lo.y = std::min(lo.y, p.y);
  lo.z = std::min(lo.z, p.z);
  hi.x = std::max(hi.x, p.x);
  hi.y = std::max(hi.y, p.y);
  hi.z = std::max(hi.z, p.z);
}

void Aabb::expand(const Aabb &o) {
  expand(o.lo);
  expand(o.hi);
}

void Aabb::pad(double r) {
  lo -= gc::Vector3{r, r, r};
  hi += gc::Vector3{r, r, r};
}

bool Aabb::overlaps(const Aabb &o) const {
  return !(o.lo.x > hi.x || o.hi.x < lo.x || o.lo.y > hi.y || o.hi.y < lo.y ||
           o.lo.z > hi.z || o.hi.z < lo.z);
}

// ─────────────────────────────────────────────────────────────────────────────
// build
// ─────────────────────────────────────────────────────────────────────────────

bool FaceBVH::needsRebuild(const gcs::ManifoldSurfaceMesh &mesh) const {
  return nodes_.empty() || faceIndex_.size() != mesh.nFaces();
}

void FaceBVH::build(gcs::ManifoldSurfaceMesh &mesh,
                    const gcs::VertexPositionGeometry &vpg) {
  const std::size_t nF = mesh.nFaces();
  nodes_.clear();
  faceIndex_.resize(nF);
  std::iota(faceIndex_.begin(), faceIndex_.end(), std::size_t{0});

  // Cache the three corner VERTEX indices per face once.  refit() runs every
  // step and walking halfedges there would dominate its cost; this turns the
  // refit inner loop into three array lookups.
  corner_.assign(3 * nF, 0);
  centroid_.assign(nF, gc::Vector3{0, 0, 0});
  for (std::size_t f = 0; f < nF; ++f) {
    gcs::Halfedge he = mesh.face(f).halfedge();
    const std::size_t a = he.vertex().getIndex();
    const std::size_t b = he.next().vertex().getIndex();
    const std::size_t c = he.next().next().vertex().getIndex();
    corner_[3 * f + 0] = a;
    corner_[3 * f + 1] = b;
    corner_[3 * f + 2] = c;
    centroid_[f] = (vpg.inputVertexPositions[mesh.vertex(a)] +
                    vpg.inputVertexPositions[mesh.vertex(b)] +
                    vpg.inputVertexPositions[mesh.vertex(c)]) /
                   3.0;
  }

  if (nF > 0) {
    nodes_.reserve(2 * (nF / kLeafSize + 1));
    buildRange(0, nF);
    refit(vpg);
  }
  centroid_.clear();
  centroid_.shrink_to_fit();
}

int FaceBVH::buildRange(std::size_t first, std::size_t count) {
  const int self = static_cast<int>(nodes_.size());
  nodes_.emplace_back();
  nodes_[self].first = first;
  nodes_[self].count = count;

  if (count <= kLeafSize) {
    nodes_[self].left = -1;
    return self;
  }

  // split on the widest axis of the centroid spread, at the median.  Median
  // rather than midpoint keeps the tree balanced when triangle density is
  // uneven -- which is the case this structure exists to handle.
  Aabb cb;
  cb.reset();
  for (std::size_t i = first; i < first + count; ++i)
    cb.expand(centroid_[faceIndex_[i]]);
  const gc::Vector3 ext = cb.hi - cb.lo;
  const int axis = (ext.x >= ext.y && ext.x >= ext.z) ? 0
                   : (ext.y >= ext.z)                 ? 1
                                                      : 2;

  const std::size_t mid = count / 2;
  auto begin = faceIndex_.begin() + static_cast<long>(first);
  std::nth_element(begin, begin + static_cast<long>(mid),
                   begin + static_cast<long>(count),
                   [this, axis](std::size_t p, std::size_t q) {
                     const gc::Vector3 &a = centroid_[p];
                     const gc::Vector3 &b = centroid_[q];
                     return (axis == 0) ? a.x < b.x
                            : (axis == 1) ? a.y < b.y
                                          : a.z < b.z;
                   });

  const int l = buildRange(first, mid);
  const int r = buildRange(first + mid, count - mid);
  nodes_[self].left = l;
  nodes_[self].right = r;
  return self;
}

// ─────────────────────────────────────────────────────────────────────────────
// refit
// ─────────────────────────────────────────────────────────────────────────────

void FaceBVH::refit(const gcs::VertexPositionGeometry &vpg) {
  if (nodes_.empty())
    return;
  refitNode(0, vpg);
}

void FaceBVH::refitNode(int node, const gcs::VertexPositionGeometry &vpg) {
  Node &n = nodes_[static_cast<std::size_t>(node)];
  n.box.reset();
  if (n.left < 0) {
    const auto &pos = vpg.inputVertexPositions;
    for (std::size_t i = n.first; i < n.first + n.count; ++i) {
      const std::size_t f = faceIndex_[i];
      for (int k = 0; k < 3; ++k)
        n.box.expand(pos[static_cast<std::size_t>(corner_[3 * f + k])]);
    }
    return;
  }
  const int l = n.left, r = n.right;
  refitNode(l, vpg);
  refitNode(r, vpg);
  Node &self = nodes_[static_cast<std::size_t>(node)];
  self.box.reset();
  self.box.expand(nodes_[static_cast<std::size_t>(l)].box);
  self.box.expand(nodes_[static_cast<std::size_t>(r)].box);
}

// ─────────────────────────────────────────────────────────────────────────────
// query
// ─────────────────────────────────────────────────────────────────────────────

void FaceBVH::query(const Aabb &box, std::vector<std::size_t> &out) const {
  if (nodes_.empty())
    return;
  /*
   Explicit stack: this query runs once per face of the other membrane, so the
   recursion overhead is worth removing.

   Depth bound: the split is at the MEDIAN, so the tree is balanced and the
   stack holds at most one entry per level -- ceil(log2(nF/kLeafSize)) + 1,
   which is 9 for 1280 faces and under 50 for anything that fits in memory.
   kMaxDepth is far above that.

   If it ever did fill, the fallback emits the node's whole face range instead
   of dropping it.  That returns a SUPERSET of the true candidates -- the narrow
   phase rejects the extras -- whereas dropping the node would lose real
   contacts, which is a wrong answer rather than a slow one.
  */
  constexpr int kMaxDepth = 64;
  int stack[kMaxDepth];
  int sp = 0;
  stack[sp++] = 0;
  while (sp > 0) {
    const Node &n = nodes_[static_cast<std::size_t>(stack[--sp])];
    if (!n.box.overlaps(box))
      continue;
    if (n.left < 0) {
      for (std::size_t i = n.first; i < n.first + n.count; ++i)
        out.push_back(faceIndex_[i]);
    } else if (sp + 2 <= kMaxDepth) {
      stack[sp++] = n.left;
      stack[sp++] = n.right;
    } else {
      for (std::size_t i = n.first; i < n.first + n.count; ++i)
        out.push_back(faceIndex_[i]);
    }
  }
}

} // namespace solver
} // namespace mem3dg
