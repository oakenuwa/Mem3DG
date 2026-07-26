// =============================================================================
// edge_face_contact.cpp
//
// Implements minimumEdgeFaceContact using Eigen::Vector3d throughout so that
// this code stays independent of the geometry-central build tree and can be
// compiled as part of the coupling layer without pulling in Mem3DG's full
// include path.
//
// The algorithms are the same as closestPointOnTriangle / closestSegmentSegment
// in face_repulsion.cpp but expressed in Eigen for consistency with the rest
// of the coupled layer.
// =============================================================================

#include "mem3dg/coupled/edge_face_contact.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace coupled {

// ─── 1  GEOMETRY HELPERS (Eigen versions) ───────────────────────────────────

// Closest point on triangle {a, b, c} to query point p.
// Returns (closest_point, {λa, λb, λc}) using Voronoi-region projection.
// Barycentric coords are non-negative and sum to 1.
// Algorithm: Ericson, "Real-Time Collision Detection", §5.1.5.
static std::pair<Eigen::Vector3d, std::array<double,3>>
closestPointOnTriangleE(const Eigen::Vector3d& a,
                        const Eigen::Vector3d& b,
                        const Eigen::Vector3d& c,
                        const Eigen::Vector3d& p)
{
    Eigen::Vector3d ab = b - a;
    Eigen::Vector3d ac = c - a;
    Eigen::Vector3d ap = p - a;

    double d1 = ab.dot(ap);
    double d2 = ac.dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0) return { a, {1.0, 0.0, 0.0} };

    Eigen::Vector3d bp = p - b;
    double d3 = ab.dot(bp);
    double d4 = ac.dot(bp);
    if (d3 >= 0.0 && d4 <= d3) return { b, {0.0, 1.0, 0.0} };

    Eigen::Vector3d cp = p - c;
    double d5 = ab.dot(cp);
    double d6 = ac.dot(cp);
    if (d6 >= 0.0 && d5 <= d6) return { c, {0.0, 0.0, 1.0} };

    // Voronoi region of edge ab
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        double v = d1 / (d1 - d3);
        return { a + v * ab, {1.0 - v, v, 0.0} };
    }

    // Voronoi region of edge ac
    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        double w = d2 / (d2 - d6);
        return { a + w * ac, {1.0 - w, 0.0, w} };
    }

    // Voronoi region of edge bc
    double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return { b + w * (c - b), {0.0, 1.0 - w, w} };
    }

    // Interior
    double denom = 1.0 / (va + vb + vc);
    double v = vb * denom;
    double w = vc * denom;
    return { a + v * ab + w * ac, {1.0 - v - w, v, w} };
}


// Closest points between segment [e0, e1] and segment [g0, g1].
// Returns (s*, t*) ∈ [0,1]^2.
// Algorithm: Ericson §5.1.9.
static std::pair<double,double>
closestSegmentSegmentE(const Eigen::Vector3d& e0, const Eigen::Vector3d& e1,
                       const Eigen::Vector3d& g0, const Eigen::Vector3d& g1)
{
    constexpr double eps = 1e-12;
    Eigen::Vector3d d1 = e1 - e0;
    Eigen::Vector3d d2 = g1 - g0;
    Eigen::Vector3d r  = e0 - g0;
    double a = d1.squaredNorm();
    double e = d2.squaredNorm();
    double f = d2.dot(r);

    double s, t;
    if (a < eps && e < eps) return {0.0, 0.0};

    if (a < eps) {
        s = 0.0;
        t = std::clamp(f / e, 0.0, 1.0);
    } else {
        double c = d1.dot(r);
        if (e < eps) {
            t = 0.0;
            s = std::clamp(-c / a, 0.0, 1.0);
        } else {
            double b     = d1.dot(d2);
            double denom = a * e - b * b;
            s = (std::abs(denom) > eps)
                    ? std::clamp((b * f - c * e) / denom, 0.0, 1.0)
                    : 0.0;
            t = (b * s + f) / e;
            if (t < 0.0) { t = 0.0; s = std::clamp(-c / a, 0.0, 1.0); }
            else if (t > 1.0) { t = 1.0; s = std::clamp((b - c) / a, 0.0, 1.0); }
        }
    }
    return {s, t};
}


// ─── 2  minimumEdgeFaceContact ──────────────────────────────────────────────

EdgeFaceContact minimumEdgeFaceContact(const Eigen::Vector3d&               p0,
                                       const Eigen::Vector3d&               p1,
                                       const std::array<Eigen::Vector3d,3>& fv)
{
    EdgeFaceContact best;   // distance = ∞

    // Helper: compare and keep the closer candidate.
    auto update = [&](double s_cand,
                      const Eigen::Vector3d& pSeg_cand,
                      const Eigen::Vector3d& pFace_cand,
                      const std::array<double,3>& w_cand)
    {
        double d = (pSeg_cand - pFace_cand).norm();
        if (d < best.distance) {
            best.s        = s_cand;
            best.pSeg     = pSeg_cand;
            best.pFace    = pFace_cand;
            best.wFace    = w_cand;
            best.distance = d;
        }
    };

    /*
     ── Sub-case 0: the segment PIERCES the face ─────────────────────────────
     The five sub-cases below (2 endpoint-vs-face + 3 segment-vs-edge) are
     complete only while the segment and the triangle do NOT intersect.  Proof
     sketch: if the closest points were interior to both, the connecting vector
     would have to be perpendicular to the segment and parallel to the face
     normal, i.e. the segment is parallel to the face -- and then its endpoints
     attain the same minimum, so sub-cases 1-2 already find it.

     When the segment passes THROUGH the triangle, however, the true distance
     is 0 while the five sub-cases return the distance to the nearest face
     EDGE.  For a filament crossing the middle of a large triangle that value
     is comparable to the triangle's size, so with any realistic cutoff the
     pair is rejected as "far apart" and the filament tunnels straight
     through -- exactly the failure this edge-face kernel exists to prevent.

     So test intersection first: if it crosses, the answer is 0 and nothing
     else can beat it.
     */
    {
        const Eigen::Vector3d e1 = fv[1] - fv[0];
        const Eigen::Vector3d e2 = fv[2] - fv[0];
        const Eigen::Vector3d n  = e1.cross(e2);
        const Eigen::Vector3d dir = p1 - p0;

        const double nd = n.dot(dir);
        if (std::abs(nd) > 1e-30)          // not parallel to the plane
        {
            // parametric crossing of the triangle's plane
            const double t = n.dot(fv[0] - p0) / nd;
            if (t >= 0.0 && t <= 1.0)
            {
                const Eigen::Vector3d P = p0 + t * dir;

                // barycentric coordinates of P within the triangle
                const Eigen::Vector3d vp = P - fv[0];
                const double d00 = e1.dot(e1), d01 = e1.dot(e2), d11 = e2.dot(e2);
                const double d20 = vp.dot(e1), d21 = vp.dot(e2);
                const double den = d00 * d11 - d01 * d01;
                if (std::abs(den) > 1e-30)
                {
                    const double v = (d11 * d20 - d01 * d21) / den;
                    const double w = (d00 * d21 - d01 * d20) / den;
                    const double u = 1.0 - v - w;
                    if (u >= 0.0 && v >= 0.0 && w >= 0.0)
                    {
                        best.s        = t;
                        best.pSeg     = P;
                        best.pFace    = P;
                        best.wFace    = {u, v, w};
                        best.distance = 0.0;
                        return best;    // nothing can be closer than 0
                    }
                }
            }
        }
    }

    // ── Sub-case 1: endpoint p0  vs  membrane face ───────────────────────────
    {
        auto [pF, wF] = closestPointOnTriangleE(fv[0], fv[1], fv[2], p0);
        // p0 has s=0; its "barycentric weight" on the segment is (1-s)=1 at p0.
        update(0.0, p0, pF, wF);
    }

    // ── Sub-case 2: endpoint p1  vs  membrane face ───────────────────────────
    {
        auto [pF, wF] = closestPointOnTriangleE(fv[0], fv[1], fv[2], p1);
        update(1.0, p1, pF, wF);
    }

    // ── Sub-cases 3,4,5: filament segment  vs  each face edge ─────────────────
    // Face edges (in halfedge order): v0–v1 (edge 0), v1–v2 (edge 1), v2–v0 (edge 2).
    constexpr int faceEdge[3][2] = {{0,1},{1,2},{2,0}};

    for (auto [ia, ib] : faceEdge) {
        auto [s, t] = closestSegmentSegmentE(p0, p1, fv[ia], fv[ib]);

        Eigen::Vector3d pSeg  = (1.0 - s) * p0 + s * p1;
        Eigen::Vector3d pFace = (1.0 - t) * fv[ia] + t * fv[ib];

        // Barycentric weight on the face: only vertices ia and ib participate;
        // the third vertex has weight 0.
        //   pFace = (1−t)·fv[ia] + t·fv[ib]
        //   ⟹  wFace[ia] = 1−t,  wFace[ib] = t,  wFace[third] = 0
        std::array<double,3> w{};
        w[ia] = 1.0 - t;
        w[ib] = t;
        // w[third] stays 0 (already zero-initialised by std::array{})

        update(s, pSeg, pFace, w);
    }

    return best;
}

} // namespace coupled
