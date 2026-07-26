#pragma once

#include "types.h"
#include <Eigen/Dense>

namespace coupled {

// Abstract interface for the filament–membrane steric interaction.
//
// Subclass this and implement compute() with the specific potential you want
// (e.g. WCA, harmonic repulsion, Lennard-Jones). The coupling layer calls
// compute() for every (segment, vertex) pair whose closest-point distance
// is below `cutoff`, then accumulates the returned forces.
//
// Example subclass skeleton:
//
//   class HarmonicKernel : public StericKernel {
//   public:
//       double k;   // spring constant
//       StericResult compute(const FilamentSegment& seg,
//                            const Eigen::Vector3d& vertex_pos,
//                            const Eigen::Vector3d& vertex_nor,
//                            double cutoff) const override
//       {
//           // ----- USER: implement your steric potential here -----
//
//           // 1. Find closest point on segment to vertex.
//           //    t=0 → p0, t=1 → p1; 'closest' is the foot of perpendicular.
//           Eigen::Vector3d seg_vec  = seg.p1 - seg.p0;
//           double len2 = seg_vec.squaredNorm();
//           double t = (len2 < 1e-14) ? 0.0
//                      : std::clamp((vertex_pos - seg.p0).dot(seg_vec) / len2,
//                                   0.0, 1.0);
//           Eigen::Vector3d closest = seg.p0 + t * seg_vec;
//
//           // 2. Separation vector: points FROM filament TOWARD membrane vertex.
//           Eigen::Vector3d delta = vertex_pos - closest;
//           double dist  = delta.norm();
//           double d0    = cutoff;             // or seg.radius + membrane_thickness
//
//           // 3. Only repel if inside cutoff and not degenerate.
//           if (dist >= d0 || dist < 1e-10) return {};   // active=false → no force
//
//           // 4. YOUR POTENTIAL HERE (example: harmonic):
//           double overlap = d0 - dist;
//           Eigen::Vector3d n_hat = delta / dist;         // unit normal, fil→vertex
//           Eigen::Vector3d F_on_vertex = k * overlap * n_hat;
//
//           // 5. Split filament-side force (Newton's 3rd) between p0 and p1.
//           StericResult res;
//           res.active          = true;
//           res.force_on_vertex = F_on_vertex;
//           res.force_on_p0     = -(1.0 - t) * F_on_vertex;
//           res.force_on_p1     = -       t  * F_on_vertex;
//           return res;
//       }
//   };
class StericKernel {
public:
    virtual ~StericKernel() = default;

    // Called for every (segment, vertex) pair closer than `cutoff`.
    //
    // Parameters
    // ----------
    // seg        : the filament segment (p0, p1, radius, owning mecable)
    // vertex_pos : 3-D position of the membrane vertex
    // vertex_nor : outward unit normal at the membrane vertex
    // cutoff     : distance threshold beyond which no force is applied
    //
    // Return a StericResult with active=false if the pair should be skipped.
    // Otherwise set all three force fields and active=true.
    virtual StericResult compute(const FilamentSegment&  seg,
                                 const Eigen::Vector3d&  vertex_pos,
                                 const Eigen::Vector3d&  vertex_nor,
                                 double                  cutoff) const = 0;
};

} // namespace coupled
