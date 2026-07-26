#pragma once

#include <Eigen/Dense>

// Forward-declare Cytosim's Mecable so this header stays independent of
// the Cytosim build tree when included from Mem3DG-only translation units.
class Mecable;

namespace coupled {

// One segment of a Cytosim filament, extracted each coupled time step.
// Segments are indexed p0 (local point k) → p1 (local point k+1).
struct FilamentSegment {
    Eigen::Vector3d p0;
    Eigen::Vector3d p1;
    double          radius;    // steric radius from FiberProp::steric_radius
    Mecable const*  mecable;   // owning Cytosim Mecable — needed by Meca::addForce
    unsigned        point0_k;  // local point index of p0 within the fiber
};

// Forces returned by StericKernel::compute() for one (segment, vertex) pair.
//
// Newton's 3rd law is enforced by the caller:
//   force_on_vertex == -(force_on_p0 + force_on_p1)
//
// force_on_p0 and force_on_p1 split the filament-side force according to
// the barycentric coordinate t of the closest point on the segment:
//   force_on_p0 = (1-t) * F_fil
//   force_on_p1 =    t  * F_fil
struct StericResult {
    Eigen::Vector3d force_on_vertex = Eigen::Vector3d::Zero();
    Eigen::Vector3d force_on_p0     = Eigen::Vector3d::Zero();
    Eigen::Vector3d force_on_p1     = Eigen::Vector3d::Zero();
    bool            active          = false;
};

} // namespace coupled
