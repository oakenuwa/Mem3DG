#pragma once

// =============================================================================
// steric_kernel_edge_face.h
//
// Abstract interface for the filament-edge to membrane-face steric interaction.
//
// Implement EdgeFaceKernel::compute() to define your repulsion potential.
// The coupling layer calls it for every (FilamentSegment, membrane face) pair
// whose minimum edge-face distance falls below the cutoff, then distributes
// the returned forces to the six involved points via DDG decomposition:
//
//   Filament side:
//     force_on_p0 = (1 − s) · result.force_on_segment
//     force_on_p1 =      s  · result.force_on_segment
//
//   Membrane side:
//     force_on_face_vertex_k = wFace[k] · result.force_on_face
//
// where s and wFace come from EdgeFaceContact (see edge_face_contact.h).
//
// NEWTON'S 3RD LAW
// ─────────────────
// For energy-conserving repulsion:
//   result.force_on_face = −result.force_on_segment
// You can relax this if the interaction is non-conservative (e.g., friction,
// active forces), but breaking it will violate momentum conservation.
// =============================================================================

#include "edge_face_contact.h"

#include <Eigen/Dense>

namespace coupled {

// Result returned by EdgeFaceKernel::compute().
struct EdgeFaceStericResult {
    // Force exerted on the filament segment at the contact point pSeg.
    // Decomposed to p0 and p1 by the coupling layer using s.
    Eigen::Vector3d force_on_segment = Eigen::Vector3d::Zero();

    // Force exerted on the membrane face at the contact point pFace.
    // Decomposed to the 3 face vertices by the coupling layer using wFace[].
    // For conservative repulsion: force_on_face = −force_on_segment.
    Eigen::Vector3d force_on_face    = Eigen::Vector3d::Zero();

    bool active = false;   // false → this pair contributes no force
};


// =============================================================================
// Abstract EdgeFaceKernel
// =============================================================================
//
// Subclass this and implement compute().  Example subclass skeleton:
//
//   class MyEdgeFaceKernel : public EdgeFaceKernel {
//   public:
//       double k;        // stiffness
//       double d0;       // contact radius (filament radius + membrane half-thickness)
//
//       EdgeFaceStericResult compute(
//           const EdgeFaceContact&  contact,
//           const Eigen::Vector3d&  face_normal,
//           double                  cutoff) const override
//       {
//           // ─── USER: fill in your steric potential here ─────────────────
//
//           double d = contact.distance;
//           if (d >= cutoff || d < 1e-10) return {};  // beyond cutoff → inactive
//
//           // 1. Unit normal pointing FROM membrane face TOWARD filament.
//           //    (contact.pSeg − contact.pFace) / d
//           Eigen::Vector3d nhat = (contact.pSeg - contact.pFace) / d;
//
//           // 2. Gap (negative when overlapping).
//           double gap = d - d0;
//
//           // 3. YOUR POTENTIAL DERIVATIVE HERE (example: harmonic).
//           //    dE/dd = k * gap   (<0 when gap<0 → repulsion pushes apart)
//           double dEdd = k * gap;
//
//           // 4. Force on the filament side at pSeg:
//           //    F_seg = −dE/dd · n̂  (positive in n̂ direction when dEdd < 0)
//           Eigen::Vector3d F = -dEdd * nhat;
//
//           EdgeFaceStericResult res;
//           res.active           = true;
//           res.force_on_segment = F;
//           res.force_on_face    = -F;   // Newton's 3rd
//           return res;
//
//           // ─────────────────────────────────────────────────────────────
//       }
//   };
// =============================================================================
class EdgeFaceKernel {
public:
    virtual ~EdgeFaceKernel() = default;

    // Called for every (filament segment, membrane face) pair whose
    // EdgeFaceContact.distance < cutoff.
    //
    // Parameters
    // ----------
    // contact     : geometry of the closest-point pair (see edge_face_contact.h)
    // face_normal : outward unit normal of the membrane face at pFace
    // cutoff      : distance threshold (same as CoupledParams::steric_cutoff)
    //
    // Return an EdgeFaceStericResult with active=false if no force should be
    // applied (e.g., the pair is actually beyond the physical interaction range
    // after a more precise check than the bounding-box pre-filter).
    virtual EdgeFaceStericResult compute(const EdgeFaceContact&  contact,
                                         const Eigen::Vector3d&  face_normal,
                                         double                   cutoff) const = 0;
};

} // namespace coupled
