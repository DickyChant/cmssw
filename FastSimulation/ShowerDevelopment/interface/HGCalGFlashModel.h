// ---------------------------------------------------------------------------
// HGCAL FastSim (GFlash revival on the Phase-2 V16 geometry).
//
//  Author : Sitian Qian
//  Date   : 2026 (implementation on CMSSW_20_1_0_pre1,
//           GeometryExtendedRun4D110)
//
//  Shower parametrization after GFlash:
//    * E. Longo, I. Sestili, Nucl. Instrum. Meth. 128 (1975) 283
//      (Gamma-function longitudinal profile);
//    * G. Grindhammer, M. Rudowicz, S. Peters,
//      Nucl. Instrum. Meth. A290 (1990) 469 (GFlash);
//    * G. Grindhammer, S. Peters, hep-ex/0001020
//      (transverse parameterization in sampling calorimeters).
//  The calo-entry design follows Jan Eysermans' HGCAL FastSim
//  demonstrator (CMSSW_11_3_0_pre3, 2021).
// ---------------------------------------------------------------------------

#ifndef FastSimulation_ShowerDevelopment_HGCalGFlashModel_H
#define FastSimulation_ShowerDevelopment_HGCalGFlashModel_H

/** \class HGCalGFlashModel
 *
 * GFlash electromagnetic shower parametrization for the HGCAL CE-E, re-derived
 * on V16 / GeometryExtendedRun4D110 FullSim samples.
 *
 * The model produces energy "spots" -- (layer, x, y, energy, time) -- in the
 * global frame. It does not touch geometry or DetIds; HGCalHitMaker turns spots
 * into PCaloHits. Keeping the two apart is what lets a different shower model
 * (e.g. a Gaussian mixture) be swapped in behind the same interface.
 *
 * Structure, following Grindhammer-Peters as recalibrated:
 *
 *  1. Longitudinal. Gamma profile in depth t (X0 along the trajectory),
 *        alpha = a0*ln(y) + a1,  T = t0*ln(y) + t1,  beta = (alpha-1)/T,
 *     with (ln alpha, ln T) sampled from a correlated bivariate Gaussian for
 *     event-to-event fluctuations. Per-layer energy is the INTEGRAL of the
 *     Gamma over that layer's X0 interval, not the density at the layer centre:
 *     the CE-E layers differ in X0 (0.47 to 1.15) and the midpoint
 *     approximation biases thin and thick layers differently.
 *
 *  2. Cassette split (w). Within a cassette (Pb | Si | CuW | Si) the two silicon
 *     layers do not share the energy in proportion to the absorber in front of
 *     them. w(x0) = a (1 - exp(b x0)) is the fraction landing on the Pb side,
 *     with a and b depending on ln(y). Without this the per-layer profile shows
 *     a characteristic +-6% alternation against FullSim.
 *
 *  3. Transverse. Radial spots about the shower axis, sampled from the
 *     two-component form  f(r) = p f(r;R_C) + (1-p) f(r;R_T),
 *     f(r;R) = 2 r R^2 / (r^2+R^2)^2, whose CDF inverts analytically to
 *     r = R sqrt(u/(1-u)). Radii grow with depth; the measured r68(tau) anchors
 *     the scale. Sampling is done PERPENDICULAR to the shower axis and then
 *     projected onto the layer plane -- at eta = 2 the incidence is ~15 degrees
 *     and the in-plane footprint is an ellipse.
 *
 * The spots carry the energy deposited in the ABSORBER, which is the quantity
 * the Gamma profile describes. Converting that to the silicon deposit a PCaloHit
 * must carry is HGCalReverseCalibration's job, applied in HGCalHitMaker once the
 * cell thickness is known from the DetId.
 */

#include "FastSimulation/CaloHitMakers/interface/HGCalShowerSpot.h"

#include <vector>

class RandomEngineAndDistribution;
class HGCalProperties;
class HGCalFastGeometry;

namespace edm {
  class ParameterSet;
}

class HGCalGFlashModel {
public:
  /// \param geom   CE-E geometry, for layers 1..kNCEE
  /// \param cehGeom CE-H silicon geometry, for the electromagnetic tail that
  ///        punches through the back of CE-E. May be null, in which case the
  ///        shower is truncated at the end of CE-E as before.
  HGCalGFlashModel(const edm::ParameterSet& ps,
                   const HGCalProperties* props,
                   const HGCalFastGeometry* geom,
                   const HGCalFastGeometry* cehGeom = nullptr);

  /// Layers 1..kNCEE are CE-E; beyond that the spot layer is global and the
  /// caller routes it to CE-H, exactly as for the hadronic model.
  static constexpr int kNCEE = 26;

  /// Generate the spots for one incident particle.
  /// \param e0        incident energy, GeV
  /// \param entry     entry point on the CE-E front face (cm)
  /// \param dir       unit direction of the particle
  /// \param isGamma   photons start after a conversion depth
  void compute(double e0,
               const double entry[3],
               const double dir[3],
               bool isGamma,
               const RandomEngineAndDistribution* random,
               std::vector<HGCalShowerSpot>& spots) const;

  /// Mean longitudinal parameters for a given y = E/E_crit (no fluctuation).
  double meanAlpha(double y) const;
  double meanT(double y) const;

  /// Event-to-event widths of (ln T, ln alpha) and their correlation, in the
  /// Grindhammer-Peters sampling-calorimeter forms 1/(s1 + s2 ln y) and
  /// r1 + r2 ln y. A sampling calorimeter fluctuates more than the homogeneous
  /// medium at low energy and the excess dies away as 1/ln y.
  double sigmaLnT(double y) const;
  double sigmaLnAlpha(double y) const;
  double rhoLnAlphaT(double y) const;

  /// Cassette split: fraction of the pair energy on the Pb side at depth x0.
  double wSplit(double x0, double y) const;

  /// Radial scale at relative depth tau = t / T_max.
  double coreRadius(double tau) const;
  double tailRadius(double tau) const;
  double coreFraction(double tau) const;

private:
  const HGCalProperties* props_;
  const HGCalFastGeometry* geom_;     ///< CE-E, for the true layer z positions
  const HGCalFastGeometry* cehGeom_;  ///< CE-H silicon, for the punch-through tail

  // longitudinal
  double a0_, a1_, t0_, t1_;
  // sigma(ln T) = 1/(sT1_ + sT2_ ln y), sigma(ln alpha) = 1/(sA1_ + sA2_ ln y),
  // rho = clamp(r1_ + r2_ ln y, 0, 0.97)
  double sT1_, sT2_, sA1_, sA2_, r1_, r2_;

  // CE-H tail: anchored attenuation instead of extrapolating the Gamma.
  // dep_k = g0 * (last CE-E layer density) * exp(-(d_k - d_exit)/lambda) * dX0_k
  double cehLambda_, cehG0_;
  bool applyCehTail_;

  // cassette split w = wa (1 - exp(wb x0)), each linear in ln(y)
  double wa0_, wa1_, wb0_, wb1_;

  // transverse: r68(tau) = rc0 + rc1 * tau, plus the core/tail decomposition
  double r68Slope_, r68Const_;
  double coreOverR68_, tailOverCore_, coreFrac0_, coreFrac1_;

  double spotEnergyGeV_;   ///< target energy per spot
  unsigned maxSpots_;
  bool applyW_;
  bool applyTransverse_;
  bool applyConversion_;  ///< off when T was fitted on photons (conversion already in it)
};

#endif
