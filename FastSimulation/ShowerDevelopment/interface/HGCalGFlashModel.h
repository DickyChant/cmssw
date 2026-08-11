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
 *  4. Per-crossing energy. Shower crossings are NOT minimum-ionising: measured
 *     on these samples even single-crossing cells give mean/MPV = 1.74 against
 *     1.34 for a true Landau. So the deposit per spot is drawn from the measured
 *     per-crossing spectrum for the relevant silicon thickness, not from a
 *     muon-derived Landau. "MIP" remains only the reco unit.
 */

#include "FastSimulation/CaloHitMakers/interface/HGCalShowerSpot.h"

#include <vector>

class RandomEngineAndDistribution;
class HGCalProperties;

namespace edm {
  class ParameterSet;
}

class HGCalGFlashModel {
public:
  HGCalGFlashModel(const edm::ParameterSet& ps, const HGCalProperties* props);

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

  /// Cassette split: fraction of the pair energy on the Pb side at depth x0.
  double wSplit(double x0, double y) const;

  /// Radial scale at relative depth tau = t / T_max.
  double coreRadius(double tau) const;
  double tailRadius(double tau) const;
  double coreFraction(double tau) const;

private:
  const HGCalProperties* props_;

  // longitudinal
  double a0_, a1_, t0_, t1_;
  double sigmaLnAlpha_, sigmaLnT_, rhoLnAlphaT_;

  // cassette split w = wa (1 - exp(wb x0)), each linear in ln(y)
  double wa0_, wa1_, wb0_, wb1_;

  // transverse: r68(tau) = rc0 + rc1 * tau, plus the core/tail decomposition
  double r68Slope_, r68Const_;
  double coreOverR68_, tailOverCore_, coreFrac0_, coreFrac1_;

  // per-crossing deposit
  std::vector<double> crossingMPV_;    ///< keV, per Si thickness type
  std::vector<double> crossingWidth_;  ///< keV
  double meanCrossingE_;               ///< GeV, mean energy per crossing

  double spotEnergyGeV_;   ///< target energy per spot
  unsigned maxSpots_;
  bool applyW_;
  bool applyTransverse_;
};

#endif
