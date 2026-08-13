#ifndef FastSimulation_ShowerDevelopment_HGCalHadronModel_H
#define FastSimulation_ShowerDevelopment_HGCalHadronModel_H

/** \class HGCalHadronModel
 *
 * Hadronic shower parametrization for HGCAL, spanning CE-E and CE-H, measured on
 * the V16/D110 single-pion FullSim samples.
 *
 * Hadronic showers are not electromagnetic ones with different numbers:
 *
 *  * they are far deeper. The measured shower maximum sits at layer 25, i.e. at
 *    the CE-E/CE-H boundary, against layer 9 for photons, and CE-H carries 20%
 *    of the visible energy at 5 GeV rising to 49% at 500 GeV. A CE-E-only model
 *    would therefore miss up to half the shower;
 *  * they are far less visible. The measured silicon sampling fraction is
 *    0.32/0.28/0.27% at 5/50/500 GeV against 0.52/0.51/0.50% for photons, i.e.
 *    a ratio of 0.62/0.56/0.55 -- the usual e/h suppression, and it is energy
 *    dependent, so it cannot be folded into a single constant;
 *  * they are much wider, spilling into the neighbouring silicon thicknesses
 *    (ten times the LD300 hits of a photon shower at the same energy).
 *
 * Rather than assume a Grindhammer-Peters hadronic form, the longitudinal shape
 * is taken directly from the data: a Gamma in LAYER units over the full 47-layer
 * stack, which the pion samples are well described by,
 *
 *     alpha_had = 0.535 * ln(y) - 1.125
 *     T_had     = 2.673 * ln(y) - 5.567     [layers]
 *
 * Layer units are used deliberately in place of interaction lengths: the CE-E
 * and CE-H stacks differ in material, so a single lambda scale would not span
 * both, whereas the layer index is exact and the geometry supplies the position.
 *
 * The per-crossing response is NOT re-measured for hadrons: the LD200 MPV agrees
 * with photons to 2% (59.9 vs 61.1 keV), as expected for a sensor property. Only
 * the visible fraction differs, and that is applied as an explicit e/h factor.
 */

#include "FastSimulation/CaloHitMakers/interface/HGCalShowerSpot.h"

#include <vector>

class RandomEngineAndDistribution;
class HGCalFastGeometry;

namespace edm {
  class ParameterSet;
}

class HGCalHadronModel {
public:
  HGCalHadronModel(const edm::ParameterSet& ps,
                   const HGCalFastGeometry* cee,
                   const HGCalFastGeometry* ceh);

  /// Generate spots across CE-E and CE-H. `layer` in the returned spots is the
  /// GLOBAL layer (1..26 CE-E, 27..47 CE-H); the caller routes them to the right
  /// hit maker.
  void compute(double e0,
               const double entry[3],
               const double dir[3],
               const RandomEngineAndDistribution* random,
               std::vector<HGCalShowerSpot>& spots) const;

  double meanAlpha(double y) const;
  double meanT(double y) const;
  /// Visible-energy suppression relative to an electromagnetic shower.
  double ehRatio(double y) const;

  static constexpr int kNCEE = 26;
  static constexpr int kNTotal = 47;

private:
  const HGCalFastGeometry* cee_;
  const HGCalFastGeometry* ceh_;

  double a0_, a1_, t0_, t1_;
  double sigmaLnAlpha_, sigmaLnT_;
  double eh0_, eh1_;
  double ecrit_;
  double r68Slope_, r68Const_, coreOverR68_, tailOverCore_, coreFraction_;
  double spotEnergyGeV_;
  unsigned maxSpots_;
};

#endif
