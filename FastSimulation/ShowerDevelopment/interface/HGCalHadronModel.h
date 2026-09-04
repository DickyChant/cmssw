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
 *    of the visible energy at 5 GeV rising to 49% at 500 GeV;
 *  * they are far less visible. The measured silicon sampling fraction gives an
 *    e/h ratio of 0.62/0.56/0.55 at 5/50/500 GeV, energy dependent, applied as
 *    an explicit factor;
 *  * they fluctuate enormously: where the shower starts (the first interaction)
 *    and how much of it is electromagnetic vary event to event far beyond what
 *    a single Gamma with lognormal (alpha, T) can carry.
 *
 * Following the structure of the 1989 three-Gamma parametrization
 * (Grindhammer-Rudowicz-Peters, SLAC-PUB-5072: several profile components with
 * FLUCTUATING fractions and a fully correlated parameter draw), the profile is
 * two compartment Gammas tied together by a correlated 5-dim Gaussian:
 *
 *   f(depth) = (1-fH) Gamma(aE, TE) [CE-E, X0 axis]  +  fH Gamma(aH, TH) [CE-H]
 *
 *   theta = (ln aE, ln TE, ln aH, ln TH, logit fH)
 *   theta = mu(ln y) + D(sigma(ln y)) . L . z,   z ~ N(0, 1)^5
 *
 * with mu and sigma linear in ln y and L the (constant, measured) Cholesky
 * factor of the correlation matrix. The wide logit-normal fH (sigma ~ 2.6!)
 * naturally produces the punch-through class of showers that start late and
 * put ~everything into CE-H -- the analog of the 1989 shower classes.
 * Population measured on 68k uniform (flat ln E, 1.5-1000 GeV) FullSim pions.
 *
 * Depth axes are those of the parameter derivation: the CE-E per-layer X0
 * table for the first compartment, and a uniform cehDepthPerLayer (3.59 X0
 * equivalent) per CE-H layer for the second. T here is the MEAN depth of the
 * compartment profile (moment convention, beta = alpha/T), unlike the
 * electromagnetic model where T is the peak.
 *
 * The per-crossing response is NOT re-measured for hadrons: the LD200 MPV agrees
 * with photons to 2% (59.9 vs 61.1 keV), as expected for a sensor property. Only
 * the visible fraction differs, and that is applied as the explicit e/h factor.
 */

#include "FastSimulation/CaloHitMakers/interface/HGCalShowerSpot.h"

#include <vector>

class RandomEngineAndDistribution;
class HGCalFastGeometry;
class HGCalProperties;

namespace edm {
  class ParameterSet;
}

class HGCalHadronModel {
public:
  HGCalHadronModel(const edm::ParameterSet& ps,
                   const HGCalFastGeometry* cee,
                   const HGCalFastGeometry* ceh,
                   const HGCalProperties* props);

  /// Generate spots across CE-E and CE-H. `layer` in the returned spots is the
  /// GLOBAL layer (1..26 CE-E, 27..47 CE-H); the caller routes them to the right
  /// hit maker.
  void compute(double e0,
               const double entry[3],
               const double dir[3],
               const RandomEngineAndDistribution* random,
               std::vector<HGCalShowerSpot>& spots) const;

  /// Visible-energy suppression relative to an electromagnetic shower.
  double ehRatio(double y) const;

  static constexpr int kNCEE = 26;
  static constexpr int kNTotal = 47;
  static constexpr int kNPar = 5;  ///< (ln aE, ln TE, ln aH, ln TH, logit fH)

private:
  const HGCalFastGeometry* cee_;
  const HGCalFastGeometry* ceh_;
  const HGCalProperties* props_;  ///< CE-E X0 table for the first compartment

  // theta_i = (muSlope_i ln y + muConst_i) + sigma_i(ln y) * (L z)_i,
  // sigma_i = sigmaSlope_i ln y + sigmaConst_i
  double muSlope_[kNPar], muConst_[kNPar];
  double sigmaSlope_[kNPar], sigmaConst_[kNPar];
  double chol_[kNPar * (kNPar + 1) / 2];  ///< row-major lower triangle

  double cehDepthPerLayer_;  ///< CE-H compartment depth axis, X0-equivalent per layer

  double eh0_, eh1_;
  double ecrit_;
  double r68Slope_, r68Const_, coreOverR68_, tailOverCore_, coreFraction_;
  double spotEnergyGeV_;
  unsigned maxSpots_;
};

#endif
