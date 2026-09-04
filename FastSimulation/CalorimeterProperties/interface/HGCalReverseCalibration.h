#ifndef FastSimulation_CalorimeterProperties_HGCalReverseCalibration_H
#define FastSimulation_CalorimeterProperties_HGCalReverseCalibration_H

/** \class HGCalReverseCalibration
 *
 * Converts the shower energy predicted by the parametrization into the energy
 * actually deposited in silicon, i.e. what a PCaloHit must carry.
 *
 * This is the inverse of what reconstruction does. HGCalRecHit computes
 *
 *    E_GeV = amplitude_MIP * dEdX_weight[layer] * 1e-3 / thicknessCorrection / cce
 *
 * so going the other way,
 *
 *    MIP    = E_GeV * 1e3 * thicknessCorrection * cce / dEdX_weight[layer]
 *    charge = MIP * fCPerMIP[thickness]                       [fC]
 *    E_Si   = charge / keV2fC                                 [keV]
 *
 * Getting this wrong is not subtle: writing the incident shower energy straight
 * into the hits overstates the silicon deposit by more than two orders of
 * magnitude (a 50 GeV shower deposits ~0.25 GeV in silicon, not 50 GeV).
 *
 * The weights used here must be the same ones reconstruction applies -- the
 * calcWeights-averaged table -- otherwise FastSim and reco disagree by
 * construction. That is the opposite of the w-parametrization, which needs the
 * raw table; the two uses are genuinely different.
 *
 * Fluctuations. The mean is only half the story: the deposit in a cell is a sum
 * over crossings, each drawn from the measured per-crossing spectrum. Shower
 * crossings are NOT minimum-ionising (measured mean/MPV = 1.74 against 1.34 for
 * a true Landau), so the spectrum comes from the shower samples themselves, not
 * from muon calibration.
 */

#include <memory>
#include <mutex>
#include <vector>

namespace edm {
  class ParameterSet;
}
class RandomEngineAndDistribution;
class LandauFluctuationGenerator;

class HGCalReverseCalibration {
public:
  explicit HGCalReverseCalibration(const edm::ParameterSet& ps);
  ~HGCalReverseCalibration();

  /// Mean silicon deposit, GeV, for a shower energy deposited in the absorber of
  /// a given CE-E layer (1-based) read out by a cell of the given thickness type
  /// (0 = HD120, 1 = LD200, 2 = LD300, 3 = HD200).
  double siEnergyGeV(double showerEnergyGeV, int layer, int thickness) const;

  /// Silicon deposit including crossing-by-crossing fluctuation: the mean is
  /// converted to a number of crossings, each sampled from the measured
  /// per-crossing spectrum. Falls back to the mean if no RNG is supplied.
  double sampleSiEnergyGeV(double showerEnergyGeV,
                           int layer,
                           int thickness,
                           const RandomEngineAndDistribution* random) const;

  /// Mean energy of one crossing in this thickness, GeV.
  double meanCrossingGeV(int thickness) const;

  bool fluctuate() const { return fluctuate_; }

  /// Reconstruction dE/dx weight for a (1-based) layer, MeV per MIP.
  /// Public so alternative shower models can close the weighted-energy
  /// sum the same way reconstruction does.
  double weight(int layer) const;

private:
  double thicknessCorrection(int thickness) const;

  std::vector<double> dEdXWeights_;           ///< MeV per MIP, indexed by layer
  std::vector<double> fCPerMIP_;              ///< fC per MIP, by thickness type
  std::vector<double> thicknessCorrection_;   ///< by thickness type (CE-E entries)
  std::vector<double> cce_;                   ///< charge collection efficiency
  double keV2fC_ = 0.044259;

  std::vector<double> crossingMPVkeV_;
  std::vector<double> crossingWidthkeV_;
  std::vector<double> crossingMeankeV_;
  bool fluctuate_ = true;
  /// Cap on a single crossing, in units of the MPV. The Landau tail is
  /// unbounded; a real crossing in thin silicon is not.
  double maxCrossingOverMPV_ = 30.;

  std::unique_ptr<LandauFluctuationGenerator> landau_;

  /// Mean of the *sampled* per-crossing shape, per thickness. This is NOT the
  /// measured crossingMean: sampling mpv + width*Landau (clamped) has its own
  /// mean, and normalising by the measured one instead biases every deposit low
  /// by a constant factor. Estimated once, on first use, from the same generator
  /// and clamp that the sampling uses.
  mutable std::vector<double> drawMeankeV_;
  mutable std::once_flag drawMeanOnce_;

  void ensureDrawMean(const RandomEngineAndDistribution* random) const;
};

#endif
