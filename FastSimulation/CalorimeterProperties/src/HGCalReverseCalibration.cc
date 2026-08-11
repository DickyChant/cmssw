#include "FastSimulation/CalorimeterProperties/interface/HGCalReverseCalibration.h"
#include "FastSimulation/Utilities/interface/LandauFluctuationGenerator.h"
#include "FastSimulation/Utilities/interface/RandomEngineAndDistribution.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <algorithm>
#include <cmath>

HGCalReverseCalibration::HGCalReverseCalibration(const edm::ParameterSet& ps) {
  dEdXWeights_ = ps.getParameter<std::vector<double> >("dEdXWeights");
  fCPerMIP_ = ps.getParameter<std::vector<double> >("fCPerMIP");
  thicknessCorrection_ = ps.getParameter<std::vector<double> >("thicknessCorrection");
  cce_ = ps.getParameter<std::vector<double> >("chargeCollectionEfficiency");
  keV2fC_ = ps.getParameter<double>("keV2fC");

  crossingMPVkeV_ = ps.getParameter<std::vector<double> >("crossingMPVkeV");
  crossingWidthkeV_ = ps.getParameter<std::vector<double> >("crossingWidthkeV");
  crossingMeankeV_ = ps.getParameter<std::vector<double> >("crossingMeankeV");
  fluctuate_ = ps.getParameter<bool>("fluctuate");

  if (dEdXWeights_.empty() || fCPerMIP_.empty())
    throw cms::Exception("HGCalReverseCalibration") << "dEdXWeights and fCPerMIP must not be empty";
  if (keV2fC_ <= 0.)
    throw cms::Exception("HGCalReverseCalibration") << "keV2fC must be positive";

  landau_ = std::make_unique<LandauFluctuationGenerator>();
}

HGCalReverseCalibration::~HGCalReverseCalibration() = default;

double HGCalReverseCalibration::weight(int layer) const {
  if (layer < 0 || layer >= static_cast<int>(dEdXWeights_.size()))
    return 0.;
  return dEdXWeights_[layer];
}

double HGCalReverseCalibration::thicknessCorrection(int thickness) const {
  if (thickness < 0 || thickness >= static_cast<int>(thicknessCorrection_.size()))
    return 1.;
  return thicknessCorrection_[thickness];
}

double HGCalReverseCalibration::siEnergyGeV(double showerEnergyGeV, int layer, int thickness) const {
  const double w = weight(layer);
  if (showerEnergyGeV <= 0. || w <= 0.)
    return 0.;

  const double cce = (thickness >= 0 && thickness < static_cast<int>(cce_.size())) ? cce_[thickness] : 1.;
  const double fcmip =
      (thickness >= 0 && thickness < static_cast<int>(fCPerMIP_.size())) ? fCPerMIP_[thickness] : fCPerMIP_.back();

  // GeV(absorber) -> MIP -> fC -> keV(silicon) -> GeV(silicon)
  const double mip = showerEnergyGeV * 1.0e3 * thicknessCorrection(thickness) * cce / w;
  const double keV = mip * fcmip / keV2fC_;
  return keV * 1.0e-6;
}

double HGCalReverseCalibration::meanCrossingGeV(int thickness) const {
  if (thickness < 0 || thickness >= static_cast<int>(crossingMeankeV_.size()))
    return 0.;
  return crossingMeankeV_[thickness] * 1.0e-6;
}

double HGCalReverseCalibration::sampleSiEnergyGeV(double showerEnergyGeV,
                                                 int layer,
                                                 int thickness,
                                                 const RandomEngineAndDistribution* random) const {
  const double mean = siEnergyGeV(showerEnergyGeV, layer, thickness);
  if (mean <= 0. || !fluctuate_ || random == nullptr)
    return mean;

  const double meanCross = meanCrossingGeV(thickness);
  if (meanCross <= 0.)
    return mean;

  // Number of crossings that make up this deposit. Poisson is a stand-in: the
  // measured multiplicity is far from Poisson (LD200: <m> = 5.6 with P(1) = 0.44,
  // against 0.02 for Poisson), so this is the known-approximate part and the
  // measured distribution should replace it.
  const double nExpected = mean / meanCross;
  unsigned n = (nExpected < 50.) ? random->poissonShoot(nExpected)
                                 : static_cast<unsigned>(std::lround(nExpected));
  if (n == 0)
    return 0.;

  const int t = std::clamp(thickness, 0, static_cast<int>(crossingMPVkeV_.size()) - 1);
  const double mpv = crossingMPVkeV_[t];
  const double width = crossingWidthkeV_[t];

  // Above a few tens of crossings the sum is Gaussian by the central limit
  // theorem; sampling each one individually would only cost time.
  double sumKeV = 0.;
  if (n <= 50) {
    for (unsigned i = 0; i < n; ++i)
      sumKeV += std::max(mpv + width * landau_->landau(random), 0.);
  } else {
    const double m = n * crossingMeankeV_[t];
    const double s = std::sqrt(static_cast<double>(n)) * width * 2.5;
    sumKeV = std::max(random->gaussShoot(m, s), 0.);
  }

  // Preserve the calibrated mean: the crossing spectrum sets the shape of the
  // fluctuation, while the reverse calibration sets the scale.
  const double expectedKeV = nExpected * crossingMeankeV_[t];
  if (expectedKeV > 0.)
    sumKeV *= (mean * 1.0e6) / expectedKeV;

  return sumKeV * 1.0e-6;
}
