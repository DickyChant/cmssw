#include "FastSimulation/CalorimeterProperties/interface/HGCalProperties.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <cmath>
#include <numeric>

HGCalProperties::HGCalProperties(const edm::ParameterSet& fastDet) {
  const edm::ParameterSet ps = fastDet.getParameter<edm::ParameterSet>("HGCalCalorimeterProperties");

  lightColl_ = ps.getParameter<double>("lightColl");
  lightCollUnif_ = ps.getParameter<double>("lightCollUnif");
  photoStatistics_ = ps.getParameter<double>("photoStatistics");
  interactionLength_ = ps.getParameter<double>("interactionLength");

  Aeff_ = ps.getParameter<double>("Aeff");
  Zeff_ = ps.getParameter<double>("Zeff");
  rho_ = ps.getParameter<double>("rho");
  radLenIngcm2_ = ps.getParameter<double>("radLenIngcm2");

  radLenIncm_ = ps.getParameter<double>("radLenIncm");
  radLenIncm_ = (radLenIncm_ < 0) ? radLenIngcm2_ / rho_ : radLenIncm_;

  criticalEnergy_ = ps.getParameter<double>("criticalEnergy");
  criticalEnergy_ =
      (criticalEnergy_ < 0) ? 2.66E-3 * std::pow(radLenIngcm2_ * Zeff_ / Aeff_, 1.1) : criticalEnergy_;

  moliereRadius_ = ps.getParameter<double>("moliereRadius");
  moliereRadius_ = (moliereRadius_ < 0) ? scaleEnergy_ / criticalEnergy_ * radLenIncm_ : moliereRadius_;

  Fs_ = ps.getParameter<double>("Fs");
  ehat_ = ps.getParameter<double>("ehat");
  resE_ = ps.getParameter<double>("resE");
  da_ = ps.getParameter<double>("da");
  dp_ = ps.getParameter<double>("dp");

  // The CE-E is a sampling calorimeter: this selects the sampling branches of
  // EMECALShowerParametrization / EMShower rather than the homogeneous ones.
  bHom_ = false;

  etaMin_ = ps.getParameter<double>("etaMin");
  etaMax_ = ps.getParameter<double>("etaMax");

  layerX0_ = ps.getParameter<std::vector<double> >("layerX0");
  if (layerX0_.empty()) {
    throw cms::Exception("HGCalProperties") << "layerX0 must not be empty: the CE-E layers are not uniform in X0 "
                                               "and cannot be approximated by a single value.";
  }

  frontX0_.assign(layerX0_.size(), 0.);
  double cum = 0.;
  for (std::size_t i = 0; i < layerX0_.size(); ++i) {
    frontX0_[i] = cum;
    cum += layerX0_[i];
  }
  totalX0_ = cum;

  // thickness_ (base class) is the CE-E depth in cm at normal incidence
  thickness_ = totalX0_ * radLenIncm_;

  edm::LogInfo("HGCalProperties") << "CE-E: " << layerX0_.size() << " layers, " << totalX0_ << " X0, "
                                  << thickness_ << " cm, X0 = " << radLenIncm_ << " cm, R_M = " << moliereRadius_
                                  << " cm, E_crit = " << criticalEnergy_ * 1000. << " MeV";
}

double HGCalProperties::thickness(double eta) const {
  const double aeta = std::abs(eta);
  if (aeta < etaMin_ || aeta > etaMax_)
    return 0.;
  // The CE-E layers are planes perpendicular to z, so the path length through
  // the stack is thickness / |cos(theta)|, and cos(theta) = tanh(eta).
  const double cth = std::tanh(aeta);
  if (cth < 1e-6)
    return thickness_;
  return thickness_ / cth;
}

double HGCalProperties::layerX0(unsigned layer) const {
  if (layer < 1 || layer > layerX0_.size())
    return 0.;
  return layerX0_[layer - 1];
}

double HGCalProperties::frontX0(unsigned layer) const {
  if (layer < 1 || layer > frontX0_.size())
    return (layer > frontX0_.size()) ? totalX0_ : 0.;
  return frontX0_[layer - 1];
}
