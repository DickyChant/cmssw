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

#ifndef FastSimulation_CalorimeterProperties_HGCalProperties_H
#define FastSimulation_CalorimeterProperties_HGCalProperties_H

/** \class HGCalProperties
 *
 * Homogenised material of the HGCAL CE-E, plus the per-layer radiation-length
 * table, for the FastSim GFlash parametrization.
 *
 * Derives from ECALProperties so that EMShower / EMECALShowerParametrization
 * remain usable, but with bHom_ = false: the CE-E is a sampling calorimeter, so
 * the sampling branches of those classes are the correct ones.
 *
 * The per-layer X0 table is NOT uniform and must not be approximated as such.
 * It is derived from the raw reco dE/dx weights (which are proportional to the
 * absorber in front of each layer) and shows the physical structure:
 *   - layer 1 is ~half the others (the half-lead first cassette),
 *   - layers then alternate 1.091 / 0.797 X0 (Pb side vs CuW side),
 *   - stepping to 1.091 / 1.149 from layer 19 where the lead thickens.
 * Using the calcWeights-averaged table instead erases that alternation, which
 * is exactly what the w-parametrization describes.
 *
 * Defaults correspond to V16 / GeometryExtendedRun4D110. Everything is settable
 * from python so the numbers can be re-derived for another geometry.
 */

#include "FastSimulation/CalorimeterProperties/interface/ECALProperties.h"

#include <vector>

namespace edm {
  class ParameterSet;
}

class HGCalProperties : public ECALProperties {
public:
  HGCalProperties(const edm::ParameterSet& fastDet);

  ~HGCalProperties() override = default;

  /// Thickness of the CE-E along the trajectory, in cm, at this eta.
  /// Zero outside the HGCAL acceptance so callers can test containment.
  double thickness(double eta) const override;

  /// Not meaningful for a silicon sampling calorimeter; kept because
  /// ECALProperties declares them pure virtual.
  double photoStatistics() const override { return photoStatistics_; }
  double lightCollectionEfficiency() const override { return lightColl_; }
  double lightCollectionUniformity() const override { return lightCollUnif_; }

  /// Number of CE-E layers.
  unsigned nLayers() const { return layerX0_.size(); }

  /// Radiation length of one CE-E layer (1-based layer number), in X0.
  double layerX0(unsigned layer) const;

  /// Cumulative X0 at the front face of a layer (1-based).
  double frontX0(unsigned layer) const;

  /// Total CE-E depth in X0 at normal incidence.
  double totalX0() const { return totalX0_; }

  /// Critical energy in GeV (y = E / E_crit).
  double criticalEnergyGeV() const { return criticalEnergy_; }

  /// eta acceptance
  double etaMin() const { return etaMin_; }
  double etaMax() const { return etaMax_; }

private:
  std::vector<double> layerX0_;
  std::vector<double> frontX0_;
  double totalX0_ = 0.;
  double etaMin_ = 1.5;
  double etaMax_ = 3.0;
};

#endif
