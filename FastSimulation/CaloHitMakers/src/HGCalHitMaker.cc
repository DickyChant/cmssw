#include "FastSimulation/CaloHitMakers/interface/HGCalHitMaker.h"
#include "FastSimulation/CaloGeometryTools/interface/HGCalFastGeometry.h"
#include "FastSimulation/CalorimeterProperties/interface/HGCalReverseCalibration.h"

#include <cmath>

HGCalHitMaker::HGCalHitMaker(const HGCalFastGeometry* geom,
                             const HGCalReverseCalibration* calib,
                             const RandomEngineAndDistribution* random,
                             int trackId)
    : geom_(geom), calib_(calib), random_(random), trackId_(trackId) {}

bool HGCalHitMaker::addSpot(const HGCalShowerSpot& spot) {
  // NaN fails every comparison, so test positively rather than for <= 0.
  if (geom_ == nullptr || !(spot.energy > 0.) || !std::isfinite(spot.x) || !std::isfinite(spot.y) ||
      !std::isfinite(spot.z) || !std::isfinite(spot.t))
    return false;

  const int zside = (spot.z >= 0.) ? 1 : -1;
  const DetId id = geom_->toDetId(spot.layer, spot.x, spot.y, zside);
  if (id == DetId()) {
    lost_ += spot.energy;
    return false;
  }

  // The spot carries the energy the shower deposited in the ABSORBER. What a
  // PCaloHit must carry is the energy deposited in SILICON, which is smaller by
  // roughly the sampling fraction. Converting here (rather than in the shower
  // model) is deliberate: the conversion depends on the cell thickness, which is
  // only known once the DetId is assigned.
  double e = spot.energy;
  if (calib_ != nullptr) {
    const int thickness = HGCalFastGeometry::siThicknessType(id);
    e = calib_->sampleSiEnergyGeV(spot.energy, spot.layer, thickness, random_);
  }
  if (!(e > 0.))
    return false;

  // FullSim (HGCalSD) scales HD120 cells by cellThickness/waferThick; without
  // the same weight the SimHit energies cannot match.
  e *= geom_->simWeight(id);

  Accum& a = cells_[id.rawId()];
  a.energy += e;
  a.eTime += e * spot.t;
  return true;
}

unsigned HGCalHitMaker::addSpots(const std::vector<HGCalShowerSpot>& spots) {
  unsigned n = 0;
  for (const auto& s : spots)
    if (addSpot(s))
      ++n;
  return n;
}

double HGCalHitMaker::depositedEnergy() const {
  double sum = 0.;
  for (const auto& kv : cells_)
    sum += kv.second.energy;
  return sum;
}

void HGCalHitMaker::fillHits(std::vector<PCaloHit>& hits) const {
  hits.reserve(hits.size() + cells_.size());
  for (const auto& kv : cells_) {
    const double e = kv.second.energy;
    if (!(e > 0.))
      continue;
    const double t = kv.second.eTime / e;
    hits.emplace_back(kv.first, static_cast<float>(e), static_cast<float>(t), trackId_);
  }
}
