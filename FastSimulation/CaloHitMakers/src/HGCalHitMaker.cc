#include "FastSimulation/CaloHitMakers/interface/HGCalHitMaker.h"
#include "FastSimulation/CaloGeometryTools/interface/HGCalFastGeometry.h"

#include <cmath>

HGCalHitMaker::HGCalHitMaker(const HGCalFastGeometry* geom, int trackId) : geom_(geom), trackId_(trackId) {}

bool HGCalHitMaker::addSpot(const HGCalShowerSpot& spot) {
  if (geom_ == nullptr || spot.energy <= 0.)
    return false;

  const int zside = (spot.z >= 0.) ? 1 : -1;
  const DetId id = geom_->toDetId(spot.layer, spot.x, spot.y, zside);
  if (id == DetId()) {
    lost_ += spot.energy;
    return false;
  }

  // FullSim (HGCalSD) scales HD120 cells by cellThickness/waferThick; without
  // the same weight the SimHit energies cannot match.
  const double e = spot.energy * geom_->simWeight(id);

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

void HGCalHitMaker::fillHits(std::vector<PCaloHit>& hits) const {
  hits.reserve(hits.size() + cells_.size());
  for (const auto& kv : cells_) {
    const double e = kv.second.energy;
    if (e <= 0.)
      continue;
    const double t = kv.second.eTime / e;
    hits.emplace_back(kv.first, static_cast<float>(e), static_cast<float>(t), trackId_);
  }
}
