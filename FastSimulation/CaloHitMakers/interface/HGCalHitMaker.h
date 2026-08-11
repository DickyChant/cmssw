#ifndef FastSimulation_CaloHitMakers_HGCalHitMaker_H
#define FastSimulation_CaloHitMakers_HGCalHitMaker_H

/** \class HGCalHitMaker
 *
 * Turns shower spots into HGCAL PCaloHits.
 *
 * Deliberately NOT derived from CaloHitMaker: that base assumes a crystal grid
 * (BaseCrystal has four-sided faces) and keys its map by hashed EB/EE indices.
 * HGCAL cells are hexagonal and their DetIds have no such dense hash, so hits
 * are keyed by rawId instead -- which is also what HGCDigitizer expects.
 *
 * Two FullSim behaviours are reproduced here because SimHit energies will not
 * match otherwise:
 *   - the HD120 thickness weight that HGCalSD applies (via HGCalFastGeometry),
 *   - aggregation of several spots landing on the same cell.
 *
 * Not yet reproduced (documented as a known systematic): the guard-ring,
 * mouse-bite and corner-mask rejections that HGCalSD applies.
 */

#include "DataFormats/DetId/interface/DetId.h"
#include "FastSimulation/CaloHitMakers/interface/HGCalShowerSpot.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include <map>
#include <vector>

class HGCalFastGeometry;

class HGCalHitMaker {
public:
  /// \param geom  fast lookup for the CE-E subdetector
  /// \param trackId  SimTrack id written into the hits
  HGCalHitMaker(const HGCalFastGeometry* geom, int trackId);

  /// Add one spot. Returns false if it did not land on a valid cell.
  bool addSpot(const HGCalShowerSpot& spot);

  /// Add a whole shower.
  unsigned addSpots(const std::vector<HGCalShowerSpot>& spots);

  /// Energy that fell outside any valid cell (leakage / dead area), GeV.
  double lostEnergy() const { return lost_; }

  /// Fill a PCaloHit container. Energies are the summed spot energies times the
  /// FullSim sim weight; time is the energy-weighted mean time of flight.
  void fillHits(std::vector<PCaloHit>& hits) const;

  std::size_t nCells() const { return cells_.size(); }

private:
  struct Accum {
    double energy = 0.;
    double eTime = 0.;  // sum of e*t, so the stored time is energy weighted
  };

  const HGCalFastGeometry* geom_;
  int trackId_;
  double lost_ = 0.;
  std::map<uint32_t, Accum> cells_;
};

#endif
