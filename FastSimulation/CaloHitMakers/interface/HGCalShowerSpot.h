#ifndef FastSimulation_CaloHitMakers_HGCalShowerSpot_H
#define FastSimulation_CaloHitMakers_HGCalShowerSpot_H

/** One energy deposit produced by an HGCAL shower model.
 *
 * Lives in CaloHitMakers rather than ShowerDevelopment because
 * ShowerDevelopment already depends on CaloHitMakers; putting it the other way
 * round would make the dependency circular. It is the only thing a shower model
 * and the hit maker need to agree on, which is what allows a different model
 * (e.g. a Gaussian mixture) to be swapped in behind the same interface.
 */
struct HGCalShowerSpot {
  int layer;       ///< CE-E layer, 1-based
  double x, y, z;  ///< global position, cm
  double energy;   ///< deposited energy, GeV
  double t;        ///< time of flight, ns
};

#endif
