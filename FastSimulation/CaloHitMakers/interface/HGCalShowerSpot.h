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
