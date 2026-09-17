#ifndef DataFormats_FlashJet_interface_FlashJetSoA_h
#define DataFormats_FlashJet_interface_FlashJetSoA_h

#include <cstdint>

#include "DataFormats/SoATemplate/interface/SoABlocks.h"
#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace flashjet {

  // One row per input particle.  The particles of all entries (independent
  // clustering problems: one event, or the constituents of one jet) are
  // stored back to back; entry b owns rows [offset[b], offset[b] + size[b]).
  //
  // inputs:   px, py, pz, e, candIdx (index of the particle in its source:
  //           the candidate View, or the daughter index within a jet)
  // outputs:  jetIdx       jet index of the particle within its entry, in
  //                        beam-merge order
  //           histP1/histP2/histChild/histD   merge history of step k
  //                        (pseudojet ids local to the entry: particles
  //                        0..n-1, new pseudojets n, n+1, ...; histP2 = -1
  //                        for a beam merge)
  //           jetPx..jetE  four-momentum of jet k, valid for k < nJets
  GENERATE_SOA_LAYOUT(FlashJetParticleLayout,
                      SOA_COLUMN(double, px),
                      SOA_COLUMN(double, py),
                      SOA_COLUMN(double, pz),
                      SOA_COLUMN(double, e),
                      SOA_COLUMN(int32_t, candIdx),
                      SOA_COLUMN(int32_t, jetIdx),
                      SOA_COLUMN(int32_t, histP1),
                      SOA_COLUMN(int32_t, histP2),
                      SOA_COLUMN(int32_t, histChild),
                      SOA_COLUMN(double, histD),
                      SOA_COLUMN(double, jetPx),
                      SOA_COLUMN(double, jetPy),
                      SOA_COLUMN(double, jetPz),
                      SOA_COLUMN(double, jetE))

  // One row per entry.
  //
  // inputs:   offset, size, source (index of the jet that was reclustered,
  //           -1 for whole-event clustering)
  // outputs:  nJets
  //           groomedPx..groomedE, zg, rg, nDropped   soft drop of the
  //                        hardest jet of the entry (zeros if disabled)
  GENERATE_SOA_LAYOUT(FlashJetEntryLayout,
                      SOA_COLUMN(int32_t, offset),
                      SOA_COLUMN(int32_t, size),
                      SOA_COLUMN(int32_t, source),
                      SOA_COLUMN(int32_t, nJets),
                      SOA_COLUMN(double, groomedPx),
                      SOA_COLUMN(double, groomedPy),
                      SOA_COLUMN(double, groomedPz),
                      SOA_COLUMN(double, groomedE),
                      SOA_COLUMN(double, zg),
                      SOA_COLUMN(double, rg),
                      SOA_COLUMN(int32_t, nDropped))

  GENERATE_SOA_BLOCKS(FlashJetBlocks,
                      SOA_BLOCK(particles, FlashJetParticleLayout),
                      SOA_BLOCK(entries, FlashJetEntryLayout))

  using FlashJetParticleSoA = FlashJetParticleLayout<>;
  using FlashJetEntrySoA = FlashJetEntryLayout<>;
  using FlashJetSoA = FlashJetBlocks<>;

}  // namespace flashjet

#endif  // DataFormats_FlashJet_interface_FlashJetSoA_h
