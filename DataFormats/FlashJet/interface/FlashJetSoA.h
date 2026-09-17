#ifndef DataFormats_FlashJet_interface_FlashJetSoA_h
#define DataFormats_FlashJet_interface_FlashJetSoA_h

#include <cstdint>

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace flashjet {

  // One row per input particle (in the order the producer selected them).
  //
  // inputs:   px, py, pz, e, candIdx (index into the source edm::View)
  // outputs:  jetIdx       jet index of the particle, in beam-merge order
  //           histP1/histP2/histChild/histD   merge history of step k
  //                        (pseudojet ids: particles are 0..n-1, new
  //                        pseudojets n, n+1, ...; histP2 = -1 for a beam
  //                        merge)
  //           jetPx..jetE  four-momentum of jet k, valid for k < nJets
  GENERATE_SOA_LAYOUT(FlashJetLayout,
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
                      SOA_COLUMN(double, jetE),
                      SOA_SCALAR(int32_t, nJets))

  using FlashJetSoA = FlashJetLayout<>;

}  // namespace flashjet

#endif  // DataFormats_FlashJet_interface_FlashJetSoA_h
