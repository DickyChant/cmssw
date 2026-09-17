#ifndef DataFormats_FlashJet_interface_FlashJetHostCollection_h
#define DataFormats_FlashJet_interface_FlashJetHostCollection_h

#include "DataFormats/FlashJet/interface/FlashJetSoA.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"

namespace flashjet {

  using FlashJetHostCollection = PortableHostCollection<FlashJetSoA>;

}  // namespace flashjet

#endif  // DataFormats_FlashJet_interface_FlashJetHostCollection_h
