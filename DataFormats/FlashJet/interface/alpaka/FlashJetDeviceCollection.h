#ifndef DataFormats_FlashJet_interface_alpaka_FlashJetDeviceCollection_h
#define DataFormats_FlashJet_interface_alpaka_FlashJetDeviceCollection_h

#include "DataFormats/FlashJet/interface/FlashJetHostCollection.h"
#include "DataFormats/FlashJet/interface/FlashJetSoA.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace flashjet {

    using namespace ::flashjet;

    using FlashJetDeviceCollection = PortableCollection<FlashJetSoA>;

  }  // namespace flashjet

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(flashjet::FlashJetDeviceCollection, flashjet::FlashJetHostCollection);

#endif  // DataFormats_FlashJet_interface_alpaka_FlashJetDeviceCollection_h
