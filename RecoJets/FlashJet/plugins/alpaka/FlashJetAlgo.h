#ifndef RecoJets_FlashJet_plugins_alpaka_FlashJetAlgo_h
#define RecoJets_FlashJet_plugins_alpaka_FlashJetAlgo_h

#include "DataFormats/FlashJet/interface/alpaka/FlashJetDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class FlashJetAlgo {
  public:
    FlashJetAlgo(double R, double p) : R_(R), p_(p) {}

    // clusters the particles held in `collection` (one event) in place
    void cluster(Queue& queue, flashjet::FlashJetDeviceCollection& collection) const;

  private:
    const double R_;
    const double p_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoJets_FlashJet_plugins_alpaka_FlashJetAlgo_h
