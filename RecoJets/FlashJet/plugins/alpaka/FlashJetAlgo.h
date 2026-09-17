#ifndef RecoJets_FlashJet_plugins_alpaka_FlashJetAlgo_h
#define RecoJets_FlashJet_plugins_alpaka_FlashJetAlgo_h

#include "DataFormats/FlashJet/interface/alpaka/FlashJetDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace edm {
  class ParameterSet;
}

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class FlashJetAlgo {
  public:
    struct SoftDrop {
      bool enable = false;
      double zcut = 0.1;
      double beta = 0.;
      double R0 = 0.8;
    };

    explicit FlashJetAlgo(edm::ParameterSet const& config);

    static void fillPSetDescription(edm::ParameterSetDescription& desc, std::string const& algorithm, double rParam);

    // Clusters every entry of `collection` in place, in one kernel call.  On
    // CPU backends one work item handles one entry; on GPU backends one block
    // does, and `maxEntrySize` (the largest entry, known to the caller that
    // filled the collection) sizes that block.
    void cluster(Queue& queue, flashjet::FlashJetDeviceCollection& collection, int32_t maxEntrySize = 0) const;

  private:
    const double R_;
    const double p_;
    SoftDrop softDrop_;
    const int32_t entriesPerThread_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoJets_FlashJet_plugins_alpaka_FlashJetAlgo_h
