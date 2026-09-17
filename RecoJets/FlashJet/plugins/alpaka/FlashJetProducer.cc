#include <cstdint>
#include <string>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/FlashJet/interface/FlashJetHostCollection.h"
#include "DataFormats/FlashJet/interface/alpaka/FlashJetDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoJets/FlashJet/interface/FlashJetInputs.h"

#include "FlashJetAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // Clusters all candidates of `src` (one entry per event) on the alpaka
  // device and puts a flashjet::FlashJetDeviceCollection.  Use
  // FlashJetRecoJetProducer to turn the host copy into reco::*Jet collections.
  class FlashJetProducer : public global::EDProducer<> {
  public:
    FlashJetProducer(edm::ParameterSet const& config)
        : EDProducer(config),
          srcToken_{consumes(config.getParameter<edm::InputTag>("src"))},
          putToken_{produces()},
          inputPtMin_{config.getParameter<double>("inputPtMin")},
          algo_{config} {}

    void produce(edm::StreamID, device::Event& event, device::EventSetup const&) const override {
      auto const& cands = event.get(srcToken_);
      auto const selected = ::flashjet::selectInputs(cands, inputPtMin_);
      const int32_t n = selected.size();

      ::flashjet::FlashJetHostCollection host{event.queue(), n, 1};
      auto particles = host.view().particles();
      for (int32_t k = 0; k < n; ++k) {
        auto const& c = cands[selected[k]];
        particles.px()[k] = c.px();
        particles.py()[k] = c.py();
        particles.pz()[k] = c.pz();
        particles.e()[k] = c.energy();
        particles.candIdx()[k] = selected[k];
      }
      auto entries = host.view().entries();
      entries.offset()[0] = 0;
      entries.size()[0] = n;
      entries.source()[0] = -1;
      entries.nJets()[0] = 0;

      flashjet::FlashJetDeviceCollection device{event.queue(), n, 1};
      alpaka::memcpy(event.queue(), device.buffer(), host.const_buffer());
      algo_.cluster(event.queue(), device, n);
      event.emplace(putToken_, std::move(device));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("src", edm::InputTag("particleFlow"));
      desc.add<double>("inputPtMin", 0.)->setComment("drop input candidates with pt below this");
      FlashJetAlgo::fillPSetDescription(desc, "AntiKt", 0.4);
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    const edm::EDGetTokenT<edm::View<reco::Candidate>> srcToken_;
    const device::EDPutToken<flashjet::FlashJetDeviceCollection> putToken_;
    const double inputPtMin_;
    const FlashJetAlgo algo_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(FlashJetProducer);
