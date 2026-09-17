#include <cmath>
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

  // Clusters the candidates of `src` into jets on the alpaka device and puts
  // a flashjet::FlashJetDeviceCollection (particle -> jet map, jet
  // four-momenta and merge history).  Use FlashJetRecoJetProducer to turn the
  // host copy into reco::*Jet collections.
  class FlashJetProducer : public global::EDProducer<> {
  public:
    FlashJetProducer(edm::ParameterSet const& config)
        : EDProducer(config),
          srcToken_{consumes(config.getParameter<edm::InputTag>("src"))},
          putToken_{produces()},
          inputPtMin_{config.getParameter<double>("inputPtMin")},
          algo_{config.getParameter<double>("rParam"),
                ::flashjet::exponentOf(config.getParameter<std::string>("jetAlgorithm"))} {}

    void produce(edm::StreamID, device::Event& event, device::EventSetup const&) const override {
      auto const& cands = event.get(srcToken_);
      auto const selected = ::flashjet::selectInputs(cands, inputPtMin_);
      const int32_t n = selected.size();

      ::flashjet::FlashJetHostCollection host{event.queue(), n};
      auto view = host.view();
      view.nJets() = 0;
      for (int32_t k = 0; k < n; ++k) {
        auto const& c = cands[selected[k]];
        view.px()[k] = c.px();
        view.py()[k] = c.py();
        view.pz()[k] = c.pz();
        view.e()[k] = c.energy();
        view.candIdx()[k] = selected[k];
      }

      flashjet::FlashJetDeviceCollection device{event.queue(), n};
      alpaka::memcpy(event.queue(), device.buffer(), host.const_buffer());
      algo_.cluster(event.queue(), device);
      event.emplace(putToken_, std::move(device));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("src", edm::InputTag("particleFlow"));
      desc.add<std::string>("jetAlgorithm", "AntiKt")->setComment("AntiKt, Kt or CambridgeAachen");
      desc.add<double>("rParam", 0.4);
      desc.add<double>("inputPtMin", 0.)->setComment("drop input candidates with pt below this");
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
