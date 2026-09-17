#include <span>
#include <string>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/FlashJet/interface/FlashJetHostCollection.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "RecoJets/FlashJet/interface/FlashJetRecoJets.h"

// Converts the host copy of the FlashJetProducer output into reco jets.
class FlashJetRecoJetProducer : public edm::global::EDProducer<> {
public:
  explicit FlashJetRecoJetProducer(edm::ParameterSet const& config)
      : srcToken_{consumes(config.getParameter<edm::InputTag>("src"))},
        clusterToken_{consumes(config.getParameter<edm::InputTag>("clusters"))},
        writer_{
            producesCollector(), config.getParameter<std::string>("jetType"), config.getParameter<double>("jetPtMin")} {
  }

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    auto const& cands = event.get(srcToken_);
    auto const& clusters = event.get(clusterToken_);
    auto const view = clusters.const_view();
    const size_t n = view.metadata().size();
    const int32_t nJets = n > 0 ? view.nJets() : 0;

    auto jets = flashjet::groupJets(cands, view.candIdx(), view.jetIdx(), nJets);
    // use the sums accumulated in merge order, as FastJet does
    for (int32_t j = 0; j < nJets; ++j)
      jets[j].p4.SetPxPyPzE(view.jetPx()[j], view.jetPy()[j], view.jetPz()[j], view.jetE()[j]);
    writer_.write(event, std::move(jets));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("particleFlow"))
        ->setComment("the candidates that were clustered (same as the FlashJetProducer src)");
    desc.add<edm::InputTag>("clusters", edm::InputTag("flashJetProducer"));
    desc.add<std::string>("jetType", "PFJet")->setComment("PFJet, GenJet or BasicJet");
    desc.add<double>("jetPtMin", 5.);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<edm::View<reco::Candidate>> srcToken_;
  const edm::EDGetTokenT<flashjet::FlashJetHostCollection> clusterToken_;
  const flashjet::RecoJetWriter writer_;
};

DEFINE_FWK_MODULE(FlashJetRecoJetProducer);
