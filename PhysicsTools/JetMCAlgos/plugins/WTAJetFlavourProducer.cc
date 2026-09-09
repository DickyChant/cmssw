#include <cmath>
#include <memory>
#include <vector>

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "PhysicsTools/JetMCAlgos/interface/WTAFlavour.h"

// The source MUST be jets made of physical, terminal shower partons.
// Hadron-level GenJets, reco jets and ghost-enriched jets are not valid inputs.
// A separate map preserves existing jet kinematics and legacy flavour labels.
class WTAJetFlavourProducer : public edm::stream::EDProducer<> {
public:
  explicit WTAJetFlavourProducer(const edm::ParameterSet& cfg)
      : jets_(consumes<edm::View<reco::Jet>>(cfg.getParameter<edm::InputTag>("src"))) {
    produces<edm::ValueMap<int>>();
    produces<edm::ValueMap<int>>("winnerIndex");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("ak4PartonJets"));
    descriptions.add("wtaJetFlavour", desc);
  }

private:
  void produce(edm::Event& event, const edm::EventSetup&) override {
    const auto jets = event.getHandle(jets_);
    std::vector<int> flavours, winners;
    flavours.reserve(jets->size());
    winners.reserve(jets->size());
    for (const auto& jet : *jets) {
      const auto constituents = jet.getJetConstituents();
      std::vector<fastjet::PseudoJet> inputs;
      inputs.reserve(constituents.size());
      for (size_t i = 0; i < constituents.size(); ++i) {
        const auto& p = constituents[i];
        if (p.isNull() || !p.isAvailable())
          throw cms::Exception("InvalidWTAInput") << "Unavailable parton constituent " << i;
        const int absId = std::abs(p->pdgId());
        if (!((absId >= 1 && absId <= 5) || absId == 21))
          throw cms::Exception("InvalidWTAInput")
              << "Expected u,d,s,c,b or gluon partons; found PDG ID " << p->pdgId()
              << ". Use jets clustered from terminal shower partons, not hadron-level GenJets.";
        inputs.emplace_back(p->px(), p->py(), p->pz(), p->energy());
        inputs.back().set_user_index(i);
      }
      int winner;
      try {
        winner = jetflavour::wtaWinner(inputs);
      } catch (const std::invalid_argument& error) {
        throw cms::Exception("InvalidWTAInput") << error.what();
      }
      winners.push_back(winner);
      flavours.push_back(winner < 0 ? 0 : constituents.at(winner)->pdgId());
    }
    auto flavourMap = std::make_unique<edm::ValueMap<int>>();
    edm::ValueMap<int>::Filler flavourFiller(*flavourMap);
    flavourFiller.insert(jets, flavours.begin(), flavours.end());
    flavourFiller.fill();
    event.put(std::move(flavourMap));
    auto winnerMap = std::make_unique<edm::ValueMap<int>>();
    edm::ValueMap<int>::Filler winnerFiller(*winnerMap);
    winnerFiller.insert(jets, winners.begin(), winners.end());
    winnerFiller.fill();
    event.put(std::move(winnerMap), "winnerIndex");
  }
  const edm::EDGetTokenT<edm::View<reco::Jet>> jets_;
};

DEFINE_FWK_MODULE(WTAJetFlavourProducer);
