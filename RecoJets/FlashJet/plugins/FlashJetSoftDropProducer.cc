#include <cmath>

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/FlashJet/interface/FlashJetHostCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "RecoJets/FlashJet/interface/FlashJetValueMaps.h"

// Soft-drop ValueMaps (mass, pt, zg, rg, nDropped) for the jets reclustered
// and groomed by FlashJetReclusterProducer@alpaka.
class FlashJetSoftDropProducer : public edm::global::EDProducer<> {
public:
  explicit FlashJetSoftDropProducer(edm::ParameterSet const& config)
      : jetsToken_{consumes(config.getParameter<edm::InputTag>("jets"))},
        clusterToken_{consumes(config.getParameter<edm::InputTag>("clusters"))},
        maps_{producesCollector()} {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    auto jets = event.getHandle(jetsToken_);
    auto const entries = event.get(clusterToken_).const_view().entries();
    auto values = flashjet::SoftDropValueMaps::make(jets->size());
    for (int32_t b = 0; b < entries.metadata().size(); ++b) {
      const int32_t j = entries.source()[b];
      const double px = entries.groomedPx()[b], py = entries.groomedPy()[b];
      const double pz = entries.groomedPz()[b], e = entries.groomedE()[b];
      const double m2 = e * e - px * px - py * py - pz * pz;
      values[flashjet::SoftDropValueMaps::kMass][j] = m2 > 0 ? std::sqrt(m2) : -std::sqrt(-m2);
      values[flashjet::SoftDropValueMaps::kPt][j] = std::hypot(px, py);
      values[flashjet::SoftDropValueMaps::kZg][j] = entries.zg()[b];
      values[flashjet::SoftDropValueMaps::kRg][j] = entries.rg()[b];
      values[flashjet::SoftDropValueMaps::kNDropped][j] = entries.nDropped()[b];
    }
    maps_.put(event, jets, values);
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("jets", edm::InputTag("ak8PFJetsPuppi"))
        ->setComment("the jets given to FlashJetReclusterProducer");
    desc.add<edm::InputTag>("clusters", edm::InputTag("flashJetReclusterProducer"));
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<edm::View<reco::Jet>> jetsToken_;
  const edm::EDGetTokenT<flashjet::FlashJetHostCollection> clusterToken_;
  const flashjet::SoftDropValueMaps maps_;
};

DEFINE_FWK_MODULE(FlashJetSoftDropProducer);
