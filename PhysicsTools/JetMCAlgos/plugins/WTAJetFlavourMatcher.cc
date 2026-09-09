#include <cmath>
#include <memory>
#include <vector>
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

// Optional stage-transfer map: an experimental association, not an IRC-safety
// statement. Require an unambiguous one-to-one geometrical match. Multiple
// parton jets are explicitly marked ambiguous, never reduced by heavy flavour
// priority. The original jets and all existing flavour products are unchanged.
class WTAJetFlavourMatcher : public edm::stream::EDProducer<> {
public:
  explicit WTAJetFlavourMatcher(const edm::ParameterSet& cfg)
      : targets_(consumes<edm::View<reco::Jet>>(cfg.getParameter<edm::InputTag>("jets"))),
        sources_(consumes<edm::View<reco::Jet>>(cfg.getParameter<edm::InputTag>("partonJets"))),
        flavours_(consumes<edm::ValueMap<int>>(cfg.getParameter<edm::InputTag>("flavour"))),
        maxDR_(cfg.getParameter<double>("maxDeltaR")) {
    if (!std::isfinite(maxDR_) || maxDR_ <= 0.) throw cms::Exception("Configuration") << "maxDeltaR must be positive and finite";
    produces<edm::ValueMap<int>>();
    produces<edm::ValueMap<int>>("matchStatus");
    produces<edm::ValueMap<int>>("sourceJetIndex");
  }
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("jets", edm::InputTag("ak4GenJets"));
    desc.add<edm::InputTag>("partonJets", edm::InputTag("ak4PartonJets"));
    desc.add<edm::InputTag>("flavour", edm::InputTag("wtaJetFlavour"));
    desc.add<double>("maxDeltaR", 0.2);
    descriptions.add("wtaJetFlavourMatcher", desc);
  }
private:
  void produce(edm::Event& event, const edm::EventSetup&) override {
    const auto targets = event.getHandle(targets_);
    const auto sources = event.getHandle(sources_);
    const auto& flavours = event.get(flavours_);
    std::vector<int> targetCounts(targets->size()), sourceCounts(sources->size());
    std::vector<int> candidate(targets->size(), -1), labels(targets->size(), 0);
    std::vector<int> indices(targets->size(), -1), status(targets->size(), 0);
    for (size_t i=0; i<targets->size(); ++i) {
      for (size_t j=0; j<sources->size(); ++j) {
        if (reco::deltaR2((*targets)[i], (*sources)[j]) < maxDR_*maxDR_) {
          ++targetCounts[i]; ++sourceCounts[j]; candidate[i] = j;
        }
      }
    }
    for (size_t i=0; i<targets->size(); ++i) {
      if (!targetCounts[i]) continue; // 0: no geometrical match
      const int j = candidate[i];
      status[i] = 2; // ambiguous (in either direction)
      if (targetCounts[i] != 1 || sourceCounts[j] != 1) continue;
      const int label = flavours[sources->refAt(j)];
      if (!label) { status[i] = 3; continue; } // empty/undefined source
      status[i] = 1;
      labels[i] = label;
      indices[i] = j;
    }
    auto put = [&](const std::vector<int>& values, const std::string& name) {
      auto map = std::make_unique<edm::ValueMap<int>>();
      edm::ValueMap<int>::Filler filler(*map);
      filler.insert(targets, values.begin(), values.end()); filler.fill();
      event.put(std::move(map), name);
    };
    put(labels, ""); put(status, "matchStatus"); put(indices, "sourceJetIndex");
  }
  const edm::EDGetTokenT<edm::View<reco::Jet>> targets_, sources_;
  const edm::EDGetTokenT<edm::ValueMap<int>> flavours_;
  const double maxDR_;
};
DEFINE_FWK_MODULE(WTAJetFlavourMatcher);
