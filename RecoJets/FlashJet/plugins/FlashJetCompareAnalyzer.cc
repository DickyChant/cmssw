#include <algorithm>
#include <atomic>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

// Compares two pt-ordered jet collections jet by jet: multiplicity, number of
// constituents and four-momentum.  Meant to validate FlashJet against
// FastjetJetProducer on the same inputs.
class FlashJetCompareAnalyzer : public edm::global::EDAnalyzer<> {
public:
  explicit FlashJetCompareAnalyzer(edm::ParameterSet const& config)
      : refToken_{consumes(config.getParameter<edm::InputTag>("reference"))},
        testToken_{consumes(config.getParameter<edm::InputTag>("test"))},
        ptMin_{config.getParameter<double>("ptMin")},
        tolerance_{config.getParameter<double>("tolerance")},
        failOnMismatch_{config.getParameter<bool>("failOnMismatch")} {}

  // the same constituents, by product and key, not merely as many
  static bool sameConstituents(reco::Jet const& a, reco::Jet const& b) {
    if (a.numberOfDaughters() != b.numberOfDaughters())
      return false;
    std::vector<std::pair<unsigned int, size_t>> ka, kb;
    ka.reserve(a.numberOfDaughters());
    kb.reserve(b.numberOfDaughters());
    for (size_t i = 0; i < a.numberOfDaughters(); ++i) {
      ka.emplace_back(a.daughterPtr(i).id().productIndex(), a.daughterPtr(i).key());
      kb.emplace_back(b.daughterPtr(i).id().productIndex(), b.daughterPtr(i).key());
    }
    std::sort(ka.begin(), ka.end());
    std::sort(kb.begin(), kb.end());
    return ka == kb;
  }

  void analyze(edm::StreamID, edm::Event const& event, edm::EventSetup const&) const override {
    auto const& ref = event.get(refToken_);
    auto const& test = event.get(testToken_);
    auto countAbove = [this](edm::View<reco::Jet> const& jets) {
      return std::count_if(jets.begin(), jets.end(), [this](reco::Jet const& j) { return j.pt() >= ptMin_; });
    };
    const size_t nRef = countAbove(ref);
    const size_t nTest = countAbove(test);
    size_t bad = (nRef != nTest) ? 1 : 0;
    const size_t n = std::min(nRef, nTest);
    for (size_t k = 0; k < n; ++k) {
      auto const& a = ref[k];
      auto const& b = test[k];
      const double scale = std::max(a.energy(), 1.);
      const double dp = std::abs(a.px() - b.px()) + std::abs(a.py() - b.py()) + std::abs(a.pz() - b.pz()) +
                        std::abs(a.energy() - b.energy());
      // NaN fails every comparison, so test for it rather than relying on one
      const bool finite =
          std::isfinite(b.px()) && std::isfinite(b.py()) && std::isfinite(b.pz()) && std::isfinite(b.energy());
      if (!finite || !(dp <= tolerance_ * scale) || !sameConstituents(a, b)) {
        ++bad;
        edm::LogWarning("FlashJetCompare")
            << "event " << event.id() << " jet " << k << ": reference pt " << a.pt() << " eta " << a.eta() << " n "
            << a.numberOfDaughters() << ", test pt " << b.pt() << " eta " << b.eta() << " n " << b.numberOfDaughters();
      }
    }
    events_ += 1;
    jets_ += nRef;
    mismatches_ += bad;
    if (bad > 0) {
      badEvents_ += 1;
      if (nRef != nTest)
        edm::LogWarning("FlashJetCompare") << "event " << event.id() << ": " << nRef << " reference jets, " << nTest
                                           << " test jets above " << ptMin_ << " GeV";
    }
  }

  void endJob() override {
    edm::LogSystem("FlashJetCompare") << "compared " << events_ << " events, " << jets_
                                      << " reference jets: " << mismatches_ << " mismatches in " << badEvents_
                                      << " events";
    if (failOnMismatch_ && mismatches_ > 0)
      throw cms::Exception("FlashJetCompare") << mismatches_ << " jet mismatches";
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("reference", edm::InputTag("ak4PFJets"));
    desc.add<edm::InputTag>("test", edm::InputTag("flashJetRecoJets"));
    desc.add<double>("ptMin", 5.);
    desc.add<double>("tolerance", 1e-9)->setComment("allowed sum |d(px,py,pz,E)| relative to the jet energy");
    desc.add<bool>("failOnMismatch", false);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<edm::View<reco::Jet>> refToken_;
  const edm::EDGetTokenT<edm::View<reco::Jet>> testToken_;
  const double ptMin_;
  const double tolerance_;
  const bool failOnMismatch_;
  mutable std::atomic<size_t> events_{0}, jets_{0}, mismatches_{0}, badEvents_{0};
};

DEFINE_FWK_MODULE(FlashJetCompareAnalyzer);
