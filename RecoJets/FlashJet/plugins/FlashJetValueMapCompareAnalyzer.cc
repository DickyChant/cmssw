#include <atomic>
#include <cmath>
#include <string>
#include <vector>

#include "DataFormats/Common/interface/ValueMap.h"
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

// Compares ValueMap<float>s of two producers (instances `names`) jet by jet.
class FlashJetValueMapCompareAnalyzer : public edm::global::EDAnalyzer<> {
public:
  explicit FlashJetValueMapCompareAnalyzer(edm::ParameterSet const& config)
      : jetsToken_{consumes(config.getParameter<edm::InputTag>("jets"))},
        tolerance_{config.getParameter<double>("tolerance")},
        absoluteTolerance_{config.getParameter<double>("absoluteTolerance")},
        failOnMismatch_{config.getParameter<bool>("failOnMismatch")} {
    const auto reference = config.getParameter<std::string>("reference");
    const auto test = config.getParameter<std::string>("test");
    for (auto const& name : config.getParameter<std::vector<std::string>>("names")) {
      names_.push_back(name);
      refTokens_.push_back(consumes<edm::ValueMap<float>>(edm::InputTag(reference, name)));
      testTokens_.push_back(consumes<edm::ValueMap<float>>(edm::InputTag(test, name)));
    }
  }

  void analyze(edm::StreamID, edm::Event const& event, edm::EventSetup const&) const override {
    auto jets = event.getHandle(jetsToken_);
    for (size_t k = 0; k < names_.size(); ++k) {
      auto const& ref = event.get(refTokens_[k]);
      auto const& test = event.get(testTokens_[k]);
      for (size_t j = 0; j < jets->size(); ++j) {
        const auto jet = jets->refAt(j);
        const float a = ref[jet];
        const float b = test[jet];
        ++values_;
        if (!std::isfinite(b) || !(std::abs(a - b) <= tolerance_ * std::max(1.f, std::abs(a)) + absoluteTolerance_)) {
          ++mismatches_;
          edm::LogWarning("FlashJetCompare")
              << "event " << event.id() << " jet " << j << " " << names_[k] << ": reference " << a << ", test " << b;
        }
      }
    }
  }

  void endJob() override {
    edm::LogSystem("FlashJetCompare") << "compared " << values_ << " values: " << mismatches_ << " mismatches";
    if (failOnMismatch_ && mismatches_ > 0)
      throw cms::Exception("FlashJetCompare") << mismatches_ << " value mismatches";
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("jets", edm::InputTag("ak8PFJetsPuppi"));
    desc.add<std::string>("reference", "fastjetSoftDrop");
    desc.add<std::string>("test", "flashJetSoftDrop");
    desc.add<std::vector<std::string>>("names", {"mass", "pt", "zg", "rg", "nDropped"});
    desc.add<double>("tolerance", 1e-5)->setComment("relative, values are stored as float");
    desc.add<double>("absoluteTolerance", 1e-4)
        ->setComment(
            "added to the tolerance: a groomed jet that is a single particle has a numerically zero mass "
            "whose sign follows the rounding");
    desc.add<bool>("failOnMismatch", false);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<edm::View<reco::Jet>> jetsToken_;
  const double tolerance_;
  const double absoluteTolerance_;
  const bool failOnMismatch_;
  std::vector<std::string> names_;
  std::vector<edm::EDGetTokenT<edm::ValueMap<float>>> refTokens_, testTokens_;
  mutable std::atomic<size_t> values_{0}, mismatches_{0};
};

DEFINE_FWK_MODULE(FlashJetValueMapCompareAnalyzer);
