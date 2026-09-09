#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include <memory>

// Extra columns only: never replace the legacy flavour table or its inputs.
class JetAlgorithmFlavourTableProducer : public edm::stream::EDProducer<> {
public:
  explicit JetAlgorithmFlavourTableProducer(const edm::ParameterSet& p)
      : jets_(consumes<reco::GenJetCollection>(p.getParameter<edm::InputTag>("src"))),
        name_(p.getParameter<std::string>("name")), cut_(p.getParameter<std::string>("cut")),
        useAlgorithms_(p.getParameter<bool>("useAlgorithms")), useWTA_(p.getParameter<bool>("useWTA")) {
    if (useAlgorithms_)
      algorithms_ = consumes<reco::JetFlavourInfoMatchingCollection>(p.getParameter<edm::InputTag>("algorithmFlavourInfos"));
    if (useWTA_) {
      wta_ = consumes<edm::ValueMap<int>>(p.getParameter<edm::InputTag>("wtaFlavour"));
      wtaStatus_ = consumes<edm::ValueMap<int>>(p.getParameter<edm::InputTag>("wtaMatchStatus"));
    }
    produces<nanoaod::FlatTable>();
  }
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription d;
    d.add<edm::InputTag>("src");
    d.add<std::string>("name", "GenJet");
    d.add<std::string>("cut", "");
    d.add<bool>("useAlgorithms", true);
    d.add<bool>("useWTA", false);
    d.add<edm::InputTag>("algorithmFlavourInfos", edm::InputTag());
    d.add<edm::InputTag>("wtaFlavour", edm::InputTag());
    d.add<edm::InputTag>("wtaMatchStatus", edm::InputTag());
    descriptions.add("jetAlgorithmFlavourTable", d);
  }
private:
  void produce(edm::Event& event, const edm::EventSetup&) override {
    const auto jets = event.getHandle(jets_);
    std::vector<const reco::JetFlavourInfo*> info(jets->size(), nullptr);
    std::vector<uint8_t> match(jets->size(), 0);
    if (useAlgorithms_) {
      for (const auto& entry : event.get(algorithms_)) {
        // Match by identity, never by first neighbour in a deltaR cone.
        if (entry.first.id() != jets.id()) continue;
        const auto key = entry.first.key();
        if (key >= jets->size()) throw cms::Exception("InvalidFlavourReference") << "Jet key outside source collection";
        if (match[key] != 0) { match[key] = 2; info[key] = nullptr; }
        else { match[key] = 1; info[key] = &entry.second; }
      }
    }
    const edm::ValueMap<int>* wta = useWTA_ ? &event.get(wta_) : nullptr;
    const edm::ValueMap<int>* wtaStatus = useWTA_ ? &event.get(wtaStatus_) : nullptr;
    std::vector<uint8_t> status, ghsValid, ifnValid;
    std::vector<uint32_t> ghsCode, ifnCode;
    std::vector<int16_t> ghsLeading, ifnLeading, wtaLabel, wtaMatch;
    for (size_t i = 0; i < jets->size(); ++i) {
      if (!cut_((*jets)[i])) continue;
      status.push_back(match[i]);
      auto append = [&](reco::FlavAlgo algo, auto& valid, auto& code, auto& leading) {
        const bool have = info[i] && info[i]->haveAlgoFlav(algo);
        valid.push_back(have);
        code.push_back(have ? info[i]->getAlgoFlavCode(algo) : 0);
        leading.push_back(have ? info[i]->getAlgoFlavLeading(algo) : 0);
      };
      if (useAlgorithms_) {
        append(reco::FlavAlgo::kGHS, ghsValid, ghsCode, ghsLeading);
        append(reco::FlavAlgo::kIFN, ifnValid, ifnCode, ifnLeading);
      }
      if (useWTA_) {
        const reco::GenJetRef ref(jets, i);
        wtaLabel.push_back((*wta)[ref]);
        wtaMatch.push_back((*wtaStatus)[ref]);
      }
    }
    auto table = std::make_unique<nanoaod::FlatTable>(status.size(), name_, false, true);
    if (useAlgorithms_) {
      table->addColumn<uint8_t>("algorithmMatchStatus", status, "0 unmatched, 1 exact jet reference, 2 duplicate association");
      table->addColumn<uint8_t>("GHSFlavValid", ghsValid, "GHS result exists; inspect before reading flavour");
      table->addColumn<uint8_t>("IFNFlavValid", ifnValid, "IFN result exists; unmatched is not flavourless");
      table->addColumn<uint32_t>("GHSFlavCode", ghsCode, "GHS flavour code; meaningful only when valid");
      table->addColumn<uint32_t>("IFNFlavCode", ifnCode, "IFN flavour code; meaningful only when valid");
      table->addColumn<int16_t>("GHSFlavLeading", ghsLeading, "Signed leading GHS flavour; meaningful only when valid");
      table->addColumn<int16_t>("IFNFlavLeading", ifnLeading, "Signed leading IFN flavour; meaningful only when valid");
    }
    if (useWTA_) {
      table->addColumn<int16_t>("WTAFlavour", wtaLabel, "Preferred WTA flavour definition; interpret with WTAMatchStatus");
      table->addColumn<int16_t>("WTAMatchStatus", wtaMatch, "0 unmatched, 1 unique, 2 ambiguous, 3 undefined source flavour");
    }
    event.put(std::move(table));
  }
  edm::EDGetTokenT<reco::GenJetCollection> jets_;
  std::string name_;
  StringCutObjectSelector<reco::GenJet> cut_;
  bool useAlgorithms_, useWTA_;
  edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> algorithms_;
  edm::EDGetTokenT<edm::ValueMap<int>> wta_, wtaStatus_;
};
DEFINE_FWK_MODULE(JetAlgorithmFlavourTableProducer);
