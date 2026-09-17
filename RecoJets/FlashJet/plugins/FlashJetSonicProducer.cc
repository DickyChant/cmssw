#include <span>
#include <string>
#include <vector>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/SonicTriton/interface/TritonEDProducer.h"
#include "RecoJets/FlashJet/interface/FlashJetInputs.h"
#include "RecoJets/FlashJet/interface/FlashJetRecoJets.h"

// Sends the candidates of one event to the FlashJet model on a Triton
// inference server and turns the returned particle -> jet map into reco jets.
//
// Model contract (RecoJets/FlashJet/data/models/flashjet/config.pbtxt):
//   inputs:  p4    TYPE_FP64 [-1, 4]  (px, py, pz, E) per particle
//            algo  TYPE_FP64 [2]      (R, p)
//   outputs: jet_idx  TYPE_INT32 [-1]  jet index per particle (beam-merge order)
//            n_jets   TYPE_INT32 [1]
class FlashJetSonicProducer : public TritonEDProducer<> {
public:
  explicit FlashJetSonicProducer(edm::ParameterSet const& config)
      : TritonEDProducer<>(config),
        srcToken_{consumes(config.getParameter<edm::InputTag>("src"))},
        rParam_{config.getParameter<double>("rParam")},
        exponent_{flashjet::exponentOf(config.getParameter<std::string>("jetAlgorithm"))},
        inputPtMin_{config.getParameter<double>("inputPtMin")},
        writer_{
            producesCollector(), config.getParameter<std::string>("jetType"), config.getParameter<double>("jetPtMin")} {
  }

  void acquire(edm::Event const& event, edm::EventSetup const&, Input& input) override {
    selected_ = flashjet::selectInputs(event.get(srcToken_), inputPtMin_);
    const int64_t n = selected_.size();
    if (n == 0) {
      // batch size 0 skips the inference call
      client_->setBatchSize(0);
      return;
    }
    client_->setBatchSize(1);
    auto const& cands = event.get(srcToken_);

    auto& p4 = input.at("p4");
    p4.setShape(0, n);
    auto p4Data = p4.allocate<double>();
    auto& flat = (*p4Data)[0];
    flat.reserve(4 * n);
    for (int32_t k : selected_) {
      auto const& c = cands[k];
      flat.push_back(c.px());
      flat.push_back(c.py());
      flat.push_back(c.pz());
      flat.push_back(c.energy());
    }
    p4.toServer(p4Data);

    auto& algo = input.at("algo");
    auto algoData = algo.allocate<double>();
    (*algoData)[0] = {rParam_, exponent_};
    algo.toServer(algoData);
  }

  void produce(edm::Event& event, edm::EventSetup const&, Output const& output) override {
    auto const& cands = event.get(srcToken_);
    std::vector<flashjet::JetInfo> jets;
    if (!selected_.empty()) {
      auto const& jetIdx = output.at("jet_idx").fromServer<int32_t>();
      auto const& nJets = output.at("n_jets").fromServer<int32_t>();
      jets = flashjet::groupJets(cands, selected_, jetIdx[0], nJets[0][0]);
    }
    writer_.write(event, std::move(jets));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    TritonClient::fillPSetDescription(desc);
    desc.add<edm::InputTag>("src", edm::InputTag("particleFlow"));
    desc.add<std::string>("jetAlgorithm", "AntiKt")->setComment("AntiKt, Kt or CambridgeAachen");
    desc.add<double>("rParam", 0.4);
    desc.add<double>("inputPtMin", 0.);
    desc.add<std::string>("jetType", "PFJet")->setComment("PFJet, GenJet or BasicJet");
    desc.add<double>("jetPtMin", 5.);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<edm::View<reco::Candidate>> srcToken_;
  const double rParam_;
  const double exponent_;
  const double inputPtMin_;
  const flashjet::RecoJetWriter writer_;
  std::vector<int32_t> selected_;
};

DEFINE_FWK_MODULE(FlashJetSonicProducer);
