#include <cmath>
#include <vector>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/contrib/SoftDrop.hh>

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "RecoJets/FlashJet/interface/FlashJetValueMaps.h"

// Reference for FlashJetReclusterProducer + FlashJetSoftDropProducer: the same
// ValueMaps computed per jet with FastJet (C/A reclustering of the
// constituents with max_allowable_R) and fastjet::contrib::SoftDrop.
class FastjetSoftDropProducer : public edm::global::EDProducer<> {
public:
  explicit FastjetSoftDropProducer(edm::ParameterSet const& config)
      : jetsToken_{consumes(config.getParameter<edm::InputTag>("jets"))},
        jetPtMin_{config.getParameter<double>("jetPtMin")},
        zcut_{config.getParameter<double>("zcut")},
        beta_{config.getParameter<double>("beta")},
        R0_{config.getParameter<double>("R0")},
        maps_{producesCollector()} {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    auto jets = event.getHandle(jetsToken_);
    auto values = flashjet::SoftDropValueMaps::make(jets->size());
    const fastjet::JetDefinition ca(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
    fastjet::contrib::SoftDrop softDrop(beta_, zcut_, R0_);
    softDrop.set_verbose_structure(true);
    std::vector<fastjet::PseudoJet> inputs;
    for (size_t j = 0; j < jets->size(); ++j) {
      auto const& jet = (*jets)[j];
      if (jet.pt() < jetPtMin_ || jet.numberOfDaughters() == 0)
        continue;
      inputs.clear();
      for (size_t k = 0; k < jet.numberOfDaughters(); ++k) {
        auto const* c = jet.daughter(k);
        inputs.emplace_back(c->px(), c->py(), c->pz(), c->energy());
      }
      fastjet::ClusterSequence cs(inputs, ca);
      auto const reclustered = fastjet::sorted_by_pt(cs.inclusive_jets(0.));
      const fastjet::PseudoJet groomed = softDrop(reclustered[0]);
      values[flashjet::SoftDropValueMaps::kMass][j] = groomed.m();
      values[flashjet::SoftDropValueMaps::kPt][j] = groomed.pt();
      const bool split = groomed.has_structure_of<fastjet::contrib::SoftDrop>() && groomed.has_pieces();
      values[flashjet::SoftDropValueMaps::kZg][j] =
          split ? groomed.structure_of<fastjet::contrib::SoftDrop>().symmetry() : 0.;
      values[flashjet::SoftDropValueMaps::kRg][j] =
          split ? groomed.structure_of<fastjet::contrib::SoftDrop>().delta_R() : 0.;
      values[flashjet::SoftDropValueMaps::kNDropped][j] =
          groomed.has_structure_of<fastjet::contrib::SoftDrop>()
              ? groomed.structure_of<fastjet::contrib::SoftDrop>().dropped_count()
              : 0;
    }
    maps_.put(event, jets, values);
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("jets", edm::InputTag("ak8PFJetsPuppi"));
    desc.add<double>("jetPtMin", 0.);
    desc.add<double>("zcut", 0.1);
    desc.add<double>("beta", 0.);
    desc.add<double>("R0", 0.8);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<edm::View<reco::Jet>> jetsToken_;
  const double jetPtMin_, zcut_, beta_, R0_;
  const flashjet::SoftDropValueMaps maps_;
};

DEFINE_FWK_MODULE(FastjetSoftDropProducer);
