#ifndef RecoJets_FlashJet_interface_FlashJetRecoJets_h
#define RecoJets_FlashJet_interface_FlashJetRecoJets_h

#include <algorithm>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ProducesCollector.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "RecoJets/JetProducers/interface/JetSpecific.h"

namespace flashjet {

  struct JetInfo {
    reco::Particle::LorentzVector p4;
    std::vector<reco::CandidatePtr> constituents;
  };

  // Groups the candidates by jet.  candIdx[k] is the View index of input k
  // and jetIdx[k] its jet; the jet four-momenta are summed from the
  // constituents (callers that have the in-merge-order sums overwrite them).
  inline std::vector<JetInfo> groupJets(edm::View<reco::Candidate> const& cands,
                                        std::span<const int32_t> candIdx,
                                        std::span<const int32_t> jetIdx,
                                        int32_t nJets) {
    std::vector<JetInfo> jets(nJets);
    for (size_t k = 0; k < jetIdx.size(); ++k) {
      const int32_t j = jetIdx[k];
      if (j < 0 || j >= nJets)
        throw cms::Exception("FlashJet") << "particle " << k << " has jet index " << j << " (nJets " << nJets << ")";
      auto ptr = cands.ptrAt(candIdx[k]);
      jets[j].p4 += ptr->p4();
      jets[j].constituents.push_back(std::move(ptr));
    }
    return jets;
  }

  // Writes reco::PFJet, reco::GenJet or reco::BasicJet collections the way
  // VirtualJetProducer does: jets above jetPtMin, sorted by decreasing pt,
  // constituents sorted by decreasing pt, reference point at the origin.
  class RecoJetWriter {
  public:
    RecoJetWriter(edm::ProducesCollector producer, std::string const& jetType, double jetPtMin) : jetPtMin_(jetPtMin) {
      if (jetType == "PFJet")
        pfToken_ = producer.produces<reco::PFJetCollection>();
      else if (jetType == "GenJet")
        genToken_ = producer.produces<reco::GenJetCollection>();
      else if (jetType == "BasicJet")
        basicToken_ = producer.produces<reco::BasicJetCollection>();
      else
        throw cms::Exception("Configuration")
            << "FlashJet: unsupported jetType '" << jetType << "' (use PFJet, GenJet or BasicJet)";
    }

    void write(edm::Event& event, std::vector<JetInfo>&& jets) const {
      if (!pfToken_.isUninitialized())
        event.emplace(pfToken_, build<reco::PFJet>(std::move(jets)));
      else if (!genToken_.isUninitialized())
        event.emplace(genToken_, build<reco::GenJet>(std::move(jets)));
      else
        event.emplace(basicToken_, build<reco::BasicJet>(std::move(jets)));
    }

  private:
    template <typename JetT>
    std::vector<JetT> build(std::vector<JetInfo>&& jets) const {
      const double ptMin2 = jetPtMin_ * jetPtMin_;
      std::erase_if(jets, [ptMin2](JetInfo const& j) { return j.p4.perp2() < ptMin2; });
      std::stable_sort(
          jets.begin(), jets.end(), [](JetInfo const& a, JetInfo const& b) { return a.p4.perp2() > b.p4.perp2(); });
      std::vector<JetT> out;
      out.reserve(jets.size());
      const reco::Particle::Point origin(0, 0, 0);
      for (auto& j : jets) {
        std::stable_sort(j.constituents.begin(),
                         j.constituents.end(),
                         [](reco::CandidatePtr const& a, reco::CandidatePtr const& b) { return a->pt() > b->pt(); });
        JetT jet;
        reco::writeSpecific(jet, j.p4, origin, j.constituents);
        out.push_back(std::move(jet));
      }
      return out;
    }

    const double jetPtMin_;
    edm::EDPutTokenT<reco::PFJetCollection> pfToken_;
    edm::EDPutTokenT<reco::GenJetCollection> genToken_;
    edm::EDPutTokenT<reco::BasicJetCollection> basicToken_;
  };

}  // namespace flashjet

#endif  // RecoJets_FlashJet_interface_FlashJetRecoJets_h
