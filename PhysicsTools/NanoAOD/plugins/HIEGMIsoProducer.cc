// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      HIEGMIsoProducer
//
/**\class HIEGMIsoProducer HIEGMIsoProducer.cc

 Heavy-ion e/gamma cone isolation ValueMaps for NanoAOD.

 Ports the core HiForest ggHiNtuplizer custom-cone PF isolation (pfIsoCalculator)
 to ValueMap<float>s that can be attached as NanoAOD extension columns to the
 standard Electron / Photon tables. For each input e/gamma it sums the pt of
 packedPFCandidates within dR cones (0.3 and 0.4) separately for charged hadrons,
 photons and neutral hadrons. Charged-hadron candidates are required to be
 compatible with the leading primary vertex (|dz| < dzMax, |dxy| < dxyMax), as in
 the forest. (The UE-subtracted / multi-cone variants are not reproduced here.)

 Output instance labels: chIso03, chIso04, phoIso03, phoIso04, nhIso03, nhIso04.
*/
//
// Original Author:  HIN NanoAOD (forest_to_nano), ported from ggHiNtuplizer/pfIsoCalculator
//

#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class HIEGMIsoProducer : public edm::stream::EDProducer<> {
public:
  explicit HIEGMIsoProducer(const edm::ParameterSet&);
  ~HIEGMIsoProducer() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  template <typename T>
  void put(edm::Event&, const edm::Handle<edm::View<reco::Candidate>>&, const std::vector<T>&, const std::string&);

  const edm::EDGetTokenT<edm::View<reco::Candidate>> srcToken_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  const double dzMax_, dxyMax_;
};

HIEGMIsoProducer::HIEGMIsoProducer(const edm::ParameterSet& iConfig)
    : srcToken_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("src"))),
      pfToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("pfCandidates"))),
      vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      dzMax_(iConfig.getParameter<double>("dzMax")),
      dxyMax_(iConfig.getParameter<double>("dxyMax")) {
  produces<edm::ValueMap<float>>("chIso03");
  produces<edm::ValueMap<float>>("chIso04");
  produces<edm::ValueMap<float>>("phoIso03");
  produces<edm::ValueMap<float>>("phoIso04");
  produces<edm::ValueMap<float>>("nhIso03");
  produces<edm::ValueMap<float>>("nhIso04");
}

template <typename T>
void HIEGMIsoProducer::put(edm::Event& iEvent,
                           const edm::Handle<edm::View<reco::Candidate>>& src,
                           const std::vector<T>& vals,
                           const std::string& label) {
  auto vm = std::make_unique<edm::ValueMap<float>>();
  typename edm::ValueMap<float>::Filler filler(*vm);
  filler.insert(src, vals.begin(), vals.end());
  filler.fill();
  iEvent.put(std::move(vm), label);
}

void HIEGMIsoProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<edm::View<reco::Candidate>> src;
  iEvent.getByToken(srcToken_, src);
  const auto& pfs = iEvent.get(pfToken_);
  const auto& vtx = iEvent.get(vtxToken_);
  const bool hasPV = !vtx.empty();
  const auto pv = hasPV ? vtx.front().position() : reco::Vertex::Point();

  const size_t n = src->size();
  std::vector<float> ch03(n, 0), ch04(n, 0), ph03(n, 0), ph04(n, 0), nh03(n, 0), nh04(n, 0);

  for (size_t i = 0; i < n; ++i) {
    const reco::Candidate& eg = (*src)[i];
    for (const auto& pf : pfs) {
      const double dr2 = reco::deltaR2(eg.eta(), eg.phi(), pf.eta(), pf.phi());
      if (dr2 > 0.16)  // outside the largest (0.4) cone
        continue;
      const int apdg = std::abs(pf.pdgId());
      const double pt = pf.pt();
      const bool in03 = (dr2 < 0.09);
      if (apdg == 211) {  // charged hadron: require compatibility with the leading PV
        if (hasPV && (std::abs(pf.dz(pv)) > dzMax_ || std::abs(pf.dxy(pv)) > dxyMax_))
          continue;
        ch04[i] += pt;
        if (in03)
          ch03[i] += pt;
      } else if (apdg == 22) {  // photon
        ph04[i] += pt;
        if (in03)
          ph03[i] += pt;
      } else if (apdg == 130) {  // neutral hadron
        nh04[i] += pt;
        if (in03)
          nh03[i] += pt;
      }
    }
  }

  put(iEvent, src, ch03, "chIso03");
  put(iEvent, src, ch04, "chIso04");
  put(iEvent, src, ph03, "phoIso03");
  put(iEvent, src, ph04, "phoIso04");
  put(iEvent, src, nh03, "nhIso03");
  put(iEvent, src, nh04, "nhIso04");
}

void HIEGMIsoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("finalElectrons"));
  desc.add<edm::InputTag>("pfCandidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<double>("dzMax", 0.2);
  desc.add<double>("dxyMax", 0.1);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(HIEGMIsoProducer);
