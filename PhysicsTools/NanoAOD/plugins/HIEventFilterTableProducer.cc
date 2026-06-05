// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      HIEventFilterTableProducer
//
/**\class HIEventFilterTableProducer HIEventFilterTableProducer.cc

 Heavy-ion event-selection flags for NanoAOD (top-level Flag_* + HF tower counts).

 Reproduces the MiniAOD-doable HiForest collision-event-selection decisions directly
 from packedPFCandidates and the primary vertices, avoiding extra filter paths:
   - the PF-based HF coincidence (HiHFFilterPF): per side (HF+ / HF-, 3<|eta|<6) count
     the neutral HF PF candidates (pdgId 1 or 2) with energy >= threshold (2/3/4/5 GeV);
     `hiHFnTowerThN` is the min of the two sides, and a coincidence flag at minnumtowers
     M is (hiHFnTowerThN >= M);
   - the primary-vertex filter (VertexSelector): leading vertex !isFake && |z|<=25 &&
     rho<=2.
 Produces a singleton FlatTable with empty name, so the columns are written as
 top-level branches (Flag_*, hiHFnTower*), as for the standard NanoAOD MET flags.
*/
//
// Original Author:  HIN NanoAOD (forest_to_nano), from HiHFFilterPF + collisionEventSelection
//

#include <array>
#include <cmath>
#include <memory>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class HIEventFilterTableProducer : public edm::stream::EDProducer<> {
public:
  explicit HIEventFilterTableProducer(const edm::ParameterSet&);
  ~HIEventFilterTableProducer() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  const double pvMaxZ_, pvMaxRho_;
  static constexpr std::array<double, 4> kThresh{{2., 3., 4., 5.}};  // GeV
};

HIEventFilterTableProducer::HIEventFilterTableProducer(const edm::ParameterSet& iConfig)
    : pfToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("pfCandidates"))),
      vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      pvMaxZ_(iConfig.getParameter<double>("pvMaxZ")),
      pvMaxRho_(iConfig.getParameter<double>("pvMaxRho")) {
  produces<nanoaod::FlatTable>();
}

void HIEventFilterTableProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  // HF coincidence tower counts per threshold (HF+ and HF- separately)
  std::array<int, 4> nPlus{{0, 0, 0, 0}}, nMinus{{0, 0, 0, 0}};
  for (const auto& pf : iEvent.get(pfToken_)) {
    const int pdg = pf.pdgId();
    if (pdg != 1 && pdg != 2)  // neutral HF PF candidates only
      continue;
    const double aeta = std::abs(pf.eta());
    if (aeta < 3. || aeta > 6.)
      continue;
    const double e = pf.energy();
    auto& side = (pf.eta() > 0) ? nPlus : nMinus;
    for (size_t t = 0; t < kThresh.size(); ++t)
      if (e >= kThresh[t])
        side[t] += 1;
  }
  std::array<int, 4> nMin{};
  for (size_t t = 0; t < kThresh.size(); ++t)
    nMin[t] = std::min(nPlus[t], nMinus[t]);

  // primary-vertex filter: pass if ANY vertex satisfies the cut (matches the forest
  // VertexSelector with filter=True), not only the leading one.
  bool pvPass = false;
  for (const auto& v : iEvent.get(vtxToken_)) {
    if (!v.isFake() && std::abs(v.z()) <= pvMaxZ_ && v.position().Rho() <= pvMaxRho_) {
      pvPass = true;
      break;
    }
  }

  auto out = std::make_unique<nanoaod::FlatTable>(1, "", /*singleton=*/true, /*extension=*/false);
  out->addColumnValue<bool>("Flag_primaryVertexFilter", pvPass, "leading PV: !isFake && |z|<=25 && rho<=2");
  out->addColumnValue<int>("hiHFnTowerTh2", nMin[0], "min(HF+,HF-) neutral PF cands with E>=2 GeV in 3<|eta|<6");
  out->addColumnValue<int>("hiHFnTowerTh3", nMin[1], "min(HF+,HF-) neutral PF cands with E>=3 GeV in 3<|eta|<6");
  out->addColumnValue<int>("hiHFnTowerTh4", nMin[2], "min(HF+,HF-) neutral PF cands with E>=4 GeV in 3<|eta|<6");
  out->addColumnValue<int>("hiHFnTowerTh5", nMin[3], "min(HF+,HF-) neutral PF cands with E>=5 GeV in 3<|eta|<6");
  // common HF-coincidence working points (PF-based)
  out->addColumnValue<bool>("Flag_phfCoincFilter1Th4", nMin[2] >= 1, "PF HF coincidence: >=1 tower E>=4 GeV each side");
  out->addColumnValue<bool>("Flag_phfCoincFilter2Th4", nMin[2] >= 2, "PF HF coincidence: >=2 towers E>=4 GeV each side");
  out->addColumnValue<bool>("Flag_phfCoincFilter3Th4", nMin[2] >= 3, "PF HF coincidence: >=3 towers E>=4 GeV each side");
  out->addColumnValue<bool>("Flag_phfCoincFilter4Th4", nMin[2] >= 4, "PF HF coincidence: >=4 towers E>=4 GeV each side");
  iEvent.put(std::move(out));
}

void HIEventFilterTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pfCandidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<double>("pvMaxZ", 25.0);
  desc.add<double>("pvMaxRho", 2.0);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(HIEventFilterTableProducer);
