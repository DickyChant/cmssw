// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      HIEvtPlaneTableProducer
//
/**\class HIEvtPlaneTableProducer HIEvtPlaneTableProducer.cc

 Heavy-ion event-plane angles for NanoAOD, computed directly from MiniAOD content.

 The standard RecoHI/HiEvtPlaneAlgos EvtPlaneProducer, when fed MiniAOD (pp_on_AA
 modifier), only fills the *tracker* event planes from packedPFCandidates — its HF /
 Castor Q-vectors require a reco::PFCandidateCollection ("particleFlow") that is absent
 in MiniAOD. Since the HF event plane is the standard reference plane for PbPb flow,
 this producer reconstructs the full official set of planes (hi::EPNames) from MiniAOD:
   - HF planes  (HFm{n}, HFp{n}, HF{n}): the neutral HF PF deposits (pdgId 1 or 2) in
     the forward acceptance, transverse-energy weighted (et = E/cosh(eta));
   - tracker planes (trackmid{n}, trackm{n}, trackp{n}): charged packedPFCandidates and
     lostTracks with track details, PV-associated (fromPV>=PVTight) and high purity, pT
     weighted. This is the MiniAOD-available approximation of the full HI track-quality
     selection (which needs dz/dxy significance, chi2/layer, Nhits not all stored here).
 For each plane the raw Q-vector is Q_n = sum_i w_i (cos n*phi_i, sin n*phi_i); the angle
 is Psi_n = atan2(Qy, Qx)/n. The eta windows, pT windows and harmonics match exactly the
 definitions in RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h.

 These are RAW, UN-FLATTENED planes: no run-by-run recentering or flattening is applied
 (that needs the HeavyIonRPRcd DB and is always done offline from the Q-vectors). The
 Qx/Qy components are stored precisely so the downstream user can recenter/flatten.

 Produces a singleton FlatTable "EvtPlane": per plane <name> the columns
 EvtPlane_psi_<name>, EvtPlane_Qx_<name>, EvtPlane_Qy_<name>, EvtPlane_sumw_<name>.
*/
//
// Original Author:  HIN NanoAOD (forest_to_nano)
//

#include <array>
#include <cmath>
#include <memory>
#include <string>

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

class HIEvtPlaneTableProducer : public edm::stream::EDProducer<> {
public:
  explicit HIEvtPlaneTableProducer(const edm::ParameterSet&);
  ~HIEvtPlaneTableProducer() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  enum Det { kHF, kTracker };
  struct EPDef {
    const char* name;
    Det det;
    int n;                            // harmonic order
    double eta1min, eta1max;          // primary eta window
    double eta2min, eta2max;          // secondary eta window (combined HF planes; empty if min==max)
    double tmin, tmax;                // transverse (et for HF, pt for tracker) acceptance
  };
  // closed eta intervals [min,max], matching the official EvtPlaneProducer passEta()
  static bool inEta(double eta, const EPDef& e) {
    if (eta >= e.eta1min && eta <= e.eta1max)
      return true;
    if (e.eta2max > e.eta2min && eta >= e.eta2min && eta <= e.eta2max)
      return true;
    return false;
  }
  // The official 12 planes, copied verbatim from
  // RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h (hi::EPNames).
  static constexpr std::array<EPDef, 12> kEP{{
      {"HFm2", kHF, 2, -5.0, -3.0, 0.0, 0.0, 0.01, 30.0},
      {"HFp2", kHF, 2, 3.0, 5.0, 0.0, 0.0, 0.01, 30.0},
      {"HF2", kHF, 2, -5.0, -3.0, 3.0, 5.0, 0.01, 30.0},
      {"trackmid2", kTracker, 2, -0.5, 0.5, 0.0, 0.0, 0.5, 3.0},
      {"trackm2", kTracker, 2, -2.0, -1.0, 0.0, 0.0, 0.5, 3.0},
      {"trackp2", kTracker, 2, 1.0, 2.0, 0.0, 0.0, 0.5, 3.0},
      {"HFm3", kHF, 3, -5.0, -3.0, 0.0, 0.0, 0.01, 30.0},
      {"HFp3", kHF, 3, 3.0, 5.0, 0.0, 0.0, 0.01, 30.0},
      {"HF3", kHF, 3, -5.0, -3.0, 3.0, 5.0, 0.01, 30.0},
      {"trackmid3", kTracker, 3, -0.5, 0.5, 0.0, 0.0, 0.5, 3.0},
      {"trackm3", kTracker, 3, -2.0, -1.0, 0.0, 0.0, 0.5, 3.0},
      {"trackp3", kTracker, 3, 1.0, 2.0, 0.0, 0.0, 0.5, 3.0}}};

  void accumulate(const edm::View<pat::PackedCandidate>& cands,
                  bool allowHF,
                  std::array<double, 12>& qx,
                  std::array<double, 12>& qy,
                  std::array<double, 12>& sumw) const;

  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfToken_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> lostToken_;
  const int minPVassoc_;  // minimum PackedCandidate::fromPV() quality for tracker planes
};

HIEvtPlaneTableProducer::HIEvtPlaneTableProducer(const edm::ParameterSet& iConfig)
    : pfToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("pfCandidates"))),
      lostToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("lostTracks"))),
      minPVassoc_(iConfig.getParameter<int>("minPVassociation")) {
  produces<nanoaod::FlatTable>();
}

void HIEvtPlaneTableProducer::accumulate(const edm::View<pat::PackedCandidate>& cands,
                                         bool allowHF,
                                         std::array<double, 12>& qx,
                                         std::array<double, 12>& qy,
                                         std::array<double, 12>& sumw) const {
  for (const auto& pf : cands) {
    const double eta = pf.eta();
    const int pdg = pf.pdgId();
    // HF deposits: neutral HF PF candidates (pdgId 1=h_HF, 2=egamma_HF) in 3<=|eta|<=5
    const bool isHF = allowHF && (pdg == 1 || pdg == 2) && std::abs(eta) >= 3.0 && std::abs(eta) <= 5.0;
    // tracker candidates: charged, with track details, PV-associated and high purity
    // (the MiniAOD-available approximation of the full HI track-quality selection)
    const bool isTrk =
        (pf.charge() != 0) && pf.hasTrackDetails() && (pf.fromPV() >= minPVassoc_) && pf.trackHighPurity();
    if (!isHF && !isTrk)
      continue;
    const double phi = pf.phi();
    const double et = pf.et();   // transverse energy (HF weight); == pt for massless HF deposits
    const double pt = pf.pt();   // tracker weight

    for (size_t i = 0; i < kEP.size(); ++i) {
      const auto& e = kEP[i];
      const bool wantHF = (e.det == kHF);
      if (wantHF ? !isHF : !isTrk)
        continue;
      const double t = wantHF ? et : pt;
      if (t < e.tmin || t > e.tmax)  // closed [tmin,tmax], matching the official cuts
        continue;
      if (!inEta(eta, e))
        continue;
      const double w = wantHF ? et : pt;
      qx[i] += w * std::cos(e.n * phi);
      qy[i] += w * std::sin(e.n * phi);
      sumw[i] += w;
    }
  }
}

void HIEvtPlaneTableProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  std::array<double, 12> qx{}, qy{}, sumw{};
  // packedPFCandidates carry both the HF deposits and the bulk of the tracks;
  // lostTracks add the charged tracks not used in PF (tracker planes only).
  accumulate(iEvent.get(pfToken_), /*allowHF=*/true, qx, qy, sumw);
  accumulate(iEvent.get(lostToken_), /*allowHF=*/false, qx, qy, sumw);

  auto out = std::make_unique<nanoaod::FlatTable>(1, "EvtPlane", /*singleton=*/true, /*extension=*/false);
  for (size_t i = 0; i < kEP.size(); ++i) {
    const auto& e = kEP[i];
    const float psi = (qx[i] != 0. || qy[i] != 0.) ? std::atan2(qy[i], qx[i]) / e.n : 0.f;
    const std::string nm = e.name;
    out->addColumnValue<float>(
        "psi_" + nm, psi, "raw event-plane angle Psi_" + std::to_string(e.n) + " for " + nm + " (0 if sumw==0)");
    out->addColumnValue<float>("Qx_" + nm, qx[i], "raw Q-vector x (sum w cos n*phi) for " + nm);
    out->addColumnValue<float>("Qy_" + nm, qy[i], "raw Q-vector y (sum w sin n*phi) for " + nm);
    out->addColumnValue<float>("sumw_" + nm, sumw[i], "sum of weights (et for HF, pt for tracker) for " + nm);
  }
  iEvent.put(std::move(out));
}

void HIEvtPlaneTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pfCandidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("lostTracks", edm::InputTag("lostTracks"));
  desc.add<int>("minPVassociation", 2)
      ->setComment("min pat::PackedCandidate::fromPV() quality (2 = PVTight) for tracker-plane charged candidates");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(HIEvtPlaneTableProducer);
