// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      HICaloJetTableProducer
//
/**\class HICaloJetTableProducer HICaloJetTableProducer.cc

 Heavy-ion calorimeter-jet table for NanoAOD, reproducing the HiForest HiCaloJetAnalyzer
 (akPu4Calo) content on MiniAOD.

 The forest reads slimmedCaloJets (the genuine pileup-subtracted akPu4Calo, which needs
 CaloTowers, is NOT available on MiniAOD) and, with doHiJetID=True, builds for each calo
 jet a PF/track "in-cone" composition by looping packedPFCandidates within dR<rParam of
 the calo-jet axis -- it does not use calo towers. This producer reproduces exactly that:
   - kinematics + jet area + the two surviving calo fractions (emEF/hadEF, from m_specific);
   - trackSum/Max/N and trackHardSum/N over good-quality packed candidates (pseudoTrack);
   - per PF type (charged h, e, mu, gamma, neutral h0): Sum/Max/N (+ hard sums where the
     forest keeps them), using reco::PFCandidate::translatePdgIdToType for the mapping.

 Produces a (non-singleton) FlatTable "CaloJet" aligned with the selected (pt>jetPtMin)
 slimmedCaloJets.
*/
//
// Original Author:  HIN NanoAOD (forest_to_nano), from HeavyIonsAnalysis HiCaloJetAnalyzer
//

#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

class HICaloJetTableProducer : public edm::stream::EDProducer<> {
public:
  explicit HICaloJetTableProducer(const edm::ParameterSet&);
  ~HICaloJetTableProducer() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<edm::View<reco::CaloJet>> jetToken_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfToken_;
  const std::string name_, doc_;
  const double jetPtMin_, jetAbsEtaMax_, r2Param_, hardPtMin_;
  const bool useQuality_;
  const reco::TrackBase::TrackQuality trackQuality_;
};

HICaloJetTableProducer::HICaloJetTableProducer(const edm::ParameterSet& iConfig)
    : jetToken_(consumes<edm::View<reco::CaloJet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      pfToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("pfCandidates"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")),
      jetPtMin_(iConfig.getParameter<double>("jetPtMin")),
      jetAbsEtaMax_(iConfig.getParameter<double>("jetAbsEtaMax")),
      r2Param_(iConfig.getParameter<double>("rParam") * iConfig.getParameter<double>("rParam")),
      hardPtMin_(iConfig.getParameter<double>("hardPtMin")),
      useQuality_(iConfig.getParameter<bool>("useQuality")),
      trackQuality_(reco::TrackBase::qualityByName(iConfig.getParameter<std::string>("trackQuality"))) {
  produces<nanoaod::FlatTable>();
}

void HICaloJetTableProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  const auto& jets = iEvent.get(jetToken_);
  const auto& pfs = iEvent.get(pfToken_);
  static const reco::PFCandidate kConverter;  // only for translatePdgIdToType

  std::vector<float> pt, eta, phi, mass, area, emEF, hadEF;
  std::vector<float> trackSum, trackMax, trackHardSum, chargedSum, chargedMax, chargedHardSum;
  std::vector<float> photonSum, photonMax, photonHardSum, neutralSum, neutralMax;
  std::vector<float> eSum, eMax, muSum, muMax;
  std::vector<int> trackN, trackHardN, chargedN, chargedHardN, photonN, photonHardN, neutralN, eN, muN;

  for (const auto& jet : jets) {
    if (jet.pt() < jetPtMin_ || std::abs(jet.eta()) > jetAbsEtaMax_)
      continue;

    float tSum = 0, tMax = 0, tHSum = 0, cSum = 0, cMax = 0, cHSum = 0;
    float gSum = 0, gMax = 0, gHSum = 0, nSum = 0, nMax = 0, eS = 0, eM = 0, mS = 0, mM = 0;
    int tN = 0, tHN = 0, cN = 0, cHN = 0, gN = 0, gHN = 0, nN = 0, eNc = 0, muNc = 0;

    for (const auto& cand : pfs) {
      // track-based block: good-quality packed candidates, using the pseudoTrack
      // direction (etaAtVtx/phiAtVtx) and pT -- exactly as the forest HiCaloJetAnalyzer
      // (deltaR2(jet, pseudoTrack) + track.pt()), not the packed-candidate 4-vector.
      if (cand.hasTrackDetails()) {
        const reco::Track& trk = cand.pseudoTrack();
        if ((!useQuality_ || trk.quality(trackQuality_)) &&
            reco::deltaR2(jet.eta(), jet.phi(), trk.eta(), trk.phi()) < r2Param_) {
          const double ptt = trk.pt();
          tSum += ptt;
          tN += 1;
          if (ptt > hardPtMin_) {
            tHSum += ptt;
            tHN += 1;
          }
          if (ptt > tMax)
            tMax = ptt;
        }
      }

      // PF-type block (charged h / e / mu / gamma / neutral h0): candidate direction + pT
      if (reco::deltaR2(jet.eta(), jet.phi(), cand.eta(), cand.phi()) >= r2Param_)
        continue;
      const double ptc = cand.pt();
      switch (kConverter.translatePdgIdToType(cand.pdgId())) {
        case reco::PFCandidate::h:
          cSum += ptc;
          cN += 1;
          if (ptc > hardPtMin_) {
            cHSum += ptc;
            cHN += 1;
          }
          if (ptc > cMax)
            cMax = ptc;
          break;
        case reco::PFCandidate::e:
          eS += ptc;
          eNc += 1;
          if (ptc > eM)
            eM = ptc;
          break;
        case reco::PFCandidate::mu:
          mS += ptc;
          muNc += 1;
          if (ptc > mM)
            mM = ptc;
          break;
        case reco::PFCandidate::gamma:
          gSum += ptc;
          gN += 1;
          if (ptc > hardPtMin_) {
            gHSum += ptc;
            gHN += 1;
          }
          if (ptc > gMax)
            gMax = ptc;
          break;
        case reco::PFCandidate::h0:
          nSum += ptc;
          nN += 1;
          if (ptc > nMax)
            nMax = ptc;
          break;
        default:
          break;
      }
    }

    pt.push_back(jet.pt());
    eta.push_back(jet.eta());
    phi.push_back(jet.phi());
    mass.push_back(jet.mass());
    area.push_back(jet.jetArea());
    emEF.push_back(jet.emEnergyFraction());
    hadEF.push_back(jet.energyFractionHadronic());
    trackSum.push_back(tSum);
    trackMax.push_back(tMax);
    trackN.push_back(tN);
    trackHardSum.push_back(tHSum);
    trackHardN.push_back(tHN);
    chargedSum.push_back(cSum);
    chargedMax.push_back(cMax);
    chargedN.push_back(cN);
    chargedHardSum.push_back(cHSum);
    chargedHardN.push_back(cHN);
    photonSum.push_back(gSum);
    photonMax.push_back(gMax);
    photonN.push_back(gN);
    photonHardSum.push_back(gHSum);
    photonHardN.push_back(gHN);
    neutralSum.push_back(nSum);
    neutralMax.push_back(nMax);
    neutralN.push_back(nN);
    eSum.push_back(eS);
    eMax.push_back(eM);
    eN.push_back(eNc);
    muSum.push_back(mS);
    muMax.push_back(mM);
    muN.push_back(muNc);
  }

  auto out = std::make_unique<nanoaod::FlatTable>(pt.size(), name_, /*singleton=*/false, /*extension=*/false);
  out->setDoc(doc_);
  out->addColumn<float>("pt", pt, "calo jet pt", 10);
  out->addColumn<float>("eta", eta, "calo jet eta", 12);
  out->addColumn<float>("phi", phi, "calo jet phi", 12);
  out->addColumn<float>("mass", mass, "calo jet mass", 10);
  out->addColumn<float>("area", area, "jet catchment area", 10);
  out->addColumn<float>("emEF", emEF, "electromagnetic energy fraction", 10);
  out->addColumn<float>("hadEF", hadEF, "hadronic energy fraction (= 1 - emEF)", 10);
  // PF/track in-cone composition (dR<rParam around the calo-jet axis), as HiForest doHiJetID
  out->addColumn<float>("trackSum", trackSum, "sum pT of good tracks in cone", 10);
  out->addColumn<float>("trackMax", trackMax, "max pT track in cone", 10);
  out->addColumn<int>("trackN", trackN, "number of good tracks in cone");
  out->addColumn<float>("trackHardSum", trackHardSum, "sum pT of tracks with pT>hardPtMin in cone", 10);
  out->addColumn<int>("trackHardN", trackHardN, "number of tracks with pT>hardPtMin in cone");
  out->addColumn<float>("chargedSum", chargedSum, "sum pT of charged-hadron PF cands in cone", 10);
  out->addColumn<float>("chargedMax", chargedMax, "max pT charged-hadron PF cand in cone", 10);
  out->addColumn<int>("chargedN", chargedN, "number of charged-hadron PF cands in cone");
  out->addColumn<float>("chargedHardSum", chargedHardSum, "sum pT of charged hadrons with pT>hardPtMin in cone", 10);
  out->addColumn<int>("chargedHardN", chargedHardN, "number of charged hadrons with pT>hardPtMin in cone");
  out->addColumn<float>("photonSum", photonSum, "sum pT of photon PF cands in cone", 10);
  out->addColumn<float>("photonMax", photonMax, "max pT photon PF cand in cone", 10);
  out->addColumn<int>("photonN", photonN, "number of photon PF cands in cone");
  out->addColumn<float>("photonHardSum", photonHardSum, "sum pT of photons with pT>hardPtMin in cone", 10);
  out->addColumn<int>("photonHardN", photonHardN, "number of photons with pT>hardPtMin in cone");
  out->addColumn<float>("neutralSum", neutralSum, "sum pT of neutral-hadron PF cands in cone", 10);
  out->addColumn<float>("neutralMax", neutralMax, "max pT neutral-hadron PF cand in cone", 10);
  out->addColumn<int>("neutralN", neutralN, "number of neutral-hadron PF cands in cone");
  out->addColumn<float>("eSum", eSum, "sum pT of electron PF cands in cone", 10);
  out->addColumn<float>("eMax", eMax, "max pT electron PF cand in cone", 10);
  out->addColumn<int>("eN", eN, "number of electron PF cands in cone");
  out->addColumn<float>("muSum", muSum, "sum pT of muon PF cands in cone", 10);
  out->addColumn<float>("muMax", muMax, "max pT muon PF cand in cone", 10);
  out->addColumn<int>("muN", muN, "number of muon PF cands in cone");
  iEvent.put(std::move(out));
}

void HICaloJetTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedCaloJets"));
  desc.add<edm::InputTag>("pfCandidates", edm::InputTag("packedPFCandidates"));
  desc.add<std::string>("name", "CaloJet");
  desc.add<std::string>("doc", "slimmed calorimeter jets (reco::CaloJet) + HiForest-style PF in-cone ID");
  desc.add<double>("jetPtMin", 15.0);
  desc.add<double>("jetAbsEtaMax", 5.1);
  desc.add<double>("rParam", 0.4);
  desc.add<double>("hardPtMin", 4.0);
  desc.add<bool>("useQuality", true);
  desc.add<std::string>("trackQuality", "highPurity");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(HICaloJetTableProducer);
