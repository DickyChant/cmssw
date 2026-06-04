// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      HiInclusiveJetTableProducer
//
/**\class HiInclusiveJetTableProducer HiInclusiveJetTableProducer.cc

 Heavy-ion inclusive-jet NanoAOD table.

 Ports the HiForest HiInclusiveJetAnalyzer content (the akCs4PF / akCs4FlowPF jet
 trees) to a nanoaod::FlatTable. For each pat::Jet it stores the kinematics + JEC
 (rawPt, pt, eta, phi, y, mass, jetArea, pileup), the PF energy fractions /
 multiplicities, the heavy-ion jet-ID composition variables computed from the
 global packedPFCandidates collection within dR<rParam of the jet axis
 (track/charged/photon/neutral/e/mu Max,Sum,N and the pt>hardPtMin "Hard" sums),
 and the Winner-Take-All reclustered axis (WTAeta/WTAphi).

 Input is the subtracted+btag-updated pat::Jet collection produced by the forest
 HI jet reco (e.g. selectedUpdatedPatJetsAK4PFBtag), so it is fully MiniAOD-doable.
*/
//
// Original Author:  HIN NanoAOD (forest_to_nano), ported from M. Nguyen's HiInclusiveJetAnalyzer
//

#include <memory>
#include <vector>
#include <string>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "fastjet/ClusterSequence.hh"

class HiInclusiveJetTableProducer : public edm::stream::EDProducer<> {
public:
  explicit HiInclusiveJetTableProducer(const edm::ParameterSet&);
  ~HiInclusiveJetTableProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const std::string name_;
  const std::string doc_;
  const int precision_;
  const edm::EDGetTokenT<edm::View<pat::Jet>> jetToken_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfToken_;
  const double jetPtMin_;
  const double jetAbsEtaMax_;
  const double rParam_;
  const double r2Param_;
  const double hardPtMin_;
  const std::string trackQuality_;
  const bool useQuality_;
  const bool doHiJetID_;
  const bool doWTARecluster_;
  const fastjet::JetDefinition wtaJetDef_;
};

HiInclusiveJetTableProducer::HiInclusiveJetTableProducer(const edm::ParameterSet& iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")),
      precision_(iConfig.getParameter<int>("precision")),
      jetToken_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      pfToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("pfCandidates"))),
      jetPtMin_(iConfig.getParameter<double>("jetPtMin")),
      jetAbsEtaMax_(iConfig.getParameter<double>("jetAbsEtaMax")),
      rParam_(iConfig.getParameter<double>("rParam")),
      r2Param_(rParam_ * rParam_),
      hardPtMin_(iConfig.getParameter<double>("hardPtMin")),
      trackQuality_(iConfig.getParameter<std::string>("trackQuality")),
      useQuality_(iConfig.getParameter<bool>("useQuality")),
      doHiJetID_(iConfig.getParameter<bool>("doHiJetID")),
      doWTARecluster_(iConfig.getParameter<bool>("doWTARecluster")),
      wtaJetDef_(fastjet::antikt_algorithm, 2, fastjet::WTA_pt_scheme) {
  produces<nanoaod::FlatTable>();
}

void HiInclusiveJetTableProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  const auto& jets = iEvent.get(jetToken_);
  edm::Handle<edm::View<pat::PackedCandidate>> pfHandle;
  iEvent.getByToken(pfToken_, pfHandle);

  // kinematics
  std::vector<float> rawPt, pt, eta, phi, y, mass, area, pileup;
  // PF fractions / multiplicities
  std::vector<float> pfCHF, pfNHF, pfCEF, pfNEF, pfMUF;
  std::vector<int> pfCHM, pfNHM, pfCEM, pfNEM, pfMUM;
  // HI jet-ID composition
  std::vector<float> trackMax, trackSum, trackHardSum, chargedMax, chargedSum, chargedHardSum;
  std::vector<float> photonMax, photonSum, photonHardSum, neutralMax, neutralSum, eMax, eSum, muMax, muSum;
  std::vector<int> trackN, trackHardN, chargedN, chargedHardN, photonN, photonHardN, neutralN, eN, muN;
  // WTA axis
  std::vector<float> wtaEta, wtaPhi;

  const auto trackQual = reco::TrackBase::qualityByName(trackQuality_);
  reco::PFCandidate pdgConverter;  // for translatePdgIdToType

  for (const auto& jet : jets) {
    if (jet.pt() < jetPtMin_ || std::abs(jet.eta()) > jetAbsEtaMax_)
      continue;

    rawPt.push_back(jet.correctedJet("Uncorrected").pt());
    pt.push_back(jet.pt());
    eta.push_back(jet.eta());
    phi.push_back(jet.phi());
    y.push_back(jet.rapidity());
    pileup.push_back(jet.pileup());
    mass.push_back(jet.mass());
    area.push_back(jet.jetArea());

    if (jet.isPFJet()) {
      pfCHF.push_back(jet.chargedHadronEnergyFraction());
      pfNHF.push_back(jet.neutralHadronEnergyFraction());
      pfCEF.push_back(jet.chargedEmEnergyFraction());
      pfNEF.push_back(jet.neutralEmEnergyFraction());
      pfMUF.push_back(jet.muonEnergyFraction());
      pfCHM.push_back(jet.chargedHadronMultiplicity());
      pfNHM.push_back(jet.neutralHadronMultiplicity());
      pfCEM.push_back(jet.electronMultiplicity());
      pfNEM.push_back(jet.photonMultiplicity());
      pfMUM.push_back(jet.muonMultiplicity());
    } else {
      pfCHF.push_back(0);
      pfNHF.push_back(0);
      pfCEF.push_back(0);
      pfNEF.push_back(0);
      pfMUF.push_back(0);
      pfCHM.push_back(0);
      pfNHM.push_back(0);
      pfCEM.push_back(0);
      pfNEM.push_back(0);
      pfMUM.push_back(0);
    }

    // Heavy-ion jet-ID composition from the global PF candidate collection
    float trkMax = 0, trkSum = 0, trkHardSum = 0;
    int trkN = 0, trkHardN = 0;
    float chMax = 0, chSum = 0, chHardSum = 0, phMax = 0, phSum = 0, phHardSum = 0;
    float neMax = 0, neSum = 0, elMax = 0, elSum = 0, mvMax = 0, mvSum = 0;
    int chN = 0, chHardN = 0, phN = 0, phHardN = 0, neN = 0, elN = 0, mvN = 0;

    if (doHiJetID_ && pfHandle.isValid()) {
      const auto& pfCandidates = *pfHandle;
      // tracks (charged candidates with track details passing quality)
      for (const auto& cand : pfCandidates) {
        if (!cand.hasTrackDetails())
          continue;
        const reco::Track& track = cand.pseudoTrack();
        if (useQuality_ && !track.quality(trackQual))
          continue;
        if (reco::deltaR2(jet, track) < r2Param_) {
          const double ptcand = track.pt();
          trkSum += ptcand;
          trkN += 1;
          if (ptcand > hardPtMin_) {
            trkHardSum += ptcand;
            trkHardN += 1;
          }
          if (ptcand > trkMax)
            trkMax = ptcand;
        }
      }
      // composition by PF type
      for (const auto& cand : pfCandidates) {
        if (reco::deltaR2(jet, cand) >= r2Param_)
          continue;
        const double ptcand = cand.pt();
        switch (pdgConverter.translatePdgIdToType(cand.pdgId())) {
          case 1:  // charged hadron
            chSum += ptcand;
            chN += 1;
            if (ptcand > hardPtMin_) {
              chHardSum += ptcand;
              chHardN += 1;
            }
            if (ptcand > chMax)
              chMax = ptcand;
            break;
          case 2:  // electron
            elSum += ptcand;
            elN += 1;
            if (ptcand > elMax)
              elMax = ptcand;
            break;
          case 3:  // muon
            mvSum += ptcand;
            mvN += 1;
            if (ptcand > mvMax)
              mvMax = ptcand;
            break;
          case 4:  // photon
            phSum += ptcand;
            phN += 1;
            if (ptcand > hardPtMin_) {
              phHardSum += ptcand;
              phHardN += 1;
            }
            if (ptcand > phMax)
              phMax = ptcand;
            break;
          case 5:  // neutral hadron
            neSum += ptcand;
            neN += 1;
            if (ptcand > neMax)
              neMax = ptcand;
            break;
          default:
            break;
        }
      }
    }
    trackMax.push_back(trkMax);
    trackSum.push_back(trkSum);
    trackN.push_back(trkN);
    trackHardSum.push_back(trkHardSum);
    trackHardN.push_back(trkHardN);
    chargedMax.push_back(chMax);
    chargedSum.push_back(chSum);
    chargedN.push_back(chN);
    chargedHardSum.push_back(chHardSum);
    chargedHardN.push_back(chHardN);
    photonMax.push_back(phMax);
    photonSum.push_back(phSum);
    photonN.push_back(phN);
    photonHardSum.push_back(phHardSum);
    photonHardN.push_back(phHardN);
    neutralMax.push_back(neMax);
    neutralSum.push_back(neSum);
    neutralN.push_back(neN);
    eMax.push_back(elMax);
    eSum.push_back(elSum);
    eN.push_back(elN);
    muMax.push_back(mvMax);
    muSum.push_back(mvSum);
    muN.push_back(mvN);

    // Winner-Take-All reclustered axis
    float weta = -999, wphi = -999;
    if (doWTARecluster_) {
      std::vector<fastjet::PseudoJet> candidates;
      for (const auto& d : jet.getJetConstituents())
        candidates.emplace_back(d->px(), d->py(), d->pz(), d->energy());
      if (!candidates.empty()) {
        fastjet::ClusterSequence cs(candidates, wtaJetDef_);
        std::vector<fastjet::PseudoJet> wtajt = fastjet::sorted_by_pt(cs.inclusive_jets(0));
        if (!wtajt.empty()) {
          weta = wtajt[0].eta();
          wphi = wtajt[0].phi_std();
        }
        // release the PseudoJets still tied to the ClusterSequence before it goes
        // out of scope, to avoid fastjet's "out of scope" warning per jet.
        wtajt.clear();
      }
    }
    wtaEta.push_back(weta);
    wtaPhi.push_back(wphi);
  }

  auto out = std::make_unique<nanoaod::FlatTable>(pt.size(), name_, /*singleton=*/false, /*extension=*/false);
  out->setDoc(doc_);
  const int p = precision_;
  out->addColumn<float>("pt", pt, "corrected jet pt", p);
  out->addColumn<float>("rawPt", rawPt, "raw (uncorrected) jet pt", p);
  out->addColumn<float>("eta", eta, "jet eta", p);
  out->addColumn<float>("phi", phi, "jet phi", p);
  out->addColumn<float>("y", y, "jet rapidity", p);
  out->addColumn<float>("mass", mass, "jet mass", p);
  out->addColumn<float>("jetArea", area, "jet catchment area", p);
  out->addColumn<float>("pileup", pileup, "pileup/UE energy subtracted from the jet", p);
  out->addColumn<float>("PfCHF", pfCHF, "charged hadron energy fraction", p);
  out->addColumn<float>("PfNHF", pfNHF, "neutral hadron energy fraction", p);
  out->addColumn<float>("PfCEF", pfCEF, "charged EM energy fraction", p);
  out->addColumn<float>("PfNEF", pfNEF, "neutral EM energy fraction", p);
  out->addColumn<float>("PfMUF", pfMUF, "muon energy fraction", p);
  out->addColumn<int>("PfCHM", pfCHM, "charged hadron multiplicity");
  out->addColumn<int>("PfNHM", pfNHM, "neutral hadron multiplicity");
  out->addColumn<int>("PfCEM", pfCEM, "electron multiplicity");
  out->addColumn<int>("PfNEM", pfNEM, "photon multiplicity");
  out->addColumn<int>("PfMUM", pfMUM, "muon multiplicity");
  out->addColumn<float>("trackMax", trackMax, "max pt track within dR<rParam", p);
  out->addColumn<float>("trackSum", trackSum, "sum pt of tracks within dR<rParam", p);
  out->addColumn<int>("trackN", trackN, "number of tracks within dR<rParam");
  out->addColumn<float>("trackHardSum", trackHardSum, "sum pt of tracks pt>hardPtMin within dR<rParam", p);
  out->addColumn<int>("trackHardN", trackHardN, "number of tracks pt>hardPtMin within dR<rParam");
  out->addColumn<float>("chargedMax", chargedMax, "max pt charged hadron within dR<rParam", p);
  out->addColumn<float>("chargedSum", chargedSum, "sum pt charged hadrons within dR<rParam", p);
  out->addColumn<int>("chargedN", chargedN, "number of charged hadrons within dR<rParam");
  out->addColumn<float>("chargedHardSum", chargedHardSum, "sum pt charged hadrons pt>hardPtMin", p);
  out->addColumn<int>("chargedHardN", chargedHardN, "number of charged hadrons pt>hardPtMin");
  out->addColumn<float>("photonMax", photonMax, "max pt photon within dR<rParam", p);
  out->addColumn<float>("photonSum", photonSum, "sum pt photons within dR<rParam", p);
  out->addColumn<int>("photonN", photonN, "number of photons within dR<rParam");
  out->addColumn<float>("photonHardSum", photonHardSum, "sum pt photons pt>hardPtMin", p);
  out->addColumn<int>("photonHardN", photonHardN, "number of photons pt>hardPtMin");
  out->addColumn<float>("neutralMax", neutralMax, "max pt neutral hadron within dR<rParam", p);
  out->addColumn<float>("neutralSum", neutralSum, "sum pt neutral hadrons within dR<rParam", p);
  out->addColumn<int>("neutralN", neutralN, "number of neutral hadrons within dR<rParam");
  out->addColumn<float>("eMax", eMax, "max pt electron within dR<rParam", p);
  out->addColumn<float>("eSum", eSum, "sum pt electrons within dR<rParam", p);
  out->addColumn<int>("eN", eN, "number of electrons within dR<rParam");
  out->addColumn<float>("muMax", muMax, "max pt muon within dR<rParam", p);
  out->addColumn<float>("muSum", muSum, "sum pt muons within dR<rParam", p);
  out->addColumn<int>("muN", muN, "number of muons within dR<rParam");
  if (doWTARecluster_) {
    out->addColumn<float>("WTAeta", wtaEta, "Winner-Take-All reclustered jet eta", p);
    out->addColumn<float>("WTAphi", wtaPhi, "Winner-Take-All reclustered jet phi", p);
  }
  iEvent.put(std::move(out));
}

void HiInclusiveJetTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("jets", edm::InputTag("selectedUpdatedPatJetsAK4PFBtag"));
  desc.add<edm::InputTag>("pfCandidates", edm::InputTag("packedPFCandidates"));
  desc.add<std::string>("name", "akCs4PFJet");
  desc.add<std::string>("doc", "Heavy-ion inclusive jets (akCs4PF)");
  desc.add<int>("precision", -1);
  desc.add<double>("jetPtMin", 15.0);
  desc.add<double>("jetAbsEtaMax", 2.5);
  desc.add<double>("rParam", 0.4);
  desc.add<double>("hardPtMin", 4.0);
  desc.add<std::string>("trackQuality", "highPurity");
  desc.add<bool>("useQuality", true);
  desc.add<bool>("doHiJetID", true);
  desc.add<bool>("doWTARecluster", true);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(HiInclusiveJetTableProducer);
