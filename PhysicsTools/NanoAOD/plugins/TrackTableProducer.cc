// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      TrackTableProducer
//
/**\class TrackTableProducer TrackTableProducer.cc

 Heavy-ion track + vertex NanoAOD tables.

 Ports the HiForest TrackAnalyzer (HeavyIonsAnalysis/TrackAnalysis) to FlatTables.
 Input is the forest TrackAndVertexUnpacker output `unpackedTracksAndVertices`
 (reco::Tracks + primary vertices rebuilt from packedPFCandidates, plus the
 index-aligned track -> pat::PackedCandidate Ptr vector). Produces:
   - "Trk": per-track kinematics + hits + PF energy + vertex-association DCAs,
   - "Vtx": per-vertex position/errors/chi2/ndof/ptSum.

 The DCA-to-associated-vertex uses the PackedCandidate (matching the forest), while
 the DCA-to-leading-vertex (FirstVtx) uses the reco::Track (also matching the forest).
*/
//
// Original Author:  HIN NanoAOD (forest_to_nano), ported from HiForest TrackAnalyzer
//

#include <cmath>
#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

namespace {
  constexpr float kBad = -999999.f;
}

class TrackTableProducer : public edm::stream::EDProducer<> {
public:
  explicit TrackTableProducer(const edm::ParameterSet&);
  ~TrackTableProducer() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const std::string name_, vtxName_, doc_;
  const int precision_;
  const double trackPtMin_, trackEtaMax_;
  const bool applyTrackSelections_;
  const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<std::vector<edm::Ptr<pat::PackedCandidate>>> track2pcToken_;
};

TrackTableProducer::TrackTableProducer(const edm::ParameterSet& iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      vtxName_(iConfig.getParameter<std::string>("vtxName")),
      doc_(iConfig.getParameter<std::string>("doc")),
      precision_(iConfig.getParameter<int>("precision")),
      trackPtMin_(iConfig.getParameter<double>("trackPtMin")),
      trackEtaMax_(iConfig.getParameter<double>("trackEtaMax")),
      applyTrackSelections_(iConfig.getParameter<bool>("applyTrackSelections")),
      trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackSrc"))),
      vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc"))),
      track2pcToken_(consumes<std::vector<edm::Ptr<pat::PackedCandidate>>>(
          iConfig.getParameter<edm::InputTag>("trackSrc"))) {
  produces<nanoaod::FlatTable>();        // Trk
  produces<nanoaod::FlatTable>("vtx");   // Vtx
}

void TrackTableProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  const auto& tracks = iEvent.get(trackToken_);
  const auto& vertices = iEvent.get(vertexToken_);
  edm::Handle<std::vector<edm::Ptr<pat::PackedCandidate>>> track2pc;
  iEvent.getByToken(track2pcToken_, track2pc);
  const bool hasPC = track2pc.isValid() && track2pc->size() == tracks.size();

  // ---- vertices ----
  std::vector<float> vx, vy, vz, vxe, vye, vze, vchi2, vndof, vptsum;
  std::vector<int> vntrk;
  std::vector<uint8_t> visfake;
  int iMaxPtSumVtx = -1;
  float maxPtSum = -1.f;
  for (const auto& v : vertices) {
    vx.push_back(v.position().x());
    vy.push_back(v.position().y());
    vz.push_back(v.position().z());
    vxe.push_back(v.xError());
    vye.push_back(v.yError());
    vze.push_back(v.zError());
    vchi2.push_back(v.chi2());
    vndof.push_back(v.ndof());
    visfake.push_back(v.isFake() ? 1 : 0);
    vntrk.push_back(v.nTracks());
    float ptsum = 0.f;
    for (auto it = v.tracks_begin(); it != v.tracks_end(); ++it)
      ptsum += (*it)->pt();
    vptsum.push_back(ptsum);
    if (ptsum > maxPtSum) {
      maxPtSum = ptsum;
      iMaxPtSumVtx = static_cast<int>(vx.size()) - 1;
    }
  }

  // ---- tracks ----
  std::vector<float> pt, ptErr, eta, phi, normChi2, pfE, pfEcal, pfHcal;
  std::vector<int> charge, nHits, nPixHits, nLayers, pdgId, assocVtxIdx, assocVtxQual, firstVtxQual;
  std::vector<uint8_t> highPurity;
  std::vector<float> dzAssoc, dzErrAssoc, dxyAssoc, dxyErrAssoc, dzFirst, dzErrFirst, dxyFirst, dxyErrFirst;

  const auto hpBit = reco::TrackBase::qualityByName("highPurity");
  for (size_t it = 0; it < tracks.size(); ++it) {
    const reco::Track& t = tracks[it];
    if (t.pt() < trackPtMin_ || std::abs(t.eta()) > trackEtaMax_)
      continue;
    const bool hp = t.quality(hpBit);
    if (applyTrackSelections_ && (!hp || t.ptError() / t.pt() > 0.1))
      continue;

    pt.push_back(t.pt());
    ptErr.push_back(t.ptError());
    eta.push_back(t.eta());
    phi.push_back(t.phi());
    charge.push_back(t.charge());
    nHits.push_back(t.numberOfValidHits());
    nPixHits.push_back(t.hitPattern().numberOfValidPixelHits());
    nLayers.push_back(t.hitPattern().trackerLayersWithMeasurement());
    normChi2.push_back(t.normalizedChi2());
    highPurity.push_back(hp ? 1 : 0);

    // PF + vertex association from the matched PackedCandidate
    int t_pdg = 0, t_avi = -1, t_avq = -999999, t_fvq = -999999;
    float t_pfE = kBad, t_pfEcal = kBad, t_pfHcal = kBad;
    float t_dzA = kBad, t_dzEA = kBad, t_dxyA = kBad, t_dxyEA = kBad;
    float t_dzF = kBad, t_dzEF = kBad, t_dxyF = kBad, t_dxyEF = kBad;
    if (hasPC) {
      const pat::PackedCandidate& c = *((*track2pc)[it]);
      t_pdg = c.pdgId();
      t_pfE = c.energy();
      t_pfEcal = c.energy() * (c.caloFraction() - c.hcalFraction());
      t_pfHcal = c.energy() * c.hcalFraction();
      const auto v = c.vertexRef();
      if (v.isNonnull() && c.hasTrackDetails()) {
        t_avi = v.key();
        t_avq = c.fromPV(v.key());
        t_dzA = c.dz(v->position());
        t_dzEA = std::sqrt(c.dzError() * c.dzError() + v->zError() * v->zError());
        t_dxyA = c.dxy(v->position());
        // NB: forest uses xError*yError here (reproduced as-is, not a sum of squares)
        t_dxyEA = std::sqrt(c.dxyError() * c.dxyError() + v->xError() * v->yError());
      }
      if (iMaxPtSumVtx >= 0) {
        const math::XYZPoint pv(vx[iMaxPtSumVtx], vy[iMaxPtSumVtx], vz[iMaxPtSumVtx]);
        t_fvq = c.fromPV(iMaxPtSumVtx);
        t_dzF = t.dz(pv);
        t_dzEF = std::sqrt(t.dzError() * t.dzError() + vze[iMaxPtSumVtx] * vze[iMaxPtSumVtx]);
        t_dxyF = t.dxy(pv);
        t_dxyEF = std::sqrt(t.dxyError() * t.dxyError() + vxe[iMaxPtSumVtx] * vye[iMaxPtSumVtx]);
      }
    }
    pdgId.push_back(t_pdg);
    pfE.push_back(t_pfE);
    pfEcal.push_back(t_pfEcal);
    pfHcal.push_back(t_pfHcal);
    assocVtxIdx.push_back(t_avi);
    assocVtxQual.push_back(t_avq);
    dzAssoc.push_back(t_dzA);
    dzErrAssoc.push_back(t_dzEA);
    dxyAssoc.push_back(t_dxyA);
    dxyErrAssoc.push_back(t_dxyEA);
    firstVtxQual.push_back(t_fvq);
    dzFirst.push_back(t_dzF);
    dzErrFirst.push_back(t_dzEF);
    dxyFirst.push_back(t_dxyF);
    dxyErrFirst.push_back(t_dxyEF);
  }

  auto trk = std::make_unique<nanoaod::FlatTable>(pt.size(), name_, false, false);
  trk->setDoc(doc_);
  const int p = precision_;
  trk->addColumn<float>("pt", pt, "track pt", p);
  trk->addColumn<float>("ptError", ptErr, "track pt error", p);
  trk->addColumn<float>("eta", eta, "track eta", p);
  trk->addColumn<float>("phi", phi, "track phi", p);
  trk->addColumn<int>("charge", charge, "track charge");
  trk->addColumn<int>("nHits", nHits, "number of valid hits");
  trk->addColumn<int>("nPixHits", nPixHits, "number of valid pixel hits");
  trk->addColumn<int>("nLayers", nLayers, "tracker layers with measurement");
  trk->addColumn<float>("normChi2", normChi2, "normalized chi2", p);
  trk->addColumn<uint8_t>("highPurity", highPurity, "passes highPurity quality");
  trk->addColumn<int>("pdgId", pdgId, "PF candidate pdgId");
  trk->addColumn<float>("pfEnergy", pfE, "associated PF candidate energy", p);
  trk->addColumn<float>("pfEcal", pfEcal, "associated PF ECAL energy", p);
  trk->addColumn<float>("pfHcal", pfHcal, "associated PF HCAL energy", p);
  trk->addColumn<int>("associatedVtxIndx", assocVtxIdx, "index of associated PV (-1 if none)");
  trk->addColumn<int>("associatedVtxQuality", assocVtxQual, "fromPV quality w.r.t. associated PV");
  trk->addColumn<float>("dzAssociatedVtx", dzAssoc, "dz w.r.t. associated PV (PackedCandidate)", p);
  trk->addColumn<float>("dzErrAssociatedVtx", dzErrAssoc, "dz error w.r.t. associated PV", p);
  trk->addColumn<float>("dxyAssociatedVtx", dxyAssoc, "dxy w.r.t. associated PV (PackedCandidate)", p);
  trk->addColumn<float>("dxyErrAssociatedVtx", dxyErrAssoc, "dxy error w.r.t. associated PV", p);
  trk->addColumn<int>("firstVtxQuality", firstVtxQual, "fromPV quality w.r.t. leading-ptSum PV");
  trk->addColumn<float>("dzFirstVtx", dzFirst, "dz w.r.t. leading-ptSum PV (reco::Track)", p);
  trk->addColumn<float>("dzErrFirstVtx", dzErrFirst, "dz error w.r.t. leading-ptSum PV", p);
  trk->addColumn<float>("dxyFirstVtx", dxyFirst, "dxy w.r.t. leading-ptSum PV (reco::Track)", p);
  trk->addColumn<float>("dxyErrFirstVtx", dxyErrFirst, "dxy error w.r.t. leading-ptSum PV", p);
  iEvent.put(std::move(trk));

  auto vtx = std::make_unique<nanoaod::FlatTable>(vx.size(), vtxName_, false, false);
  vtx->setDoc("Heavy-ion primary vertices (unpacked)");
  vtx->addColumn<float>("x", vx, "vertex x", p);
  vtx->addColumn<float>("y", vy, "vertex y", p);
  vtx->addColumn<float>("z", vz, "vertex z", p);
  vtx->addColumn<float>("xError", vxe, "vertex x error", p);
  vtx->addColumn<float>("yError", vye, "vertex y error", p);
  vtx->addColumn<float>("zError", vze, "vertex z error", p);
  vtx->addColumn<float>("chi2", vchi2, "vertex chi2", p);
  vtx->addColumn<float>("ndof", vndof, "vertex ndof", p);
  vtx->addColumn<int>("nTracks", vntrk, "number of tracks with weight>0.5");
  vtx->addColumn<float>("ptSum", vptsum, "sum pt of vertex tracks", p);
  vtx->addColumn<uint8_t>("isFake", visfake, "vertex isFake");
  iEvent.put(std::move(vtx), "vtx");
}

void TrackTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackSrc", edm::InputTag("unpackedTracksAndVertices"));
  desc.add<edm::InputTag>("vertexSrc", edm::InputTag("unpackedTracksAndVertices"));
  desc.add<std::string>("name", "Trk");
  desc.add<std::string>("vtxName", "Vtx");
  desc.add<std::string>("doc", "Heavy-ion tracks (unpacked from packedPFCandidates)");
  desc.add<int>("precision", -1);
  desc.add<double>("trackPtMin", 0.01);
  desc.add<double>("trackEtaMax", 4.0);
  desc.add<bool>("applyTrackSelections", false);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(TrackTableProducer);
