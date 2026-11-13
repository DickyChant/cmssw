/////////////////////////// HFCandidateBuilder ////////////////////////////////
// Generic N-track (2..5) heavy-flavour hadron finder for NanoAOD: a port of the
// HI "Dfinder" (boundino/Bfinder, Dfinder_14XX_miniAOD) combinatorics into the
// BPHNano composite-candidate style, without any lepton requirement.
//
// One module instance = one decay channel, e.g.
//   D0 -> K- pi+            daughters = [K-, pi+]
//   D+ -> K- pi+ pi+        daughters = [K-, pi+, pi+]
//   Ds+ -> phi(K+ K-) pi+   daughters = [K+, K-, pi+],  resonance = [0,1] (phi)
//   D*+ -> D0(K- pi+) pi+   daughters = [K-, pi+, pi+], resonance = [0,1] (D0)
//   Lc+ -> p K- pi+         daughters = [p, K-, pi+]
// With chargeConjugate=True the charge-conjugate mode (all daughter charges flipped)
// is built as well; userInt("flavour") = +1 (as configured) / -1 (conjugate).
// Identical daughters (same mass and charge, same resonance group) are combined
// once (ordered track indices). The resonance daughters must be the first ones.
//
// Inputs: the BPHTrackMerger track collection + its transient tracks (with
// requireDilepton=False for hadronic use), the beam spot and the primary vertices.
// Output: pat::CompositeCandidateCollection with the vertex-fitted candidates.

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "KinVtxFitter.h"
#include "helper.h"

namespace {
  struct Daughter {
    double mass;
    int charge;  // +1 / -1 for the nominal (flavour = +1) mode
    StringCutObjectSelector<pat::CompositeCandidate> selection;
  };

  struct FitSummary {
    bool ok = false;
    reco::Vertex vtx;
    math::PtEtaPhiMLorentzVector p4;
    float massErr = -1.f, chi2 = -1.f, ndof = -1.f, prob = 0.f;
  };

  // reco::Vertex from a KinVtxFitter result (for VertexDistance3D/XY)
  reco::Vertex toVertex(const KinVtxFitter &f, unsigned nTrk) {
    const GlobalPoint p = f.fitted_vtx();
    return reco::Vertex(reco::Vertex::Point(p.x(), p.y(), p.z()),
                        f.fitted_vtx_uncertainty().matrix(),
                        f.chi2(),
                        f.dof(),
                        nTrk);
  }

  // angle between a flight direction and a momentum
  float angle3D(const reco::Vertex &from, const reco::Vertex &to, const math::PtEtaPhiMLorentzVector &p4) {
    const double dx = to.x() - from.x(), dy = to.y() - from.y(), dz = to.z() - from.z();
    const double d = std::sqrt(dx * dx + dy * dy + dz * dz), p = p4.P();
    if (d <= 0 || p <= 0)
      return -1.f;
    const double c = std::clamp((dx * p4.Px() + dy * p4.Py() + dz * p4.Pz()) / (d * p), -1.0, 1.0);
    return std::acos(c);
  }
  float angle2D(const reco::Vertex &from, const reco::Vertex &to, const math::PtEtaPhiMLorentzVector &p4) {
    const double dx = to.x() - from.x(), dy = to.y() - from.y();
    const double d = std::hypot(dx, dy), p = p4.Pt();
    if (d <= 0 || p <= 0)
      return -1.f;
    const double c = std::clamp((dx * p4.Px() + dy * p4.Py()) / (d * p), -1.0, 1.0);
    return std::acos(c);
  }
}  // namespace

class HFCandidateBuilder : public edm::global::EDProducer<> {
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit HFCandidateBuilder(const edm::ParameterSet &cfg)
      : tracksToken_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
        ttracksToken_(consumes<TransientTrackCollection>(cfg.getParameter<edm::InputTag>("transientTracks"))),
        beamSpotToken_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))),
        pvToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
        chargeConjugate_(cfg.getParameter<bool>("chargeConjugate")),
        pdgId_(cfg.getParameter<int>("pdgId")),
        preMassMin_(cfg.getParameter<double>("preMassMin")),
        preMassMax_(cfg.getParameter<double>("preMassMax")),
        massMin_(cfg.getParameter<double>("massMin")),
        massMax_(cfg.getParameter<double>("massMax")),
        ptMin_(cfg.getParameter<double>("ptMin")),
        absYMax_(cfg.getParameter<double>("absYMax")),
        vtxProbMin_(cfg.getParameter<double>("vtxProbMin")),
        svpvSigMin_(cfg.getParameter<double>("svpvSigMin")),
        alphaMax_(cfg.getParameter<double>("alphaMax")),
        maxDocaMax_(cfg.getParameter<double>("maxDocaMax")),
        resMass_(cfg.getParameter<double>("resMass")),
        resMassWindow_(cfg.getParameter<double>("resMassWindow")),
        resVtxProbMin_(cfg.getParameter<double>("resVtxProbMin")),
        maxTracks_(cfg.getParameter<unsigned>("maxTracks")),
        maxCandidates_(cfg.getParameter<unsigned>("maxCandidates")),
        postSelection_(cfg.getParameter<std::string>("postSelection")) {
    for (const auto &d : cfg.getParameter<std::vector<edm::ParameterSet>>("daughters")) {
      daughters_.push_back(Daughter{d.getParameter<double>("mass"),
                                    d.getParameter<int>("charge"),
                                    StringCutObjectSelector<pat::CompositeCandidate>(d.getParameter<std::string>("selection"))});
    }
    for (auto i : cfg.getParameter<std::vector<unsigned>>("resonance"))
      resIdx_.push_back(i);
    if (daughters_.size() < 2 || daughters_.size() > 5)
      throw cms::Exception("Configuration") << "HFCandidateBuilder: 2 to 5 daughters supported";
    for (unsigned k = 0; k < resIdx_.size(); ++k)
      if (resIdx_[k] != k)
        throw cms::Exception("Configuration") << "HFCandidateBuilder: resonance daughters must be the first ones";
    if (resIdx_.size() == 1 || resIdx_.size() >= daughters_.size())
      throw cms::Exception("Configuration") << "HFCandidateBuilder: resonance must have >= 2 and < N daughters";
    // identical daughters: same mass and charge, both inside or both outside the resonance
    const unsigned n = daughters_.size();
    orderAfter_.assign(n, -1);
    for (unsigned j = 1; j < n; ++j)
      for (int i = int(j) - 1; i >= 0; --i)
        if (daughters_[i].mass == daughters_[j].mass && daughters_[i].charge == daughters_[j].charge &&
            inRes(i) == inRes(j)) {
          orderAfter_[j] = i;
          break;
        }
    produces<pat::CompositeCandidateCollection>();
  }

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracks", edm::InputTag("hfTracks", "SelectedTracks"));
    desc.add<edm::InputTag>("transientTracks", edm::InputTag("hfTracks", "SelectedTransientTracks"));
    desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
    desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
    edm::ParameterSetDescription dau;
    dau.add<double>("mass", 0.13957039);
    dau.add<int>("charge", 1);
    dau.add<std::string>("selection", "");
    desc.addVPSet("daughters", dau, std::vector<edm::ParameterSet>());
    desc.add<std::vector<unsigned>>("resonance", std::vector<unsigned>());
    desc.add<bool>("chargeConjugate", true);
    desc.add<int>("pdgId", 421);
    desc.add<double>("preMassMin", 1.4);
    desc.add<double>("preMassMax", 2.3);
    desc.add<double>("massMin", 1.6);
    desc.add<double>("massMax", 2.1);
    desc.add<double>("ptMin", 0.);
    desc.add<double>("absYMax", 99.);
    desc.add<double>("vtxProbMin", 0.);
    desc.add<double>("svpvSigMin", -1.);
    desc.add<double>("alphaMax", 99.);
    desc.add<double>("maxDocaMax", -1.);
    desc.add<double>("resMass", 0.);
    desc.add<double>("resMassWindow", -1.);
    desc.add<double>("resVtxProbMin", 0.);
    desc.add<unsigned>("maxTracks", 300);
    desc.add<unsigned>("maxCandidates", 2000);
    desc.add<std::string>("postSelection", "");
    descriptions.add("hfCandidateBuilder", desc);
  }

  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;

private:
  bool inRes(unsigned i) const { return std::find(resIdx_.begin(), resIdx_.end(), i) != resIdx_.end(); }

  FitSummary fit(const std::vector<reco::TransientTrack> &tt, const std::vector<double> &masses, KinVtxFitter &f) const;

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> tracksToken_;
  const edm::EDGetTokenT<TransientTrackCollection> ttracksToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  const bool chargeConjugate_;
  const int pdgId_;
  const double preMassMin_, preMassMax_, massMin_, massMax_, ptMin_, absYMax_;
  const double vtxProbMin_, svpvSigMin_, alphaMax_, maxDocaMax_;
  const double resMass_, resMassWindow_, resVtxProbMin_;
  const unsigned maxTracks_, maxCandidates_;
  const StringCutObjectSelector<pat::CompositeCandidate> postSelection_;
  std::vector<Daughter> daughters_;
  std::vector<unsigned> resIdx_;
  std::vector<int> orderAfter_;  // daughter j must use a larger track index than daughter orderAfter_[j]
};

FitSummary HFCandidateBuilder::fit(const std::vector<reco::TransientTrack> &tt,
                                   const std::vector<double> &masses,
                                   KinVtxFitter &f) const {
  FitSummary s;
  std::vector<float> sigmas(masses.size(), bph::PI_SIGMA);
  try {
    f = KinVtxFitter(tt, masses, sigmas);
  } catch (const std::exception &e) {
    return s;
  }
  if (!f.success())
    return s;
  s.ok = true;
  s.vtx = toVertex(f, tt.size());
  s.p4 = f.fitted_p4();
  const auto &err = f.fitted_candidate().kinematicParametersError().matrix();
  s.massErr = err(6, 6) > 0 ? std::sqrt(err(6, 6)) : -1.f;
  s.chi2 = f.chi2();
  s.ndof = f.dof();
  s.prob = f.prob();
  return s;
}

void HFCandidateBuilder::produce(edm::StreamID, edm::Event &evt, const edm::EventSetup &) const {
  auto out = std::make_unique<pat::CompositeCandidateCollection>();

  edm::Handle<pat::CompositeCandidateCollection> tracks;
  evt.getByToken(tracksToken_, tracks);
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracksToken_, ttracks);
  const reco::BeamSpot &bs = evt.get(beamSpotToken_);
  const auto &pvs = evt.get(pvToken_);

  const unsigned nTrk = tracks->size();
  if (nTrk != ttracks->size())
    throw cms::Exception("LogicError") << "HFCandidateBuilder: track / transient-track size mismatch";
  if (nTrk > maxTracks_) {
    edm::LogWarning("HFCandidateBuilder") << nTrk << " tracks > maxTracks=" << maxTracks_ << ", event skipped";
    evt.put(std::move(out));
    return;
  }
  // primary vertex: leading PV, or the beam spot if there is none
  const reco::Vertex pv = pvs.empty() ? reco::Vertex(bs.position(), bs.covariance3D()) : pvs.front();

  const unsigned nD = daughters_.size();
  // per-daughter list of usable tracks (selection only; the charge is checked per flavour)
  std::vector<std::vector<bool>> passSel(nD, std::vector<bool>(nTrk, false));
  for (unsigned d = 0; d < nD; ++d)
    for (unsigned t = 0; t < nTrk; ++t)
      passSel[d][t] = daughters_[d].selection((*tracks)[t]);

  // minimal remaining daughter mass sum after position k (for invariant-mass pruning)
  std::vector<double> massAfter(nD + 1, 0.);
  for (int k = int(nD) - 1; k >= 0; --k)
    massAfter[k] = massAfter[k + 1] + daughters_[k].mass;

  bool capped = false;
  std::vector<int> flavours = {+1};
  if (chargeConjugate_)
    flavours.push_back(-1);

  std::vector<unsigned> idx(nD, 0);
  for (int flav : flavours) {
    // depth-first enumeration of ordered track tuples
    std::vector<math::PtEtaPhiMLorentzVector> partial(nD + 1);
    std::function<void(unsigned)> recurse = [&](unsigned k) {
      if (capped)
        return;
      if (k == nD) {
        // ---- full tuple: pre-fit selection
        const auto &raw = partial[nD];
        if (raw.Pt() < ptMin_ || raw.M() < preMassMin_ || raw.M() > preMassMax_)
          return;
        std::vector<reco::TransientTrack> tt;
        std::vector<double> masses;
        for (unsigned d = 0; d < nD; ++d) {
          tt.push_back((*ttracks)[idx[d]]);
          masses.push_back(daughters_[d].mass);
        }
        // max distance of closest approach between any two daughters
        float maxDoca = -1.f;
        for (unsigned a = 0; a < nD; ++a)
          for (unsigned b = a + 1; b < nD; ++b) {
            TwoTrackMinimumDistance md;
            if (md.calculate(tt[a].initialFreeState(), tt[b].initialFreeState()))
              maxDoca = std::max(maxDoca, float(md.distance()));
          }
        if (maxDocaMax_ > 0 && maxDoca > maxDocaMax_)
          return;
        // ---- vertex fit of all daughters
        KinVtxFitter fitter;
        const FitSummary s = fit(tt, masses, fitter);
        if (!s.ok || s.prob < vtxProbMin_)
          return;
        if (s.p4.M() < massMin_ || s.p4.M() > massMax_)
          return;
        if (std::abs(s.p4.Rapidity()) > absYMax_)
          return;
        const Measurement1D d3 = VertexDistance3D().distance(pv, s.vtx);
        const Measurement1D d2 = VertexDistanceXY().distance(pv, s.vtx);
        const float sig3 = d3.error() > 0 ? d3.value() / d3.error() : -1.f;
        if (svpvSigMin_ > 0 && sig3 < svpvSigMin_)
          return;
        const float alpha = angle3D(pv, s.vtx, s.p4);
        if (alpha > alphaMax_)
          return;

        pat::CompositeCandidate cand;
        cand.setP4(s.p4);
        cand.setCharge(0);
        int q = 0;
        for (unsigned d = 0; d < nD; ++d) {
          const auto &trk = (*tracks)[idx[d]];
          q += trk.charge();
          const std::string pre = "trk" + std::to_string(d + 1) + "_";
          cand.addUserInt(pre + "idx", idx[d]);
          cand.addUserFloat(pre + "massHypo", daughters_[d].mass);
          cand.addUserFloat(pre + "pt", trk.pt());
          cand.addUserFloat(pre + "eta", trk.eta());
          cand.addUserFloat(pre + "phi", trk.phi());
          cand.addUserFloat(pre + "dxyS", trk.hasUserFloat("dxyS") ? trk.userFloat("dxyS") : -999.f);
          cand.addUserFloat(pre + "dEdx", trk.hasUserFloat("dEdx") ? trk.userFloat("dEdx") : -1.f);
        }
        cand.setCharge(q);
        cand.setPdgId(flav * pdgId_);
        cand.setVertex(reco::Candidate::Point(s.vtx.x(), s.vtx.y(), s.vtx.z()));
        cand.addUserInt("flavour", flav);
        cand.addUserFloat("unfitted_mass", raw.M());
        cand.addUserFloat("unfitted_pt", raw.Pt());
        cand.addUserFloat("fitted_mass", s.p4.M());
        cand.addUserFloat("fitted_massErr", s.massErr);
        cand.addUserFloat("fitted_pt", s.p4.Pt());
        cand.addUserFloat("fitted_eta", s.p4.Eta());
        cand.addUserFloat("fitted_phi", s.p4.Phi());
        cand.addUserFloat("fitted_y", s.p4.Rapidity());
        cand.addUserFloat("sv_chi2", s.chi2);
        cand.addUserFloat("sv_ndof", s.ndof);
        cand.addUserFloat("sv_prob", s.prob);
        cand.addUserFloat("vtx_x", s.vtx.x());
        cand.addUserFloat("vtx_y", s.vtx.y());
        cand.addUserFloat("vtx_z", s.vtx.z());
        cand.addUserFloat("vtx_cxx", s.vtx.covariance(0, 0));
        cand.addUserFloat("vtx_cyy", s.vtx.covariance(1, 1));
        cand.addUserFloat("vtx_czz", s.vtx.covariance(2, 2));
        cand.addUserFloat("vtx_cyx", s.vtx.covariance(1, 0));
        cand.addUserFloat("vtx_czx", s.vtx.covariance(2, 0));
        cand.addUserFloat("vtx_czy", s.vtx.covariance(2, 1));
        // Dfinder: svpvDistance(_2D), svpvDisErr(_2D), alpha, ip3d, MaxDoca
        cand.addUserFloat("svpvDistance", d3.value());
        cand.addUserFloat("svpvDisErr", d3.error());
        cand.addUserFloat("svpvDistance_2D", d2.value());
        cand.addUserFloat("svpvDisErr_2D", d2.error());
        cand.addUserFloat("alpha", alpha);
        cand.addUserFloat("alpha_2D", angle2D(pv, s.vtx, s.p4));
        cand.addUserFloat("cos_theta_2D_bs", bph::cos_theta_2D(fitter, bs, s.p4));
        float ip3d = -999.f, ip3dErr = -1.f;
        try {
          const auto ip = IPTools::absoluteImpactParameter3D(fitter.fitted_candidate_ttrk(), pv);
          if (ip.first) {
            ip3d = ip.second.value();
            ip3dErr = ip.second.error();
          }
        } catch (const std::exception &) {
        }
        cand.addUserFloat("ip3d", ip3d);
        cand.addUserFloat("ip3dErr", ip3dErr);
        cand.addUserFloat("maxDoca", maxDoca);

        // ---- resonance (first nRes daughters): own vertex fit, Dfinder tktkRes_* block
        if (!resIdx_.empty()) {
          std::vector<reco::TransientTrack> rtt(tt.begin(), tt.begin() + resIdx_.size());
          std::vector<double> rm(masses.begin(), masses.begin() + resIdx_.size());
          KinVtxFitter rfit;
          const FitSummary r = fit(rtt, rm, rfit);
          if (!r.ok || r.prob < resVtxProbMin_)
            return;
          if (resMassWindow_ > 0 && std::abs(r.p4.M() - resMass_) > resMassWindow_)
            return;
          const Measurement1D rd = VertexDistance3D().distance(pv, r.vtx);
          const Measurement1D rdSV = VertexDistance3D().distance(s.vtx, r.vtx);
          cand.addUserFloat("res_unfitted_mass", partial[resIdx_.size()].M());
          cand.addUserFloat("res_mass", r.p4.M());
          cand.addUserFloat("res_massErr", r.massErr);
          cand.addUserFloat("res_pt", r.p4.Pt());
          cand.addUserFloat("res_eta", r.p4.Eta());
          cand.addUserFloat("res_phi", r.p4.Phi());
          cand.addUserFloat("res_sv_prob", r.prob);
          cand.addUserFloat("res_vtx_x", r.vtx.x());
          cand.addUserFloat("res_vtx_y", r.vtx.y());
          cand.addUserFloat("res_vtx_z", r.vtx.z());
          cand.addUserFloat("res_svpvDistance", rd.value());
          cand.addUserFloat("res_svpvDisErr", rd.error());
          cand.addUserFloat("res_svpvDistanceToSV", rdSV.value());
          cand.addUserFloat("res_svpvDisErrToSV", rdSV.error());
          cand.addUserFloat("res_alpha", angle3D(pv, r.vtx, r.p4));
          cand.addUserFloat("res_alphaToSV", angle3D(s.vtx, r.vtx, r.p4));
          cand.addUserFloat("deltaM", s.p4.M() - r.p4.M());
          cand.addUserFloat("unfitted_deltaM", raw.M() - partial[resIdx_.size()].M());
        }
        if (!postSelection_(cand))
          return;
        if (out->size() >= maxCandidates_) {
          capped = true;
          return;
        }
        out->push_back(cand);
        return;
      }
      // ---- choose the track for daughter k
      const int q = flav * daughters_[k].charge;
      const unsigned tStart = orderAfter_[k] >= 0 ? idx[orderAfter_[k]] + 1 : 0;
      for (unsigned t = tStart; t < nTrk; ++t) {
        if ((*tracks)[t].charge() != q || !passSel[k][t])
          continue;
        bool used = false;
        for (unsigned j = 0; j < k; ++j)
          used |= (idx[j] == t);
        if (used)
          continue;
        const auto &trk = (*tracks)[t];
        math::PtEtaPhiMLorentzVector p4(trk.pt(), trk.eta(), trk.phi(), daughters_[k].mass);
        partial[k + 1] = (k == 0 ? p4 : partial[k] + p4);
        // invariant-mass pruning: the full mass is at least the partial mass plus the remaining masses
        if (partial[k + 1].M() + massAfter[k + 1] > preMassMax_)
          continue;
        // resonance pre-window on the raw resonance mass once it is complete
        if (!resIdx_.empty() && k + 1 == resIdx_.size() && resMassWindow_ > 0 &&
            std::abs(partial[k + 1].M() - resMass_) > 2. * resMassWindow_)
          continue;
        idx[k] = t;
        recurse(k + 1);
      }
    };
    recurse(0);
  }
  if (capped)
    edm::LogWarning("HFCandidateBuilder") << "candidate cap maxCandidates=" << maxCandidates_ << " reached";
  evt.put(std::move(out));
}

DEFINE_FWK_MODULE(HFCandidateBuilder);
