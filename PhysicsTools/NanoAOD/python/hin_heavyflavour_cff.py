"""HI heavy-flavour finders for the HIN NanoAOD flavours (HINUPC).

NanoAOD port of the HI forest add-on finders, built on the BPHNano composite-candidate machinery:
  * D finder  (boundino/Bfinder "Dfinder"): HFCandidateBuilder, one instance per channel, no lepton requirement
                D0->K-pi+, D*+->D0(K-pi+)pi+, D+->K-pi+pi+, Ds+->phi(K+K-)pi+, Lc+->pK-pi+ (+ charge conjugates)
  * onia      (HiOnia / HiOnia2MuMuPAT):  low-pT HI muon selection + dimuon vertex fit (OS and SS)
  * B finder  (boundino/Bfinder "Bfinder"): BPHNano B->J/psi K and B->J/psi hh (K*0, phi) builders on the HI inputs
All modules use their own labels (hf*/hiOnia*), so they do not interfere with the central @BPH flavour.

Tables: HFTrack, D0, Dstar, Dplus, Ds, LambdaC, OniaMuon, Onia, BToJpsiK, BToJpsiTrkTrk (+ HFGenPart and genPartIdx in MC).
Track indices (trkN_idx / k_idx / trk1_idx, ...) point to HFTrack rows; muon indices (l1_idx, l2_idx) to OniaMuon rows.
"""
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.NanoAOD.simplePATMuonFlatTableProducer_cfi import simplePATMuonFlatTableProducer
from PhysicsTools.NanoAOD.simpleGenParticleFlatTableProducer_cfi import simpleGenParticleFlatTableProducer
from PhysicsTools.NanoAOD.genparticles_cff import finalGenParticles, genParticleTable

PI_MASS, K_MASS, P_MASS = 0.13957039, 0.493677, 0.938272
D0_MASS, DPLUS_MASS, DS_MASS, DSTAR_MASS, LC_MASS, PHI_MASS = 1.86484, 1.86966, 1.96835, 2.01026, 2.28646, 1.019461


def _ufloat(expr, doc, precision=-1):
    return Var("userFloat('%s')" % expr, float, doc=doc, precision=precision)


def _uint(expr, doc):
    return Var("userInt('%s')" % expr, int, doc=doc)


# ---------------------------------------------------------------------------
#  Tracks: BPHTrackMerger without the dilepton requirement (+ dE/dx, full helix covariance)
# ---------------------------------------------------------------------------
hfTracks = cms.EDProducer("BPHTrackMerger",
    beamSpot        = cms.InputTag("offlineBeamSpot"),
    dileptons       = cms.InputTag("hiOniaMuMu", "SelectedDiLeptons"),   # only used for dzTrg
    tracks          = cms.InputTag("packedPFCandidates"),
    lostTracks      = cms.InputTag("lostTracks"),
    trackSelection  = cms.string("pt > 0.3 && abs(eta) < 2.4"),
    muons           = cms.InputTag("slimmedMuons"),
    electrons       = cms.InputTag("slimmedElectrons"),
    pvSrc           = cms.InputTag("offlineSlimmedPrimaryVertices"),
    maxDzDilep      = cms.double(-1.),
    dcaSig          = cms.double(-100000.),
    requireDilepton = cms.bool(False),
    dedx            = cms.InputTag("dedxEstimator", "dedxAllLikelihood"),
)

_covNames = ["QopQop", "LamLam", "PhiPhi", "QopLam", "QopPhi", "LamPhi", "QopDxy", "QopDsz", "LamDxy", "LamDsz",
             "PhiDxy", "PhiDsz", "DxyDxy", "DxyDsz", "DszDsz"]
hfTrackTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer",
    src  = cms.InputTag("hfTracks", "SelectedTracks"),
    cut  = cms.string(""),
    name = cms.string("HFTrack"),
    doc  = cms.string("charged tracks (packed PF candidates with track details + lost tracks, |pdgId|=211, pt>0.3) "
                      "used by the HF finders; no lepton requirement"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        CandVars,
        vx = Var("vx()", float, doc="x of the reference point [cm]", precision=10),
        vy = Var("vy()", float, doc="y of the reference point [cm]", precision=10),
        vz = Var("vz()", float, doc="z of the reference point [cm]", precision=10),
        dxy = _ufloat("dxy", "dxy wrt PV[0] [cm]", 10),
        dxyS = _ufloat("dxyS", "dxy/err wrt PV[0]", 10),
        dz = _ufloat("dz", "dz wrt PV[0] [cm]", 10),
        dzS = _ufloat("dzS", "dz/err", 10),
        DCASig = _ufloat("DCASig", "transverse DCA significance wrt the beam spot", 10),
        ptErr = _ufloat("ptErr", "pt uncertainty", 10),
        normChi2 = _ufloat("normChi2", "track normalized chi2", 10),
        dEdx = _ufloat("dEdx", "dE/dx (dedxAllLikelihood) [MeV/cm], -1 if not available", 10),
        dEdxNMeas = _uint("dEdxNMeas", "number of dE/dx measurements, -1 if not available"),
        nValidHits = _uint("nValidHits", "number of valid hits"),
        nValidPixelHits = _uint("nValidPixelHits", "number of valid pixel hits"),
        highPurity = _uint("highPurity", "highPurity track"),
        isLostTrk = _uint("isLostTrk", "1 if from lostTracks"),
        keyPacked = _uint("keyPacked", "index in packedPFCandidates (= PFCand row), or nPacked + index in lostTracks (LostTrack row)"),
        isMatchedToMuon = _uint("isMatchedToMuon", "track is a muon's source candidate"),
        isMatchedToEle = _uint("isMatchedToEle", "track is an electron's source candidate"),
        **{"cov" + n: _ufloat("cov" + n, "helix covariance " + n, 10) for n in _covNames}
    ),
)


# ---------------------------------------------------------------------------
#  D finder channels
# ---------------------------------------------------------------------------
def _dau(mass, charge, selection=""):
    return cms.PSet(mass=cms.double(mass), charge=cms.int32(charge), selection=cms.string(selection))


_trkSel = "pt > 0.5"
_hfBase = dict(
    tracks          = cms.InputTag("hfTracks", "SelectedTracks"),
    transientTracks = cms.InputTag("hfTracks", "SelectedTransientTracks"),
    beamSpot        = cms.InputTag("offlineBeamSpot"),
    vertices        = cms.InputTag("offlineSlimmedPrimaryVertices"),
    chargeConjugate = cms.bool(True),
    ptMin           = cms.double(0.),
    absYMax         = cms.double(99.),
    vtxProbMin      = cms.double(0.01),    # keeps 98% of matched D0 (UPC non-prompt D0 MC)
    svpvSigMin      = cms.double(-1.),
    alphaMax        = cms.double(99.),
    maxDocaMax      = cms.double(-1.),
    resMass         = cms.double(0.),
    resMassWindow   = cms.double(-1.),
    resVtxProbMin   = cms.double(0.),
    resonance       = cms.vuint32(),
    maxTracks       = cms.uint32(300),
    maxCandidates   = cms.uint32(2000),
    postSelection   = cms.string(""),
)


def _hfBuilder(**kw):
    p = dict(_hfBase)
    p.update(kw)
    return cms.EDProducer("HFCandidateBuilder", **p)


# window conventions follow Dfinder: D0 +-0.4 GeV (pre-fit), D+/Ds/Lc +-0.2 GeV, D* +-0.2 GeV
hfD0 = _hfBuilder(pdgId=cms.int32(421),
    daughters=cms.VPSet(_dau(K_MASS, -1, _trkSel), _dau(PI_MASS, +1, _trkSel)),
    preMassMin=cms.double(D0_MASS - 0.4), preMassMax=cms.double(D0_MASS + 0.4),
    massMin=cms.double(D0_MASS - 0.3), massMax=cms.double(D0_MASS + 0.3))
hfDstar = _hfBuilder(pdgId=cms.int32(413),
    daughters=cms.VPSet(_dau(K_MASS, -1, _trkSel), _dau(PI_MASS, +1, _trkSel), _dau(PI_MASS, +1)),   # soft pion: any HFTrack
    resonance=cms.vuint32(0, 1), resMass=cms.double(D0_MASS), resMassWindow=cms.double(0.1),
    preMassMin=cms.double(DSTAR_MASS - 0.3), preMassMax=cms.double(DSTAR_MASS + 0.3),
    massMin=cms.double(DSTAR_MASS - 0.2), massMax=cms.double(DSTAR_MASS + 0.2),
    postSelection=cms.string("userFloat('deltaM') < 0.2"))
hfDplus = _hfBuilder(pdgId=cms.int32(411),
    daughters=cms.VPSet(_dau(K_MASS, -1, _trkSel), _dau(PI_MASS, +1, _trkSel), _dau(PI_MASS, +1, _trkSel)),
    preMassMin=cms.double(DPLUS_MASS - 0.25), preMassMax=cms.double(DPLUS_MASS + 0.25),
    massMin=cms.double(DPLUS_MASS - 0.2), massMax=cms.double(DPLUS_MASS + 0.2))
hfDs = _hfBuilder(pdgId=cms.int32(431),
    daughters=cms.VPSet(_dau(K_MASS, +1, _trkSel), _dau(K_MASS, -1, _trkSel), _dau(PI_MASS, +1, _trkSel)),
    resonance=cms.vuint32(0, 1), resMass=cms.double(PHI_MASS), resMassWindow=cms.double(0.05),
    preMassMin=cms.double(DS_MASS - 0.35), preMassMax=cms.double(DS_MASS + 0.25),   # covers D+ -> phi pi too
    massMin=cms.double(DS_MASS - 0.3), massMax=cms.double(DS_MASS + 0.2))
hfLambdaC = _hfBuilder(pdgId=cms.int32(4122),
    daughters=cms.VPSet(_dau(P_MASS, +1, _trkSel), _dau(K_MASS, -1, _trkSel), _dau(PI_MASS, +1, _trkSel)),
    preMassMin=cms.double(LC_MASS - 0.25), preMassMax=cms.double(LC_MASS + 0.25),
    massMin=cms.double(LC_MASS - 0.2), massMax=cms.double(LC_MASS + 0.2))


def _hfTable(src, name, doc, nDau, hasRes=False, resName="res"):
    v = cms.PSet(
        CandVars,
        flavour = _uint("flavour", "+1: decay as named, -1: charge conjugate"),
        unfitted_mass = _ufloat("unfitted_mass", "pre-fit mass"),
        unfitted_pt = _ufloat("unfitted_pt", "pre-fit pt"),
        massErr = _ufloat("fitted_massErr", "mass uncertainty from the vertex fit"),
        y = _ufloat("fitted_y", "rapidity"),
        svprob = _ufloat("sv_prob", "vertex-fit probability"),
        svchi2 = _ufloat("sv_chi2", "vertex-fit chi2"),
        svndof = _ufloat("sv_ndof", "vertex-fit ndof"),
        vtx_x = _ufloat("vtx_x", "SV x [cm]", 10), vtx_y = _ufloat("vtx_y", "SV y [cm]", 10), vtx_z = _ufloat("vtx_z", "SV z [cm]", 10),
        vtx_cxx = _ufloat("vtx_cxx", "SV cov xx", 10), vtx_cyy = _ufloat("vtx_cyy", "SV cov yy", 10), vtx_czz = _ufloat("vtx_czz", "SV cov zz", 10),
        vtx_cyx = _ufloat("vtx_cyx", "SV cov yx", 10), vtx_czx = _ufloat("vtx_czx", "SV cov zx", 10), vtx_czy = _ufloat("vtx_czy", "SV cov zy", 10),
        svpvDistance = _ufloat("svpvDistance", "3D SV-PV distance [cm]"),
        svpvDisErr = _ufloat("svpvDisErr", "its uncertainty (SV and PV errors)"),
        svpvDistance_2D = _ufloat("svpvDistance_2D", "transverse SV-PV distance [cm]"),
        svpvDisErr_2D = _ufloat("svpvDisErr_2D", "its uncertainty"),
        alpha = _ufloat("alpha", "3D angle between SV-PV and the candidate momentum [rad]"),
        alpha_2D = _ufloat("alpha_2D", "transverse pointing angle [rad]"),
        cosTheta2D_bs = _ufloat("cos_theta_2D_bs", "cos of the transverse pointing angle wrt the beam spot"),
        ip3d = _ufloat("ip3d", "3D impact parameter of the candidate wrt PV[0] [cm]"),
        ip3dErr = _ufloat("ip3dErr", "its uncertainty"),
        maxDoca = _ufloat("maxDoca", "max distance of closest approach between two daughters [cm]"),
    )
    for i in range(1, nDau + 1):
        p = "trk%d_" % i
        setattr(v, p + "idx", _uint(p + "idx", "daughter %d: index in HFTrack" % i))
        setattr(v, p + "massHypo", _ufloat(p + "massHypo", "daughter %d mass hypothesis" % i))
        setattr(v, p + "pt", _ufloat(p + "pt", "daughter %d track pt" % i))
        setattr(v, p + "eta", _ufloat(p + "eta", "daughter %d track eta" % i))
        setattr(v, p + "phi", _ufloat(p + "phi", "daughter %d track phi" % i))
        setattr(v, p + "dxyS", _ufloat(p + "dxyS", "daughter %d dxy significance wrt PV" % i))
        setattr(v, p + "dEdx", _ufloat(p + "dEdx", "daughter %d dE/dx" % i))
    if hasRes:
        for n, d in [("unfitted_mass", "pre-fit mass"), ("mass", "fitted mass"), ("massErr", "mass uncertainty"),
                     ("pt", "pt"), ("eta", "eta"), ("phi", "phi"), ("sv_prob", "vertex probability"),
                     ("vtx_x", "SV x"), ("vtx_y", "SV y"), ("vtx_z", "SV z"),
                     ("svpvDistance", "3D distance to PV"), ("svpvDisErr", "its uncertainty"),
                     ("svpvDistanceToSV", "3D distance to the mother SV"), ("svpvDisErrToSV", "its uncertainty"),
                     ("alpha", "pointing angle wrt PV"), ("alphaToSV", "pointing angle wrt the mother SV")]:
            setattr(v, resName + "_" + n, _ufloat("res_" + n, resName + " (resonance) " + d))
        v.deltaM = _ufloat("deltaM", "fitted mass - fitted resonance mass")
        v.unfitted_deltaM = _ufloat("unfitted_deltaM", "pre-fit mass difference")
    return cms.EDProducer("SimpleCompositeCandidateFlatTableProducer",
        src=cms.InputTag(src), cut=cms.string(""), name=cms.string(name), doc=cms.string(doc),
        singleton=cms.bool(False), extension=cms.bool(False), variables=v)


hfD0Table = _hfTable("hfD0", "D0", "D0 -> K- pi+ (+cc) candidates (Dfinder channels 0/1)", 2)
hfDstarTable = _hfTable("hfDstar", "Dstar", "D*+ -> D0(K- pi+) pi+ (+cc) candidates (Dfinder channels 8/9)", 3, True, "D0")
hfDplusTable = _hfTable("hfDplus", "Dplus", "D+ -> K- pi+ pi+ (+cc) candidates (Dfinder channels 2/3)", 3)
hfDsTable = _hfTable("hfDs", "Ds", "Ds+ (and D+) -> phi(K+ K-) pi+ (+cc) candidates (Dfinder channels 6/7)", 3, True, "phi")
hfLambdaCTable = _hfTable("hfLambdaC", "LambdaC", "Lc+ -> p K- pi+ (+cc) candidates (Dfinder channels 14/15)", 3)

# ---------------------------------------------------------------------------
#  Onia (HiOnia-like): low-pT HI muons + dimuon vertex fit, OS and SS
# ---------------------------------------------------------------------------
hiOniaMuons = cms.EDProducer("MuonTriggerSelector",
    muonCollection = cms.InputTag("slimmedMuons"),
    bits           = cms.InputTag("TriggerResults", "", "HLT"),
    prescales      = cms.InputTag("patTrigger"),
    objects        = cms.InputTag("slimmedPatTrigger"),
    maxdR_matching = cms.double(0.3),
    # HiOnia2MuMuPAT lowerPuritySelection, without the gen-match branch (the flags are stored instead)
    muonSelection  = cms.string("(isTrackerMuon || isGlobalMuon) && innerTrack.isNonnull && abs(eta) < 2.4"
                                " && abs(innerTrack.dxy) < 4 && abs(innerTrack.dz) < 35"),
    HLTPaths       = cms.vstring(),
)
hiOniaMuMu = cms.EDProducer("DiMuonBuilder",
    src                = cms.InputTag("hiOniaMuons", "SelectedMuons"),
    transientTracksSrc = cms.InputTag("hiOniaMuons", "SelectedTransientMuons"),
    lep1Selection      = cms.string(""),
    lep2Selection      = cms.string(""),
    beamSpot           = cms.InputTag("offlineBeamSpot"),
    preVtxSelection    = cms.string('mass() < 15.0 && userFloat("lep_deltaR") > 0.001'),   # OS and SS (HiOnia keeps both)
    postVtxSelection   = cms.string('userFloat("sv_prob") > 0.'),
)
hiOniaMuonTable = simplePATMuonFlatTableProducer.clone(
    src  = cms.InputTag("hiOniaMuons", "SelectedMuons"),
    cut  = cms.string(""),
    name = cms.string("OniaMuon"),
    doc  = cms.string("HI onia muons (tracker or global, |eta|<2.4, no pt cut; HiOnia2MuMuPAT-like)"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        CandVars,
        isGlobal = Var("isGlobalMuon", bool, doc="global muon"),
        isTracker = Var("isTrackerMuon", bool, doc="tracker muon"),
        isPF = Var("isPFMuon", bool, doc="PF muon"),
        softId = Var("passed('SoftCutBasedId')", bool, doc="soft cut-based ID"),
        looseId = Var("passed('CutBasedIdLoose')", bool, doc="loose cut-based ID"),
        TMOneStationTight = Var("muonID('TMOneStationTight')", bool, doc="TMOneStationTight"),
        nTrkLayers = Var("?innerTrack.isNonnull?innerTrack().hitPattern().trackerLayersWithMeasurement():-1", int, doc="tracker layers with measurement"),
        nPixLayers = Var("?innerTrack.isNonnull?innerTrack().hitPattern().pixelLayersWithMeasurement():-1", int, doc="pixel layers with measurement"),
        nPixHits = Var("?innerTrack.isNonnull?innerTrack().hitPattern().numberOfValidPixelHits():-1", int, doc="valid pixel hits"),
        dxy = Var("dB('PV2D')", float, doc="dxy wrt PV [cm]", precision=10),
        dz = Var("dB('PVDZ')", float, doc="dz wrt PV [cm]", precision=10),
        innerDxy = Var("?innerTrack.isNonnull?innerTrack().dxy():-999", float, doc="inner-track dxy wrt (0,0) [cm] (HiOnia cut)", precision=10),
        innerDz = Var("?innerTrack.isNonnull?innerTrack().dz():-999", float, doc="inner-track dz wrt (0,0) [cm] (HiOnia cut)", precision=10),
        globalNormChi2 = Var("?isGlobalMuon?globalTrack().normalizedChi2():-1", float, doc="global-track normalized chi2", precision=10),
        nMuonValidHits = Var("?isGlobalMuon?globalTrack().hitPattern().numberOfValidMuonHits():-1", int, doc="valid muon hits (global track)"),
        nMatchedStations = Var("numberOfMatchedStations()", int, doc="matched stations"),
        isTriggering = Var("userInt('isTriggering')", int, doc="matched to one of the configured HLT paths"),
    ),
)
hiOniaTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer",
    src  = cms.InputTag("hiOniaMuMu", "SelectedDiLeptons"),
    cut  = cms.string(""),
    name = cms.string("Onia"),
    doc  = cms.string("HI onia dimuon candidates (OS and SS, m < 15 GeV), HiOnia-like"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        CandVars,
        y = Var("rapidity()", float, doc="rapidity (pre-fit)"),
        fitted_mass = _ufloat("fitted_mass", "vertex-fitted mass"),
        svprob = _ufloat("sv_prob", "vertex probability"),
        svchi2 = _ufloat("sv_chi2", "vertex chi2"),
        svndof = _ufloat("sv_ndof", "vertex ndof"),
        vtx_x = _ufloat("vtx_x", "SV x [cm]", 10), vtx_y = _ufloat("vtx_y", "SV y [cm]", 10), vtx_z = _ufloat("vtx_z", "SV z [cm]", 10),
        l_xy = _ufloat("l_xy", "transverse decay length wrt the beam spot [cm]"),
        l_xy_unc = _ufloat("l_xy_unc", "its uncertainty"),
        cosTheta2D = _ufloat("cos_theta_2D", "cos of the transverse pointing angle (beam spot)"),
        ctau2D = Var("userFloat('l_xy') * mass() / pt()", float, doc="pseudo-proper decay length l_xy*m/pt [cm]"),
        lep_deltaR = _ufloat("lep_deltaR", "deltaR between the muons"),
        l1_idx = _uint("l1_idx", "muon 1 index in OniaMuon"),
        l2_idx = _uint("l2_idx", "muon 2 index in OniaMuon"),
    ),
)

# ---------------------------------------------------------------------------
#  B finder on the HI onia + HF tracks (BPHNano builders, Bfinder-like loose cuts)
# ---------------------------------------------------------------------------
hfBToJpsiK = cms.EDProducer("BToTrkLLBuilder",
    dileptons             = cms.InputTag("hiOniaMuMu", "SelectedDiLeptons"),
    leptonTransientTracks = cms.InputTag("hiOniaMuons", "SelectedTransientMuons"),
    dileptonKinVtxs       = cms.InputTag("hiOniaMuMu", "SelectedDiLeptonKinVtxs"),
    kaons                 = cms.InputTag("hfTracks", "SelectedTracks"),
    kaonsTransientTracks  = cms.InputTag("hfTracks", "SelectedTransientTracks"),
    trackMass             = cms.double(493.677),
    beamSpot              = cms.InputTag("offlineBeamSpot"),
    PUtracks              = cms.InputTag("hfTracks", "SelectedTracks"),
    preVtxSelection       = cms.string('4.5 < mass && mass < 6.5 && userFloat("min_dr") > 0.01'),
    postVtxSelection      = cms.string('4.7 < userFloat("fitted_mass") && userFloat("fitted_mass") < 6.0 && userFloat("sv_prob") > 0.001'),
    dileptonMassContraint = cms.bool(True),
)
hfDiTrack = cms.EDProducer("DiTrackBuilder",
    tracks          = cms.InputTag("hfTracks", "SelectedTracks"),
    transientTracks = cms.InputTag("hfTracks", "SelectedTransientTracks"),
    trk1Selection   = cms.string("pt > 0.5"),
    trk2Selection   = cms.string("pt > 0.5"),
    trk1Mass        = cms.double(0.139),
    trk2Mass        = cms.double(0.494),
    beamSpot        = cms.InputTag("offlineBeamSpot"),
    preVtxSelection = cms.string('((userFloat("unfitted_mass_KK")>0.95 && userFloat("unfitted_mass_KK")<1.12) || (userFloat("unfitted_mass_Kpi")>0.6 && userFloat("unfitted_mass_Kpi")<1.2) || (userFloat("unfitted_mass_piK")>0.6 && userFloat("unfitted_mass_piK")<1.2)) && charge() == 0'),
    postVtxSelection = cms.string('((userFloat("fitted_mass_KK")>0.95 && userFloat("fitted_mass_KK")<1.12) || (userFloat("fitted_mass_Kpi")>0.6 && userFloat("fitted_mass_Kpi")<1.2) || (userFloat("fitted_mass_piK")>0.6 && userFloat("fitted_mass_piK")<1.2)) && userFloat("sv_prob") > 0.001'),
)
hfBToJpsiTrkTrk = cms.EDProducer("BToTrkTrkLLBuilder",
    dileptons             = cms.InputTag("hiOniaMuMu", "SelectedDiLeptons"),
    leptonTransientTracks = cms.InputTag("hiOniaMuons", "SelectedTransientMuons"),
    ditracks              = cms.InputTag("hfDiTrack"),
    transientTracks       = cms.InputTag("hfTracks", "SelectedTransientTracks"),
    PUtracks              = cms.InputTag("hfTracks", "SelectedTracks"),
    beamSpot              = cms.InputTag("offlineBeamSpot"),
    preVtxSelection       = cms.string('userFloat("min_dr") > 0.01 && ((4.5<userFloat("unfitted_B_mass_KK") && userFloat("unfitted_B_mass_KK")<6.5) || (4.5<userFloat("unfitted_B_mass_Kpi") && userFloat("unfitted_B_mass_Kpi")<6.5) || (4.5<userFloat("unfitted_B_mass_piK") && userFloat("unfitted_B_mass_piK")<6.5))'),
    postVtxSelection      = cms.string('userFloat("sv_prob") > 1.e-3 && ((4.6<userFloat("fitted_mass_KK") && userFloat("fitted_mass_KK")<6.0) || (4.6<userFloat("fitted_mass_Kpi") && userFloat("fitted_mass_Kpi")<6.0) || (4.6<userFloat("fitted_mass_piK") && userFloat("fitted_mass_piK")<6.0))'),
    dileptonMassContraint = cms.bool(True),
)


def _bTables():
    from PhysicsTools.BPHNano.BToKLL_cff import BToKmumuTable
    from PhysicsTools.BPHNano.BToTrkTrkLL_cff import BToTrkTrkMuMuTable
    from PhysicsTools.BPHNano.DiTrack_cff import DiTrack  # noqa: F401 (ensures the package python is importable)
    t1 = BToKmumuTable.clone(src=cms.InputTag("hfBToJpsiK"), name=cms.string("BToJpsiK"),
                             doc=cms.string("B+ -> mu mu K+ candidates on HI onia muons + HFTrack (Bfinder channel 0); "
                                            "l*_idx -> OniaMuon, k_idx -> HFTrack"))
    t2 = BToTrkTrkMuMuTable.clone(src=cms.InputTag("hfBToJpsiTrkTrk", "SelectedBToTrkTrkMuMu"), name=cms.string("BToJpsiTrkTrk"),
                                  doc=cms.string("B -> mu mu h h (K*0 -> K pi, phi -> K K) candidates on HI onia muons + HFTrack "
                                                 "(Bfinder channels 3-5); l*_idx -> OniaMuon, trk*_idx -> HFTrack"))
    return t1, t2


# ---------------------------------------------------------------------------
#  MC: heavy-flavour gen particles + candidate truth matching
# ---------------------------------------------------------------------------
hfGenParticles = finalGenParticles.clone(
    src = cms.InputTag("mergedGenParticles"),
    select = cms.vstring(
        "drop *",
        # charm / beauty hadrons (mesons 4xx/5xx incl. onia, baryons 4xxx/5xxx) and their (grand)daughters
        "keep++ (abs(pdgId) >= 400 && abs(pdgId) < 600) || (abs(pdgId) >= 4000 && abs(pdgId) < 6000)"
        " || abs(pdgId) == 100443 || abs(pdgId) == 100553 || abs(pdgId) == 200553 || abs(pdgId) == 20443",
    ),
)
hfGenParticleTable = simpleGenParticleFlatTableProducer.clone(
    src = cms.InputTag("hfGenParticles"),
    name = cms.string("HFGenPart"),
    doc = cms.string("charm/beauty hadrons and quarkonia with their daughters and granddaughters"),
    variables = cms.PSet(
        genParticleTable.variables,
        vx = Var("vx", float, doc="production vertex x [cm]", precision=12),
        vy = Var("vy", float, doc="production vertex y [cm]", precision=12),
        vz = Var("vz", float, doc="production vertex z [cm]", precision=12),
    ),
)


def _mcMatch(src, pdgIds, status=None, dR=0.05, dPtRel=0.5, resolveAmbiguities=True):
    return cms.EDProducer("MCMatcher",
        src=cms.InputTag(src) if isinstance(src, str) else src,
        matched=cms.InputTag("hfGenParticles"),
        mcPdgId=cms.vint32(*pdgIds), checkCharge=cms.bool(False),
        mcStatus=cms.vint32(*(status or [])),
        maxDeltaR=cms.double(dR), maxDPtRel=cms.double(dPtRel),
        resolveAmbiguities=cms.bool(resolveAmbiguities), resolveByMatchQuality=cms.bool(True))


def _mcTable(src, mapLabel, objName, doc, objType="Other"):
    return cms.EDProducer("CandMCMatchTableProducer",
        src=cms.InputTag(src) if isinstance(src, str) else src, mcMap=cms.InputTag(mapLabel),
        objName=cms.string(objName), objType=cms.string(objType),   # Muon/Electron/Tau/Photon/Track/Other
        branchName=cms.string("genPart"), docString=cms.string(doc))


def _isMC(process):
    return (hasattr(process, "nanoSequenceMC") and getattr(process, "schedule", None) is not None
            and process.schedule.contains(process.nanoSequenceMC))


def addHIHeavyFlavour(process, doCharm=True, doOnia=True, doBeauty=True, isMC=None, preset="loose"):
    """Add HFTrack + D finder + onia + B finder tables (and their MC truth) to a NANO process.

    preset = "loose" (default): vertex prob > 0.01, no displacement cut (prompt + non-prompt charm).
    preset = "tight": Dfinder-production-like for higher multiplicity: vertex prob > 0.05,
             3D decay-length significance > 2.5, daughter pt > 0.7 (~10x fewer candidates at ~30 tracks/event,
             ~50% efficiency for non-prompt D0).
    """
    if isMC is None:
        isMC = _isMC(process)
    prods, tables = [], []

    def _add(label, mod, lst):
        setattr(process, label, mod.clone() if hasattr(mod, "clone") else mod)
        lst.append(getattr(process, label))

    _add("hiOniaMuons", hiOniaMuons, prods)   # also needed for dzTrg of the tracks
    _add("hiOniaMuMu", hiOniaMuMu, prods)
    _add("hfTracks", hfTracks, prods)
    _add("hfTrackTable", hfTrackTable, tables)
    charm = [("hfD0", hfD0, hfD0Table, [421]), ("hfDstar", hfDstar, hfDstarTable, [413]),
             ("hfDplus", hfDplus, hfDplusTable, [411]), ("hfDs", hfDs, hfDsTable, [431, 411]),
             ("hfLambdaC", hfLambdaC, hfLambdaCTable, [4122])] if doCharm else []
    for label, mod, tab, _ in charm:
        _add(label, mod, prods)
        if preset == "tight":
            m = getattr(process, label)
            m.vtxProbMin = 0.05
            m.svpvSigMin = 2.5
            for d in m.daughters:
                if d.selection.value():
                    d.selection = "pt > 0.7"
        elif preset != "loose":
            raise ValueError("addHIHeavyFlavour: unknown preset '%s'" % preset)
        _add(label + "Table", tab, tables)
    if doOnia:
        _add("hiOniaMuonTable", hiOniaMuonTable, tables)
        _add("hiOniaTable", hiOniaTable, tables)
    if doBeauty:
        t1, t2 = _bTables()
        _add("hfBToJpsiK", hfBToJpsiK, prods)
        _add("hfDiTrack", hfDiTrack, prods)
        _add("hfBToJpsiTrkTrk", hfBToJpsiTrkTrk, prods)
        _add("hfBToJpsiKTable", t1, tables)
        _add("hfBToJpsiTrkTrkTable", t2, tables)
    if isMC:
        if not hasattr(process, "mergedGenParticles"):
            from PhysicsTools.NanoAOD.particlelevel_cff import mergedGenParticles
            process.mergedGenParticles = mergedGenParticles.clone()
            prods.append(process.mergedGenParticles)
        _add("hfGenParticles", hfGenParticles, prods)
        _add("hfGenParticleTable", hfGenParticleTable, tables)
        _add("hfTracksMCMatch", _mcMatch(cms.InputTag("hfTracks", "SelectedTracks"), [211, 321, 2212, 13, 11], [1]), prods)
        _add("hfTrackMCTable", _mcTable(cms.InputTag("hfTracks", "SelectedTracks"), "hfTracksMCMatch", "HFTrack",
                                        "MC matching to status-1 pi/K/p/mu/e in HFGenPart"), tables)
        for label, _, tab, pdg in charm:
            # both K/pi hypotheses of one track pair may match the same gen hadron (flavour tells them apart);
            # exact truth: trkN_idx -> HFTrack_genPartIdx -> common HFGenPart mother
            _add(label + "MCMatch", _mcMatch(label, pdg, dR=0.05, dPtRel=0.3, resolveAmbiguities=False), prods)
            _add(label + "MCTable", _mcTable(label, label + "MCMatch", tab.name.value(),
                                             "MC matching to %s in HFGenPart (deltaR/pt)" % "/".join(map(str, pdg))), tables)
        if doOnia:
            _add("hiOniaMuonsMCMatch", _mcMatch(cms.InputTag("hiOniaMuons", "SelectedMuons"), [13], [1]), prods)
            _add("hiOniaMuonMCTable", _mcTable(cms.InputTag("hiOniaMuons", "SelectedMuons"), "hiOniaMuonsMCMatch",
                                               "OniaMuon", "MC matching to status-1 muons in HFGenPart", objType="Muon"), tables)
    process.hiHeavyFlavourTask = cms.Task(*prods)
    process.hiHeavyFlavourTableTask = cms.Task(*tables)
    if hasattr(process, "schedule") and process.schedule is not None:
        process.schedule.associate(process.hiHeavyFlavourTask)
        process.schedule.associate(process.hiHeavyFlavourTableTask)
    return process
