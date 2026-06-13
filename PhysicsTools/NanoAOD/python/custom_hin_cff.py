"""
Heavy-ion (HIN) NanoAOD flavours: HINUPC and HINHAD.

Two customisation entry points, registered in autoNANO.py and selectable with
`cmsDriver.py ... --nano HINUPC` / `--nano HINHAD`:

  * HINUPCCustomNanoAOD  -- ultra-peripheral collisions (UPC). Adds the PF-candidate
    and (lost-)track tables (BTVNano-like, i.e. PF candidates carrying their track
    details), the ZDC rechit ntuple, and the HF / event-activity table.  The
    centrality *bin* is not defined in UPC, so it is left out; the HF and ZDC sums
    (the important UPC observables) are still written.

  * HINHADCustomNanoAOD  -- hadronic heavy-ion events (e.g. PbPb).  Same PF/track +
    ZDC content, plus the full centrality table including the integer hiBin.

Both build on the standard PHYS NanoAOD content; they only *add* tables, mirroring
the HiForest HiEvtAnalyzer (event activity) and ZDC ntuple.

The new FlatTables are produced by the C++ plugins CentralityTableProducer and
ZDCTableProducer (PhysicsTools/NanoAOD/plugins).  The ZDC rechits are reconstructed
on the fly with the central RecoLocalCalo ZdcHitReconstructor_Run3 from the
hcalDigis:ZDC digis kept in the heavy-ion MiniAOD; the centrality (reco::Centrality
'hiCentrality' and the int 'centralityBin:HFtowers') is taken from the input, as in
the forest.
"""

import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import Var, CandVars, ExtVar
from Configuration.ProcessModifiers.run3_nanoAOD_HIN_cff import run3_nanoAOD_HIN

# ---------------------------------------------------------------------------
#  PF candidates / tracks (BTVNano style: PF candidates with their track details)
# ---------------------------------------------------------------------------
# Variables carried by a packed candidate, including the embedded (pseudo-)track.
_pfCandVars = cms.PSet(
    CandVars,
    puppiWeight=Var("puppiWeight()", float, doc="Puppi weight", precision=10),
    puppiWeightNoLep=Var("puppiWeightNoLep()", float, doc="Puppi weight removing leptons", precision=10),
    vtxChi2=Var("?hasTrackDetails()?vertexChi2():-1", float, doc="vertex chi2", precision=10),
    trkChi2=Var("?hasTrackDetails()?pseudoTrack().normalizedChi2():-1", float, doc="normalized track chi2", precision=10),
    dz=Var("?hasTrackDetails()?dz():-1", float, doc="dz w.r.t. PV", precision=10),
    dzErr=Var("?hasTrackDetails()?dzError():-1", float, doc="dz error", precision=10),
    d0=Var("?hasTrackDetails()?dxy():-1", float, doc="dxy w.r.t. PV", precision=10),
    d0Err=Var("?hasTrackDetails()?dxyError():-1", float, doc="dxy error", precision=10),
    pvAssocQuality=Var("pvAssociationQuality()", int, doc="primary-vertex association quality"),
    lostInnerHits=Var("lostInnerHits()", int, doc="lost inner hits"),
    numberOfHits=Var("numberOfHits()", int, doc="number of tracker hits"),
    numberOfPixelHits=Var("numberOfPixelHits()", int, doc="number of pixel hits"),
    trkHighPurity=Var("?hasTrackDetails()?pseudoTrack().quality('highPurity'):0", bool, doc="track is high purity"),
    trkAlgo=Var("?hasTrackDetails()?pseudoTrack().algo():-1", int, doc="track algorithm"),
    trkPt=Var("?hasTrackDetails()?pseudoTrack().pt():-1", float, doc="track pt", precision=-1),
    trkEta=Var("?hasTrackDetails()?pseudoTrack().eta():-1", float, doc="track eta", precision=12),
    trkPhi=Var("?hasTrackDetails()?pseudoTrack().phi():-1", float, doc="track phi", precision=12),
)

hinPFCandTable = cms.EDProducer(
    "SimplePATCandidateFlatTableProducer",
    src=cms.InputTag("packedPFCandidates"),
    cut=cms.string(""),
    name=cms.string("PFCand"),
    doc=cms.string("PF candidates (packedPFCandidates) with track details"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=_pfCandVars,
)

hinLostTrackTable = cms.EDProducer(
    "SimplePATCandidateFlatTableProducer",
    src=cms.InputTag("lostTracks"),
    cut=cms.string(""),
    name=cms.string("LostTrack"),
    doc=cms.string("Lost tracks (charged PF candidates not used in PF) with track details"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=_pfCandVars,
)

# ---------------------------------------------------------------------------
#  ZDC : reconstruct Run-3 ZDC rechits and tabulate them (HIN UPC ZDC ntuple)
# ---------------------------------------------------------------------------
zdcrecoRun3 = cms.EDProducer(
    "ZdcHitReconstructor_Run3",
    digiLabelQIE10ZDC=cms.InputTag("hcalDigis", "ZDC"),
    Subdetector=cms.string("ZDC"),
    dropZSmarkedPassed=cms.bool(True),
    skipRPD=cms.bool(True),
    recoMethod=cms.int32(1),
    correctionMethodEM=cms.int32(1),
    correctionMethodHAD=cms.int32(1),
    correctionMethodRPD=cms.int32(0),
    ootpuRatioEM=cms.double(3.0),
    ootpuRatioHAD=cms.double(3.0),
    ootpuRatioRPD=cms.double(-1.0),
    ootpuFracEM=cms.double(0.4),
    ootpuFracHAD=cms.double(0.4),
    ootpuFracRPD=cms.double(1.0),
    chargeRatiosEM=cms.vdouble(1.0, 0.23157, 0.10477, 0.06312),
    chargeRatiosHAD=cms.vdouble(1.0, 0.23157, 0.10477, 0.06312),
    chargeRatiosRPD=cms.vdouble(1.0, 0.23157, 0.10477, 0.06312),
    bxTs=cms.vuint32(0, 2, 4),
    nTs=cms.int32(6),
    forceSOI=cms.bool(False),
    signalSOI=cms.vuint32(2),
    noiseSOI=cms.vuint32(1),
    setSaturationFlags=cms.bool(False),
    saturationParameters=cms.PSet(maxADCvalue=cms.int32(255)),
)

zdcTable = cms.EDProducer(
    "ZDCTableProducer",
    src=cms.InputTag("zdcrecoRun3"),
    name=cms.string("ZDC"),
    doc=cms.string("ZDC rechits (forest ZDC ntuple ported to NanoAOD)"),
    precision=cms.int32(-1),
)

# ---------------------------------------------------------------------------
#  Centrality / HF event activity (HiEvtAnalyzer content)
# ---------------------------------------------------------------------------
# "GO" = heavy-ion Global Observables (event-level): centrality bin + HF/ECAL/ZDC
# sums and pixel/track multiplicities (the HiEvtAnalyzer event-level content).
centralityTable = cms.EDProducer(
    "CentralityTableProducer",
    src=cms.InputTag("hiCentrality"),
    srcBin=cms.InputTag("centralityBin", "HFtowers"),
    name=cms.string("GO"),
    doc=cms.string("Heavy-ion global observables: centrality bin, HF/ECAL/ZDC sums, multiplicities"),
    precision=cms.int32(10),
)

# UPC variant: HF / ZDC sums but no centrality bin.
hfTable = centralityTable.clone(
    srcBin=cms.InputTag(""),
    doc=cms.string("Heavy-ion global observables (UPC): HF / ECAL / ZDC sums (no centrality bin)"),
)


# ---------------------------------------------------------------------------
#  Helpers that schedule the producers above into the NANO schedule
# ---------------------------------------------------------------------------
def _associate(process, task):
    """Attach an (unscheduled) Task to the process schedule if there is one."""
    if hasattr(process, "schedule") and process.schedule is not None:
        process.schedule.associate(task)


def addHIPFCands(process):
    """All packed PF candidates (with track details) + lost tracks.

    Central Run-3 NanoAOD already ships a 'PFCand' table built from the jet-constituent
    merge (finalPFCandidates). Heavy-ion analyses want ALL packed PF candidates, so we
    repoint that existing table (and the jet-constituent index table, to keep its indices
    consistent) to packedPFCandidates with our track-detail variables, instead of adding a
    second, colliding 'PFCand' table.
    """
    if hasattr(process, "pfCandidatesTable"):
        process.pfCandidatesTable.src = cms.InputTag("packedPFCandidates")
        process.pfCandidatesTable.variables = _pfCandVars
        for t in ("finalJetsAK8ConstituentsTable", "finalJetsConstituentsTable",
                  "customAK4ConstituentsTable", "customAK8ConstituentsTable"):
            if hasattr(process, t):
                getattr(process, t).candidates = cms.InputTag("packedPFCandidates")
    else:
        process.hinPFCandTable = hinPFCandTable.clone()
        process.hinPFCandsTask = cms.Task(process.hinPFCandTable)
        _associate(process, process.hinPFCandsTask)
    process.hinLostTrackTable = hinLostTrackTable.clone()
    process.hinLostTrackTask = cms.Task(process.hinLostTrackTable)
    _associate(process, process.hinLostTrackTask)
    return process


def addZDCTable(process):
    """Reconstruct Run-3 ZDC rechits and add the ZDC rechit + per-side-sum tables."""
    # ZdcHitReconstructor_Run3 needs the HcalSeverityLevelComputerRcd record
    # (same ES producer the forest loads before its ZDC sequence).
    process.load("RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi")
    process.zdcrecoRun3 = zdcrecoRun3.clone()
    process.zdcTable = zdcTable.clone()
    process.zdcTableTask = cms.Task(process.zdcrecoRun3, process.zdcTable)
    _associate(process, process.zdcTableTask)
    return process


def addCentralityTable(process, addBin=True):
    """Add the HF / centrality event-activity table.

    addBin=True  -> full centrality table including the integer hiBin (HINHAD).
    addBin=False -> HF / ZDC sums only, no centrality bin (HINUPC).

    For HINHAD the hiBin is recomputed in-process by CentralityTableProducer
    itself from the 2025 recalibrated table (hiCentralityBinTable2025_cff), the
    same EtHFtowerSum lookup as the forest's HICentralityBinProducer, instead of
    the PromptReco bin in the file - the two differ by ~11 bins (~5.5%
    centrality) on 2025 data. So GO_hiBin matches the forest hiBin, with no
    separate producer and no HeavyIonsAnalysis dependency.
    """
    process.centralityTable = (centralityTable if addBin else hfTable).clone()
    if addBin:
        from PhysicsTools.NanoAOD.hiCentralityBinTable2025_cff import hiCentralityBinTable2025
        process.centralityTable.table = hiCentralityBinTable2025
    process.centralityTableTask = cms.Task(process.centralityTable)
    _associate(process, process.centralityTableTask)
    return process


def addHIJets(process, labels=("4", "4Flow"), doBtagging=False, jetPtMin=15.0):
    """Schedule the heavy-ion jet reconstruction (constituent-subtracted akCs PF jets
    and flow-subtracted akCsFlow PF jets) from the forest setup, and tabulate the
    resulting pat::Jets as NanoAOD FlatTables via HiInclusiveJetTableProducer.

    Reuses PhysicsTools.NanoAOD.hi_setupJets_PbPb_cff.candidateBtaggingMiniAOD
    (the same reco the forest config runs). On a plain CMSSW area without the forest
    package this is a no-op (so HINHAD still imports on the non-forest branch).
    """
    try:
        from PhysicsTools.NanoAOD.hi_setupJets_PbPb_cff import candidateBtaggingMiniAOD
    except ImportError:
        print("[custom_hin_cff] HeavyIonsAnalysis jet setup not found -> skipping HI jets. "
              "This should not happen on the self-contained nano branch.")
        return process

    # data vs MC: is the MC nano sequence actually scheduled? (same check BTVNano uses)
    isMC = (hasattr(process, "nanoSequenceMC") and getattr(process, "schedule", None) is not None
            and process.schedule.contains(process.nanoSequenceMC))
    jetCorrLevels = ["L2Relative", "L3Absolute"] if isMC else ["L2Relative", "L2L3Residual"]

    # candidateBtaggingMiniAOD associates its reco tasks to process.forest; provide a stub
    # path (added to the schedule) so the unscheduled reco modules actually run.
    if not hasattr(process, "forest"):
        process.forest = cms.Path()
        if getattr(process, "schedule", None) is not None:
            process.schedule.append(process.forest)

    if not hasattr(process, "hiJetTableTask"):
        process.hiJetTableTask = cms.Task()
        _associate(process, process.hiJetTableTask)

    for label in labels:
        candidateBtaggingMiniAOD(process, isMC=isMC, jetPtMin=jetPtMin,
                                 jetCorrLevels=jetCorrLevels, doBtagging=doBtagging, labelR=label)
        coll = "selectedUpdatedPatJetsAK" + label + "PFBtag"
        tabname = "akCs" + label + "PFJet"  # akCs4PFJet, akCs4FlowPFJet
        rParam = 0.4 if label == "0" else float(label.replace("Flow", "")) * 0.1
        # tabulate the embedded UParT outputs (PbPb training; un-postfixed names)
        discTags, discNames = [], []
        if doBtagging:
            _uparT = "pfUnifiedParticleTransformerAK4JetTags"
            for f in ("probb", "probbb", "problepb", "probc", "probg", "probu", "probd", "probs",
                      "probtaup1h0p", "probtaup1h1p", "probtaup1h2p", "probtaup3h0p", "probtaup3h1p",
                      "probtaum1h0p", "probtaum1h1p", "probtaum1h2p", "probtaum3h0p", "probtaum3h1p",
                      "probele", "probmu", "ptcorr", "ptnu"):
                discTags.append(_uparT + ":" + f)
                discNames.append("UParT_" + f)
                discTags.append("pfNegative" + _uparT[2:] + ":" + f)
                discNames.append("UParTNeg_" + f)
        setattr(process, tabname + "Table", cms.EDProducer(
            "HiInclusiveJetTableProducer",
            discriminatorTags=cms.vstring(discTags),
            discriminatorNames=cms.vstring(discNames),
            jets=cms.InputTag(coll),
            pfCandidates=cms.InputTag("packedPFCandidates"),
            name=cms.string(tabname),
            doc=cms.string("Heavy-ion inclusive jets (" + tabname + ")"),
            precision=cms.int32(-1),
            jetPtMin=cms.double(jetPtMin),
            jetAbsEtaMax=cms.double(2.5),
            rParam=cms.double(rParam),
            hardPtMin=cms.double(4.0),
            trackQuality=cms.string("highPurity"),
            useQuality=cms.bool(True),
            doHiJetID=cms.bool(True),
            doWTARecluster=cms.bool(True),
        ))
        process.hiJetTableTask.add(getattr(process, tabname + "Table"))
    return process


def addHITracks(process):
    """Unpack reco::Tracks + vertices from packedPFCandidates (forest
    TrackAndVertexUnpacker) and tabulate them as the Trk + Vtx FlatTables.
    No-op if the forest TrackAnalysis package is absent."""
    try:
        from PhysicsTools.NanoAOD.unpackedTracksAndVertices_cfi import unpackedTracksAndVertices
    except ImportError:
        print("[custom_hin_cff] HeavyIonsAnalysis TrackAnalysis not found -> skipping HI tracks.")
        return process
    process.unpackedTracksAndVertices = unpackedTracksAndVertices.clone()
    process.trkTable = cms.EDProducer(
        "TrackTableProducer",
        trackSrc=cms.InputTag("unpackedTracksAndVertices"),
        vertexSrc=cms.InputTag("unpackedTracksAndVertices"),
        name=cms.string("Trk"),
        vtxName=cms.string("Vtx"),
        doc=cms.string("Heavy-ion tracks (unpacked from packedPFCandidates)"),
        precision=cms.int32(-1),
        trackPtMin=cms.double(0.01),
        trackEtaMax=cms.double(4.0),
        applyTrackSelections=cms.bool(False),
    )
    process.hiTrackTask = cms.Task(process.unpackedTracksAndVertices, process.trkTable)
    _associate(process, process.hiTrackTask)
    return process


def _removeNanoModules(process, names):
    """Remove the named modules from every task/sequence/path/endpath so they are not
    scheduled (and their tables are not written)."""
    for n in names:
        if not hasattr(process, n):
            continue
        mod = getattr(process, n)
        for attr in ("tasks", "sequences", "paths", "endpaths"):
            d = getattr(process, attr, None)
            if d is None:
                continue
            try:
                for cont in d.values():
                    try:
                        cont.remove(mod)
                    except Exception:
                        pass
            except AttributeError:
                pass
    return process


def tolerateMissingProducts(process):
    """Safety net: let the job continue and still write output if some producer's input
    product is absent (heavy-ion MiniAOD content varies). HI table producers consume
    in-process products and are unaffected."""
    if not hasattr(process, "options"):
        process.options = cms.untracked.PSet()
    process.options.TryToContinue = cms.untracked.vstring("ProductNotFound")
    return process


def stripPPonlyContent(process):
    """Hadronic heavy-ion MiniAOD (e.g. HINPbPbWinter25) drops several pp-only
    collections. Remove the standard NanoAOD tables that read the absent ones, and give
    the cross-linked low-pT-electron chain an empty input so it runs harmlessly."""
    # AK8 soft-drop subjets (slimmedJetsAK8PFPuppiSoftDropPacked:SubJets) and Puppi MET
    # (slimmedMETsPuppi) are absent -> just drop those standalone tables.
    _removeNanoModules(process, ["subJetTable", "subjetMCTable", "puppiMetTable", "rawPuppiMetTable"])
    # low-pT electrons (slimmedLowPtElectrons) are absent but cross-linked via
    # linkedObjects; feed the chain an always-empty stand-in collection.
    if not hasattr(process, "slimmedLowPtElectrons"):
        process.slimmedLowPtElectrons = cms.EDFilter(
            "PATElectronSelector",
            src=cms.InputTag("slimmedElectrons"),
            cut=cms.string("pt < 0"),  # always empty
        )
        process.hiEmptyLowPtEleTask = cms.Task(process.slimmedLowPtElectrons)
        _associate(process, process.hiEmptyLowPtEleTask)
    return process


def addHIMuons(process):
    """Add an HI-specific extension to the standard NanoAOD Muon table: inner/global-track
    detail and R03 track isolation (the forest MuonAnalyzer content not in standard nano).

    Note: the forest also re-unpacks muons (unpackedMuons, recovering HI muons from
    packed/lost tracks). Routing the standard muon chain through it requires giving
    slimmedMuonsWithUserData's value-map embedder the right parentSrcs; left as a
    follow-up. This extension uses the standard (slimmedMuons-based) Muon collection."""
    # HI extension aligned with the standard Muon table
    if hasattr(process, "muonTable"):
        process.hiMuonExtTable = cms.EDProducer(
            "SimplePATMuonFlatTableProducer",
            src=process.muonTable.src,
            cut=cms.string(""),
            name=cms.string("Muon"),
            doc=cms.string("HI muon extension (unpacked source; inner/global track, R03 iso)"),
            singleton=cms.bool(False),
            extension=cms.bool(True),
            variables=cms.PSet(
                isoTrkR03=Var("isolationR03().sumPt", float, doc="sum track pt in R03 cone", precision=10),
                innerPt=Var("?innerTrack.isNonnull()?innerTrack().pt():-1", float, doc="inner-track pt", precision=-1),
                innerPtErr=Var("?innerTrack.isNonnull()?innerTrack().ptError():-1", float, doc="inner-track pt error", precision=10),
                innerEta=Var("?innerTrack.isNonnull()?innerTrack().eta():-99", float, doc="inner-track eta", precision=12),
                innerNTrkHits=Var("?innerTrack.isNonnull()?innerTrack().hitPattern().numberOfValidTrackerHits():-1", int, doc="inner-track valid tracker hits"),
                innerNPixHits=Var("?innerTrack.isNonnull()?innerTrack().hitPattern().numberOfValidPixelHits():-1", int, doc="inner-track valid pixel hits"),
                innerHighPurity=Var("?innerTrack.isNonnull()?innerTrack().quality('highPurity'):0", bool, doc="inner-track high purity"),
                globalPt=Var("?globalTrack.isNonnull()?globalTrack().pt():-1", float, doc="global-track pt", precision=-1),
                globalNormChi2=Var("?globalTrack.isNonnull()?globalTrack().normalizedChi2():-1", float, doc="global-track norm chi2", precision=10),
                globalNMuonHits=Var("?globalTrack.isNonnull()?globalTrack().hitPattern().numberOfValidMuonHits():-1", int, doc="global-track valid muon hits"),
            ),
        )
        process.hiMuonExtTask = cms.Task(process.hiMuonExtTable)
        _associate(process, process.hiMuonExtTask)
    return process


def addHIEGM(process):
    """Add HI custom-cone PF isolation (charged/photon/neutral, dR 0.3 & 0.4) as
    extension columns to the standard Electron and Photon tables, ported from the
    ggHiNtuplizer pfIsoCalculator. The standard nano e/gamma tables already cover the
    rest of the ggHiNtuplizer content (kinematics, SC, shower shapes, IDs, std iso)."""
    specs = [
        ("electron", "electronTable", "hiElectronIso", "SimplePATElectronFlatTableProducer"),
        ("photon", "photonTable", "hiPhotonIso", "SimplePATPhotonFlatTableProducer"),
    ]
    for obj, table, isoMod, prod in specs:
        if not hasattr(process, table):
            continue
        src = getattr(process, table).src
        setattr(process, isoMod, cms.EDProducer(
            "HIEGMIsoProducer",
            src=src,
            pfCandidates=cms.InputTag("packedPFCandidates"),
            vertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
            dzMax=cms.double(0.2),
            dxyMax=cms.double(0.1),
        ))
        setattr(process, isoMod + "ExtTable", cms.EDProducer(
            prod,
            src=src,
            cut=cms.string(""),
            name=getattr(process, table).name,
            doc=cms.string("HI cone isolation"),
            singleton=cms.bool(False),
            extension=cms.bool(True),
            variables=cms.PSet(),
            externalVariables=cms.PSet(
                hiPFChIso03=ExtVar(cms.InputTag(isoMod, "chIso03"), float, doc="HI charged-hadron iso, dR<0.3", precision=10),
                hiPFChIso04=ExtVar(cms.InputTag(isoMod, "chIso04"), float, doc="HI charged-hadron iso, dR<0.4", precision=10),
                hiPFPhoIso03=ExtVar(cms.InputTag(isoMod, "phoIso03"), float, doc="HI photon iso, dR<0.3", precision=10),
                hiPFPhoIso04=ExtVar(cms.InputTag(isoMod, "phoIso04"), float, doc="HI photon iso, dR<0.4", precision=10),
                hiPFNeuIso03=ExtVar(cms.InputTag(isoMod, "nhIso03"), float, doc="HI neutral-hadron iso, dR<0.3", precision=10),
                hiPFNeuIso04=ExtVar(cms.InputTag(isoMod, "nhIso04"), float, doc="HI neutral-hadron iso, dR<0.4", precision=10),
            ),
        ))
        t = cms.Task(getattr(process, isoMod), getattr(process, isoMod + "ExtTable"))
        setattr(process, obj + "HIIsoTask", t)
        _associate(process, t)
    return process


def addFullHadSkim(process,
                   nJets=5, jetPtMin=25.0, jetAbsEta=2.1,
                   requireBTags=False, nBJets=2, bDiscrMin=0.15,
                   usePreCount=True):
    """Skim the NanoAOD for the fully hadronic ttbar / weak-supervision (CWoLa)
    PbPb analysis: keep events with >= nJets akCs3PF jets with UParT-regressed
    pT (= rawPt * UParT ptcorr) >= jetPtMin and |eta| <= jetAbsEta (the offline
    analysis acceptance: eta is not recalibrated, so unlike pT no margin is needed).

    Stays strictly looser than the offline selection (>=6 jets, pT > 30,
    |eta| < 2.1 signal region; nJet == 5 sideband), and applies NO b-tag
    requirement so the full b-tag spectrum (0B/1B mixtures, negative tags)
    survives for the data-driven region building; requireBTags=True
    re-enables a >= nBJets cut on the normalized UParT b discriminant
    ((b+bb+lepb)/(b+bb+lepb+c+s+u+d+g) >= bDiscrMin).

    Implementation mirrors the validated HiForest full-had skim
    (forest_miniAOD_ParticleTransformer_run3_SKIM_FULLHAD_DATA.py): an
    unbiased >= nJets precount on the bare patJetsAKCs3PF (raw pT at the
    15 GeV clustering threshold) keeps the b-tag inference off for most
    events. NB embedded discriminators carry the UN-postfixed tagger names
    (pfUnifiedParticleTransformerAK4JetTags:*). Measured on 2025 PbPb data
    (HiForest twin, run 400059): ~8e-5 pass rate, all 0-5% central.

    The skim gates the scheduled nano tables (filters prepended to
    nanoAOD_step) and the NanoAOD output (SelectEvents on nanoAOD_step);
    task-resident HI tables are gated automatically through the output
    module prefetch.
    """
    # make sure the akCs3PF chain with b-tagging (UParT incl. ptcorr) and its
    # table exist; addHIJets is a no-op without HeavyIonsAnalysis in the area
    if not hasattr(process, "selectedUpdatedPatJetsAK3PFBtag"):
        addHIJets(process, labels=("3",), doBtagging=True)
    if not hasattr(process, "selectedUpdatedPatJetsAK3PFBtag"):
        raise RuntimeError(
            "addFullHadSkim: akCs3 jet chain unavailable - check out "
            "HeavyIonsAnalysis/JetAnalysis (forest branch) in this area.")

    # pin the skim to the forest's 2024 PbPb UParT training so scores (and the
    # measured skim rate) are identical to the HiForest full-had skim twin; the
    # central cms-data 2023 model regresses jets ~35% higher at low pT
    # (ptcorr median 0.80 vs 0.58 on 2025 data), which would loosen the pT cut
    uparTModel = "PhysicsTools/NanoAOD/data/PbPb_AK3_2024_v6.onnx"
    for m in ("pfUnifiedParticleTransformerAK4JetTagsAK3PFBtag",
              "pfNegativeUnifiedParticleTransformerAK4JetTagsAK3PFBtag"):
        if hasattr(process, m):
            getattr(process, m).model_path = uparTModel

    uparT = "pfUnifiedParticleTransformerAK4JetTags"
    uparTB = "+".join("bDiscriminator('%s:%s')" % (uparT, x)
                      for x in ("probb", "probbb", "problepb"))
    uparTAll = uparTB + "+" + "+".join("bDiscriminator('%s:%s')" % (uparT, x)
                                       for x in ("probc", "probs", "probu", "probd", "probg"))
    uparTPt = "correctedJet('Uncorrected').pt * bDiscriminator('%s:ptcorr')" % uparT

    # cheap unbiased precount on the bare subtracted jets (threshold = clustering pT min)
    process.fullHadPreJets = cms.EDFilter("PATJetSelector",
        src = cms.InputTag("patJetsAKCs3PF"),
        cut = cms.string("correctedJet('Uncorrected').pt >= 15.0 && abs(eta) <= 2.3"))
    process.fullHadPreJetCount = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag("fullHadPreJets"), minNumber = cms.uint32(nJets))

    # jet counting on the UParT-regressed pT (offline: jtupartpt = jtrawPt * ptcorr)
    process.fullHadJets = cms.EDFilter("PATJetSelector",
        src = cms.InputTag("selectedUpdatedPatJetsAK3PFBtag"),
        cut = cms.string("(%s) >= %g && abs(eta) <= %g" % (uparTPt, jetPtMin, jetAbsEta)))
    process.fullHadJetCount = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag("fullHadJets"), minNumber = cms.uint32(nJets))

    _seq = process.fullHadJets * process.fullHadJetCount
    if requireBTags:
        process.fullHadBJets = cms.EDFilter("PATJetSelector",
            src = cms.InputTag("fullHadJets"),
            cut = cms.string("(%s) > 0 && (%s)/(%s) >= %g" % (uparTB, uparTB, uparTAll, bDiscrMin)))
        process.fullHadBJetCount = cms.EDFilter("CandViewCountFilter",
            src = cms.InputTag("fullHadBJets"), minNumber = cms.uint32(nBJets))
        _seq = _seq * process.fullHadBJets * process.fullHadBJetCount
    if usePreCount:
        _seq = process.fullHadPreJets * process.fullHadPreJetCount * _seq
    process.fullHadSkimSequence = cms.Sequence(_seq)

    # gate the scheduled nano tables and skim the output
    if hasattr(process, "nanoAOD_step"):
        process.nanoAOD_step._seq = process.fullHadSkimSequence * process.nanoAOD_step._seq
    for out in process.outputModules_().values():
        if out.type_() == "NanoAODOutputModule":
            out.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("nanoAOD_step"))
    return process


# ---------------------------------------------------------------------------
#  Flavour entry points (referenced from autoNANO.py)
# ---------------------------------------------------------------------------
def HINUPCCustomNanoAOD(process):
    addHIPFCands(process)
    addZDCTable(process)
    addCentralityTable(process, addBin=False)  # HF info, no centrality bin (UPC)
    tolerateMissingProducts(process)
    return process


def HINHADCustomNanoAOD(process):
    addHIPFCands(process)
    addZDCTable(process)
    addCentralityTable(process, addBin=True)  # full centrality incl. hiBin
    addHIJets(process)                        # akCs4PF + akCs4FlowPF HI jets (forest reco)
    addHITracks(process)                      # unpacked reco::Tracks + vertices
    addHIMuons(process)                       # unpacked muons -> Muon table + HI extension
    addHIEGM(process)                         # HI cone isolation on Electron/Photon tables
    stripPPonlyContent(process)               # drop pp-only tables absent from HI MiniAOD
    tolerateMissingProducts(process)
    return process


# When the run3_nanoAOD_HIN process modifier is active (cmsDriver --procModifiers
# run3_nanoAOD_HIN) apply a soft pT threshold to the PF-candidate / lost-track tables:
# heavy-ion events have very high multiplicity, so this keeps the NanoAOD size under
# control without affecting the high-pT physics objects.
run3_nanoAOD_HIN.toModify(hinPFCandTable, cut="pt > 0.5")
run3_nanoAOD_HIN.toModify(hinLostTrackTable, cut="pt > 0.5")
