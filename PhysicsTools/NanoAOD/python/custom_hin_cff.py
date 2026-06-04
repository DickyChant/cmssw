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

from PhysicsTools.NanoAOD.common_cff import Var, CandVars
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
centralityTable = cms.EDProducer(
    "CentralityTableProducer",
    src=cms.InputTag("hiCentrality"),
    srcBin=cms.InputTag("centralityBin", "HFtowers"),
    name=cms.string("Cent"),
    doc=cms.string("Heavy-ion event activity: HF / ECAL / ZDC sums, multiplicities and centrality bin"),
    precision=cms.int32(10),
)

# UPC variant: HF / ZDC sums but no centrality bin.
hfTable = centralityTable.clone(
    srcBin=cms.InputTag(""),
    doc=cms.string("Heavy-ion event activity for UPC: HF / ECAL / ZDC sums (no centrality bin)"),
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
    """
    process.centralityTable = (centralityTable if addBin else hfTable).clone()
    process.centralityTableTask = cms.Task(process.centralityTable)
    _associate(process, process.centralityTableTask)
    return process


# ---------------------------------------------------------------------------
#  Flavour entry points (referenced from autoNANO.py)
# ---------------------------------------------------------------------------
def HINUPCCustomNanoAOD(process):
    addHIPFCands(process)
    addZDCTable(process)
    addCentralityTable(process, addBin=False)  # HF info, no centrality bin (UPC)
    return process


def HINHADCustomNanoAOD(process):
    addHIPFCands(process)
    addZDCTable(process)
    addCentralityTable(process, addBin=True)  # full centrality incl. hiBin
    return process


# When the run3_nanoAOD_HIN process modifier is active (cmsDriver --procModifiers
# run3_nanoAOD_HIN) apply a soft pT threshold to the PF-candidate / lost-track tables:
# heavy-ion events have very high multiplicity, so this keeps the NanoAOD size under
# control without affecting the high-pT physics objects.
run3_nanoAOD_HIN.toModify(hinPFCandTable, cut="pt > 0.5")
run3_nanoAOD_HIN.toModify(hinLostTrackTable, cut="pt > 0.5")
