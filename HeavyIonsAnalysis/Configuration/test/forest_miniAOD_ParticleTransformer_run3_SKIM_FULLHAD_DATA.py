### HiForest Configuration
# Input: miniAOD
# Type: data
#
# Skim: full-hadronic ttbar multijet selection
#       (>= 5 akCs3PF jets with UParT-regressed pT >= 25 GeV, |eta| <= 2.1;
#        no b-tag requirement - optional via requireBTags)
#
# Derived from forest_miniAOD_ParticleTransformer_run3_SKIM_DATA.py
# (stahlleiton/cmssw @ HIForest_TTBAR_Run3_2025_PbPb). Only the
# "Event Selection" block at the bottom differs from the lepton-skim config:
# the lepton/photon gate is replaced by the full-hadronic multijet gate,
# kept strictly looser than the offline hin2023 full-hadronic selection.

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_pp_on_PbPb_2025_cff import Run3_pp_on_PbPb_2025
process = cms.Process('HiForest', Run3_pp_on_PbPb_2025)

###############################################################################

# HiForest info
process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 151X, data, fullhad skim")

###############################################################################

# input files
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/hidata/HIRun2025A/HIPhysicsRawPrime0/MINIAOD/PbPbEW-PromptReco-v1/000/400/414/00000/ff20df76-fe62-433b-8543-710276680707.root'),
)

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

###############################################################################

# load Global Tag, geometry, etc.
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '151X_dataRun3_Prompt_v1', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag

###############################################################################

# Define centrality binning (2025 recalibrated table, shared cfi - also used by
# the NanoAOD HIN customization so GO_hiBin matches the forest hiBin)
from HeavyIonsAnalysis.EventAnalysis.hiCentralityBin2025_cfi import centralityBin2025
process.centralityBin = centralityBin2025.clone()

#process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
#process.GlobalTag.toGet.extend([
#    cms.PSet(
#        record = cms.string("HeavyIonRcd"),
#        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v140x01_offline_Nominal"),
#        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
#        label = cms.untracked.string("HFtowers")
#    ),
#])

###############################################################################

# root output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("HiForestMiniAOD.root"))

###############################################################################

# event analysis
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')
process.metFilters = process.skimanalysis.clone(hltresults = "TriggerResults::RECO")

from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_data_2025_skimmed
process.hltobject.triggerNames = trigger_list_data_2025_skimmed

process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
################################
# electrons, photons, muons
process.load('HeavyIonsAnalysis.EGMAnalysis.ggHiNtuplizer_cfi')
process.load('HeavyIonsAnalysis.EGMAnalysis.hiElectrons_cfi')
process.hiElectrons.file_idModel = "HeavyIonsAnalysis/EGMAnalysis/data/Run3_2024_PbPb/eleid_BDT.ubj"
process.hiElectrons.file_isoModel = "HeavyIonsAnalysis/EGMAnalysis/data/Run3_2024_PbPb/eleiso_BDT.ubj"
process.hiElectrons.file_corr = "HeavyIonsAnalysis/Configuration/data/lepton_spectra_train_weights_Run3_2024_PbPb.json.gz"
process.hiElectrons.era = "Run3_2024_PbPb"
process.ggHiNtuplizer.electronSrc = "hiElectrons"
process.egammaSequence = cms.Sequence(process.hiElectrons * process.ggHiNtuplizer)
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
################################
# jet reco sequence
process.load('HeavyIonsAnalysis.JetAnalysis.akCs4PFJetSequence_pponPbPb_data_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.hiFJRhoAnalyzer_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.hiFJSoftKillerAnalyzer_cff')
process.load('HeavyIonsAnalysis.JetAnalysis.hiFlowRhoAnalyzer_cff')
################################
# tracks
process.load("HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzers_cff")
# muons
process.load("HeavyIonsAnalysis.MuonAnalysis.unpackedMuons_cfi")
process.load('HeavyIonsAnalysis.MuonAnalysis.hiMuons_cfi')
process.hiMuons.muon_minPt = 10
process.hiMuons.file_isoModel = "HeavyIonsAnalysis/MuonAnalysis/data/muiso_BDT.ubj"
process.hiMuons.file_isoCorr = "HeavyIonsAnalysis/Configuration/data/lepton_spectra_train_weights_Run3_2024_PbPb.json.gz"
process.hiMuons.era = "Run3_2024_PbPb"
process.unpackedMuons.muons = "hiMuons"
process.muonSequence = cms.Sequence(process.hiMuons * process.unpackedMuons)
process.load("HeavyIonsAnalysis.MuonAnalysis.muonAnalyzer_cfi")
###############################################################################

#########################
# ZDC RecHit Producer && Analyzer
#########################
# to prevent crash related to HcalSeverityLevelComputerRcd record
process.load("RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi")
process.load('HeavyIonsAnalysis.ZDCAnalysis.ZDCAnalyzersPbPb_cff')
process.zdcanalyzer.doZdcDigis = False

###############################################################################
# main forest sequence
process.forest = cms.Path(
    process.HiForestInfo +
    process.centralityBin +
    process.hiEvtAnalyzer +
    process.hltanalysis +
    process.hltobject +
    process.l1object +
    process.unpackedTracksAndVertices +
    process.particleFlowAnalyser +
    process.rhoSequence +
    process.muonSequence +
    process.egammaSequence +
    process.hiFJSoftKillerAnalyzer +
    process.rhoFlowDataSequence +
    process.metFilters +
    process.zdcSequencePbPb
    )

#customisation
process.particleFlowAnalyser.ptMin = 0.0
process.ggHiNtuplizer.muonPtMin = 0.0

# Select the types of jets filled
jetPtMin = 15
jetAbsEtaMax = 2.5

# Choose which additional information is added to jet trees
doHIJetID = True             # Fill jet ID and composition information branches
doWTARecluster = True        # Add jet phi and eta for WTA axis

# add candidate tagging
for jetR, doFlow in zip([0.3, 0.3, 0.4], [False, True, True]):
    R = str(int(jetR*10))
    from HeavyIonsAnalysis.JetAnalysis.deepNtupleSettings_cff import candidateBtaggingMiniAOD
    candidateBtaggingMiniAOD(process, isMC = False, jetPtMin = jetPtMin, jetR = jetR, jetCorrLevels = ['L2Relative', 'L2L3Residual'], doFlow = doFlow)

    # setup jet analyzer
    jL = f"Cs{R}Flow" if doFlow else f"Cs{R}"
    setattr(process,f'ak{jL}PFJetAnalyzer', process.akCs4PFJetAnalyzer.clone())
    getattr(process,f'ak{jL}PFJetAnalyzer').jetTag = f'selectedUpdatedPatJetsAK{jL}DeepFlavour'
    getattr(process,f'ak{jL}PFJetAnalyzer').jetName = f'ak{jL}PF'
    getattr(process,f'ak{jL}PFJetAnalyzer').rParam = jetR
    getattr(process,f'ak{jL}PFJetAnalyzer').doHiJetID = doHIJetID
    getattr(process,f'ak{jL}PFJetAnalyzer').doWTARecluster = doWTARecluster
    getattr(process,f'ak{jL}PFJetAnalyzer').jetPtMin = jetPtMin
    getattr(process,f'ak{jL}PFJetAnalyzer').useRawPt = True
    getattr(process,f'ak{jL}PFJetAnalyzer').jetAbsEtaMax = cms.untracked.double(jetAbsEtaMax)
    getattr(process,f'ak{jL}PFJetAnalyzer').pfJetProbabilityBJetTag = cms.untracked.string(f"pfJetProbabilityBJetTagsAK{jL}DeepFlavour")
    getattr(process,f'ak{jL}PFJetAnalyzer').pfNegativeOnlyJetProbabilityBJetTag = cms.untracked.string(f"pfNegativeOnlyJetProbabilityBJetTagsAK{jL}DeepFlavour")
    getattr(process,f'ak{jL}PFJetAnalyzer').pfDeepCSVJetTags = cms.untracked.string(f"pfDeepCSVJetTagsAK{jL}DeepFlavour")
    getattr(process,f'ak{jL}PFJetAnalyzer').pfDeepFlavourJetTags = cms.untracked.string(f"pfDeepFlavourJetTagsAK{jL}DeepFlavour")
    getattr(process,f'ak{jL}PFJetAnalyzer').pfParticleTransformerAK4JetTags = cms.untracked.string(f"pfParticleTransformerAK4JetTagsAK{jL}DeepFlavour")
    getattr(process,f'ak{jL}PFJetAnalyzer').pfUnifiedParticleTransformerAK4JetTags = cms.untracked.string(f"pfUnifiedParticleTransformerAK4JetTagsAK{jL}DeepFlavour")
    getattr(process,f'ak{jL}PFJetAnalyzer').pfNegativeUnifiedParticleTransformerAK4JetTags = cms.untracked.string(f"pfNegativeUnifiedParticleTransformerAK4JetTagsAK{jL}DeepFlavour")
    process.forest += getattr(process,f'ak{jL}PFJetAnalyzer')


#########################
# Event Selection -> add the needed filters here
#########################

process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.load('HeavyIonsAnalysis.EventAnalysis.hffilter_cfi')
process.pphfCoincFilter4Th2 = cms.Path(process.phfCoincFilter4Th2)
process.pphfCoincFilter1Th3 = cms.Path(process.phfCoincFilter1Th3)
process.pphfCoincFilter2Th3 = cms.Path(process.phfCoincFilter2Th3)
process.pphfCoincFilter3Th3 = cms.Path(process.phfCoincFilter3Th3)
process.pphfCoincFilter4Th3 = cms.Path(process.phfCoincFilter4Th3)
process.pphfCoincFilter5Th3 = cms.Path(process.phfCoincFilter5Th3)
process.pphfCoincFilter1Th4 = cms.Path(process.phfCoincFilter1Th4)
process.pphfCoincFilter2Th4 = cms.Path(process.phfCoincFilter2Th4)
process.pphfCoincFilter3Th4 = cms.Path(process.phfCoincFilter3Th4)
process.pphfCoincFilter4Th4 = cms.Path(process.phfCoincFilter4Th4)
process.pphfCoincFilter5Th4 = cms.Path(process.phfCoincFilter5Th4)
process.pphfCoincFilter1Th5 = cms.Path(process.phfCoincFilter1Th5)
process.pphfCoincFilter2Th5 = cms.Path(process.phfCoincFilter2Th5)
process.pphfCoincFilter3Th5 = cms.Path(process.phfCoincFilter3Th5)
process.pphfCoincFilter4Th5 = cms.Path(process.phfCoincFilter4Th5)
process.pphfCoincFilter5Th5 = cms.Path(process.phfCoincFilter5Th5)
process.pAna = cms.EndPath(process.skimanalysis)

#########################
# The full-had skim gate depends (via the Cs3 jet chain in patAlgosToolsTask)
# on PackedPFTowers and hiPuRho, which extraFlowJetsData schedules on the
# forest path AFTER the prepended gate -> unrunnable schedule (deadlock).
# Drop them from the explicit forest scheduling; they remain in
# patAlgosToolsTask and run on demand (for every event, as inputs of the
# ungated fullHadSkimPath).
#########################
for _mod in [process.PackedPFTowers, process.hiPuRho]:
    while process.forest.remove(_mod):
        pass

#########################
# Skim: full-hadronic ttbar multijet selection (akCs3PF, UParT)
#
# Offline (hin2023 new_full_had) selection this must stay LOOSER than:
#   jets:   UParT-regressed pT = rawPt * UParT(ptcorr), pT > 30, |eta| < 2.1
#   filter: nJet >= 6 && nJetBLoose >= 2 (+JEC-varied counterparts)
#   BLoose: discr_btag = UParT(b+bb+lepb)/(b+bb+lepb+c+s+u+d+g) > wp(cen),
#           with wp(cen) in [0.21 (0-5%), 0.45 (peripheral)]  (2024-derived WPs)
#########################
fullHadNJets     = 5      # min number of selected jets (offline: >= 6 signal region,
                          # nJet == 5 kept as weak-supervision sideband -> skim at >= 5)
fullHadJetPtMin  = 25.0   # GeV, on UParT-regressed pT (offline: 30; margin for JEC variations/residuals)
fullHadJetAbsEta = 2.1    # = offline cut; eta is not recalibrated, so no margin is
                          # needed (unlike pT) - an offline jet at |eta|<2.1 is
                          # the same jet at skim level
requireBTags     = False  # data-driven (CWoLa/anomaly-detection) analysis needs the FULL
                          # b-tag spectrum in the forest (0B/1B depleted mixtures, negative
                          # tags); measured on 2025 data the b-cut also gives ~no rate
                          # reduction (the jet-count cut is the bottleneck) -> keep OFF
fullHadNBJets    = 2      # (only if requireBTags) min loose b-tagged jets (offline: >= 2 loose)
fullHadBDiscrMin = 0.15   # (only if requireBTags) normalized UParT b discr (offline loose WP >= 0.21)
usePreCountFilter = True  # cheap >=6-jet precount on untagged PAT jets, so the
                          # b-tagging chain only runs for multijet candidates;
                          # threshold = clustering jetPtMin, so it counts every
                          # clustered jet in |eta|<=2.6 -> structurally unbiased

# NB: discriminators embedded by updateJetCollection keep the UN-postfixed
# tagger names (verified on 2025 data): 'pfUnifiedParticleTransformerAK4JetTags:*',
# not 'pfUnifiedParticleTransformerAK4JetTagsAKCs3DeepFlavour:*'
_uparT  = 'pfUnifiedParticleTransformerAK4JetTags'
_uparTB   = '+'.join(f"bDiscriminator('{_uparT}:{p}')" for p in ['probb','probbb','problepb'])
_uparTAll = _uparTB + '+' + '+'.join(f"bDiscriminator('{_uparT}:{p}')" for p in ['probc','probs','probu','probd','probg'])
_uparTPt  = f"correctedJet('Uncorrected').pt * bDiscriminator('{_uparT}:ptcorr')"

# cheap precount on the bare (untagged) subtracted PAT jets: raw pT cut at the
# clustering threshold (15), so this only rejects events with < 6 Cs3 jets
process.fullHadPreJets = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("patJetsAKCs3PF"),
    cut = cms.string("correctedJet('Uncorrected').pt >= 15.0 && abs(eta) <= 2.3")
)
process.fullHadPreJetCount = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("fullHadPreJets"),
    minNumber = cms.uint32(fullHadNJets)
)

# jet counting on UParT-regressed pT (mirrors offline jtupartpt = jtrawPt*ptcorr)
process.fullHadJets = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("selectedUpdatedPatJetsAKCs3DeepFlavour"),
    cut = cms.string(f"({_uparTPt}) >= {fullHadJetPtMin} && abs(eta) <= {fullHadJetAbsEta}")
)
process.fullHadJetCount = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("fullHadJets"),
    minNumber = cms.uint32(fullHadNJets)
)

# loose b-tagged jets among the selected jets, normalized UParT b discriminant
# (mirrors offline discr_btag = (b+bb+lepb)/(b+bb+lepb+c+s+u+d+g))
process.fullHadBJets = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("fullHadJets"),
    cut = cms.string(f"({_uparTB}) > 0 && ({_uparTB})/({_uparTAll}) >= {fullHadBDiscrMin}")
)
process.fullHadBJetCount = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("fullHadBJets"),
    minNumber = cms.uint32(fullHadNBJets)
)

_fullHadSeq = process.fullHadJets * process.fullHadJetCount
if requireBTags:
    _fullHadSeq = _fullHadSeq * process.fullHadBJets * process.fullHadBJetCount
if usePreCountFilter:
    _fullHadSeq = process.fullHadPreJets * process.fullHadPreJetCount * _fullHadSeq
process.fullHadSelection = cms.Sequence(_fullHadSeq)

# the b-tagging chain must run before the skim decision, so it lives in its
# own path (with the PAT tasks) whose status gates everything else
process.fullHadSkimPath = cms.Path(process.fullHadSelection, process.patAlgosToolsTask, process.svTask)

process.fullHadSkimGate = cms.EDFilter("PathStatusFilter",
    logicalExpression = cms.string("fullHadSkimPath")
)

process.filterSequence = cms.Sequence(
    process.primaryVertexFilter *
    process.fullHadSkimGate
)

process.superFilterPath = cms.Path(process.filterSequence)
process.skimanalysis.superFilters = cms.vstring("superFilterPath")

selectionPaths = ["superFilterPath", "fullHadSkimPath"]
for path in process.paths:
    if path not in selectionPaths:
        getattr(process, path)._seq = process.filterSequence * getattr(process,path)._seq
