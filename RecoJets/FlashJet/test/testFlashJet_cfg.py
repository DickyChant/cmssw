# Validates FlashJet against FastjetJetProducer on generator-level jets.
#
#   cmsDriver.py TTbar_14TeV_TuneCP5_cfi -s GEN --conditions auto:phase1_2024_realistic \
#       --era Run3_2024 --beamspot Realistic25ns13p6TeVEarly2023Collision \
#       --datatier GEN --eventcontent RAWSIM -n 20 --fileout file:gen.root
#   cmsRun testFlashJet_cfg.py --inputFiles file:gen.root --alpaka serial_sync --sonic
#
# --alpaka BACKEND  run FlashJetProducer@alpaka (serial_sync, cuda_async, rocm_async, or "" to skip)
# --sonic           also run FlashJetSonicProducer (starts a local fallback Triton server
#                   unless --address is given; run test/setupFlashJetModel.sh first)
import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.SonicTriton.customize import getDefaultClientPSet, getParser, getOptions, applyOptions, applyClientOptions

parser = getParser()
parser.add_argument("--inputFiles", default=["file:gen.root"], nargs="+", type=str)
parser.add_argument("--alpaka", default="serial_sync", type=str, help="alpaka backend, empty to skip")
parser.add_argument("--sonic", default=False, action="store_true", help="run the SONIC producer")
parser.add_argument("--mode", default="PseudoAsync", choices=["Async", "PseudoAsync", "Sync"], help="SONIC client mode")
parser.add_argument("--jetAlgorithm", default="AntiKt", choices=["AntiKt", "Kt", "CambridgeAachen"])
parser.add_argument("--rParam", default=0.4, type=float)
parser.add_argument("--failOnMismatch", default=False, action="store_true")
options = getOptions(parser, verbose=True)

modifiers = []
if options.sonic:
    from Configuration.ProcessModifiers.enableSonicTriton_cff import enableSonicTriton
    modifiers.append(enableSonicTriton)
process = cms.Process("FLASHJET", *modifiers)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.load("Configuration.StandardSequences.Accelerators_cff")
process.load("HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi")

process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))

process.task = cms.Task()
process.path = cms.Path(process.task)

# reference: FastJet
process.load("RecoJets.Configuration.GenJetParticles_cff")
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets

process.fastjetGenJets = ak4GenJets.clone(
    src="genParticlesForJetsNoNu", jetAlgorithm=options.jetAlgorithm, rParam=options.rParam
)
process.task.add(process.genParticlesForJetsNoNu, process.fastjetGenJets)

common = dict(
    src=cms.InputTag("genParticlesForJetsNoNu"),
    jetAlgorithm=cms.string(options.jetAlgorithm),
    rParam=cms.double(options.rParam),
)
jetParams = dict(jetType=cms.string("GenJet"), jetPtMin=process.fastjetGenJets.jetPtMin)

if options.alpaka:
    process.flashJetClusters = cms.EDProducer(
        "FlashJetProducer@alpaka",
        inputPtMin=cms.double(0.0),
        alpaka=cms.untracked.PSet(backend=cms.untracked.string(options.alpaka)),
        **common,
    )
    process.flashJetGenJets = cms.EDProducer(
        "FlashJetRecoJetProducer",
        src=common["src"],
        clusters=cms.InputTag("flashJetClusters"),
        **jetParams,
    )
    process.compareAlpaka = cms.EDAnalyzer(
        "FlashJetCompareAnalyzer",
        reference=cms.InputTag("fastjetGenJets"),
        test=cms.InputTag("flashJetGenJets"),
        ptMin=process.fastjetGenJets.jetPtMin,
        tolerance=cms.double(1e-9),
        failOnMismatch=cms.bool(options.failOnMismatch),
    )
    process.task.add(process.flashJetClusters, process.flashJetGenJets)
    process.path += process.compareAlpaka

if options.sonic:
    process.load("HeterogeneousCore.SonicTriton.TritonService_cff")
    process.flashJetSonicGenJets = cms.EDProducer(
        "FlashJetSonicProducer",
        Client=applyClientOptions(getDefaultClientPSet().clone(), options).clone(
            mode=cms.string(options.mode),
            preferredServer=cms.untracked.string(""),
            modelName=cms.string("flashjet"),
            modelVersion=cms.string(""),
            modelConfigPath=cms.FileInPath("RecoJets/FlashJet/data/models/flashjet/config.pbtxt"),
        ),
        inputPtMin=cms.double(0.0),
        **common,
        **jetParams,
    )
    process.compareSonic = cms.EDAnalyzer(
        "FlashJetCompareAnalyzer",
        reference=cms.InputTag("fastjetGenJets"),
        test=cms.InputTag("flashJetSonicGenJets"),
        ptMin=process.fastjetGenJets.jetPtMin,
        # the SONIC path sums constituents in input order, not merge order
        tolerance=cms.double(1e-9),
        failOnMismatch=cms.bool(options.failOnMismatch),
    )
    process.task.add(process.flashJetSonicGenJets)
    process.path += process.compareSonic

process = applyOptions(process, options)
