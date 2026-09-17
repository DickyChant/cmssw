# Throughput benchmark: FastJet vs FlashJet (alpaka / SONIC) on synthetic events.
#
#   cmsRun benchmarkFlashJet_cfg.py --workflow ak4 --impl alpaka --backend cuda_async \
#       --nSoft 3000 --maxEvents 2000 --threads 8 --json ak4_alpaka_cuda.json
#
# workflows:
#   ak4       anti-kt R=0.4 of all particles of the event (one clustering problem per event)
#   softdrop  C/A reclustering + soft drop of every AK8 jet (pt > 100) of the event
#             (FlashJet: all jets of the event in one batched kernel call)
# impls: fastjet, alpaka (--backend serial_sync|cuda_async|rocm_async), sonic (ak4 only)
# Timings are written by FastTimerService to --json.
import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.SonicTriton.customize import getDefaultClientPSet, getParser, getOptions, applyOptions, applyClientOptions

parser = getParser()
parser.add_argument("--workflow", default="ak4", choices=["ak4", "softdrop"])
parser.add_argument("--impl", default="fastjet", choices=["fastjet", "alpaka", "sonic"])
parser.add_argument("--backend", default="serial_sync", help="alpaka backend")
parser.add_argument("--mode", default="Async", choices=["Async", "PseudoAsync", "Sync"], help="SONIC client mode")
parser.add_argument("--nSoft", default=1500, type=int)
parser.add_argument("--nJets", default=6, type=int)
parser.add_argument("--nPerJet", default=80, type=int)
parser.add_argument("--entriesPerThread", default=1, type=int)
parser.add_argument("--json", default="benchmark.json", type=str)
options = getOptions(parser, verbose=True)
if options.impl == "sonic" and options.workflow != "ak4":
    raise ValueError("the SONIC producer implements the ak4 workflow only")

modifiers = []
if options.impl == "sonic":
    from Configuration.ProcessModifiers.enableSonicTriton_cff import enableSonicTriton
    modifiers.append(enableSonicTriton)
process = cms.Process("BENCH", *modifiers)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.load("Configuration.StandardSequences.Accelerators_cff")
process.load("HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi")
process.source = cms.Source("EmptySource")

process.load("HLTrigger.Timer.FastTimerService_cfi")
process.FastTimerService.writeJSONSummary = True
process.FastTimerService.jsonFileName = options.json
process.FastTimerService.enableDQM = False
process.ThroughputService = cms.Service(
    "ThroughputService",
    enableDQM=cms.untracked.bool(False),
    printEventSummary=cms.untracked.bool(True),
    eventRange=cms.untracked.uint32(1000000),
    eventResolution=cms.untracked.uint32(100),
)

process.particles = cms.EDProducer(
    "FlashJetRandomCandidateProducer", nSoft=cms.int32(options.nSoft), nJets=cms.int32(options.nJets),
    nPerJet=cms.int32(options.nPerJet)
)
process.path = cms.Path(process.particles)

alpakaPSet = cms.untracked.PSet(backend=cms.untracked.string(options.backend))
sdParams = dict(zcut=cms.double(0.1), beta=cms.double(0.0), R0=cms.double(0.8))

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets as _fastjet

if options.workflow == "ak4":
    params = dict(src=cms.InputTag("particles"), jetAlgorithm=cms.string("AntiKt"), rParam=cms.double(0.4))
    if options.impl == "fastjet":
        process.jets = _fastjet.clone(src="particles", jetType="BasicJet", jetPtMin=5.0)
        process.path += process.jets
    elif options.impl == "alpaka":
        process.clusters = cms.EDProducer("FlashJetProducer@alpaka", entriesPerThread=cms.int32(options.entriesPerThread),
                                          alpaka=alpakaPSet, **params)
        process.jets = cms.EDProducer("FlashJetRecoJetProducer", src=params["src"], clusters=cms.InputTag("clusters"),
                                      jetType=cms.string("BasicJet"), jetPtMin=cms.double(5.0))
        process.path += process.clusters + process.jets
    else:
        process.load("HeterogeneousCore.SonicTriton.TritonService_cff")
        process.jets = cms.EDProducer(
            "FlashJetSonicProducer",
            Client=applyClientOptions(getDefaultClientPSet().clone(), options).clone(
                mode=cms.string(options.mode),
                preferredServer=cms.untracked.string(""),
                modelName=cms.string("flashjet"),
                modelVersion=cms.string(""),
                modelConfigPath=cms.FileInPath("RecoJets/FlashJet/data/models/flashjet/config.pbtxt"),
            ),
            jetType=cms.string("BasicJet"), jetPtMin=cms.double(5.0), inputPtMin=cms.double(0.0), **params,
        )
        process.path += process.jets
else:
    # AK8 jets are an input here, not part of the measurement
    process.ak8 = _fastjet.clone(src="particles", jetType="BasicJet", rParam=0.8, jetPtMin=100.0)
    process.path += process.ak8
    if options.impl == "fastjet":
        process.softdrop = cms.EDProducer("FastjetSoftDropProducer", jets=cms.InputTag("ak8"),
                                          jetPtMin=cms.double(0.0), **sdParams)
        process.path += process.softdrop
    else:
        process.clusters = cms.EDProducer(
            "FlashJetReclusterProducer@alpaka", src=cms.InputTag("ak8"), jetPtMin=cms.double(0.0),
            jetAlgorithm=cms.string("CambridgeAachen"), rParam=cms.double(1000.0),
            entriesPerThread=cms.int32(options.entriesPerThread),
            softDrop=cms.PSet(enable=cms.bool(True), **sdParams), alpaka=alpakaPSet,
        )
        process.softdrop = cms.EDProducer("FlashJetSoftDropProducer", jets=cms.InputTag("ak8"),
                                          clusters=cms.InputTag("clusters"))
        process.path += process.clusters + process.softdrop

process = applyOptions(process, options)
