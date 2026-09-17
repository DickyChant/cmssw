# Throughput benchmark: FastJet vs FlashJet (alpaka / SONIC) on synthetic or real events.
#
#   cmsRun benchmarkFlashJet_cfg.py --workflow ak4 --impl alpaka --backend cuda_async \
#       --nSoft 3000 --maxEvents 2000 --threads 8 --json ak4_alpaka_cuda.json
#
# workflows:
#   ak4       anti-kt R=0.4 of all particles of the event (one clustering problem per event)
#   softdrop  C/A reclustering + soft drop of every AK8 jet (pt > 100) of the event
#             (FlashJet: all jets of the event in one batched kernel call)
# impls: fastjet, alpaka (--backend serial_sync|cuda_async|rocm_async), sonic (ak4 only)
# inputs:
#   synthetic  FlashJetRandomCandidateProducer (--nSoft, --nJets, --nPerJet)
#   miniaod    packedPFCandidates; softdrop on slimmedJetsAK8 (--inputFiles)
#   scouting   Run-3 PF scouting particles (hltScoutingPFPacker) converted to reco::PFCandidate;
#              softdrop on FastJet AK8 jets built from them (--inputFiles)
# --compare also runs the FastJet reference and checks FlashJet against it (validation, not timing)
# Timings are written by FastTimerService to --json.
import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.SonicTriton.customize import getDefaultClientPSet, getParser, getOptions, applyOptions, applyClientOptions

parser = getParser()
parser.add_argument("--workflow", default="ak4", choices=["ak4", "softdrop"])
parser.add_argument("--impl", default="fastjet", choices=["fastjet", "alpaka", "sonic"])
parser.add_argument("--backend", default="serial_sync", help="alpaka backend")
parser.add_argument("--mode", default="Async", choices=["Async", "PseudoAsync", "Sync"], help="SONIC client mode")
parser.add_argument("--input", default="synthetic", choices=["synthetic", "miniaod", "scouting"])
parser.add_argument("--inputFiles", default=[], nargs="+", type=str)
parser.add_argument("--compare", default=False, action="store_true")
parser.add_argument("--nSoft", default=1500, type=int)
parser.add_argument("--nJets", default=6, type=int)
parser.add_argument("--nPerJet", default=80, type=int)
parser.add_argument("--entriesPerThread", default=1, type=int)
parser.add_argument("--json", default="benchmark.json", type=str)
parser.add_argument("--warmup", default=-1, type=int,
                    help="events excluded from the throughput measurement (default: 10%% of --maxEvents)")
parser.add_argument("--resolution", default=0, type=int,
                    help="ThroughputService sampling in events (default: about 50 samples)")
options = getOptions(parser, verbose=True)
if options.input != "synthetic" and not options.inputFiles:
    raise ValueError("--input miniaod/scouting needs --inputFiles")
if options.compare and options.impl == "fastjet":
    raise ValueError("--compare needs a FlashJet implementation (alpaka or sonic)")
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
if options.input == "synthetic":
    process.source = cms.Source("EmptySource")
else:
    process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
    if options.input == "scouting":
        # read only the scouting PF candidates: raw HLTSCOUT files carry products (e.g. CTPPS digis)
        # whose classes may be gone from newer releases
        process.source.inputCommands = cms.untracked.vstring("drop *", "keep *_hltScoutingPFPacker_*_*")
        process.source.dropDescendantsOfDroppedBranches = cms.untracked.bool(False)

process.load("HLTrigger.Timer.FastTimerService_cfi")
process.FastTimerService.writeJSONSummary = True
process.FastTimerService.jsonFileName = options.json
process.FastTimerService.enableDQM = False
# throughput: wall-clock event rate after a warm-up (GPU context, JIT, caching allocators)
nEvents = max(options.maxEvents, 1)
warmup = options.warmup if options.warmup >= 0 else nEvents // 10
resolution = options.resolution if options.resolution > 0 else max(1, (nEvents - warmup) // 50)
process.ThroughputService = cms.Service(
    "ThroughputService",
    enableDQM=cms.untracked.bool(False),
    printEventSummary=cms.untracked.bool(False),
    eventRange=cms.untracked.uint32(nEvents + 1),
    eventResolution=cms.untracked.uint32(resolution),
    eventSkip=cms.untracked.uint32(warmup),
)
# the summary ("Average throughput: ...") is a LogInfo
process.MessageLogger.cerr.ThroughputService = cms.untracked.PSet(limit=cms.untracked.int32(100))
process.MessageLogger.ThroughputService = dict()

# the input particles ("particles") and AK8 input jets ("ak8") are not part of the measurement
if options.input == "synthetic":
    process.particles = cms.EDProducer(
        "FlashJetRandomCandidateProducer", nSoft=cms.int32(options.nSoft), nJets=cms.int32(options.nJets),
        nPerJet=cms.int32(options.nPerJet)
    )
    process.path = cms.Path(process.particles)
    particles = "particles"
elif options.input == "miniaod":
    process.path = cms.Path()
    particles = "packedPFCandidates"
else:
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")  # the converter looks up particle masses
    process.particles = cms.EDProducer(
        "Run3ScoutingParticleToRecoPFCandidateProducer",
        scoutingparticle=cms.InputTag("hltScoutingPFPacker"), softKiller=cms.bool(False), CHS=cms.bool(False)
    )
    process.path = cms.Path(process.particles)
    particles = "particles"

alpakaPSet = cms.untracked.PSet(backend=cms.untracked.string(options.backend))
sdParams = dict(zcut=cms.double(0.1), beta=cms.double(0.0), R0=cms.double(0.8))

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets as _fastjet

if options.workflow == "ak4":
    params = dict(src=cms.InputTag(particles), jetAlgorithm=cms.string("AntiKt"), rParam=cms.double(0.4))
    if options.impl == "fastjet" or options.compare:
        process.fastjetJets = _fastjet.clone(src=particles, jetType="BasicJet", jetPtMin=5.0)
        process.path += process.fastjetJets
    if options.impl == "fastjet":
        pass
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
    if options.input == "miniaod":
        ak8 = "slimmedJetsAK8"
    else:
        process.ak8 = _fastjet.clone(src=particles, jetType="BasicJet", rParam=0.8, jetPtMin=100.0)
        process.path += process.ak8
        ak8 = "ak8"
    if options.impl == "fastjet" or options.compare:
        process.fastjetSoftdrop = cms.EDProducer("FastjetSoftDropProducer", jets=cms.InputTag(ak8),
                                                 jetPtMin=cms.double(0.0), **sdParams)
        process.path += process.fastjetSoftdrop
    if options.impl != "fastjet":
        process.clusters = cms.EDProducer(
            "FlashJetReclusterProducer@alpaka", src=cms.InputTag(ak8), jetPtMin=cms.double(0.0),
            jetAlgorithm=cms.string("CambridgeAachen"), rParam=cms.double(1000.0),
            entriesPerThread=cms.int32(options.entriesPerThread),
            softDrop=cms.PSet(enable=cms.bool(True), **sdParams), alpaka=alpakaPSet,
        )
        process.softdrop = cms.EDProducer("FlashJetSoftDropProducer", jets=cms.InputTag(ak8),
                                          clusters=cms.InputTag("clusters"))
        process.path += process.clusters + process.softdrop
        if options.compare:
            process.compare = cms.EDAnalyzer("FlashJetValueMapCompareAnalyzer", jets=cms.InputTag(ak8),
                                             reference=cms.string("fastjetSoftdrop"), test=cms.string("softdrop"))
            process.path += process.compare

if options.compare and options.workflow == "ak4":
    process.compare = cms.EDAnalyzer("FlashJetCompareAnalyzer", reference=cms.InputTag("fastjetJets"),
                                     test=cms.InputTag("jets"), ptMin=cms.double(5.0), tolerance=cms.double(1e-9))
    process.path += process.compare

process = applyOptions(process, options)

# metadata for summarizeBenchmarks.py, next to the timing JSON
import json as _json
with open(options.json.removesuffix(".json") + ".meta.json", "w") as _meta:
    _json.dump(dict(workflow=options.workflow, impl=options.impl, input=options.input,
                    backend=options.backend if options.impl == "alpaka" else
                    (f"sonic-{options.mode}" if options.impl == "sonic" else options.impl),
                    nSoft=options.nSoft, nJets=options.nJets, nPerJet=options.nPerJet,
                    threads=process.options.numberOfThreads.value(),
                    streams=process.options.numberOfStreams.value() or process.options.numberOfThreads.value(),
                    events=options.maxEvents, warmup=warmup), _meta, indent=1)
