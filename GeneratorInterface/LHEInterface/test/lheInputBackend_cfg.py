import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing()
options.register("inputFiles", [], VarParsing.multiplicity.list, VarParsing.varType.string, "Input files")
options.register("outputFile", "lhe.root", VarParsing.multiplicity.singleton, VarParsing.varType.string, "Output file")
options.register("maxEvents", -1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Event limit")
options.register("encoding", "xml", VarParsing.multiplicity.singleton, VarParsing.varType.string, "xml or hdf5")
options.register("skip", 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Events to skip")
options.register("allowLoss", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Allow metadata loss")
options.register("shower", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Run a Pythia smoke test")
options.parseArguments()

process = cms.Process("LHEINPUTTEST")
process.source = cms.Source(
    "LHESource",
    fileNames=cms.untracked.vstring(options.inputFiles),
    inputFormat=cms.untracked.string(options.encoding),
    skipEvents=cms.untracked.uint32(options.skip),
    hdf5AllowUnsupportedMetadata=cms.untracked.bool(options.allowLoss),
)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.options = cms.untracked.PSet(numberOfThreads=cms.untracked.uint32(1), numberOfStreams=cms.untracked.uint32(1))
process.output = cms.OutputModule("PoolOutputModule", fileName=cms.untracked.string(options.outputFile))
process.out = cms.EndPath(process.output)
if options.shower:
    process.RandomNumberGeneratorService = cms.Service(
        "RandomNumberGeneratorService", generator=cms.PSet(initialSeed=cms.untracked.uint32(12345))
    )
    process.generator = cms.EDFilter(
        "Pythia8HadronizerFilter",
        pythiaHepMCVerbosity=cms.untracked.bool(False),
        pythiaPylistVerbosity=cms.untracked.int32(0),
        maxEventsToPrint=cms.untracked.int32(0),
        filterEfficiency=cms.untracked.double(1.0),
        comEnergy=cms.double(91.0),
        LHEInputTag=cms.InputTag("source"),
        PythiaParameters=cms.PSet(
            processParameters=cms.vstring("PartonLevel:MPI = off"),
            parameterSets=cms.vstring("processParameters"),
        ),
    )
    process.generation = cms.Path(process.generator)
