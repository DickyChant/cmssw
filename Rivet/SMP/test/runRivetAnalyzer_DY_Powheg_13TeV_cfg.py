import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('standard')
options.register('runOnly', '', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Run only specified analysis")
options.register('yodafile', 'test.yoda', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Name of yoda output file")
options.setDefault('maxEvents', 10000)
if(hasattr(sys, "argv")):
    options.parseArguments()
print(options)

process = cms.Process("runRivetAnalysis")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource")
process.source.fileNames = cms.untracked.vstring(
    # '/store/mc/RunIISummer16MiniAODv3/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/240000/FA67CC2B-1190-E911-B11F-34E6D7E3879B.root'
    '/store/mc/RunIISummer16MiniAODv3/ZToEE_NNPDF30_13TeV-powheg_M_50_120/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/270000/6AF659B2-0070-E911-8D06-3464A9B95FF0.root'
)

process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
    inputPruned = cms.InputTag("prunedGenParticles"),
    inputPacked = cms.InputTag("packedGenParticles"),
)
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.generator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("mergedGenParticles"),
    genEventInfo = cms.InputTag("generator", "", "SIM"),
    signalParticlePdgIds = cms.vint32()
)

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

if options.runOnly:
    process.rivetAnalyzer.AnalysisNames = cms.vstring(options.runOnly)
else:
    process.rivetAnalyzer.AnalysisNames = cms.vstring(
        'CMS_2018_I1711625', 'CMS_2019_I1753680:LMODE=EL', 'CMS_2018_I1667854:LMODE=EL', 'MC_XS',
    )
process.rivetAnalyzer.OutputFile      = options.yodafile
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.CrossSection = 6077.22/3

process.p = cms.Path(process.mergedGenParticles*process.generator*process.rivetAnalyzer)


