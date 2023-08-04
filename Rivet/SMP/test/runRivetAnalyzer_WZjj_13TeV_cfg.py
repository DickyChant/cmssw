import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('standard')
options.register('runOnly', '', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Run only specified analysis")
options.register('yodafile', 'wzjj.yoda', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Name of yoda output file")
options.setDefault('maxEvents', 100000)
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
    '/store/mc/RunIIFall17MiniAODv2/WLLJJ_WToLNu_EWK_TuneCP5_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/00000/3C6BB476-3450-E811-876A-E0071B7A8570.root',
    '/store/mc/RunIIFall17MiniAODv2/WLLJJ_WToLNu_EWK_TuneCP5_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/00000/DA97AC81-764F-E811-BA47-0025905A607E.root',
    '/store/mc/RunIIFall17MiniAODv2/WLLJJ_WToLNu_EWK_TuneCP5_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/00000/E42AF3F8-C04F-E811-8C77-0025905A60CE.root',
    '/store/mc/RunIIFall17MiniAODv2/WLLJJ_WToLNu_EWK_TuneCP5_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/00000/4095EC31-FA4F-E811-9333-0CC47A4D7606.root',
    '/store/mc/RunIIFall17MiniAODv2/WLLJJ_WToLNu_EWK_TuneCP5_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/00000/B2C057B5-6350-E811-A96C-0CC47A7C340E.root',
    '/store/mc/RunIIFall17MiniAODv2/WLLJJ_WToLNu_EWK_TuneCP5_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/00000/36E8520C-7A50-E811-A8CF-0CC47A7C35D2.root',
    '/store/mc/RunIIFall17MiniAODv2/WLLJJ_WToLNu_EWK_TuneCP5_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/00000/58FE69C8-CD50-E811-A39D-0CC47A7C3420.root',
    '/store/mc/RunIIFall17MiniAODv2/WLLJJ_WToLNu_EWK_TuneCP5_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/00000/7ACE9AA5-3151-E811-836F-0CC47A4C8E26.root',
    '/store/mc/RunIIFall17MiniAODv2/WLLJJ_WToLNu_EWK_TuneCP5_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/00000/F88822FF-9451-E811-8A17-00248C55CC9D.root',
    '/store/mc/RunIIFall17MiniAODv2/WLLJJ_WToLNu_EWK_TuneCP5_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/00000/02A0739C-D155-E811-AF6D-00000086FE80.root',
    '/store/mc/RunIIFall17MiniAODv2/WLLJJ_WToLNu_EWK_TuneCP5_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/00000/6C6110BB-D858-E811-A93F-A0369FE2C0E2.root',
    '/store/mc/RunIIFall17MiniAODv2/WLLJJ_WToLNu_EWK_TuneCP5_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/00000/843FC8A2-D958-E811-9102-7845C4FBB764.root',
    '/store/mc/RunIIFall17MiniAODv2/WLLJJ_WToLNu_EWK_TuneCP5_13TeV_madgraph-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/00000/5A7A1622-D858-E811-AC19-C4346BC08440.root'
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
        'CMS_2020_I1794169:LMODE=WZ', 'MC_WEIGHTS', 'MC_XS',
    )
process.rivetAnalyzer.OutputFile      = options.yodafile
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.CrossSection = 0.0176340

process.p = cms.Path(process.mergedGenParticles*process.generator*process.rivetAnalyzer)


