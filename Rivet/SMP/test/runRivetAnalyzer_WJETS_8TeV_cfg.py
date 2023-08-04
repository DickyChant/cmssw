import FWCore.ParameterSet.Config as cms

process = cms.Process("runRivetAnalysis")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")



process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
									  #'file:gen.root'
									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/0022F532-C2CE-E111-B33F-003048D476A6.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00929F05-20CF-E111-AC61-003048D45FEC.root'
									  
)
)

#process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cfi")

process.generator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("genParticles"),
    genEventInfo = cms.InputTag("generator", "", "SIM"),
)

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2016_I1491953')

process.rivetAnalyzer.OutputFile= cms.string('wjets.yoda')

process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")

process.p = cms.Path(process.generator*process.rivetAnalyzer)

#process.p = cms.Path(process.rivetAnalyzer)


