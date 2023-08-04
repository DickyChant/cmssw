import FWCore.ParameterSet.Config as cms
process = cms.Process("runRivetAnalysis")
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
		'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/90000/FE8A7852-66E4-E611-B5D0-002590E7E01A.root',
		'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/3A0E5113-28E3-E611-8143-002481DE4668.root',
		'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/60000/F4AD65E7-68F1-E611-B967-02163E012034.root',
		
		#'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/90000/FE6DFC4B-91E5-E611-9EA2-0025905A6076.root',
		#'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/3A080C98-8AF1-E611-9F59-02163E014479.root',
		#'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/60000/FA52C61E-5DF1-E611-B54B-02163E011E3F.root',
		
		#'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/90000/FACC175D-A5E5-E611-940F-70106F4A44AC.root',
        #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/38901FD9-22E3-E611-ABE6-001517F807D4.root',
        #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/60000/F60874D5-25E5-E611-A69C-008CFA06488C.root',
        
        #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/90000/FA68D6E8-7DE5-E611-9158-0CC47A4C8EC6.root',
        #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/387B3E03-A8E3-E611-8B10-FA163EE71BDA.root',
        #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/60000/F4ECBA3C-46E4-E611-8301-FA163E448773.root',
        
        #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/90000/FA48B1F1-0DF1-E611-A5C1-02163E011CDB.root'
        #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/38732585-19E4-E611-88EF-FA163EB96B44.root'
        #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/60000/FE373B97-72F1-E611-A4B9-02163E01A312.root'
        #'/store/mc/RunIISummer20UL16MiniAODAPV/DYJetsToEE_M-50_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v8-v2/00000/00CCAEA8-C03A-C148-82C3-425B3AD2EFC5.root'
        
        #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/90000/F8E67C19-C7E4-E611-9333-F832E4CC4D39.root',
        #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/90000/F81F42AB-9FE5-E611-8BB1-0025905B8572.root',
        #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/90000/F6F73F10-13E4-E611-A8A2-0025901AC12A.root',
        #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/90000/F6E03A52-4CE4-E611-AA96-BC305B390A9A.root',
        #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/90000/F6D062B0-14F1-E611-B7BA-02163E0141C2.root'
        #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/90000/FE8A7852-66E4-E611-B5D0-002590E7E01A.root'
        
    )
)
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.generatorMix = cms.EDProducer("MergedGenParticleProducer",
    inputPacked = cms.InputTag("packedGenParticles"),
    inputPruned = cms.InputTag("prunedGenParticles")
)

process.generator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("generatorMix"),
    genEventInfo = cms.InputTag("generator"),
    lheSrc = cms.InputTag("externalLHEProducer"),
    signalParticlePdgIds = cms.vint32()
)

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")
process.rivetAnalyzer.HepMCCollection = cms.InputTag('generator','unsmeared')
process.rivetAnalyzer.CrossSection = cms.double(5930.9)
process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2022_I2079374')
process.p = cms.Path(process.generatorMix*
                     process.generator*process.rivetAnalyzer)
process.rivetAnalyzer.useLHEweights = cms.bool(True)
process.rivetAnalyzer.matchWeightNames=cms.string("lhapdf=306.*")
process.rivetAnalyzer.OutputFile = cms.string('CMS_2022_I2079374.out.yoda')
