import FWCore.ParameterSet.Config as cms
import sys

# if len(sys.argv) < 1:
# 	print "usage: python DYJetsToLL_13TeV_cfg.py <mass range>"
# 	print "ex> python DYJetsToLL_13TeV_cfg.py M100to200"
# 	sys.exit()

def FindCrossSectionByFileName( fileName ):
    ####################################################################
    # -- automatically set a cross section for a sample by fileName -- #
    ####################################################################
    SF = 1.033504 # -- to resolve a discontinuity near M=100 GeV found in generator level of aMC@NLO sample
                  # -- details: https://indico.cern.ch/event/555648/contributions/2240001/attachments/1308706/1957318/2016.07.13_ARC_DYDiffXsec_Muon_KPLee_v3_1.pdf
    # -- x-section of each mass-binned samples
    Dic_xSec = {
        "DYJetsToLL_M-10to50":    18610.0       * (1.0/3.0),
        "DYJetsToLL_M-50_":        5869.58346   * (1.0/3.0), # -- 50 < M < 100 GeV cross section (above 100 GeV will be truncated)
        "DYJetsToLL_M-100to200":   226     * SF * (1.0/3.0),
        "DYJetsToLL_M-200to400":   7.67    * SF * (1.0/3.0),
        "DYJetsToLL_M-400to500":   0.423   * SF * (1.0/3.0),
        "DYJetsToLL_M-500to700":   0.24    * SF * (1.0/3.0),
        "DYJetsToLL_M-700to800":   0.035   * SF * (1.0/3.0),
        "DYJetsToLL_M-800to1000":  0.03    * SF * (1.0/3.0),
        "DYJetsToLL_M-1000to1500": 0.016   * SF * (1.0/3.0),
        "DYJetsToLL_M-1500to2000": 0.002   * SF * (1.0/3.0),
        "DYJetsToLL_M-2000to3000": 0.00054 * SF * (1.0/3.0)
    }

    xSecForThis = 0
    for massRange in Dic_xSec.keys():
        if massRange in fileName:
            xSecForThis = Dic_xSec[massRange]
            print "massRange: %s -> x-section: %lf" % ( massRange, xSecForThis )
            break

    if xSecForThis == 0:
        print "no corresponding x-section is found ... need to be checked!"
        sys.exit()

    return xSecForThis;


process = cms.Process("runRivetAnalysis")

# -- https://github.com/cms-sw/cmssw/blob/master/FWCore/MessageService/python/MessageLogger_cfi.py
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    # input = cms.untracked.int32(1000)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

process.source.fileNames = cms.untracked.vstring(
    # -- M10-50 sample
    "file:/pnfs/knu.ac.kr/data/cms/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0005EA25-8CB8-E511-A910-00266CF85DA0.root",
    "file:/pnfs/knu.ac.kr/data/cms/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/004DA7F9-8AB8-E511-9857-848F69FD0BD5.root",
    # "file:/pnfs/knu.ac.kr/data/cms/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/088BCFC6-8BB8-E511-8E24-848F69FD0BD5.root",
    # "file:/pnfs/knu.ac.kr/data/cms/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/08D83F38-8FB8-E511-A3BD-00266CF85DA0.root",

    # -- M50 sample
    # 'file:/pnfs/knu.ac.kr/data/cms/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/06532BBC-05C8-E511-A60A-F46D043B3CE5.root'

    # -- M100to200 sample
    # "file:/pnfs/knu.ac.kr/data/cms/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/60000/027D7153-29BF-E511-A2BC-0025B3E025B6.root",
    # "file:/pnfs/knu.ac.kr/data/cms/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/60000/06663B51-29BF-E511-8B11-002481DE4A28.root",
)

process.genDYTauTauFilter = cms.EDFilter("GenDYTauTauFilter",
    GenParticle = cms.untracked.InputTag("prunedGenParticles")
    )

# -- only for M50 sample
process.genDileptonMassSelector = cms.EDFilter( "GenDiLeptonMassSelector", 
    lowerEdge = cms.untracked.double(50.0),
    upperEdge = cms.untracked.double(100.0),
    GenParticle = cms.untracked.InputTag("prunedGenParticles"),
    )

process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
    inputPruned = cms.InputTag("prunedGenParticles"),
    inputPacked = cms.InputTag("packedGenParticles"),
)

process.generator = cms.EDProducer("GenParticles2HepMCConverter",
    # genParticles = cms.InputTag("genParticles"),
    genParticles = cms.InputTag("mergedGenParticles"),
    genEventInfo = cms.InputTag("generator", "", "SIM"),
    signalParticlePdgIds = cms.vint32(),
)

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

xSec = FindCrossSectionByFileName( process.source.fileNames[0] )

process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2018_I1711625')
process.rivetAnalyzer.CrossSection = cms.double(xSec)
process.rivetAnalyzer.OutputFile = cms.string('CMS_2018_I1711625_MC.yoda')
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")

common = cms.Sequence(
    process.genDYTauTauFilter *
    process.mergedGenParticles * 
    process.generator * 
    process.rivetAnalyzer)

# common = cms.Sequence(
#     process.mergedGenParticles * 
#     process.generator * 
#     process.rivetAnalyzer)

process.finalPath = cms.Path(common)
if "DYJetsToLL_M-50_" in process.source.fileNames[0]:
    process.finalPath = cms.Path(process.genDileptonMassSelector * common)
    
print "final path: ", process.finalPath
