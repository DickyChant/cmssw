import FWCore.ParameterSet.Config as cms

process = cms.Process("Gen")

process.load("Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff")
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Generator_cff')
#process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic25ns13TeV2016Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8aMCatNLOSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *
#process.source = cms.Source("EmptySource")
'''process.source = cms.Source("LHESource",
        fileNames = cms.untracked.vstring(
               '/store/mc/RunIIWinter15wmLHE/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/LHE/MCRUN2_71_V1-v1/50000/E678D56B-65DC-E411-8033-3417EBE644B3.root'

         )
)'''

process.source = cms.Source( "PoolSource",
                             fileNames = cms.untracked.vstring(
			     '/store/mc/RunIIWinter15wmLHE/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/LHE/MCRUN2_71_V1_ext1-v1/00000/0024C5EC-3391-E511-BC7F-02163E013E39.root'
#                             'rootFileInput'
                             )
)

process.options = cms.untracked.PSet(

)

#from genAnalysis.analysis.CFIFILE import *
#from GeneratorInterface.ExternalDecays.TauolaSettings_cff import *


#-------Produce charged jets-----------
from SimGeneral.HepPDTESSource.pythiapdt_cfi import *
from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJetsNoMuNoNu
#from RecoJets.JetProducers.GenJetParameters_cfi import *

#----------------------------------------------

process.generator = cms.EDFilter("Pythia8HadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        pythia8aMCatNLOSettingsBlock,
        processParameters = cms.vstring(
            'JetMatching:setMad = off',
            'JetMatching:scheme = 1',
            'JetMatching:merge = on',
            'JetMatching:jetAlgorithm = 2',
            'JetMatching:etaJetMax = 999.',
            'JetMatching:coneRadius = 1.',
            'JetMatching:slowJetPower = 1',
            'JetMatching:qCut = 30.', #this is the actual merging scale
            'JetMatching:doFxFx = on',
            'JetMatching:qCutME = 10.',#this must match the ptj cut in the lhe generation step
            'JetMatching:nQmatch = 5', #4 corresponds to 4-flavour scheme (no matching of b-quarks), 5 for 5-flavour scheme
            'JetMatching:nJetMax = 2', #number of partons in born matrix element for highest multiplicity
            'TimeShower:mMaxGamma = 4.0',
        ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CP5Settings',
                                    'pythia8aMCatNLOSettings',
                                    'processParameters',
                                    )
    )
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

process.rivetAnalyzer_zmm = cms.EDAnalyzer('RivetAnalyzer',
  AnalysisNames = cms.vstring('CMS_2021_I1866118'),
  HepMCCollection = cms.InputTag('generatorSmeared'),
  genLumiInfo = cms.InputTag("generator"),
  weightCap = cms.double(0.),
  NLOSmearing = cms.double(0.),
  skipMultiWeights = cms.bool(False),
  selectMultiWeights = cms.string(''),
  deselectMultiWeights = cms.string(''),
  setNominalWeightName = cms.string(''),
  UseExternalWeight = cms.bool(False),
  GenEventInfoCollection = cms.InputTag('generator'),
  useLHEweights = cms.bool(False),
  LHECollection = cms.InputTag('externalLHEProducer'),
  CrossSection = cms.double(-1),
  DoFinalize = cms.bool(True),
  ProduceDQMOutput = cms.bool(False),
  OutputFile = cms.string('test1.yoda')
)

'''process.rivetAnalyzer_zee = cms.EDAnalyzer('RivetAnalyzer',
  AnalysisNames = cms.vstring('CMS_2020_PAS_SMP_20_YYYY'),
  HepMCCollection = cms.InputTag('generatorSmeared'),
  UseExternalWeight = cms.bool(False),
  GenEventInfoCollection = cms.InputTag('generator'),
  CrossSection = cms.double(-1),
  DoFinalize = cms.bool(True),
  ProduceDQMOutput = cms.bool(False),
  OutputFile = cms.string('ee_yodaFileOutput')
)'''



process.load("Configuration.EventContent.EventContent_cff")
process.p = cms.Path(process.generator * process.pgen * process.genParticlesForJetsNoMuNoNu * process.ak4GenJetsNoMuNoNu * process.rivetAnalyzer_zmm)
