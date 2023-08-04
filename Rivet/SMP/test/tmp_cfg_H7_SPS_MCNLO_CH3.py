import FWCore.ParameterSet.Config as cms

process = cms.Process("Gen")

process.load("Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff")
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')

# Input source
process.source = cms.Source("EmptySource")
# Production Info
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Configuration/GenProduction/python/SMP-RunIIFall17wmLHEGS-00175-fragment.py nevts:16'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

process.externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    args = cms.vstring('/cvmfs/cms.cern.ch/phys_generator/gridpacks/slc6_amd64_gcc481/13TeV/madgraph/V5_2.2.2/dyellell012j_5f_NLO_FXFX/v1/dyellell012j_5f_NLO_FXFX_tarball.tar.xz'),
    nEvents = cms.untracked.uint32(10000),
    numberOfParameters = cms.uint32(1),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')
)

process.generator = cms.EDFilter("Herwig7GeneratorFilter",
    configFiles = cms.vstring(),
    crossSection = cms.untracked.double(-1),
    dataLocation = cms.string('${HERWIGPATH:-6}'),
    eventHandlers = cms.string('/Herwig/EventHandlers'),
    filterEfficiency = cms.untracked.double(1.0),
    generatorModule = cms.string('/Herwig/Generators/EventGenerator'),
    herwig7CH3AlphaS = cms.vstring('cd /Herwig/Shower',
        'set AlphaQCD:AlphaMZ 0.118',
        'cd /'),
    herwig7CH3MPISettings = cms.vstring('read snippets/SoftModel.in',
        'set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 0.4712',
        'set /Herwig/UnderlyingEvent/MPIHandler:pTmin0 3.04',
        'set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 1.284',
        'set /Herwig/UnderlyingEvent/MPIHandler:Power 0.1362',
        'set /Herwig/Partons/RemnantDecayer:ladderPower -0.08',
        'set /Herwig/Partons/RemnantDecayer:ladderNorm 0.95'),
    herwig7CH3PDF = cms.vstring('cd /Herwig/Partons',
        'create ThePEG::LHAPDF PDFSet_nnlo ThePEGLHAPDF.so',
        'set PDFSet_nnlo:PDFName NNPDF31_nnlo_as_0118.LHgrid',
        'set PDFSet_nnlo:RemnantHandler HadronRemnants',
        'set /Herwig/Particles/p+:PDF PDFSet_nnlo',
        'set /Herwig/Particles/pbar-:PDF PDFSet_nnlo',
        'set /Herwig/Partons/PPExtractor:FirstPDF  PDFSet_nnlo',
        'set /Herwig/Partons/PPExtractor:SecondPDF PDFSet_nnlo',
        'set /Herwig/Shower/ShowerHandler:PDFA PDFSet_nnlo',
        'set /Herwig/Shower/ShowerHandler:PDFB PDFSet_nnlo',
        'create ThePEG::LHAPDF PDFSet_lo ThePEGLHAPDF.so',
        'set PDFSet_lo:PDFName NNPDF31_lo_as_0130.LHgrid',
        'set PDFSet_lo:RemnantHandler HadronRemnants',
        'set /Herwig/Shower/ShowerHandler:PDFARemnant PDFSet_lo',
        'set /Herwig/Shower/ShowerHandler:PDFBRemnant PDFSet_lo',
        'set /Herwig/Partons/MPIExtractor:FirstPDF PDFSet_lo',
        'set /Herwig/Partons/MPIExtractor:SecondPDF PDFSet_lo',
        'cd /'),
    herwig7StableParticlesForDetector = cms.vstring('set /Herwig/Decays/DecayHandler:MaxLifeTime 10*mm',
        'set /Herwig/Decays/DecayHandler:LifeTimeOption Average'),
    hw_PSWeights_settings = cms.vstring('cd /',
        'cd /Herwig/Shower',
        'do ShowerHandler:AddVariation RedHighAll 1.141 1.141  All',
        'do ShowerHandler:AddVariation RedLowAll 0.707 0.707 All',
        'do ShowerHandler:AddVariation DefHighAll 2 2 All',
        'do ShowerHandler:AddVariation DefLowAll 0.5 0.5 All',
        'do ShowerHandler:AddVariation ConHighAll 4 4 All',
        'do ShowerHandler:AddVariation ConLowAll 0.25 0.25 All',
        'do ShowerHandler:AddVariation RedHighHard 1.141 1.141  Hard',
        'do ShowerHandler:AddVariation RedLowHard 0.707 0.707 Hard',
        'do ShowerHandler:AddVariation DefHighHard 2 2 Hard',
        'do ShowerHandler:AddVariation DefLowHard 0.5 0.5 Hard',
        'do ShowerHandler:AddVariation ConHighHard 4 4 Hard',
        'do ShowerHandler:AddVariation ConLowHard 0.25 0.25 Hard',
        'do ShowerHandler:AddVariation RedHighSecondary 1.141 1.141  Secondary', 
        'do ShowerHandler:AddVariation RedLowSecondary 0.707 0.707 Secondary',
        'do ShowerHandler:AddVariation DefHighSecondary 2 2 Secondary',
        'do ShowerHandler:AddVariation DefLowSecondary 0.5 0.5 Secondary',
        'do ShowerHandler:AddVariation ConHighSecondary 4 4 Secondary',
        'do ShowerHandler:AddVariation ConLowSecondary 0.25 0.25 Secondary',
        'set SplittingGenerator:Detuning 2.0',
        'cd /'),
    hw_common_merging_settings = cms.vstring('read snippets/PPCollider.in',
        'cd /Herwig/EventHandlers',
        'library FxFx.so',
        'create Herwig::FxFxEventHandler theLesHouchesHandler',
        'cd /Herwig/EventHandlers',
        'library FxFx.so',
        'create Herwig::FxFxFileReader theLHReader',
        'cd /Herwig/Shower',
        'library FxFxHandler.so',
        'create Herwig::FxFxFileReader theLHReader',
        'cd /Herwig/Shower',
        'library FxFxHandler.so',
        'create Herwig::FxFxHandler FxFxHandler',
        'set /Herwig/Shower/FxFxHandler:ShowerModel /Herwig/Shower/ShowerModel', 
        'set /Herwig/Shower/FxFxHandler:SplittingGenerator /Herwig/Shower/SplittingGenerator',
        'cd /Herwig/EventHandlers',
        'create ThePEG::Cuts   /Herwig/Cuts/NoCuts',
        'cd /Herwig/EventHandlers',
        'insert theLesHouchesHandler:FxFxReaders[0] theLHReader',
        'set theLesHouchesHandler:WeightOption VarNegWeight',
        'set theLesHouchesHandler:PartonExtractor /Herwig/Partons/PPExtractor',
        'set theLesHouchesHandler:CascadeHandler /Herwig/Shower/FxFxHandler',
        'set theLesHouchesHandler:HadronizationHandler /Herwig/Hadronization/ClusterHadHandler',
        'set theLesHouchesHandler:DecayHandler /Herwig/Decays/DecayHandler',
        'set /Herwig/Shower/ShowerHandler:MaxPtIsMuF Yes',
        'set /Herwig/Shower/ShowerHandler:RestrictPhasespace Yes',
        'set /Herwig/Shower/PartnerFinder:PartnerMethod Random',
        'set /Herwig/Shower/PartnerFinder:ScaleChoice Partner',
        'cd /Herwig/EventHandlers',
        'set theLHReader:AllowedToReOpen No',
        'set theLHReader:WeightWarnings    false',
        'set theLHReader:FileName cmsgrid_final.lhe',
        'set theLHReader:MomentumTreatment      RescaleEnergy',
        'set theLHReader:Cuts  /Herwig/Cuts/NoCuts',
        'cd /Herwig/Generators',
        'set EventGenerator:EventHandler  /Herwig/EventHandlers/theLesHouchesHandler',
        'set EventGenerator:PrintEvent     1',
        'set EventGenerator:MaxErrors      10000',
        'cd /Herwig/Shower',
        'set /Herwig/Shower/FxFxHandler:MPIHandler  /Herwig/UnderlyingEvent/MPIHandler',
        'set /Herwig/Shower/FxFxHandler:RemDecayer  /Herwig/Partons/RemnantDecayer',
        'set /Herwig/Shower/FxFxHandler:ShowerAlpha  AlphaQCD',
        'set FxFxHandler:HeavyQVeto Yes',
        'set FxFxHandler:HardProcessDetection Automatic',
        'set FxFxHandler:ihrd        3',
        'set FxFxHandler:njetsmax      2',
        'set FxFxHandler:drjmin      0',
        'cd /Herwig/Shower',
        'set FxFxHandler:VetoIsTurnedOff VetoingIsOn',
        'set FxFxHandler:MergeMode FxFx',
        'set FxFxHandler:ETClus 20*GeV',
        'set FxFxHandler:RClus 1.0',
        'set FxFxHandler:EtaClusMax 5',
        'set FxFxHandler:RClusFactor 1.5'),
    hw_user_settings = cms.vstring('cd /Herwig/EventHandlers',
        'set EventHandler:LuminosityFunction:Energy 13000*GeV',
        'cd /Herwig/Shower',
        'set FxFxHandler:njetsmax 2',
        'set FxFxHandler:MergeMode FxFx',
        'cd /',
        'set /Herwig/Particles/h0:NominalMass 125.0'),
    parameterSets = cms.vstring('hw_common_merging_settings',
        'herwig7CH3PDF',
        'herwig7CH3AlphaS',
        'herwig7CH3MPISettings',
        'herwig7StableParticlesForDetector',
        'hw_PSWeights_settings',
        'hw_user_settings'),
    repository = cms.string('${HERWIGPATH}/HerwigDefaults.rpo'),
    run = cms.string('InterfaceMatchboxTest'),
    runModeList = cms.untracked.string('read,run')
)


#-------Produce charged jets-----------
from SimGeneral.HepPDTESSource.pythiapdt_cfi import *
#----------------------------------------------



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
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
  OutputFile = cms.string('ee_yodaFileOutput.yoda')
)'''



process.load("Configuration.EventContent.EventContent_cff")
process.p = cms.Path(process.externalLHEProducer * process.generator * process.pgen * process.rivetAnalyzer_zmm)
