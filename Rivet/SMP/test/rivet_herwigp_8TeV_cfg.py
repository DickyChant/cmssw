# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/2760GeV/MinBias_TuneEE3C_2760GeV_herwigpp_cff.py -s GEN --datatier=GEN-SIM-RAW --conditions auto:mc --eventcontent RAWSIM --no_exec -n 500000 --python_filename=rivet_herwigpp_TuneEE3C_cfg.py --customise=Configuration/GenProduction/rivet_customize.py
import FWCore.ParameterSet.Config as cms
import os

import random
from random import SystemRandom
_MAXINT = 900000000
_inst = SystemRandom()
SEED = _inst.randint(1, _MAXINT)
#SEED=RPLME_SEED
print 'SEED = ',str(SEED)

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
)

# Input source
process.source = cms.Source("EmptySource")

process.MessageLogger.cerr.FwkReport.reportEvery = 100000

process.options = cms.untracked.PSet(

)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('MinBias_TuneEE5C_2760GeV_herwigpp_cff_py_GEN.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)


# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.RandomNumberGeneratorService.generator.initialSeed = 25554

process.generator = cms.EDFilter("ThePEGGeneratorFilter",

    dataLocation = cms.string('${HERWIGPATH}'),
    run = cms.string('LHC'),
    repository = cms.string('HerwigDefaults.rpo'),
    generatorModule = cms.string('/Herwig/Generators/LHCGenerator'),
    configFiles = cms.vstring(),
    crossSection = cms.untracked.double(-1),
    eventHandlers = cms.string('/Herwig/EventHandlers'),

    ##################################################
    # Center-of-mass energy
    ##################################################
    set8TeV  = cms.vstring(
        'cd /Herwig/Generators',
        'set LHCGenerator:EventHandler:LuminosityFunction:Energy 8000.0',
        'cd /'),
    set13TeV  = cms.vstring(
        'cd /Herwig/Generators',
        'set LHCGenerator:EventHandler:LuminosityFunction:Energy 13000.0',
        'cd /'),
    set1960GeV  = cms.vstring(
        'cd /Herwig/Generators',
        'set LHCGenerator:EventHandler:LuminosityFunction:Energy 1960.0',
        'cd /'),
    set900GeV  = cms.vstring(
        'cd /Herwig/Generators',
        'set LHCGenerator:EventHandler:LuminosityFunction:Energy 900.0',
        'cd /'),
    set300GeV  = cms.vstring(
        'cd /Herwig/Generators',
        'set LHCGenerator:EventHandler:LuminosityFunction:Energy 300.0',
        'cd /'),

    ##################################################
    # Technical parameters for this run
    ##################################################
    basicSetup = cms.vstring(
        'cd /Herwig/Generators',
        'create ThePEG::RandomEngineGlue /Herwig/RandomGlue', 
        'set LHCGenerator:RandomNumberGenerator /Herwig/RandomGlue', 
        'set LHCGenerator:NumberOfEvents 10000000',
        'set LHCGenerator:PrintEvent 10',
        'set LHCGenerator:MaxErrors 1000000',
        'set LHCGenerator:DebugLevel 0',
        'set LHCGenerator:DumpPeriod -1',
        'set LHCGenerator:DebugEvent 0',
        'cd /Herwig/Particles', 
        'set p+:PDF /Herwig/Partons/cmsPDFSet', 
        'set pbar-:PDF /Herwig/Partons/cmsPDFSet', 
        'set K0:Width 1e300*GeV', 
        'set Kbar0:Width 1e300*GeV', 
        'cd /'),
         
    ################################################# 
    # Stable Particles (CMS Config)
    #################################################
    setParticlesStableForDetector = cms.vstring(
        'cd /Herwig/Particles', 
        'set mu-:Stable Stable', 
        'set mu+:Stable Stable', 
        'set Sigma-:Stable Stable', 
        'set Sigmabar+:Stable Stable', 
        'set Lambda0:Stable Stable', 
        'set Lambdabar0:Stable Stable', 
        'set Sigma+:Stable Stable', 
        'set Sigmabar-:Stable Stable', 
        'set Xi-:Stable Stable', 
        'set Xibar+:Stable Stable', 
        'set Xi0:Stable Stable', 
        'set Xibar0:Stable Stable', 
        'set Omega-:Stable Stable', 
        'set Omegabar+:Stable Stable', 
        'set pi+:Stable Stable', 
        'set pi-:Stable Stable', 
        'set K+:Stable Stable', 
        'set K-:Stable Stable', 
        'set K_S0:Stable Stable', 
        'set K_L0:Stable Stable', 
        'cd /'),


     ##################################################
     # Running of alpha_s
     ##################################################
     setAlphaS = cms.vstring(
        'cd /Herwig',
        'create Herwig::O2AlphaS O2AlphaS',
        'set Model:QCD/RunningAlphaS O2AlphaS',
        'cd /'),


     ##################################################
     # PDF SET : CTEQ6L1 
     ##################################################

     #pdfCTEQ6L1 = cms.vstring(
     #   'cd /Herwig/Partons',
     #   'create ThePEG::LHAPDF myPDFset ThePEGLHAPDF.so',
     #   'set myPDFset:PDFName cteq6ll.LHpdf',
     #   'set myPDFset:RemnantHandler HadronRemnants',
     #   'set /Herwig/Particles/p+:PDF myPDFset',
     #   'set /Herwig/Particles/pbar-:PDF myPDFset',
     #   'cd /'),

     pdfCTEQ6L1 = cms.vstring(
        'cd /Herwig/Partons', 
        'create ThePEG::LHAPDF CTEQ6L1 ThePEGLHAPDF.so', 
        'set CTEQ6L1:PDFName cteq6ll.LHpdf', 
        'set CTEQ6L1:RemnantHandler HadronRemnants', 
        'cp CTEQ6L1 cmsPDFSet', 
        'cd /'),


     ##################################################
     # MPI energy extrapolation settings
     ##################################################
     EE4CEnergyExtrapol =  cms.vstring(
        'set /Herwig/UnderlyingEvent/MPIHandler:EnergyExtrapolation Power',
        'set /Herwig/UnderlyingEvent/MPIHandler:ReferenceScale 7000.*GeV',
        'set /Herwig/UnderlyingEvent/MPIHandler:Power 0.24',
        #'set /Herwig/UnderlyingEvent/MPIHandler:pTmin0 2.81*GeV',
	'set /Herwig/UnderlyingEvent/MPIHandler:pTmin0 3.91*GeV',
        ),

     EE5CEnergyExtrapol =  cms.vstring(
        'set /Herwig/UnderlyingEvent/MPIHandler:EnergyExtrapolation Power',
        'set /Herwig/UnderlyingEvent/MPIHandler:ReferenceScale 7000.*GeV',
        'set /Herwig/UnderlyingEvent/MPIHandler:Power 3.705288e-01',
	#'set /Herwig/UnderlyingEvent/MPIHandler:Power 0.33',
        'set /Herwig/UnderlyingEvent/MPIHandler:pTmin0 3.91*GeV',
        ),

     ##################################################
     # Energy-independent MPI parameters
     ##################################################
     EE4C =  cms.vstring('+EE4CEnergyExtrapol',
        # Colour reconnection settings
        'set /Herwig/Hadronization/ColourReconnector:ColourReconnection Yes',
        'set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 0.61',
        # Colour Disrupt settings
        'set /Herwig/Partons/RemnantDecayer:colourDisrupt 0.75',
        # inverse hadron radius
        'set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 1.35',
        # MPI model settings
        'set /Herwig/UnderlyingEvent/MPIHandler:softInt Yes',
        'set /Herwig/UnderlyingEvent/MPIHandler:twoComp Yes',
        'set /Herwig/UnderlyingEvent/MPIHandler:DLmode 2',
        ),

      EE5C =  cms.vstring('+EE5CEnergyExtrapol',
        # Colour reconnection settings
        'set /Herwig/Hadronization/ColourReconnector:ColourReconnection Yes',
        #'set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 0.49',
			  'set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 5.278926e-01',
        # Colour Disrupt settings
        #'set /Herwig/Partons/RemnantDecayer:colourDisrupt 0.80', 			  
			  'set /Herwig/Partons/RemnantDecayer:colourDisrupt 6.284222e-01',
        # inverse hadron radius
        #'set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 2.30',
			  'set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 2.254998e+00',	 		  
        # MPI model settings
        'set /Herwig/UnderlyingEvent/MPIHandler:softInt Yes',
        'set /Herwig/UnderlyingEvent/MPIHandler:twoComp Yes',
        'set /Herwig/UnderlyingEvent/MPIHandler:DLmode 2',
        ),


      #####################################
      # FINAL PARAMETERS CHOICE
      #####################################

      MinBias = cms.vstring(

        # Matrix element
        'cd /Herwig/MatrixElements/',
        'insert SimpleQCD:MatrixElements[0] MEMinBias',
        # tell the MPI handler that the primary subprocess is QCD jet production
        'set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0',

        # Cuts for soft-inclusive events
        'cd /Herwig/Cuts',
	'set JetKtCut:MinKT 0.0*GeV', 
	'set QCDCuts:MHatMin 0.0*GeV', 
        'set QCDCuts:X1Min 0.01',
        'set QCDCuts:X2Min 0.01', 
	),


      HardFlat = cms.vstring(
		'mkdir /Herwig/Weights',
		'cd /Herwig/Weights',
		'create ThePEG::ReweightMinPT reweightMinPT ReweightMinPT.so',
		'cd /Herwig/MatrixElements/',
		'insert SimpleQCD:MatrixElements[0] MEQCD2to2',
		'insert SimpleQCD:Preweights[0] /Herwig/Weights/reweightMinPT',
		
		'cd /',
		'set /Herwig/Weights/reweightMinPT:Power 4.5',
		'set /Herwig/Weights/reweightMinPT:Scale 10*GeV',
		'set /Herwig/Cuts/JetKtCut:MaxKT 7000*GeV',
		'set /Herwig/Cuts/QCDCuts:MHatMin 0.0*GeV',
		'set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0',
		#'set /Herwig/EventHandlers/LHCHandler:HadronizationHandler NULL',
		#'set /Herwig/EventHandlers/LHCHandler:DecayHandler NULL',
		#'set /Herwig/EventHandlers/LHCHandler:CascadeHandler NULL',
		#'set /Herwig/Shower/ShowerHandler:MPIHandler NULL',
        ),

      #####################################
      # FINAL PARAMETERS CHOICE
      #####################################
      parameterSets = cms.vstring(

        #### Choose Energy (Uncomment one of them):
        'set8TeV',
        #'set13TeV',
        #'set2760GeV',
        #'set900GeV',

        #### AlphaS & PDF
        'setAlphaS',
        'pdfCTEQ6L1', 

        #### UE Tune
        #'EE4C',
        'EE5C', 

        #### Process
        #'MinBias',
        'HardFlat',

        #### Basic Setup
        'basicSetup',
        'setParticlesStableForDetector',
      ),

    filterEfficiency = cms.untracked.double(1)
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from Configuration.GenProduction.rivet_customize
#from Configuration.GenProduction.rivet_customize import customise 

def customise(process):
	process.load('GeneratorInterface.RivetInterface.rivetAnalyzer_cfi')
	process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_SMP_16_011')
        #	process.rivetAnalyzer.CrossSection = cms.double(0.469e+09)
	process.rivetAnalyzer.OutputFile = cms.string('output_Klaus_noxsec.yoda')
	process.generation_step+=process.rivetAnalyzer
        process.schedule.remove(process.RAWSIMoutput_step)
	return(process)


#call to customisation function customise imported from Configuration.GenProduction.rivet_customize
process = customise(process)

# End of customisation functions
#process.rivetAnalyzer.OutputFile = cms.string('Rivet_HPPEE5C_HardFlat_'+str(SEED)+'.aida')

