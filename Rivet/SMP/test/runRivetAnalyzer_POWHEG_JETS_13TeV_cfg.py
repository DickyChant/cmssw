# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/MinBias_TuneZ2_7TeV_pythia6_cff.py -s GEN --datatier=GEN-SIM-RAW --conditions auto:mc --eventcontent RAWSIM --no_exec -n 10000 --python_filename=rivet_cfg.py --customise=Configuration/GenProduction/rivet_customize.py

import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *

import os
import sys
process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic7TeV2011Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')


inputfile_list = open( 'files_powheg.txt' ).readlines()
dirstr = 'root://cmsxrootd.fnal.gov//store/user/rappocc/POWHEG_DIJET_LHE/'
#infilestr = '/store/user/rappocc/POWHEG_DIJET_LHE/Dijet_NNPDF30_born10_LHESET/pwgevents-%04d.lhe'
infiles = [ ]
for ifile in inputfile_list : 
    infiles.append( dirstr + '/' + ( ifile.rstrip() ) )
    print 'Added ', infiles[ len(infiles)-1]

process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring(infiles),
    eventsToPrint = cms.untracked.int32(0)
)
process.options = cms.untracked.PSet(

)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    #input = cms.untracked.int32(200)
)

print 'process.source: '
print process.source


#process.RandomNumberGeneratorService.generator.initialSeed = SEED
#print 'Seed = ', int (SEED)


# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.3 $'),
    annotation = cms.untracked.string('PYTHIA6-MinBias TuneZ2 at 7TeV'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/CMSSW/Configuration/GenProduction/python/Attic/MinBias_TuneZ2_7TeV_pythia6_cff.py,v $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('out.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# Other statements
#process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'

process.genstepfilter.triggerConditions=cms.vstring("generation_step")

process.generator = cms.EDFilter("Pythia8HadronizerFilter",
				 comEnergy = cms.double(13000.0),
				 crossSection = cms.untracked.double(3.278e+09),
				 filterEfficiency = cms.untracked.double(1),
				 maxEventsToPrint = cms.untracked.int32(0),
				 pythiaHepMCVerbosity = cms.untracked.bool(False),
				 pythiaPylistVerbosity = cms.untracked.int32(0),
				 PythiaParameters = cms.PSet(
				     pythia8CUEP8M1SettingsBlock,
				     pythia8PowhegEmissionVetoSettingsBlock,
	processParameters = cms.vstring(
#	'PartonLevel:MPI = off',
#	'HadronLevel:all = off',
#	'SpaceShower:QCDshower = off',
#	'TimeShower:QCDshower  = off',
#	'PartonLevel:ISR=off',
#       'PartonLevel:FSR=off',
	'Main:timesAllowErrors    = 10000',
	'ParticleDecays:limitTau0 = on',
	'ParticleDecays:tauMax = 10',
	'POWHEG:pThard=2', #userhook = 2
	),
	parameterSets = cms.vstring('pythia8CUEP8M1Settings','pythia8PowhegEmissionVetoSettings','processParameters')
	)
				 )

process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

# customisation of the process.
process.load('GeneratorInterface.RivetInterface.rivetAnalyzer_cfi')

# Automatic addition of the customisation function from Configuration.GenProduction.rivet_customize
#from Configuration.GenProduction.rivet_customize import customise 

def customise(process):
	process.load('GeneratorInterface.RivetInterface.rivetAnalyzer_cfi')
	process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2018_I1682495')
	process.rivetAnalyzer.CrossSection = cms.double(9757000000)
	#process.rivetAnalyzer.OutputFile = cms.string('CMS_2018_I1682495.yoda')
	process.rivetAnalyzer.OutputFile = cms.string ('powheg.yoda')
	process.generation_step+=process.rivetAnalyzer
	process.schedule.remove(process.RAWSIMoutput_step)
        return(process)      
#call to customisation function customise imported from Configuration.GenProduction.rivet_customize
process = customise(process)

# End of customisation functions
