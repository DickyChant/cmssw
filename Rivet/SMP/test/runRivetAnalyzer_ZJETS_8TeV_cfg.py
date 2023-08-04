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

process.source = cms.Source("EmptySource")

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.generator = cms.EDFilter("Pythia8GeneratorFilter",
     comEnergy = cms.double(8000.0),
     filterEfficiency = cms.untracked.double(1),
     maxEventsToPrint = cms.untracked.int32(0),
     pythiaHepMCVerbosity = cms.untracked.bool(False),
     pythiaPylistVerbosity = cms.untracked.int32(1),
     PythiaParameters = cms.PSet(
         processParameters = cms.vstring(
             'Main:timesAllowErrors = 10000',
             'ParticleDecays:limitTau0 = on',
             'ParticleDecays:tauMax = 10',
             'Tune:ee 7',
             'Tune:pp 14',      # Monash tune
             '# Switch on gamma*/Z',
             'WeakSingleBoson:ffbar2gmZ = On',
             '# Switch off all Z decays',
             '23:onMode = off',
             '# Switch back on those to electrons and muons',
             '23:onIfAny = 11 13',
         ),
         parameterSets = cms.vstring('processParameters')
     )
)

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

if options.runOnly:
    process.rivetAnalyzer.AnalysisNames = cms.vstring(options.runOnly)
else:
    process.rivetAnalyzer.AnalysisNames = cms.vstring(
        'CMS_2017_I1497519'
    )

process.rivetAnalyzer.OutputFile      = options.yodafile
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")

process.p = cms.Path(process.generator*process.rivetAnalyzer)


