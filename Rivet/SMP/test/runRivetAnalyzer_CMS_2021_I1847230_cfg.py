import FWCore.ParameterSet.Config as cms

import os
#mode = "QCD13TeV"
#mode = "QCD8TeV"
#mode = "ZJet"
mode = os.getenv("RIVETMODE")
model = "PYTHIA"
process = cms.Process("runRivetAnalysis")

process.options.numberOfThreads=cms.untracked.uint32(24)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    generator = cms.PSet(
        initialSeed = cms.untracked.uint32(5161993),
        engineName = cms.untracked.string('HepJamesRandom')
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000000)
)
process.source = cms.Source("EmptySource")

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *

if mode.startswith("QCD"):
  if mode.endswith("8TeV"): QCDCM = 8000.0
  else: QCDCM = 13000.0
  process.generator = cms.EDFilter("Pythia8GeneratorFilter",
                           pythiaPylistVerbosity = cms.untracked.int32(0),
                           filterEfficiency = cms.untracked.double(1.0),
                           pythiaHepMCVerbosity = cms.untracked.bool(False),
                           comEnergy = cms.double(QCDCM),
                           maxEventsToPrint = cms.untracked.int32(0),
                           PythiaParameters = cms.PSet(
          pythia8CommonSettingsBlock,
          pythia8CUEP8M1SettingsBlock,
          processParameters = cms.vstring(
              'HardQCD:all = on',
              'PhaseSpace:pTHatMin = 300.0'
              ),
          parameterSets = cms.vstring('pythia8CommonSettings',
                                      'pythia8CUEP8M1Settings',
                                      'processParameters',
                                      )
          )
                           )
 
else:
  process.generator = cms.EDFilter("Pythia8GeneratorFilter",
                           pythiaPylistVerbosity = cms.untracked.int32(1),
                           filterEfficiency = cms.untracked.double(1.0),
                           pythiaHepMCVerbosity = cms.untracked.bool(False),
                           comEnergy = cms.double(8000.0),
                           maxEventsToPrint = cms.untracked.int32(0),
                           PythiaParameters = cms.PSet(
          pythia8CommonSettingsBlock,
          pythia8CUEP8M1SettingsBlock,
          processParameters = cms.vstring(
              'WeakSingleBoson:ffbar2gmZ = on',
              '23:onMode = off',
              '23:onIfAny = 13',
              '23:mMin = 50.'
              ),
          parameterSets = cms.vstring('pythia8CommonSettings',
                                      'pythia8CUEP8M1Settings',
                                      'processParameters',
                                      )
          )
                           )
process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2021_I1847230:MODE={}'.format(mode))
outPutName = "{}.yoda".format(mode)
if model == "PYTHIA":
  if mode == "QCD8TeV": outPutName = "PYTHIA LO 2jets+PS (8TeV).yoda"
  if mode == "QCD13TeV": outPutName = "PYTHIA LO 2jets+PS (13TeV).yoda"
  if mode == "ZJet": outPutName = "PYTHIA LO Z+1j+PS (8TeV).yoda"
process.rivetAnalyzer.OutputFile = cms.string(outPutName)
process.p = cms.Path(process.generator*process.rivetAnalyzer)

