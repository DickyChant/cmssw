"""
First HGCAL FastSim job: single photons at eta = 2.0 into the CE-E.

Checks, in order of what would break first:
  1. the job runs at all (Phase-2 calorimetry used to segfault),
  2. fastSimProducer:HGCHitsEE is non-empty,
  3. the total silicon energy is ~0.5% of the incident energy, not equal to it,
  4. the per-layer profile has the shape we fitted.

Run:  cmsRun hgcal_photon_fastsim_cfg.py [maxEvents=N] [energy=50]
"""

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Phase2C17I13M9_FastSim_cff import Phase2C17I13M9_FastSim

opts = VarParsing.VarParsing('analysis')
opts.register('energy', 50.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'photon energy [GeV]')
opts.register('seed', 12345, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'RNG seed (must differ per batch job)')
opts.setDefault('maxEvents', 20)
opts.parseArguments()

process = cms.Process('HGCALFS', Phase2C17I13M9_FastSim)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.maxEvents))
process.source = cms.Source('EmptySource')

process.RandomNumberGeneratorService = cms.Service(
    'RandomNumberGeneratorService',
    generator=cms.PSet(initialSeed=cms.untracked.uint32(opts.seed)),
    VtxSmeared=cms.PSet(initialSeed=cms.untracked.uint32(opts.seed + 1)),
    fastSimProducer=cms.PSet(initialSeed=cms.untracked.uint32(opts.seed + 2)),
)

# Fixed eta = 2.0, matching the samples the parametrization was derived from.
process.generator = cms.EDProducer(
    'FlatRandomEGunProducer',
    PGunParameters=cms.PSet(
        PartID=cms.vint32(22),
        MinEta=cms.double(1.999), MaxEta=cms.double(2.001),
        MinPhi=cms.double(-3.14159265359), MaxPhi=cms.double(3.14159265359),
        MinE=cms.double(opts.energy), MaxE=cms.double(opts.energy),
    ),
    AddAntiParticle=cms.bool(False),
    firstRun=cms.untracked.uint32(1),
    psethack=cms.string('single photon E %.0f eta 2' % opts.energy),
)

process.load('Configuration.StandardSequences.VtxSmearedNoSmear_cff')
# fastSimProducer reads "generatorSmeared"; the vertex-smearing module itself is
# called VtxSmeared, so the renaming producer has to be declared explicitly.
process.generatorSmeared = cms.EDProducer('GeneratorSmearedProducer')
process.load('FastSimulation.SimplifiedGeometryPropagator.fastSimProducer_cff')

# FastSim normally propagates through a deliberately misaligned tracker geometry.
# For a calorimetry smoke test the nominal one is enough, so point the propagator
# at the unlabelled GeometricSearchTracker and load the standard producer.
process.load('RecoTracker.GeometryESProducer.TrackerRecoGeometryESProducer_cfi')
process.load('Geometry.CaloEventSetup.CaloTopology_cfi')
for _blk in (process.fastSimProducer.trackerDefinition, process.fastSimProducer.caloDefinition):
    if hasattr(_blk, 'trackerAlignmentLabel'):
        _blk.trackerAlignmentLabel = cms.untracked.string('')

process.MessageLogger.cerr.FwkReport.reportEvery = 5
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit=cms.untracked.int32(-1))
process.MessageLogger.CalorimetryManager = cms.untracked.PSet(limit=cms.untracked.int32(-1))
process.MessageLogger.HGCalProperties = cms.untracked.PSet(limit=cms.untracked.int32(-1))

process.out = cms.OutputModule(
    'PoolOutputModule',
    fileName=cms.untracked.string('hgcal_photon_fastsim.root'),
    outputCommands=cms.untracked.vstring(
        'drop *',
        'keep *_fastSimProducer_*_*',
        'keep *_generatorSmeared_*_*',
        'keep *_generator_*_*',
    ),
)

process.gen = cms.Path(process.generator * process.VtxSmeared * process.generatorSmeared)
process.sim = cms.Path(process.fastSimProducer)
process.outpath = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.gen, process.sim, process.outpath)
