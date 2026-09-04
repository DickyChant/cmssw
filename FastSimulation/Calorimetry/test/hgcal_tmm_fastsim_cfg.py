"""Smoke test for the TMM splatting EM model: identical to
hgcal_photon_fastsim_cfg.py but with simulateTMMSplat = True and the
parameter file produced by test/export_tmm_params.py.

Run:  cmsRun hgcal_tmm_fastsim_cfg.py [maxEvents=N] [energy=50] \
          [tmmParams=/path/to/tmm_splat_photon.params]
Pass: job completes and fastSimProducer:HGCHitsEE is non-empty.
"""
import importlib.util as _ilu
import os
import sys

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# reuse the photon config wholesale
_here = os.path.dirname(os.path.abspath(__file__))
opts = VarParsing.VarParsing('analysis')
opts.register('energy', 50.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'gun energy [GeV]')
opts.register('seed', 12345, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'random seed')
opts.register('tmmParams',
              '/afs/cern.ch/work/s/sqian/hgcal_build/tmm_splat_photon.params',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'TMM parameter file (from export_tmm_params.py)')
opts.setDefault('maxEvents', 20)
opts.parseArguments()

_spec = _ilu.spec_from_file_location(
    'hgcal_photon_base', os.path.join(_here, 'hgcal_photon_fastsim_cfg.py'))
# the base config parses sys.argv itself; hand it only known options
sys.argv = [sys.argv[0], f'maxEvents={opts.maxEvents}',
            f'energy={opts.energy}', f'seed={opts.seed}']
_base = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_base)
process = _base.process

process.fastSimProducer.Calorimetry.HGCal.simulateTMMSplat = cms.bool(True)
process.fastSimProducer.Calorimetry.HGCal.tmmParamsFile = \
    cms.string(opts.tmmParams)
