# Validate and benchmark HGCalFastGeometry on any Phase-2 geometry.
#
# Usage:
#   cmsRun hgcalFastGeometryTester_cfg.py --geom D110 --era Phase2C17I13M9
#   cmsRun hgcalFastGeometryTester_cfg.py --geom D122 --era Phase2C26I13M9 \
#       --sense HGCalHESiliconSensitive
import argparse
import importlib
import sys

import FWCore.ParameterSet.Config as cms

parser = argparse.ArgumentParser()
parser.add_argument("--geom", default="D110", help="GeometryExtendedRun4<geom>Reco_cff")
parser.add_argument("--era", default="Phase2C17I13M9", help="era matching the geometry")
parser.add_argument("--sense", default="HGCalEESensitive",
                    help="HGCalEESensitive | HGCalHESiliconSensitive | HGCalHEScintillatorSensitive")
parser.add_argument("--maxCells", type=int, default=200000)
args = parser.parse_args(sys.argv[2:] if sys.argv[0].endswith("cmsRun") else sys.argv[1:])

era = getattr(importlib.import_module(f"Configuration.Eras.Era_{args.era}_cff"), args.era)
process = cms.Process("GEOTEST", era)
process.load(f"Configuration.Geometry.GeometryExtendedRun4{args.geom}Reco_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "INFO"
process.MessageLogger.HGCalFastGeom = cms.untracked.PSet(limit=cms.untracked.int32(-1))

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

process.test = cms.EDAnalyzer("HGCalFastGeometryTester",
                              nameSense=cms.string(args.sense),
                              benchmark=cms.bool(True),
                              maxCells=cms.int32(args.maxCells))
process.p = cms.Path(process.test)
