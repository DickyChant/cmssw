import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process("GEOTEST", Phase2C17I13M9)
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit=cms.untracked.int32(-1))
process.MessageLogger.cerr.default = cms.untracked.PSet(limit=cms.untracked.int32(-1))
process.MessageLogger.HGCalFastGeom = cms.untracked.PSet(limit=cms.untracked.int32(-1))

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

process.load('FastSimulation.CaloGeometryTools.hgcalFastGeometryTester_cfi') \
    if False else None
process.test = cms.EDAnalyzer("HGCalFastGeometryTester",
    nameSense = cms.string("HGCalEESensitive"),
    benchmark = cms.bool(True),
    maxCells  = cms.int32(200000),
)
process.p = cms.Path(process.test)
