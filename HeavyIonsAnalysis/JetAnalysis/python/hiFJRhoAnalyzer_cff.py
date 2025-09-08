import FWCore.ParameterSet.Config as cms
from RecoHI.HiJetAlgos.HiRecoPFJets_cff import kt4PFJetsForRho
from RecoHI.HiJetAlgos.hiFJRhoProducer import hiFJRhoProducerFinerBins
from RecoHI.HiJetAlgos.hiFJGridEmptyAreaCalculator_cff import hiFJGridEmptyAreaCalculatorFinerBins

hiFJRhoAnalyzer = cms.EDAnalyzer(
    'HiFJRhoAnalyzer',
    etaMap        = cms.InputTag('hiFJRhoProducer','mapEtaEdges','HiForest'),
    rho           = cms.InputTag('hiFJRhoProducer','mapToRho'),
    rhom          = cms.InputTag('hiFJRhoProducer','mapToRhoM'),
    rhoCorr       = cms.InputTag('hiFJGridEmptyAreaCalculator','mapToRhoCorr'),
    rhomCorr      = cms.InputTag('hiFJGridEmptyAreaCalculator','mapToRhoMCorr'),
    rhoCorr1Bin   = cms.InputTag('hiFJGridEmptyAreaCalculator','mapToRhoCorr1Bin'),
    rhomCorr1Bin  = cms.InputTag('hiFJGridEmptyAreaCalculator','mapToRhoMCorr1Bin'),
    rhoGrid       = cms.InputTag('hiFJGridEmptyAreaCalculator','mapRhoVsEtaGrid'),
    meanRhoGrid   = cms.InputTag('hiFJGridEmptyAreaCalculator','mapMeanRhoVsEtaGrid'),
    etaMaxRhoGrid = cms.InputTag('hiFJGridEmptyAreaCalculator','mapEtaMaxGrid'),
    etaMinRhoGrid = cms.InputTag('hiFJGridEmptyAreaCalculator','mapEtaMinGrid'),
    ptJets        = cms.InputTag('hiFJRhoProducer','ptJets'),
    etaJets       = cms.InputTag('hiFJRhoProducer','etaJets'),
    areaJets      = cms.InputTag('hiFJRhoProducer','areaJets'),
    useModulatedRho = cms.bool(False),
)

hiFJRhoAnalyzerFinerBins = cms.EDAnalyzer(
    'HiFJRhoAnalyzer',
    etaMap        = cms.InputTag('hiFJRhoProducerFinerBins','mapEtaEdges','HiForest'),
    rho           = cms.InputTag('hiFJRhoProducerFinerBins','mapToRho'),
    rhom          = cms.InputTag('hiFJRhoProducerFinerBins','mapToRhoM'),
    rhoCorr       = cms.InputTag('hiFJGridEmptyAreaCalculatorFinerBins','mapToRhoCorr'),
    rhomCorr      = cms.InputTag('hiFJGridEmptyAreaCalculatorFinerBins','mapToRhoMCorr'),
    rhoCorr1Bin   = cms.InputTag('hiFJGridEmptyAreaCalculatorFinerBins','mapToRhoCorr1Bin'),
    rhomCorr1Bin  = cms.InputTag('hiFJGridEmptyAreaCalculatorFinerBins','mapToRhoMCorr1Bin'),
    rhoGrid       = cms.InputTag('hiFJGridEmptyAreaCalculatorFinerBins','mapRhoVsEtaGrid'),
    meanRhoGrid   = cms.InputTag('hiFJGridEmptyAreaCalculatorFinerBins','mapMeanRhoVsEtaGrid'),
    etaMaxRhoGrid = cms.InputTag('hiFJGridEmptyAreaCalculatorFinerBins','mapEtaMaxGrid'),
    etaMinRhoGrid = cms.InputTag('hiFJGridEmptyAreaCalculatorFinerBins','mapEtaMinGrid'),
    ptJets        = cms.InputTag('hiFJRhoProducerFinerBins','ptJets'),
    etaJets       = cms.InputTag('hiFJRhoProducerFinerBins','etaJets'),
    areaJets      = cms.InputTag('hiFJRhoProducerFinerBins','areaJets'),
    useModulatedRho = cms.bool(False),
)

hiPuRhoR3Analyzer = hiFJRhoAnalyzer.clone(
    etaMap = cms.InputTag('hiPuRhoR3Producer','mapEtaEdges','HiForest'),
    rho = cms.InputTag('hiPuRhoR3Producer','mapToRho'),
    rhoExtra = cms.InputTag('hiPuRhoR3Producer','mapToRhoExtra'),
    rhom = cms.InputTag('hiPuRhoR3Producer','mapToRhoM'),
    rhoCorr = cms.InputTag('hiPuRhoR3Producer','mapToRhoMedian'),
    rhomCorr = cms.InputTag('hiPuRhoR3Producer','mapToRhoM'),
    rhoCorr1Bin = cms.InputTag('hiPuRhoR3Producer','mapToRho'),
    rhomCorr1Bin = cms.InputTag('hiPuRhoR3Producer','mapToRhoM'),
    nTow = cms.InputTag('hiPuRhoR3Producer','mapToNTow'),
    towExcludePt = cms.InputTag('hiPuRhoR3Producer','mapToTowExcludePt'),
    towExcludePhi = cms.InputTag('hiPuRhoR3Producer','mapToTowExcludePhi'),
    towExcludeEta = cms.InputTag('hiPuRhoR3Producer','mapToTowExcludeEta'),
    rhoGrid = cms.InputTag('hiFJGridEmptyAreaCalculator','mapRhoVsEtaGrid'),
    meanRhoGrid = cms.InputTag('hiFJGridEmptyAreaCalculator','mapMeanRhoVsEtaGrid'),
    etaMaxRhoGrid = cms.InputTag('hiFJGridEmptyAreaCalculator','mapEtaMaxGrid'),
    etaMinRhoGrid = cms.InputTag('hiFJGridEmptyAreaCalculator','mapEtaMinGrid'),
    rhoFlowFitParams = cms.InputTag('hiFJRhoFlowModulationProducer','rhoFlowFitParams'),
    ptJets = cms.InputTag('hiPuRhoR3Producer', 'ptJets'),
    etaJets = cms.InputTag('hiPuRhoR3Producer', 'etaJets'),
    areaJets = cms.InputTag('hiPuRhoR3Producer', 'areaJets'),
    useModulatedRho = cms.bool(True),
)

pfFilter = cms.EDFilter('CandViewSelector',
    src = cms.InputTag('packedPFCandidates'),
    cut = cms.string("")
)

# Add rho estimator
def addRhoSequence(process, pid = 0):
    tag = f'pid{pid}' if pid > 0 else ''
    setattr(process, f'kt4PFJetsForRho{tag}', kt4PFJetsForRho.clone(src = f'packedPFCandidates{tag}'))
    setattr(process, f'hiFJRhoProducerFinerBins{tag}', hiFJRhoProducerFinerBins.clone(jetSource = f'kt4PFJetsForRho{tag}'))
    setattr(process, f'hiFJGridEmptyAreaCalculatorFinerBins{tag}', hiFJGridEmptyAreaCalculatorFinerBins.clone(
        mapEtaEdges = f'hiFJRhoProducerFinerBins{tag}:mapEtaEdges',
        mapToRho = f'hiFJRhoProducerFinerBins{tag}:mapToRho',
        mapToRhoM = f'hiFJRhoProducerFinerBins{tag}:mapToRhoM',
        pfCandSource = f'packedPFCandidates{tag}',
        jetSource = f'kt4PFJetsForRho{tag}'
    ))
    setattr(process, f'hiFJRhoAnalyzerFinerBins{tag}', hiFJRhoAnalyzerFinerBins.clone(
        etaMap        = f'hiFJRhoProducerFinerBins{tag}:mapEtaEdges',
        rho           = f'hiFJRhoProducerFinerBins{tag}:mapToRho',
        rhom          = f'hiFJRhoProducerFinerBins{tag}:mapToRhoM',
        rhoCorr       = f'hiFJGridEmptyAreaCalculatorFinerBins{tag}:mapToRhoCorr',
        rhomCorr      = f'hiFJGridEmptyAreaCalculatorFinerBins{tag}:mapToRhoMCorr',
        rhoCorr1Bin   = f'hiFJGridEmptyAreaCalculatorFinerBins{tag}:mapToRhoCorr1Bin',
        rhomCorr1Bin  = f'hiFJGridEmptyAreaCalculatorFinerBins{tag}:mapToRhoMCorr1Bin',
        rhoGrid       = f'hiFJGridEmptyAreaCalculatorFinerBins{tag}:mapRhoVsEtaGrid',
        meanRhoGrid   = f'hiFJGridEmptyAreaCalculatorFinerBins{tag}:mapMeanRhoVsEtaGrid',
        etaMaxRhoGrid = f'hiFJGridEmptyAreaCalculatorFinerBins{tag}:mapEtaMaxGrid',
        etaMinRhoGrid = f'hiFJGridEmptyAreaCalculatorFinerBins{tag}:mapEtaMinGrid',
        ptJets        = f'hiFJRhoProducerFinerBins{tag}:ptJets',
        etaJets       = f'hiFJRhoProducerFinerBins{tag}:etaJets',
        areaJets      = f'hiFJRhoProducerFinerBins{tag}:areaJets'
    ))
    setattr(process, f'rhoSequence{tag}', cms.Sequence(getattr(process, f'kt4PFJetsForRho{tag}') + getattr(process, f'hiFJRhoProducerFinerBins{tag}') + getattr(process, f'hiFJGridEmptyAreaCalculatorFinerBins{tag}') + getattr(process, f'hiFJRhoAnalyzerFinerBins{tag}')))
    if pid > 0:
        if pid == 1:
            setattr(process, f'packedPFCandidates{tag}', pfFilter.clone(cut = f'pdgId() == 211 || pdgId() == -211'))
        elif pid == 2:
            setattr(process, f'packedPFCandidates{tag}', pfFilter.clone(cut = f'pdgId() == 130 || pdgId() == 1'))
        elif pid == 3:
            setattr(process, f'packedPFCandidates{tag}', pfFilter.clone(cut = f'pdgId() == 22 || pdgId() == 2'))
        else:
            raise RuntimeError(f'Incorrect particle ID {pid} for setting rho analyzer')
        getattr(process, f'rhoSequence{tag}').insert(0, getattr(process, f'packedPFCandidates{tag}'))
    if not hasattr(process, 'rhoSequences'):
        process.rhoSequences = cms.Sequence()
    process.rhoSequences += getattr(process, f'rhoSequence{tag}')
