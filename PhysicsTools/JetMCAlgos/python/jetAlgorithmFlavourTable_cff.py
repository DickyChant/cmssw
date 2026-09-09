import FWCore.ParameterSet.Config as cms

def addJetAlgorithmFlavourTable(process, baseTable, task, label,
                               wtaFlavour, wtaMatchStatus, algorithmFlavourInfos=None):
    """Append columns to an existing Nano table without changing legacy labels.

    The supplied WTA ValueMaps and optional GHS/IFN association must be keyed
    to baseTable.src. WTA maps must come from an explicitly defined parton-level
    input and matching step; this helper never infers partons from MiniAOD.
    """
    if hasattr(process, label):
        raise ValueError("Refusing to replace existing module " + label)
    def tag(value):
        return value if isinstance(value, cms.InputTag) else cms.InputTag(value)
    table = cms.EDProducer("JetAlgorithmFlavourTableProducer",
        src=baseTable.src, name=baseTable.name, cut=baseTable.cut,
        useWTA=cms.bool(True), useAlgorithms=cms.bool(algorithmFlavourInfos is not None),
        wtaFlavour=tag(wtaFlavour), wtaMatchStatus=tag(wtaMatchStatus),
        algorithmFlavourInfos=tag(algorithmFlavourInfos if algorithmFlavourInfos is not None else ""))
    setattr(process, label, table)
    task.add(getattr(process, label))
    return getattr(process, label)
