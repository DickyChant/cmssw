import FWCore.ParameterSet.Config as cms

def candidateBtaggingMiniAOD(process, isMC = True, jetPtMin = 15, jetR = 0.4, jetCorrLevels = ['L2Relative', 'L3Absolute'], doFlow = False):
    # DeepNtuple settings
    R = str(int(jetR*10))
    jetCorrectionsAK = ('AK4PF', jetCorrLevels, 'None')

    bTagInfos = [
        'pfDeepCSVTagInfos',
        'pfDeepFlavourTagInfos',
        'pfImpactParameterTagInfos',
        'pfInclusiveSecondaryVertexFinderTagInfos',
        'pfParticleTransformerAK4TagInfos',
        'pfUnifiedParticleTransformerAK4TagInfos'
    ]

    bTagDiscriminators = [
        'pfJetProbabilityBJetTags',
        'pfNegativeOnlyJetProbabilityBJetTags',
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probbb',
        'pfDeepCSVJetTags:probc',
        'pfDeepCSVJetTags:probudsg',
        'pfDeepFlavourJetTags:probb',
        'pfDeepFlavourJetTags:probbb',
        'pfDeepFlavourJetTags:probc',
        'pfDeepFlavourJetTags:probg',
        'pfDeepFlavourJetTags:problepb',
        'pfDeepFlavourJetTags:probuds',
        'pfParticleTransformerAK4JetTags:probb',
        'pfParticleTransformerAK4JetTags:probbb',
        'pfParticleTransformerAK4JetTags:probc',
        'pfParticleTransformerAK4JetTags:probg',
        'pfParticleTransformerAK4JetTags:problepb',
        'pfParticleTransformerAK4JetTags:probuds',
    ]
    from RecoBTag.ONNXRuntime.pfUnifiedParticleTransformerAK4JetTags_cfi import pfUnifiedParticleTransformerAK4JetTags as _UParTJetTags
    bTagDiscriminators += ['pfUnifiedParticleTransformerAK4JetTags:' + f for f in _UParTJetTags.flav_names]
    bTagDiscriminators += ['pfNegativeUnifiedParticleTransformerAK4JetTags:' + f for f in _UParTJetTags.flav_names]

    # Create gen-level information
    if isMC:
        from RecoHI.HiJetAlgos.hiSignalParticleProducer_cfi import hiSignalParticleProducer as hiSignalGenParticles
        process.hiSignalGenParticles = hiSignalGenParticles.clone(
            src = "prunedGenParticles"
        )
        from PhysicsTools.PatAlgos.producersHeavyIons.heavyIonJets_cff import allPartons
        process.allPartons = allPartons.clone(
            src = 'hiSignalGenParticles'
        )
        process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector",
            src = cms.InputTag("packedGenParticlesSignal"),
            cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")
        )
        from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
        setattr(process,f'ak{R}GenJetsRecluster', ak4GenJets.clone(
            src = 'packedGenParticlesForJetsNoNu',
            rParam = jetR
        ))
        setattr(process,f'genAK{R}Task', cms.Task(process.hiSignalGenParticles, process.allPartons, process.packedGenParticlesForJetsNoNu, getattr(process,f'ak{R}GenJetsRecluster')))

    # Add negative secondary vertices
    from RecoVertex.AdaptiveVertexFinder.inclusiveNegativeVertexing_cff import inclusiveCandidateNegativeVertexFinder, candidateNegativeVertexMerger, candidateNegativeVertexArbitrator, inclusiveCandidateNegativeSecondaryVertices
    process.inclusiveCandidateNegativeVertexFinder = inclusiveCandidateNegativeVertexFinder.clone(
        tracks = "packedPFCandidates",
        primaryVertices = "offlineSlimmedPrimaryVertices",
    )
    process.candidateNegativeVertexMerger = candidateNegativeVertexMerger.clone()
    process.candidateNegativeVertexArbitrator = candidateNegativeVertexArbitrator.clone(
        tracks = "packedPFCandidates",
        primaryVertices = "offlineSlimmedPrimaryVertices"
    )
    process.inclusiveCandidateNegativeSecondaryVertices = inclusiveCandidateNegativeSecondaryVertices.clone()
    process.svTask = cms.Task(process.inclusiveCandidateNegativeVertexFinder, process.candidateNegativeVertexMerger, process.candidateNegativeVertexArbitrator, process.inclusiveCandidateNegativeSecondaryVertices)

    # Create unsubtracted reco jets
    from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
    addJetCollection(
        process,
        postfix            = "UnsubJets",
        labelName          = f"AK{R}PF",
        jetSource          = cms.InputTag(f"ak{R}PFUnsubJets"),
        algo               = "ak", #name of algo must be in this format
        rParam             = jetR,
        pvSource           = cms.InputTag("offlineSlimmedPrimaryVertices"),
        pfCandidates       = cms.InputTag("packedPFCandidates"),
        svSource           = cms.InputTag('slimmedSecondaryVertices'),
        muSource           = cms.InputTag("slimmedMuons"),
        elSource           = cms.InputTag("slimmedElectrons"),
        getJetMCFlavour    = isMC,
        genJetCollection   = cms.InputTag(f"ak{R}GenJetsRecluster" if isMC else ""),
        genParticles       = cms.InputTag("hiSignalGenParticles" if isMC else ""),
        jetCorrections     = jetCorrectionsAK,
    )
    getattr(process,f'patJetsAK{R}PFUnsubJets').useLegacyJetMCFlavour = False
    getattr(process,f'patJetPartonMatchAK{R}PFUnsubJets').maxDeltaR = jetR

    from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff import ak4PFJets
    setattr(process,f'ak{R}PFUnsubJets', ak4PFJets.clone(
        src = 'packedPFCandidates',
        rParam = jetR,
        jetPtMin = jetPtMin
    ))
    process.patAlgosToolsTask.add(getattr(process,f'ak{R}PFUnsubJets'))

    # Create HIN subtracted reco jets
    jL = f"Cs{R}Flow" if doFlow else f"Cs{R}"
    from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
    addJetCollection(
        process,
        postfix            = "",
        labelName          = f"AK{jL}PF",
        jetSource          = cms.InputTag(f"ak{jL}PFJets"),
        algo               = "ak", #name of algo must be in this format
        rParam             = jetR,
        pvSource           = cms.InputTag("offlineSlimmedPrimaryVertices"),
        pfCandidates       = cms.InputTag("packedPFCandidates"),
        svSource           = cms.InputTag('slimmedSecondaryVertices'),
        muSource           = cms.InputTag("slimmedMuons"),
        elSource           = cms.InputTag("slimmedElectrons"),
        getJetMCFlavour    = isMC,
        genJetCollection   = cms.InputTag(f"ak{R}GenJetsRecluster" if isMC else ""),
        genParticles       = cms.InputTag("hiSignalGenParticles" if isMC else ""),
        jetCorrections     = jetCorrectionsAK,
    )
    getattr(process,f'patJetsAK{jL}PF').embedPFCandidates = True
    getattr(process,f'patJetPartonMatchAK{jL}PF').maxDeltaR = jetR

    if not isMC:
        for label in [f"patJetsAK{R}PFUnsubJets", f"patJetsAK{jL}PF"]:
            getattr(process, label).addGenJetMatch = False
            getattr(process, label).addGenPartonMatch = False
            getattr(process, label).embedGenJetMatch = False
            getattr(process, label).embedGenPartonMatch = False
            getattr(process, label).genJetMatch = ""
            getattr(process, label).genPartonMatch = ""
    else:
        getattr(process,f'patJetPartonAssociationLegacyAK{R}PFUnsubJets').coneSizeToAssociate = min(jetR, 0.3)
        getattr(process,f'patJetPartonAssociationLegacyAK{jL}PF').coneSizeToAssociate = min(jetR, 0.3)

    from PhysicsTools.PatAlgos.producersHeavyIons.heavyIonJets_cff import PackedPFTowers, hiPuRho
    process.PackedPFTowers = PackedPFTowers.clone()
    process.hiPuRho = hiPuRho.clone(
        src = 'PackedPFTowers'
    )
    if doFlow:
        from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff import akCs4PFJets, ak4PFJetsForFlow, hiFJRhoFlowModulation
        process.ak4PFJetsForFlow = ak4PFJetsForFlow.clone(
            src = "PackedPFTowers"
        )
        process.hiFJRhoFlowModulation = hiFJRhoFlowModulation.clone(
            jetTag = "ak4PFJetsForFlow",
            pfCandidateEtaCut = 2,
            minPfCandidatesPerEvent = 60
        )
        setattr(process,f'akCs{R}PFJetsForFlow', akCs4PFJets.clone(
            src = 'packedPFCandidates',
            rParam = jetR,
            jetPtMin = 40,
            useModulatedRho = True,
            rhoFlowFitParams = "hiFJRhoFlowModulation:rhoFlowFitParams"
        ))
        setattr(process,f'hiAKCs{R}FJRhoFlowModulationIter', hiFJRhoFlowModulation.clone(
            jetTag = f"akCs{R}PFJetsForFlow",
            doJettyExclusion = True,
            pfCandidateEtaCut = 2,
            minPfCandidatesPerEvent = 60
        ))
        setattr(process,f'ak{jL}PFJets', akCs4PFJets.clone(
            src = 'packedPFCandidates',
            rParam = jetR,
            jetPtMin = jetPtMin,
            useModulatedRho = True,
            rhoFlowFitParams = f"hiAKCs{R}FJRhoFlowModulationIter:rhoFlowFitParams",
            minFlowChi2Prob = cms.double(0),
            maxFlowChi2Prob = cms.double(1)
        ))
        for mod in ["ak4PFJetsForFlow", "hiFJRhoFlowModulation", f"akCs{R}PFJetsForFlow", f"hiAKCs{R}FJRhoFlowModulationIter"]:
            process.patAlgosToolsTask.add(getattr(process, mod))
    else:
        from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff import akCs4PFJets
        setattr(process,f'ak{jL}PFJets', akCs4PFJets.clone(
            src = 'packedPFCandidates',
            rParam = jetR,
            jetPtMin = jetPtMin
        ))
    for mod in ["PackedPFTowers", "hiPuRho", f"ak{jL}PFJets"]:
        process.patAlgosToolsTask.add(getattr(process, mod))

    # Create b-tagging sequence ----------------
    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    updateJetCollection(
        process,
        labelName = f"AK{jL}DeepFlavour",
        jetSource = cms.InputTag(f'patJetsAK{jL}PF'),
        jetCorrections = jetCorrectionsAK,
        pfCandidates = cms.InputTag('packedPFCandidates'),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        explicitJTA = False
    )

    setattr(process,f'unsubUpdatedPatJetsAK{jL}DeepFlavour', cms.EDProducer("JetMatcherDR",
        source = cms.InputTag(f"updatedPatJetsAK{jL}DeepFlavour"),
        matched = cms.InputTag(f"patJetsAK{R}PFUnsubJets")
    ))
    process.patAlgosToolsTask.add(getattr(process,f'unsubUpdatedPatJetsAK{jL}DeepFlavour'))

    for mod in [f'pfUnifiedParticleTransformerAK4', f'pfNegativeUnifiedParticleTransformerAK4']:
        getattr(process, f'{mod}JetTagsAK{jL}DeepFlavour').model_path = 'HeavyIonsAnalysis/Configuration/data/'+(f'PbPb_AK3_2024_v6.onnx' if jetR==0.3 else f'UParTAK{R}_PbPb_2024.onnx')
        getattr(process, f'{mod}TagInfosAK{jL}DeepFlavour').sort_cand_by_pt = True
        getattr(process, f'{mod}TagInfosAK{jL}DeepFlavour').fix_lt_sorting = True
        getattr(process, f'{mod}TagInfosAK{jL}DeepFlavour').secondary_vertices = 'inclusiveCandidateNegativeSecondaryVertices' if 'Negative' in mod else 'slimmedSecondaryVertices'

    getattr(process,f'pfImpactParameterTagInfosAK{jL}DeepFlavour').maxDeltaR = jetR
    for taginfo in [f"pfDeepFlavourTagInfosAK{jL}DeepFlavour", f"pfParticleTransformerAK4TagInfosAK{jL}DeepFlavour", f"pfUnifiedParticleTransformerAK4TagInfosAK{jL}DeepFlavour", f"pfNegativeUnifiedParticleTransformerAK4TagInfosAK{jL}DeepFlavour"]:
        getattr(process, taginfo).jet_radius = jetR

    if hasattr(process,f'updatedPatJetsTransientCorrectedAK{jL}DeepFlavour'):
        getattr(process,f'updatedPatJetsTransientCorrectedAK{jL}DeepFlavour').addTagInfos = True
        getattr(process,f'updatedPatJetsTransientCorrectedAK{jL}DeepFlavour').addBTagInfo = True
    else:
        raise ValueError(f'I could not find updatedPatJetsTransientCorrectedAK{jL}DeepFlavour to embed the tagInfos, please check the cfg')

    # Remove PUPPI
    process.patAlgosToolsTask.remove(process.packedpuppi)
    process.patAlgosToolsTask.remove(process.packedpuppiNoLep)
    getattr(process,f'pfInclusiveSecondaryVertexFinderTagInfosAK{jL}DeepFlavour').weights = ""
    for taginfo in [f"pfDeepFlavourTagInfosAK{jL}DeepFlavour", f"pfParticleTransformerAK4TagInfosAK{jL}DeepFlavour", f"pfUnifiedParticleTransformerAK4TagInfosAK{jL}DeepFlavour", f"pfNegativeUnifiedParticleTransformerAK4TagInfosAK{jL}DeepFlavour"]:
        getattr(process, taginfo).fallback_puppi_weight = True
        getattr(process, taginfo).fallback_vertex_association = True
        getattr(process, taginfo).unsubjet_map = f"unsubUpdatedPatJetsAK{jL}DeepFlavour"
        getattr(process, taginfo).puppi_value_map = ""

    # Match with unsubtracted jets
    setattr(process,f'unsubAK{jL}JetMap', getattr(process,f'unsubUpdatedPatJetsAK{jL}DeepFlavour').clone(
        source = f"selectedUpdatedPatJetsAK{jL}DeepFlavour"
    ))
    process.patAlgosToolsTask.add(getattr(process,f'unsubAK{jL}JetMap'))

    # Add extra b tagging algos
    from RecoBTag.ImpactParameter.pfJetProbabilityBJetTags_cfi import pfJetProbabilityBJetTags
    setattr(process,f'pfJetProbabilityBJetTagsAK{jL}DeepFlavour', pfJetProbabilityBJetTags.clone(tagInfos = [f"pfImpactParameterTagInfosAK{jL}DeepFlavour"]))
    process.patAlgosToolsTask.add(getattr(process,f'pfJetProbabilityBJetTagsAK{jL}DeepFlavour'))

    # Associate to forest sequence
    if isMC:
        process.forest.associate(getattr(process,f'genAK{R}Task'))
    process.forest.associate(process.svTask)
    process.forest.associate(process.patAlgosToolsTask)
