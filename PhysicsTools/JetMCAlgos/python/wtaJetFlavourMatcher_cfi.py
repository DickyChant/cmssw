import FWCore.ParameterSet.Config as cms

wtaJetFlavourMatcher = cms.EDProducer("WTAJetFlavourMatcher",
    jets=cms.InputTag("ak4GenJets"),
    partonJets=cms.InputTag("ak4PartonJets"),
    flavour=cms.InputTag("wtaJetFlavour"),
    maxDeltaR=cms.double(0.2),
)
