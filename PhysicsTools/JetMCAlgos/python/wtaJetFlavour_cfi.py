import FWCore.ParameterSet.Config as cms

# Terminal-parton flavour at the shower cutoff. Do not point this at ak4GenJets
# from a hadronized event: those constituents no longer carry parton flavour.
wtaJetFlavour = cms.EDProducer(
    "WTAJetFlavourProducer",
    src=cms.InputTag("ak4PartonJets"),
)
