import FWCore.ParameterSet.Config as cms

genDecayFilter = cms.EDFilter("GenDecayFilter",
    sourceTag = cms.InputTag("genParticles"),
    useGenParticles = cms.untracked.bool(True),
    filterExpression = cms.string(""),
    veto = cms.bool(False)
)
