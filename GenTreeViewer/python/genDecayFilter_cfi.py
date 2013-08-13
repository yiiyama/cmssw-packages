import FWCore.ParameterSet.Config as cms

genDecayFilter = cms.EDFilter("GenDecayFilter",
    genParticlesTag = cms.untracked.InputTag("genParticles"),
    filterExpression = cms.untracked.string("")
)
