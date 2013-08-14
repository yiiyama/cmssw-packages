import FWCore.ParameterSet.Config as cms

genDecayFilter = cms.EDFilter("GenDecayFilter",
    genParticlesTag = cms.InputTag("genParticles"),
    filterExpression = cms.string(""),
    veto = cms.bool(False)
)
