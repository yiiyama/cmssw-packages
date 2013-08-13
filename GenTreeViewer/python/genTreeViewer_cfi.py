import FWCore.ParameterSet.Config as cms

genTreeViewer = cms.EDAnalyzer("GenTreeViewer",
    genParticlesTag = cms.untracked.InputTag("genParticles"),
    cleaningMode = cms.untracked.int32(1),
    minPt = cms.untracked.double(2.)
)
