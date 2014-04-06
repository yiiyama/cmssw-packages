import FWCore.ParameterSet.Config as cms

genTreeViewer = cms.EDAnalyzer("GenTreeViewer",
    genParticlesTag = cms.untracked.InputTag("genParticles"),
    showMomentum = cms.untracked.bool(True),
    showMass = cms.untracked.bool(True),
    usePtEtaPhi = cms.untracked.bool(True),
    cleaningMode = cms.untracked.int32(1), #0 -> only full tree, 1 -> only cleaned tree, 2 -> both
    minPt = cms.untracked.double(2.)
)
