import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/data/yiiyama/AODSIM/GJet_Pt40_doubleEMEnriched_pythia6_Summer12.root'
    )
)

process.demo = cms.EDAnalyzer('GenTreeViewer',
    genParticlesTag = cms.InputTag("genParticles"),
    cleaningMode = cms.int32(2), # 0 -> show only full tree, 1 -> show only cleaned tree, 2 -> show both
    minPt = cms.double(2.) # minimum Pt required for final state particles
)


process.p = cms.Path(process.demo)
