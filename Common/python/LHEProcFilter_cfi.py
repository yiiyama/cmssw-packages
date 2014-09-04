import FWCore.ParameterSet.Config as cms

lheProcFilter = cms.EDFilter("LHEProcFilter",
                             lheEventsTag = cms.untracked.InputTag("source"),
                             allowedProcs = cms.untracked.vint32(0, 1, 2)
)
