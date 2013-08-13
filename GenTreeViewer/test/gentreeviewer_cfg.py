import FWCore.ParameterSet.Config as cms

process = cms.Process("View")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(30) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://xrootd.unl.edu//store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/00001/FEAC1583-031B-E211-AF54-00215E2283FA.root'
    )
)

process.load("Toolset.GenTreeViewer.genTreeViewer_cfi")

process.p = cms.Path(process.genTreeViewer)
