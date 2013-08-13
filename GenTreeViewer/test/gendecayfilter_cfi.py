import FWCore.ParameterSet.Config as cms

process = cms.Process("Filter")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(30) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://xrootd.unl.edu//store/mc/Summer12_DR53X/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_v2/AODSIM/PU_S10_START53_V7A-v1/00000/02A9AA94-F720-E211-AF5C-001E67396E3C.root'
    )
)

process.load("Toolset.GenTreeViewer.genDecayFilter_cfi")
process.genDecayFilter.filterExpression = '15>11 || 15>13'

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("output.root"),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )
)

process.p = cms.Path(process.genDecayFilter)
process.ep = cms.EndPath(process.out)
