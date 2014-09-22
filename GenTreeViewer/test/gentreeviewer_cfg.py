import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.parseArguments()

process = cms.Process('VIEW')

process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.load('Toolset.GenTreeViewer.genTreeViewer_cfi')

process.view = cms.Path(process.genTreeViewer)
