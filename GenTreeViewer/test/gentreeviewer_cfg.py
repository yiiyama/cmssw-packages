import FWCore.ParameterSet.Config as cms

process = cms.Process("View")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(30) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'root://xrootd.unl.edu//store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/00001/FEAC1583-031B-E211-AF54-00215E2283FA.root'
    'root://xrootd.ba.infn.it//store/user/htholen/LHE2EDM_WHIZARD_2to5_ttA/20130610ReRecoMC_smp_WHIZ/196f8124f9470d4cf6c3b680d1da9040/out_patTuple_100_1_Wh9.root'
    ),
    eventsToProcess = cms.untracked.VEventRange('1:356:586', '1:356:600', '1:747:186', '1:747:216', '1:747:220')
)

process.load("Toolset.GenTreeViewer.genTreeViewer_cfi")

process.p = cms.Path(process.genTreeViewer)
