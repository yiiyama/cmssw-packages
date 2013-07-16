import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/data/yiiyama/AOD/195397/DoubleElectron.root'
#        'file:/data/yiiyama/AOD/195658/MuEG.root'
#        'file:/data/yiiyama/AOD/195658/MinimumBias.root'
#        'file:/data/yiiyama/AOD/194317/DoublePhoton.root'
    )
)

process.demo = cms.EDAnalyzer('TriggerEventAnalyzer',
    pathNames = cms.vstring(
        "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50",
        "HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50",
        "HLT_Photon36_R9Id85_Photon22_R9Id85"        
    ),
    accept = cms.vint32(1, 1, 0),
    maxDisplay = cms.untracked.int32(2)
)


process.p = cms.Path(process.demo)
