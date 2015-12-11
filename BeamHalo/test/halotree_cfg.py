import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.parseArguments()

process = cms.Process('Analysis')

process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('beamhalo.root')
)

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '74X_dataRun2_v2'

process.load("RecoEgamma.PhotonIdentification.PhotonIDValueMapProducer_cfi")

process.ntuples = cms.EDAnalyzer('BeamHaloTreeMaker',
    photonCollectionTag = cms.InputTag('gedPhotons'),
    metCollectionTag = cms.InputTag('pfMet'),
    muonCollectionTag = cms.InputTag('muons'),
    genParticleCollectionTag = cms.InputTag('genParticles'),
    clusterCollectionTag = cms.InputTag('particleFlowEGamma:EBEEClusters'),
    ebRecHitCollectionTag = cms.InputTag('reducedEcalRecHitsEB'),
    eeRecHitCollectionTag = cms.InputTag('reducedEcalRecHitsEE'),
    sieieMapTag   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
    chIsoMapTag = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
    nhIsoMapTag = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
    phIsoMapTag = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
    rhoTag = cms.InputTag("fixedGridRhoFastjetAll")
)

process.path = cms.Path(process.photonIDValueMapProducer + process.ntuples)
