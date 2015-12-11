#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "../interface/TreeEntries_halotree.h"

#include "TVector2.h"

#include <cmath>
#include <limits>

class BeamHaloTreeMaker : public edm::EDAnalyzer {
public:
  explicit BeamHaloTreeMaker(edm::ParameterSet const&);
  ~BeamHaloTreeMaker() {}

private:
  void analyze(edm::Event const&, edm::EventSetup const&) override;

  edm::EDGetTokenT<edm::View<reco::Photon>> photonCollectionToken_;
  edm::EDGetTokenT<reco::PFMETCollection> metCollectionToken_;
  edm::EDGetTokenT<reco::MuonCollection> muonCollectionToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollectionToken_;
  edm::EDGetTokenT<reco::CaloClusterCollection> clusterCollectionToken_;
  edm::EDGetTokenT<EcalRecHitCollection> ebRecHitCollectionToken_;
  edm::EDGetTokenT<EcalRecHitCollection> eeRecHitCollectionToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> sieieMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> chIsoMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> nhIsoMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> phIsoMapToken_;
  edm::EDGetTokenT<double> rhoToken_;

  TTree* output_{0};
  halotree::Event event_;
};

BeamHaloTreeMaker::BeamHaloTreeMaker(edm::ParameterSet const& _ps) :
  photonCollectionToken_(consumes<edm::View<reco::Photon>>(_ps.getParameter<edm::InputTag>("photonCollectionTag"))),
  metCollectionToken_(consumes<reco::PFMETCollection>(_ps.getParameter<edm::InputTag>("metCollectionTag"))),
  muonCollectionToken_(consumes<reco::MuonCollection>(_ps.getParameter<edm::InputTag>("muonCollectionTag"))),
  genParticleCollectionToken_(consumes<reco::GenParticleCollection>(_ps.getParameter<edm::InputTag>("genParticleCollectionTag"))),
  clusterCollectionToken_(consumes<reco::CaloClusterCollection>(_ps.getParameter<edm::InputTag>("clusterCollectionTag"))),
  ebRecHitCollectionToken_(consumes<EcalRecHitCollection>(_ps.getParameter<edm::InputTag>("ebRecHitCollectionTag"))),
  eeRecHitCollectionToken_(consumes<EcalRecHitCollection>(_ps.getParameter<edm::InputTag>("eeRecHitCollectionTag"))),
  sieieMapToken_(consumes<edm::ValueMap<float>>(_ps.getParameter<edm::InputTag>("sieieMapTag"))),
  chIsoMapToken_(consumes<edm::ValueMap<float>>(_ps.getParameter<edm::InputTag>("chIsoMapTag"))),
  nhIsoMapToken_(consumes<edm::ValueMap<float>>(_ps.getParameter<edm::InputTag>("nhIsoMapTag"))),
                 phIsoMapToken_(consumes<edm::ValueMap<float>>(_ps.getParameter<edm::InputTag>("phIsoMapTag"))),
                                rhoToken_(consumes<double>(_ps.getParameter<edm::InputTag>("rhoTag")))
{
  edm::Service<TFileService> fileService;
  output_ = fileService->make<TTree>("halotree", "halo tree");
  event_.book(*output_);
}

void
BeamHaloTreeMaker::analyze(edm::Event const& _event, edm::EventSetup const& _eventSetup)
{
  edm::Handle<edm::View<reco::Photon>> photonHndl;
  if (!_event.getByToken(photonCollectionToken_, photonHndl))
    throw cms::Exception("ProductNotFound") << "photons";
  auto& photons(*photonHndl);

  edm::Handle<reco::PFMETCollection> metHndl;
  if (!_event.getByToken(metCollectionToken_, metHndl))
    throw cms::Exception("ProductNotFound") << "met";
  auto& met(metHndl->at(0));

  edm::Handle<reco::MuonCollection> muonHndl;
  if (!_event.getByToken(muonCollectionToken_, muonHndl))
    throw cms::Exception("ProductNotFound") << "muons";
  auto& muons(*muonHndl);

  edm::Handle<reco::GenParticleCollection> genParticleHndl;
  if (!_event.getByToken(genParticleCollectionToken_, genParticleHndl))
    throw cms::Exception("ProductNotFound") << "genParticles";
  auto& genParticles(*genParticleHndl);

  if (genParticles.size() != 1)
    throw cms::Exception("RuntimeError") << "Gen particle collection size != 1";

  edm::Handle<reco::CaloClusterCollection> clusterHndl;
  if (!_event.getByToken(clusterCollectionToken_, clusterHndl))
    throw cms::Exception("ProductNotFound") << "clusters";
  auto& clusters(*clusterHndl);
  
  auto& genHalo(genParticles.at(0));

  edm::Handle<edm::ValueMap<float>> sieieHndl;
  if (!_event.getByToken(sieieMapToken_, sieieHndl))
    throw cms::Exception("ProductNotFound") << "sieie";
  auto& sieieMap(*sieieHndl);

  edm::Handle<edm::ValueMap<float>> chIsoHndl;
  if (!_event.getByToken(chIsoMapToken_, chIsoHndl))
    throw cms::Exception("ProductNotFound") << "chIso";
  auto& chIsoMap(*chIsoHndl);

  edm::Handle<edm::ValueMap<float>> nhIsoHndl;
  if (!_event.getByToken(nhIsoMapToken_, nhIsoHndl))
    throw cms::Exception("ProductNotFound") << "nhIso";
  auto& nhIsoMap(*nhIsoHndl);

  edm::Handle<edm::ValueMap<float>> phIsoHndl;
  if (!_event.getByToken(phIsoMapToken_, phIsoHndl))
    throw cms::Exception("ProductNotFound") << "phIso";
  auto& phIsoMap(*phIsoHndl);

  edm::Handle<double> rhoHndl;
  if (!_event.getByToken(rhoToken_, rhoHndl))
    throw cms::Exception("ProductNotFound") << "rho";
  double rho(*rhoHndl);

  EcalClusterLazyTools lazyTools(_event, _eventSetup, ebRecHitCollectionToken_, eeRecHitCollectionToken_);

  event_.mip.energy = genHalo.energy();
  event_.mip.phi = std::atan2(genHalo.vy(), genHalo.vx());
  event_.mip.r = std::sqrt(genHalo.vx() * genHalo.vx() + genHalo.vy() * genHalo.vy());

  event_.mip.nclus = 0;
  for (auto&& cluster : clusters) {
    if (std::abs(TVector2::Phi_mpi_pi(cluster.position().phi() - event_.mip.phi)) < 0.1)
      event_.mip.nclus += 1;
  }

  event_.muons.clear();
  for (auto&& muon : muons) {
    if (std::abs(TVector2::Phi_mpi_pi(muon.phi() - event_.mip.phi)) < 0.1) {
      event_.muons.resize(event_.muons.size() + 1);
      auto& outMuon(event_.muons.back());
      outMuon.pt = muon.pt();
      outMuon.eta = muon.eta();
      outMuon.phi = muon.phi();
    }
  }

  event_.met.met = met.pt();
  event_.met.phi = met.phi();
  event_.met.sumEt = met.sumEt();

  double etaBinning[] = {0., 1., 1.479, 2., 2.2, 2.3, 2.4, std::numeric_limits<double>::max()};
  double chAreas[] = {0.0234, 0.0189, 0.0171, 0.0129, 0.011, 0.0074, 0.0034};
  double nhAreas[] = {0.0053, 0.0103, 0.0057, 0.007, 0.0152, 0.0232, 0.1709};
  double phAreas[] = {0.078, 0.0629, 0.0264, 0.0462, 0.074, 0.0924, 0.1484};

  event_.photons.clear();
  for (unsigned iPh(0); iPh != photons.size(); ++iPh) {
    auto& photon(photons.at(iPh));

    if (std::abs(TVector2::Phi_mpi_pi(photon.phi() - event_.mip.phi)) < 0.1) {
      auto&& ptr(photons.ptrAt(iPh));

      event_.photons.resize(event_.photons.size() + 1);
      auto& outPhoton(event_.photons.back());
      outPhoton.pt = photon.pt();
      outPhoton.eta = photon.eta();
      outPhoton.phi = photon.phi();

      auto&& superCluster(*photon.superCluster());
      double scEta(std::abs(superCluster.position().eta()));
      outPhoton.isEB = scEta < 1.444;

      outPhoton.sieie = sieieMap[ptr];
      outPhoton.hOverE = photon.hadTowOverEm();
      outPhoton.pixelVeto = !photon.hasPixelSeed();

      unsigned etaBin = 0;
      while (scEta >= etaBinning[etaBin + 1])
        ++etaBin;

      outPhoton.chIso = chIsoMap[ptr] - rho * chAreas[etaBin];
      outPhoton.nhIso = nhIsoMap[ptr] - rho * nhAreas[etaBin];
      outPhoton.phIso = phIsoMap[ptr] - rho * phAreas[etaBin];
      if (outPhoton.isEB) {
        outPhoton.nhIso -= std::exp(0.0028 * outPhoton.pt + 0.5408);
        outPhoton.phIso -= 0.0014 * outPhoton.pt;
      }
      else {
        outPhoton.nhIso -= 0.01725 * outPhoton.pt;
        outPhoton.phIso -= 0.0091 * outPhoton.pt;
      }

      outPhoton.seedTime = lazyTools.SuperClusterTime(superCluster, _event);
    }
  }

  output_->Fill();
}

DEFINE_FWK_MODULE(BeamHaloTreeMaker);
