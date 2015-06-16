#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

class CollectionMerger : public edm::EDProducer {
public:
  explicit CollectionMerger(edm::ParameterSet const&);
  ~CollectionMerger();

private:
  typedef edm::View<reco::Candidate> CandidateView;
  
  void produce(edm::Event&, edm::EventSetup const&);

  std::vector<edm::EDGetTokenT<CandidateView> > viewTokens_;
  std::string colType_;
  //  bool clone_;
  unsigned maxOutputSize_;
};

CollectionMerger::CollectionMerger(edm::ParameterSet const& iConfig) :
  colType_(iConfig.getParameter<std::string>("type")),
  //  clone_(false),
  maxOutputSize_(-1)
{
  for (auto&& tag : iConfig.getParameter<std::vector<edm::InputTag> >("src"))
    viewTokens_.push_back(consumes<CandidateView>(tag));

  //   if(iConfig.existsAs<bool>("clone"))
  //     clone_ = iConfig.getParameter<bool>("clone");

  if(iConfig.existsAs<int>("maxOutputSize"))
    maxOutputSize_ = iConfig.getParameter<int>("maxOutputSize");

  if(colType_ == "Photon"){
    //     if(clone_)
    produces<reco::PhotonCollection>();
    //     else
    //       produces<edm::RefToBaseVector<reco::Photon> >();
  }
  else if(colType_ == "Electron"){
    //     if(clone_)
    produces<reco::GsfElectronCollection>();
    //     else
    //       produces<reco::GsfElectronRefVector>();
  }
  else if(colType_ == "Muon"){
    //     if(clone_)
    produces<reco::MuonCollection>();
    //     else
    //       produces<reco::MuonRefVector>();
  }
  else
    throw cms::Exception("Configuration") << colType_;
}


CollectionMerger::~CollectionMerger()
{
}

void
CollectionMerger::produce(edm::Event& iEvent, edm::EventSetup const&)
{
  using namespace edm;

  std::vector<reco::CandidatePtr> merged;

  for (auto&& token : viewTokens_) {
    std::vector<reco::CandidatePtr> vSrc;

    Handle<CandidateView> srcHndl;
    if(iEvent.getByToken(token, srcHndl))
      vSrc = srcHndl->ptrs();
    else{
      // Handle<RefToBaseVector<reco::Candidate> > bsrcHndl;
      // if(iEvent.getByLabel(src_[iS], bsrcHndl))
      //   vSrc = bsrcHndl.product();
      // else
      //   throw cms::Exception("ProductNotFound") << src_[iS];
      throw cms::Exception("ProductNotFound");
    }

    for (auto&& sPtr : vSrc) {
      auto mItr(merged.begin());
      for (; mItr != merged.end(); ++mItr) {
        if (mItr->get() == sPtr.get())
          break;
      }
      if (mItr != merged.end())
        continue;

      merged.push_back(sPtr);
    }
  }

  if(colType_ == "Photon"){
    //     if(clone_){
    std::auto_ptr<reco::PhotonCollection> output(new reco::PhotonCollection);
    unsigned iOut(0);
    for (auto&& mPtr : merged) {
      auto* photon(dynamic_cast<reco::Photon const*>(mPtr.get()));
      if (photon)
        output->push_back(*photon);

      if (++iOut == maxOutputSize_)
        break;
    }
    iEvent.put(output);
    //    }
    //    else{
    //      std::auto_ptr<RefToBaseVector<reco::Photon> > output(new RefToBaseVector<>);
    //      unsigned iOut(0);
    //      for(reco::PhotonRefVector::const_iterator eItr(merged.begin()); eItr != merged.end() && iOut != maxOutputSize_; ++eItr, iOut++)
    //        output->push_back(*eItr);
    //      iEvent.put(output);

    //       std::auto_
    //       produces<reco::PhotonCollection>();
    //     else
    //       produces<edm::RefToBaseVector<reco::Photon> >();
  }
  else if(colType_ == "Electron"){
    std::auto_ptr<reco::GsfElectronCollection> output(new reco::GsfElectronCollection);
    unsigned iOut(0);
    for (auto&& mPtr : merged) {
      auto* electron(dynamic_cast<reco::GsfElectron const*>(mPtr.get()));
      if (electron)
        output->push_back(*electron);

      if (++iOut == maxOutputSize_)
        break;
    }
    iEvent.put(output);
  }
  else if(colType_ == "Muon"){
    std::auto_ptr<reco::MuonCollection> output(new reco::MuonCollection);
    unsigned iOut(0);
    for (auto&& mPtr : merged) {
      auto* muon(dynamic_cast<reco::Muon const*>(mPtr.get()));
      if (muon)
        output->push_back(*muon);

      if (++iOut == maxOutputSize_)
        break;
    }
    iEvent.put(output);
  }
}

DEFINE_FWK_MODULE(CollectionMerger);
