// -*- C++ -*-
//
// Package:    CollectionMerger
// Class:      CollectionMerger
// 
/**\class CollectionMerger CollectionMerger.cc Toolset/CollectionMerger/src/CollectionMerger.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yutaro Iiyama,512 1-005,+41227670489,
//         Created:  Thu Jun 21 16:48:38 CEST 2012
// $Id: CollectionMerger.cc,v 1.1 2012/10/06 09:09:21 yiiyama Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBaseVector.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

//
// class declaration
//

class CollectionMerger : public edm::EDProducer {
   public:
      explicit CollectionMerger(const edm::ParameterSet&);
      ~CollectionMerger();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

  std::vector<edm::InputTag> src_;
  std::string colType_;
  //  bool clone_;
  unsigned maxOutputSize_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
CollectionMerger::CollectionMerger(const edm::ParameterSet& iConfig) :
  src_(iConfig.getParameter<std::vector<edm::InputTag> >("src")),
  colType_(iConfig.getParameter<std::string>("type")),
  //  clone_(false),
  maxOutputSize_(-1)
{

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


//
// member functions
//

// ------------ method called to produce the data  ------------
void
CollectionMerger::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  RefToBaseVector<reco::Candidate> merged;

  for(unsigned iS(0); iS < src_.size(); iS++){
    RefToBaseVector<reco::Candidate> const* vSrc(0);

    Handle<View<reco::Candidate> > srcHndl;
    if(iEvent.getByLabel(src_[iS], srcHndl))
      vSrc = &(srcHndl->refVector());
    else{
      Handle<RefToBaseVector<reco::Candidate> > bsrcHndl;
      if(iEvent.getByLabel(src_[iS], bsrcHndl))
        vSrc = bsrcHndl.product();
      else
        throw cms::Exception("ProductNotFound") << src_[iS];
    }

    for(RefToBaseVector<reco::Candidate>::const_iterator cItr(vSrc->begin()); cItr != vSrc->end(); ++cItr){
      bool counted(false);
      for(RefToBaseVector<reco::Candidate>::const_iterator it(merged.begin()); it != merged.end(); ++it){
        if(it->get() == cItr->get()){
          counted = true;
          break;
        }
      }
      if(counted) continue;

      merged.push_back(*cItr);
    }
  }

  if(colType_ == "Photon"){
//     if(clone_){
    std::auto_ptr<reco::PhotonCollection> output(new reco::PhotonCollection);
    unsigned iOut(0);
    for(RefToBaseVector<reco::Candidate>::const_iterator eItr(merged.begin()); eItr != merged.end() && iOut != maxOutputSize_; ++eItr, iOut++)
      output->push_back(*(static_cast<reco::Photon const*>(eItr->get())));
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
    for(RefToBaseVector<reco::Candidate>::const_iterator eItr(merged.begin()); eItr != merged.end() && iOut != maxOutputSize_; ++eItr, iOut++)
      output->push_back(*(static_cast<reco::GsfElectron const*>(eItr->get())));
    iEvent.put(output);
  }
  else if(colType_ == "Muon"){
    std::auto_ptr<reco::MuonCollection> output(new reco::MuonCollection);
    unsigned iOut(0);
    for(RefToBaseVector<reco::Candidate>::const_iterator eItr(merged.begin()); eItr != merged.end() && iOut != maxOutputSize_; ++eItr, iOut++)
      output->push_back(*(static_cast<reco::Muon const*>(eItr->get())));
    iEvent.put(output);
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
CollectionMerger::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CollectionMerger::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
CollectionMerger::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CollectionMerger::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CollectionMerger::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CollectionMerger::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CollectionMerger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CollectionMerger);
