// -*- C++ -*-
//
// Package:    ObjectCounter
// Class:      ObjectCounter
// 
/**\class ObjectCounter ObjectCounter.cc CommonTools/ObjectCounter/src/ObjectCounter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yutaro Iiyama,512 1-005,+41227670489,
//         Created:  Fri Jun  1 11:20:59 CEST 2012
// $Id: ObjectCounter.cc,v 1.2 2012/09/19 11:28:22 yiiyama Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

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

class ObjectCounter : public edm::EDFilter {
   public:
      explicit ObjectCounter(const edm::ParameterSet&);
      ~ObjectCounter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

  edm::InputTag srcTag_;
  edm::InputTag matchTag_;
  unsigned threshold_;
  std::string object_;
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
ObjectCounter::ObjectCounter(const edm::ParameterSet& iConfig) :
  srcTag_(iConfig.getParameter<edm::InputTag>("src")),
  threshold_(iConfig.getUntrackedParameter<int>("threshold", 1)),
  object_(iConfig.getParameter<std::string>("object"))
{
  if(iConfig.existsAs<edm::InputTag>("matchTag"))
    matchTag_ = iConfig.getParameter<edm::InputTag>("matchTag");

}


ObjectCounter::~ObjectCounter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ObjectCounter::filter(edm::Event& iEvent, const edm::EventSetup&)
{
   using namespace edm;

   if(object_ == "Photon"){
     Handle<reco::PhotonCollection> collection;
     if(!iEvent.getByLabel(srcTag_, collection))
       throw cms::Exception("ProductNotFound") << srcTag_;

     return collection->size() >= threshold_;
   }
   else if(object_ == "Electron"){
     Handle<reco::GsfElectronCollection> collection;
     if(!iEvent.getByLabel(srcTag_, collection))
       throw cms::Exception("ProductNotFound") << srcTag_;

     return collection->size() >= threshold_;
   }
   else if(object_ == "Muon"){
     Handle<reco::MuonCollection> collection;
     if(!iEvent.getByLabel(srcTag_, collection))
       throw cms::Exception("ProductNotFound") << srcTag_;

     return collection->size() >= threshold_;
   }
   else if(object_ == "PhotonRef"){
     Handle<RefToBaseVector<reco::Photon> > collection;
     if(!iEvent.getByLabel(srcTag_, collection))
       throw cms::Exception("ProductNotFound") << srcTag_;

     return collection->size() >= threshold_;
   }
   else if(object_ == "ElectronRef"){
     Handle<reco::GsfElectronRefVector> collection;
     if(iEvent.getByLabel(srcTag_, collection))
       return collection->size() >= threshold_;
     else{
       Handle<RefToBaseVector<reco::GsfElectron> > collection2;
       if(iEvent.getByLabel(srcTag_, collection2))
         return collection2->size() >= threshold_;
       else
         throw cms::Exception("ProductNotFound") << srcTag_;
     }
   }
   else if(object_ == "MuonRef"){
     Handle<reco::MuonRefVector> collection;
     if(iEvent.getByLabel(srcTag_, collection))
       return collection->size() >= threshold_;
     else{
       Handle<RefToBaseVector<reco::Muon> > collection2;
       if(iEvent.getByLabel(srcTag_, collection2))
         return collection2->size() >= threshold_;
       else
         throw cms::Exception("ProductNotFound") << srcTag_;
     }
   }
   else if(object_ == "Candidate"){
     Handle<reco::CandidateView> collection;
     if(!iEvent.getByLabel(srcTag_, collection)){
       throw cms::Exception("ProductNotFound") << srcTag_;
     }

     return collection->size() >= threshold_;
   }
   else
     return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
ObjectCounter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ObjectCounter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
ObjectCounter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
ObjectCounter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
ObjectCounter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
ObjectCounter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ObjectCounter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ObjectCounter);
