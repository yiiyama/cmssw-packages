// -*- C++ -*-
//
// Package:    TriggerEventAnalyzer
// Class:      TriggerEventAnalyzer
// 
/**\class TriggerEventAnalyzer TriggerEventAnalyzer.cc SusyAnalysis/TriggerEventAnalyzer/src/TriggerEventAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yutaro Iiyama,512 1-005,+41227670489,
//         Created:  Sat Jun 16 12:57:28 CEST 2012
// $Id: TriggerEventAnalyzer.cc,v 1.1 2012/09/19 12:30:07 yiiyama Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h" 
//
// class declaration
//

class TriggerEventAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TriggerEventAnalyzer(const edm::ParameterSet&);
      ~TriggerEventAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

  std::vector<std::string> pathNames_;
  std::vector<int> accept_;
  int iDisplay_;
  int maxDisplay_;
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
TriggerEventAnalyzer::TriggerEventAnalyzer(const edm::ParameterSet& iConfig) :
  pathNames_(iConfig.getParameter<std::vector<std::string> >("pathNames")),
  accept_(iConfig.getParameter<std::vector<int> >("accept")),
  iDisplay_(0),
  maxDisplay_(iConfig.getUntrackedParameter<int>("maxDisplay", 1))
{
   //now do what ever initialization is needed

}


TriggerEventAnalyzer::~TriggerEventAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggerEventAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   if(iDisplay_ == maxDisplay_) return;

   Handle<trigger::TriggerEvent> teHndl;
   if(!iEvent.getByLabel("hltTriggerSummaryAOD", teHndl)) return;

   Handle<TriggerResults> trHndl;
   if(!iEvent.getByLabel(InputTag("TriggerResults", "", "HLT"), trHndl)) return;

   vector<bool> passed(pathNames_.size());
   for(unsigned iP(0); iP < passed.size(); iP++)
     passed[iP] = !bool(accept_.at(iP));

   edm::TriggerNames const& triggerNames(iEvent.triggerNames(*trHndl));

   for(unsigned iP(0); iP < pathNames_.size(); iP++){
     for(unsigned iT(0); iT < triggerNames.size(); iT++){
       if(triggerNames.triggerName(iT).find(pathNames_[iP]) != string::npos)
         passed[iP] = (trHndl->accept(iT) ^ passed[iP]);
     }
   }

   bool process(true);
   for(unsigned iP(0); iP < passed.size(); iP++) process &= passed[iP];

   if(!process) return;

   iDisplay_ += 1;

   cout << "TriggerEvent" << endl;

   cout << "COLLECTIONTAGS:" << endl;
   vector<string> const& tags(teHndl->collectionTags());
   for(unsigned iT = 0; iT < tags.size(); iT++)
     cout << " " << iT << ". " << tags[iT] << endl;

   cout << endl << "FILTERS: " << endl;
   unsigned nF(teHndl->sizeFilters());
   for(unsigned iF = 0; iF < nF; iF++){
     cout << " " << iF << ". " << teHndl->filterTag(iF) << ": ";
     trigger::Vids const& vids(teHndl->filterIds(iF));
     trigger::Keys const& keys(teHndl->filterKeys(iF));
     for(unsigned iI = 0; iI < vids.size(); iI++){
       cout << keys.at(iI) << "(" << vids.at(iI) << ") ";
     }
     cout << endl;
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
TriggerEventAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerEventAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TriggerEventAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TriggerEventAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TriggerEventAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TriggerEventAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerEventAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerEventAnalyzer);
