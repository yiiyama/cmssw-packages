// -*- C++ -*-
//
// Package:    GenTreeViewer
// Class:      GenTreeViewer
// 
/**\class GenTreeViewer GenTreeViewer.cc DataInspection/GenTreeViewer/src/GenTreeViewer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yutaro Iiyama,512 1-005,+41227670489,
//         Created:  Sun Jul  8 21:52:04 CEST 2012
// $Id: GenTreeViewer.cc,v 1.2 2012/10/08 12:31:54 yiiyama Exp $
//
//


// system include files
#include <memory>
#include <algorithm>
#include <iomanip>
#include <sstream>

#include "TString.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
//
// class declaration
//

struct PNode {
  PNode() : mother(0), daughters(0), pdgId(0), status(0), mass(0.) {}
  PNode& operator=(PNode const& _rhs)
  {
    mother = _rhs.mother;
    daughters = _rhs.daughters;
    pdgId = _rhs.pdgId;
    status = _rhs.status;
    mass = _rhs.mass;
    pt = _rhs.pt;
    return *this;
  }
  PNode* mother;
  std::vector<PNode*> daughters;
  int pdgId;
  int status;
  double mass;
  double pt;
  std::string print(std::vector<bool>&);
  bool hasMother() { return mother != 0; }
};

std::string
PNode::print(std::vector<bool>& last)
{
  using namespace std;

  stringstream ss;
  ss << "+" << setw(7) << pdgId << " (" << setw(5) << fixed << setprecision(1) << pt << ") " << status << " ";
  for(unsigned iD(0); iD < daughters.size(); iD++){
    if(iD > 0){
      ss << endl;
      for(unsigned i(0); i < last.size(); i++)
        ss << (last[i] ? " " : "|") << "                  ";
    }
    last.push_back(iD == daughters.size() - 1);
    ss << daughters[iD]->print(last);
    last.pop_back();
  }

  return ss.str();
}

PNode*
setDaughters(reco::GenParticle const* gen, std::map<reco::GenParticle const*, PNode*>& nodeMap, float _minPt)
{
  if(gen->status() == 1 && gen->pt() < _minPt) return 0;
  unsigned nD(gen->numberOfDaughters());
  if(nD == 0 && gen->status() != 1) return 0;

  PNode* node(new PNode);
  node->pdgId = gen->pdgId();
  node->status = gen->status();
  node->mass = gen->mass();
  node->pt = gen->pt();
  nodeMap[gen] = node;


  for(unsigned iD(0); iD < nD; iD++){
    reco::GenParticle const* daughter(static_cast<reco::GenParticle const*>(gen->daughter(iD)));
    if(nodeMap.find(daughter) != nodeMap.end()){
      PNode* dnode(nodeMap[daughter]);
      PNode* mother(dnode->mother);
      if(mother){
        if(mother == node) continue;

        int mpdg(std::abs(mother->pdgId));
        bool mhad((mpdg / 100) % 10 != 0 || mpdg == 21 || (mpdg > 80 && mpdg < 101));
        int pdg(std::abs(node->pdgId));
        bool had((pdg / 100) % 10 != 0 || pdg == 21 || (pdg > 80 && pdg < 101));
        bool takeAway(false);
        if(had && mhad)
          takeAway = node->pt > mother->pt;
        else if(!had && mhad)
          takeAway = true;

        if(takeAway){
          dnode->mother = node;
          node->daughters.push_back(dnode);
          std::vector<PNode*>::iterator dItr(std::find(mother->daughters.begin(), mother->daughters.end(), dnode));
          mother->daughters.erase(dItr);
        }
      }
      else{
        node->daughters.push_back(dnode);
        dnode->mother = node;
      }
    }
    else{
      PNode* dnode(setDaughters(daughter, nodeMap, _minPt));
      if(dnode){
        dnode->mother = node;
        node->daughters.push_back(dnode);
      }
    }
  }

  return node;
}

void
cleanDaughters(PNode* node)
{
  std::vector<PNode*> daughters(node->daughters);
  int motherPdg(std::abs(node->pdgId));
  for(unsigned iD(0); iD < daughters.size(); iD++){
    PNode* dnode(daughters[iD]);
    cleanDaughters(dnode);

    unsigned nGD(dnode->daughters.size());
    bool intermediateTerminal(nGD == 0 && dnode->status != 1);
    bool noDecay(nGD == 1 && dnode->pdgId == dnode->daughters.front()->pdgId);
    int pdg(std::abs(dnode->pdgId));
    bool hadronicIntermediate(dnode->status != 1 &&
                              ((pdg / 100) % 10 != 0 || pdg == 21 || (pdg > 80 && pdg < 101)));
    bool firstHeavyHadron((motherPdg / 1000) % 10 < 4 && (motherPdg / 100) % 10 < 4 &&
                          ((pdg / 1000) % 10 >= 4 || (pdg / 100) % 10 >= 4));
    bool lightDecayingToLight(pdg < 4);
    for(unsigned iGD(0); iGD < nGD; iGD++){
      if(std::abs(dnode->daughters[iGD]->pdgId) > 3){
        lightDecayingToLight = false;
        break;
      }
    }
    if(intermediateTerminal || noDecay || (hadronicIntermediate && !firstHeavyHadron) || lightDecayingToLight){
      node->daughters[iD] = 0;
      for(unsigned iGD(0); iGD < nGD; iGD++){
        node->daughters.push_back(dnode->daughters[iGD]);
        dnode->daughters[iGD]->mother = node;
        daughters.push_back(dnode->daughters[iGD]);
      }
    }
  }

  for(unsigned iD(0); iD < node->daughters.size(); iD++){
    if(!node->daughters[iD]){
      node->daughters.erase(node->daughters.begin() + iD);
      iD--;
    }
  }
}


class GenTreeViewer : public edm::EDAnalyzer {
   public:
      explicit GenTreeViewer(const edm::ParameterSet&);
      ~GenTreeViewer();

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

  edm::InputTag genParticlesTag_;
  int cleaningMode_;
  float minPt_;
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
GenTreeViewer::GenTreeViewer(const edm::ParameterSet& iConfig) :
  genParticlesTag_(iConfig.getParameter<edm::InputTag>("genParticlesTag")),
  cleaningMode_(iConfig.getParameter<int>("cleaningMode")),
  minPt_(iConfig.getParameter<double>("minPt"))
{
   //now do what ever initialization is needed

}


GenTreeViewer::~GenTreeViewer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenTreeViewer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::GenParticleCollection> gpHndl;
   if(!iEvent.getByLabel(genParticlesTag_, gpHndl)) return;

   std::map<reco::GenParticle const*, PNode*> nodeMap;
   std::vector<PNode*> rootNodes;

   for(reco::GenParticleCollection::const_iterator genItr(gpHndl->begin()); genItr != gpHndl->end(); ++genItr){
     reco::GenParticle const& gen(*genItr);

     if(gen.numberOfMothers() == 0){
       PNode* node(setDaughters(&gen, nodeMap, minPt_));
       if(node) rootNodes.push_back(node);
     }
   }

   if(cleaningMode_ == 0 || cleaningMode_ == 2){
     std::cout << "=== FULL DECAY TREE ===" << std::endl << std::endl;
     for(unsigned iN(0); iN < rootNodes.size(); iN++){
       std::vector<bool> last(1, true);
       std::cout << rootNodes[iN]->print(last);
       std::cout << std::endl;
     }
   }

   if(cleaningMode_ == 1 || cleaningMode_ == 2){
     std::cout << "--- CLEANED DECAY TREE ---" << std::endl << std::endl;
     for(unsigned iN(0); iN < rootNodes.size(); iN++){
       cleanDaughters(rootNodes[iN]);
       std::vector<bool> last(1, true);
       std::cout << rootNodes[iN]->print(last);
       std::cout << std::endl;
     }
   }

   for(std::map<reco::GenParticle const*, PNode*>::iterator nItr(nodeMap.begin()); nItr != nodeMap.end(); ++nItr)
     delete nItr->second;

   std::cout << std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void 
GenTreeViewer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenTreeViewer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
GenTreeViewer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
GenTreeViewer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GenTreeViewer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GenTreeViewer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenTreeViewer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenTreeViewer);
