#include "../interface/GenTreeViewer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "../interface/Utilities.h"

#include <iostream>
#include <stdexcept>

GenTreeViewer::GenTreeViewer(edm::ParameterSet const& _ps) :
  genParticlesTag_(_ps.getUntrackedParameter<edm::InputTag>("genParticlesTag", edm::InputTag("genParticles"))),
  pMode_(PNode::nMomentumDispModes),
  mMode_(PNode::nMassDispModes),
  usePtEtaPhi_(_ps.getUntrackedParameter<bool>("usePtEtaPhi", true)),
  cleaningMode_(_ps.getUntrackedParameter<int>("cleaningMode", 1)),
  minPt_(_ps.getUntrackedParameter<double>("minPt", 2.))
{
  if(_ps.existsAs<bool>("showMomentum", false)){
    if(_ps.getUntrackedParameter<bool>("showMomentum"))
      pMode_ = PNode::kShowAllP;
  }
  else if(_ps.existsAs<std::string>("showMomentum", false)){
    std::string showMomentum(_ps.getUntrackedParameter<std::string>("showMomentum"));
    if(showMomentum == "All")
      pMode_ = PNode::kShowAllP;
    else if(showMomentum == "Final")
      pMode_ = PNode::kShowFinalP;
    else if(showMomentum == "None")
      pMode_ = PNode::kNoP;
    else{
      std::cerr << "Unrecognized showMomentum value " << showMomentum << std::endl;
      throw std::runtime_error("Configuration");
    }
  }

  if(_ps.existsAs<bool>("showMass", false)){
    if(_ps.getUntrackedParameter<bool>("showMass"))
      pMode_ = PNode::kShowAllP;
  }
  else if(_ps.existsAs<std::string>("showMass", false)){
    std::string showMass(_ps.getUntrackedParameter<std::string>("showMass"));
    if(showMass == "All")
      pMode_ = PNode::kShowAllP;
    else if(showMass == "HardScat")
      pMode_ = PNode::kShowFinalP;
    else if(showMass == "None")
      pMode_ = PNode::kNoP;
    else{
      std::cerr << "Unrecognized showMass value " << showMass << std::endl;
      throw std::runtime_error("Configuration");
    }
  }
}

void
GenTreeViewer::analyze(edm::Event const& _event, edm::EventSetup const&)
{
  edm::Handle<reco::GenParticleCollection> gpHndl;
  if(!_event.getByLabel(genParticlesTag_, gpHndl)) return;

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
      rootNodes[iN]->generateInfo(pMode_, mMode_, usePtEtaPhi_);
      std::cout << rootNodes[iN]->print() << std::endl;
    }
  }

  if(cleaningMode_ == 1 || cleaningMode_ == 2){
    std::cout << "--- CLEANED DECAY TREE ---" << std::endl << std::endl;
    for(unsigned iN(0); iN < rootNodes.size(); iN++){
      rootNodes[iN]->cleanDaughters();
      rootNodes[iN]->generateInfo(pMode_, mMode_, usePtEtaPhi_);
      std::cout << rootNodes[iN]->print() << std::endl;
    }
  }
  std::cout << std::endl;

  for(unsigned iN(0); iN != rootNodes.size(); ++iN)
    delete rootNodes[iN];
}

void
GenTreeViewer::fillDescriptions(edm::ConfigurationDescriptions& _descriptions)
{
  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("genParticlesTag", edm::InputTag("genParticles"));
  desc.addUntracked<std::string>("showMomentum", "All");
  desc.addUntracked<std::string>("showMass", "All");
  desc.addUntracked<bool>("usePtEtaPhi", true);
  desc.addUntracked<int>("cleaningMode", 1);
  desc.addUntracked<double>("minPt", 2.);

  _descriptions.add("genTreeViewer", desc);
}

DEFINE_FWK_MODULE(GenTreeViewer);
