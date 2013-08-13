#include "../interface/GenTreeViewer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "../interface/PNode.h"
#include "../interface/Utilities.h"

GenTreeViewer::GenTreeViewer(edm::ParameterSet const& _ps) :
  genParticlesTag_(_ps.getUntrackedParameter<edm::InputTag>("genParticlesTag", edm::InputTag("genParticles"))),
  cleaningMode_(_ps.getUntrackedParameter<int>("cleaningMode", 1)),
  minPt_(_ps.getUntrackedParameter<double>("minPt", 2.))
{
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
      std::cout << rootNodes[iN]->print(true);
      std::cout << std::endl;
    }
  }

  if(cleaningMode_ == 1 || cleaningMode_ == 2){
    std::cout << "--- CLEANED DECAY TREE ---" << std::endl << std::endl;
    for(unsigned iN(0); iN < rootNodes.size(); iN++){
      cleanDaughters(rootNodes[iN]);
      std::cout << rootNodes[iN]->print(true);
      std::cout << std::endl;
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
  desc.addUntracked<int>("cleaningMode", 1);
  desc.addUntracked<double>("minPt", 2.);

  _descriptions.add("genTreeViewer", desc);
}

DEFINE_FWK_MODULE(GenTreeViewer);
