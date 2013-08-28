// system include files
#include <vector>
#include <iostream>

#include "TString.h"

// user include files
#include "SusyEvent.h"

#include "../interface/PNode.h"

void
viewGenTreeRA3(susy::Event const& _event, bool _showMass = false)
{
  susy::ParticleCollection const& particles(_event.genParticles);

  std::vector<PNode> allNodes(particles.size());
  std::vector<PNode*> rootNodes;

  for(unsigned iP(0); iP != particles.size(); ++iP){
    susy::Particle const& gen(particles[iP]);

    PNode& node(allNodes[iP]);
    node.pdgId = gen.pdgId;
    node.status = gen.status;
    node.mass = gen.momentum.M();
    node.pt = gen.momentum.Pt();
    node.eta = node.pt > 0. ? gen.momentum.Eta() : (gen.momentum.Z() > 0. ? 10000. : -10000.);
    node.phi = gen.momentum.Phi();
  }

  for(unsigned iP(0); iP != particles.size(); ++iP){
    susy::Particle const& gen(particles[iP]);

    PNode& node(allNodes[iP]);

    if(gen.motherIndex == -1) rootNodes.push_back(&node);
    else{
      node.mother = &allNodes[gen.motherIndex];
      allNodes[gen.motherIndex].daughters.push_back(&node);
    }
  }

  std::cout << "--- CLEANED DECAY TREE ---" << std::endl << std::endl;
  for(unsigned iN(0); iN < rootNodes.size(); iN++){
    std::cout << rootNodes[iN]->print(_showMass);
    std::cout << std::endl;
  }
  std::cout << std::endl;
}
