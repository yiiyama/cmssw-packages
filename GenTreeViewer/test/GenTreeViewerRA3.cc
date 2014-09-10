// system include files
#include <vector>
#include <iostream>
#include <map>

#include "TString.h"

// user include files
#include "SusyEvent.h"

#include "../interface/PNode.h"

double PNode::matchEta = 0.;
double PNode::matchPhi = 0.;
double PNode::matchDR = 0.;

void
viewGenTreeRA3(susy::Event const& _event, double _ptThreshold = 0., bool _showMass = false)
{
  susy::ParticleCollection const& particles(_event.genParticles);

  std::vector<PNode> allNodes;
  std::vector<PNode*> rootNodes;
  std::map<unsigned, unsigned> indexMap;

  for(unsigned iP(0); iP != particles.size(); ++iP){
    susy::Particle const& gen(particles[iP]);

    if(gen.status == 1 && gen.momentum.Pt() < _ptThreshold) continue;

    PNode node;
    node.pdgId = gen.pdgId;
    node.status = gen.status;
    node.mass = gen.momentum.M();
    node.pt = gen.momentum.Pt();
    node.eta = node.pt > 0. ? gen.momentum.Eta() : (gen.momentum.Z() > 0. ? 10000. : -10000.);
    node.phi = gen.momentum.Phi();

    allNodes.push_back(node);
    indexMap[iP] = allNodes.size() - 1;
  }

  for(unsigned iP(0); iP != particles.size(); ++iP){
    susy::Particle const& gen(particles[iP]);

    if(gen.status == 1 && gen.momentum.Pt() < _ptThreshold) continue;

    unsigned iN(indexMap[iP]);

    PNode& node(allNodes[iN]);

    if(gen.motherIndex == -1) rootNodes.push_back(&node);
    else{
      unsigned iM(indexMap[gen.motherIndex]);
      node.mother = &allNodes[iM];
      allNodes[iM].daughters.push_back(&node);
    }
  }

  std::cout << "--- CLEANED DECAY TREE ---" << std::endl << std::endl;
  for(unsigned iN(0); iN < rootNodes.size(); iN++){
    std::cout << rootNodes[iN]->print(true, _showMass);
    std::cout << std::endl;
  }
  std::cout << std::endl;
}
