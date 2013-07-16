// system include files
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <vector>
#include <iostream>

#include "TString.h"

// user include files
#include "SusyAnalysis/SusyNtuplizer/src/SusyEvent.h"

#include "PNode.h"

void
viewGenTreeRA3(susy::Event& _event, bool _showMass = false)
{
   std::vector<susy::Particle>& particles(_event.genParticles);

   std::vector<PNode> allNodes(particles.size());
   std::vector<PNode*> rootNodes;

   for(unsigned iP(0); iP != particles.size(); ++iP){
     susy::Particle& gen(particles[iP]);

     PNode& node(allNodes[iP]);
     node.pdgId = gen.pdgId;
     node.status = gen.status;
     node.mass = gen.momentum.M();
     node.pt = gen.momentum.Pt();
     node.eta = node.pt > 0. ? gen.momentum.Eta() : (gen.momentum.Z() > 0. ? 10000. : -10000.);
     node.phi = gen.momentum.Phi();

     if(gen.motherIndex == -1) rootNodes.push_back(&node);
     else allNodes[gen.motherIndex].daughters.push_back(&node);
   }

   std::cout << "--- CLEANED DECAY TREE ---" << std::endl << std::endl;
   for(unsigned iN(0); iN < rootNodes.size(); iN++){
     std::vector<bool> isLastChildAt(1, true);
     std::cout << rootNodes[iN]->print(isLastChildAt, _showMass);
     std::cout << std::endl;
   }

   std::cout << std::endl;
}
