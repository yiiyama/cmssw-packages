#include <vector>
#include <map>
#include <iostream>

#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "../interface/PNode.h"

double PNode::matchEta(0.);
double PNode::matchPhi(0.);
double PNode::matchDR(0.);

void
GenTreeViewerBambu(mithep::MCParticleCol const* _particles, double _ptThreshold = 0., PNode::MomentumDispMode _showP = PNode::kShowFinalP, PNode::MassDispMode _showM = PNode::kShowHardScatM, bool _cleanDaughters = true)
{
  std::vector<PNode> allNodes;
  std::vector<PNode*> rootNodes;
  std::map<mithep::MCParticle const*, unsigned> indexMap;

  for (unsigned iP(0); iP != _particles->GetEntries(); ++iP) {
    auto& part(*_particles->At(iP));

    if (part.Status() == 1 && part.Pt() < _ptThreshold)
      continue;

    allNodes.emplace_back();
    auto& node(allNodes.back());
    node.pdgId = part.PdgId();
    node.status = part.Status();
    for (unsigned iB(0); iB != PNode::nStatusBits; ++iB)
      node.statusBits.set(iB, part.StatusFlag(iB));
    node.mass = part.Mass();
    node.pt = part.Pt();
    node.eta = node.pt > 0. ? part.Eta() : (part.Pz() > 0. ? 10000. : -10000.);
    node.phi = part.Phi();

    indexMap[&part] = allNodes.size() - 1;
  }

  for(unsigned iP(0); iP != _particles->GetEntries(); ++iP){
    auto& part(*_particles->At(iP));

    if(part.Status() == 1 && part.Pt() < _ptThreshold) continue;

    unsigned iN(indexMap[&part]);

    PNode& node(allNodes[iN]);

    if(part.Mother() == 0)
      rootNodes.push_back(&node);
    else{
      unsigned iM(indexMap[part.Mother()]);
      node.mother = &allNodes[iM];
      allNodes[iM].daughters.push_back(&node);
    }
  }
  
  if (_cleanDaughters)
    std::cout << "--- CLEANED DECAY TREE ---" << std::endl << std::endl;
  else
    std::cout << "--- FULL DECAY TREE ---" << std::endl << std::endl;

  for (auto* rootNode : rootNodes) {
    if (_cleanDaughters)
      rootNode->cleanDaughters(false);
    rootNode->generateInfo(_showP, _showM, true);
    std::cout << rootNode->print();
  }
  std::cout << std::endl;
}
