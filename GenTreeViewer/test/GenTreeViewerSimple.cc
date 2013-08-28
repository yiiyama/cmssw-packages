#include <vector>
#include <iostream>

#include "../interface/PNode.h"

void
viewGenTreeSimple(unsigned _size, unsigned short const* _status, short const* _motherIndex, int const* _pdgId, float const* _pt, float const* _eta, float const* _phi, float const* _mass, bool _showMass = false)
{
  std::vector<PNode*> rootNodes;
  std::vector<PNode*> allNodes(_size);

  for(unsigned iP(0); iP != _size; ++iP){
    PNode* node(new PNode);
    node->pdgId = _pdgId[iP];
    node->status = _status[iP];
    node->mass = _mass[iP];
    node->pt = _pt[iP];
    node->eta = _eta[iP];
    node->phi = _phi[iP];

    allNodes[iP] = node;
  }

  for(unsigned iP(0); iP != _size; ++iP){
    PNode* node(allNodes[iP]);

    if(_motherIndex[iP] == -1) rootNodes.push_back(node);
    else{
      node->mother = allNodes[_motherIndex[iP]];
      allNodes[_motherIndex[iP]]->daughters.push_back(node);
    }
  }

  std::cout << "--- CLEANED DECAY TREE ---" << std::endl << std::endl;
  for(unsigned iN(0); iN < rootNodes.size(); iN++){
    std::cout << rootNodes[iN]->print(_showMass);
    std::cout << std::endl;
  }
  std::cout << std::endl;

  for(unsigned iN(0); iN < rootNodes.size(); iN++)
    delete rootNodes[iN];
}
