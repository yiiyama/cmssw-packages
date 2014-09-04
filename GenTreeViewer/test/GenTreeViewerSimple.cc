#include <vector>
#include <iostream>
#include <sstream>

#include "TTree.h"

#include "../interface/PNode.h"

double PNode::matchEta = 0.;
double PNode::matchPhi = 0.;
double PNode::matchDR = 0.;

bool gViewGenTreeStringOutput(false);

std::vector<PNode*>
formDecayTree(unsigned _size, unsigned short const* _status, short const* _motherIndex, int const* _pdgId, float const* _pt, float const* _eta, float const* _phi, float const* _mass, std::vector<PNode>& _allNodes)
{
  std::vector<PNode*> rootNodes;

  _allNodes.resize(_size);

  for(unsigned iP(0); iP != _size; ++iP){
    PNode& node(_allNodes[iP]);
    node.pdgId = _pdgId[iP];
    node.status = _status[iP];
    if(_mass) node.mass = _mass[iP];
    node.pt = _pt[iP];
    node.eta = _eta[iP];
    node.phi = _phi[iP];
  }

  for(unsigned iP(0); iP != _size; ++iP){
    PNode& node(_allNodes[iP]);

    if(_motherIndex[iP] == -1) rootNodes.push_back(&node);
    else{
      node.mother = &_allNodes[_motherIndex[iP]];
      _allNodes[_motherIndex[iP]].daughters.push_back(&node);
    }
  }

  return rootNodes;
}

char const*
viewGenTreeSimple(unsigned _size, unsigned short const* _status, short const* _motherIndex, int const* _pdgId, float const* _pt = 0, float const* _eta = 0, float const* _phi = 0, float const* _mass = 0)
{
  std::vector<PNode> allNodes;
  std::vector<PNode*> rootNodes(formDecayTree(_size, _status, _motherIndex, _pdgId, _pt, _eta, _phi, _mass, allNodes));

  std::ostream* out(&std::cout);
  std::stringstream ss;
  if(gViewGenTreeStringOutput) out = &ss;

  (*out) << "--- CLEANED DECAY TREE ---" << std::endl << std::endl;
  for(unsigned iN(0); iN < rootNodes.size(); iN++){
    (*out) << rootNodes[iN]->print(_pt != 0 && _eta != 0 && _phi != 0, _mass != 0);
    (*out) << std::endl;
  }
  (*out) << std::endl;

  return ss.str().c_str();
}

char const*
viewGenTreeSimple(TTree* _eventVars, long _iEntry = 0, bool _showMass = false, double _eta = -100., double _phi = 0.)
{
  if(!_eventVars) return "";

  if(_eta > -5.){
    PNode::matchDR = 0.3;
    PNode::matchEta = _eta;
    PNode::matchPhi = _phi;
  }

  unsigned size;
  unsigned short status[2048];
  short motherIndex[2048];
  int pdgId[2048];
  float pt[2048];
  float eta[2048];
  float phi[2048];
  float mass[2048];

  _eventVars->ResetBranchAddresses();

  _eventVars->SetBranchAddress("gen.size", &size);
  _eventVars->SetBranchAddress("gen.status", status);
  _eventVars->SetBranchAddress("gen.motherIndex", motherIndex);
  _eventVars->SetBranchAddress("gen.pdgId", pdgId);
  _eventVars->SetBranchAddress("gen.pt", pt);
  _eventVars->SetBranchAddress("gen.eta", eta);
  _eventVars->SetBranchAddress("gen.phi", phi);
  _eventVars->SetBranchAddress("gen.mass", mass);

  _eventVars->GetEntry(_iEntry);

  char const* result(viewGenTreeSimple(size, status, motherIndex, pdgId, pt, eta, phi, _showMass ? mass : 0));

  _eventVars->ResetBranchAddresses();

  return result;
}
