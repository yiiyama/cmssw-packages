// system include files
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <vector>
#include <iostream>

#include "TString.h"

// user include files
#include "SusyAnalysis/SusyNtuplizer/src/SusyEvent.h"

struct PNode {
  PNode() : mother(0), daughters(0), pdgId(0), status(0), mass(0.) {}
  ~PNode();
  PNode& operator=(PNode const& _rhs)
  {
    mother = _rhs.mother;
    daughters = _rhs.daughters;
    pdgId = _rhs.pdgId;
    status = _rhs.status;
    mass = _rhs.mass;
    pt = _rhs.pt;
    eta = _rhs.eta;
    phi = _rhs.phi;
    return *this;
  }
  PNode* mother;
  std::vector<PNode*> daughters;
  int pdgId;
  int status;
  double mass;
  double pt;
  double eta;
  double phi;
  std::string print(std::vector<bool>&, bool = false);
  bool hasMother() { return mother != 0; }
};

PNode::~PNode()
{
  for(unsigned iD(0); iD < daughters.size(); ++iD){
    delete daughters[iD];
  }
}

std::string
PNode::print(std::vector<bool>& _isLastChildAt, bool _showMass/* = false*/)
{
  using namespace std;

  stringstream ss;
  ss << "+" << setw(8) << pdgId;
  if(_showMass)
    ss << " [" << setw(6) << fixed << setprecision(2) << mass << "]";
  ss << " (" << setw(5) << fixed << setprecision(1) << pt << ",";
  if(eta > 3.5)
    ss << "  inf,";
  else if(eta < -3.5)
    ss << " -inf,";
  else
    ss << setw(5) << fixed << setprecision(1) << eta << ",";
  ss << setw(5) << fixed << setprecision(1) << phi << ") ";
  ss << status << " ";
  for(unsigned iD(0); iD < daughters.size(); iD++){
    if(iD > 0){
      ss << endl;
      for(unsigned i(0); i < _isLastChildAt.size(); i++){
        ss << (_isLastChildAt[i] ? " " : "|") << "                               ";
        if(_showMass)
          ss << "         ";
      }
    }
    _isLastChildAt.push_back(iD == daughters.size() - 1);
    ss << daughters[iD]->print(_isLastChildAt, _showMass);
    _isLastChildAt.pop_back();
  }

  return ss.str();
}

PNode*
setDaughters(short _index, std::vector<susy::Particle>& _particles, float _minPt)
{
  susy::Particle& gen(_particles[_index]);
  if(gen.status == 1 && gen.momentum.Pt() < _minPt) return 0;

  PNode* node(new PNode);
  node->pdgId = gen.pdgId;
  node->status = gen.status;
  node->mass = gen.momentum.M();
  node->pt = gen.momentum.Pt();
  node->eta = node->pt > 0. ? gen.momentum.Eta() : (gen.momentum.Z() > 0. ? 10000. : -10000.);
  node->phi = gen.momentum.Phi();

  for(unsigned iP(0); iP < _particles.size(); ++iP){
    if(_particles[iP].motherIndex == _index){
      PNode* daughter(setDaughters(iP, _particles, _minPt));
      daughter->mother = node;
      node->daughters.push_back(daughter);
    }
  }

  return node;
}

void
viewGenTreeRA3(susy::Event& _event, float _minPt, bool _showMass = false)
{
   std::vector<PNode*> rootNodes;

   std::vector<susy::Particle>& particles(_event.genParticles);

   for(short iP(0); iP < short(particles.size()); ++iP){
     susy::Particle& gen(particles[iP]);

     if(gen.motherIndex == -1){
       PNode* node(setDaughters(iP, particles, _minPt));
       rootNodes.push_back(node);
     }
   }

   std::cout << "--- CLEANED DECAY TREE ---" << std::endl << std::endl;
   for(unsigned iN(0); iN < rootNodes.size(); iN++){
     std::vector<bool> isLastChildAt(1, true);
     std::cout << rootNodes[iN]->print(isLastChildAt, _showMass);
     std::cout << std::endl;

     delete rootNodes[iN];
   }

   std::cout << std::endl;
}
