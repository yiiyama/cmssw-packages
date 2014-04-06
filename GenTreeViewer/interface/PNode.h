#ifndef ToolsetGenTreeViewerPNode_h
#define ToolsetGenTreeViewerPNode_h

#include <vector>
#include <sstream>
#include <string>
#include <iomanip>
#include <map>
#include <algorithm>
#include <cmath>

struct PNode {
  PNode* mother;
  std::vector<PNode*> daughters;
  int pdgId;
  int status;
  double mass;
  double pt;
  double eta;
  double phi;
  double px;
  double py;
  double pz;
  double energy;
  bool ownDaughters;

  PNode() : mother(0), daughters(0), pdgId(0), status(0),
            mass(0.), pt(0.), eta(0.), phi(0.),
            px(0.), py(0.), pz(0.), energy(0.),
            ownDaughters(false) {}
  ~PNode()
  {
    if(ownDaughters){
      for(unsigned iD(0); iD != daughters.size(); ++iD)
        delete daughters[iD];
    }
  }
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
    px = _rhs.px;
    py = _rhs.py;
    pz = _rhs.pz;
    energy = _rhs.energy;
    return *this;
  }
  std::string print(bool _showMomentum = true, bool _showMass = false, bool _usePtEtaPhi = true)
  {
    using namespace std;

    stringstream ss;
    ss << "+" << setw(9) << pdgId;

    if(_showMass)
      ss << " [" << setw(6) << fixed << setprecision(2) << mass << "]";

    if(_showMomentum){
      if(_usePtEtaPhi){
        ss << " (";
        ss << setw(5) << fixed << setprecision(pt < 1000. ? 1 : 0) << pt;
        ss << ",";
        if(eta > 10.)
          ss << "  inf";
        else if(eta < -10.)
          ss << " -inf";
        else
          ss << setw(5) << fixed << setprecision(1) << eta;
        ss << ",";
        ss << setw(5) << fixed << setprecision(1) << phi;
        ss << ")";
      } 
      else{
        double p4[4] = {px, py, pz, energy};
        ss << " (";
        for(unsigned iP(0); iP != 4; ++iP){
          ss << setw(5) << fixed << setprecision(p4[iP] < 1000. ? 1 : 0) << p4[iP];
          if(iP != 3) ss << ",";
        }
        ss << ")";
      }
    }

    ss << " " << status << " ";

    string spaceHolder(ss.str().length() - 1, ' ');

    for(unsigned iD(0); iD < daughters.size(); iD++){
      string line;
      if(iD > 0){
        ss << endl;
        PNode* node(this);
        while(node){
          PNode* motherNode(node->mother);
          bool isLastChild(!motherNode);
          if(motherNode){
            std::vector<PNode*>& siblings(motherNode->daughters);
            std::vector<PNode*>::iterator nItr(std::find(siblings.begin(), siblings.end(), node));
            isLastChild = (++nItr == siblings.end());
          }
          line = (isLastChild ? " " : "|") + spaceHolder + line;

          node = motherNode;
        }
      }
      ss << line << daughters[iD]->print(_showMomentum, _showMass, _usePtEtaPhi);
    }

    return ss.str();
  }
  bool hasDescendant(int _pdgId, bool _signed = false){
    for(unsigned iD(0); iD != daughters.size(); ++iD){
      if(_signed){
        if(daughters[iD]->pdgId == _pdgId) return true;
      }
      else{
        if(std::abs(daughters[iD]->pdgId) == _pdgId) return true;
      }
      if(daughters[iD]->hasDescendant(_pdgId, _signed)) return true;
    }
    return false;
  }
  bool hasAncestor(int _pdgId, bool _signed = false){
    if(!mother) return false;

    if(_signed){
      if(mother->pdgId == _pdgId) return true;
    }
    else{
      if(std::abs(mother->pdgId) == _pdgId) return true;
    }
    return mother->hasAncestor(_pdgId, _signed);
  }
};

#endif
