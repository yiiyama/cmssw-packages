#ifndef ToolsetGenTreeViewerPNode_h
#define ToolsetGenTreeViewerPNode_h

#include <vector>
#include <sstream>
#include <string>
#include <iomanip>
#include <map>
#include <algorithm>

struct PNode {
  PNode* mother;
  std::vector<PNode*> daughters;
  int pdgId;
  int status;
  double mass;
  double pt;
  double eta;
  double phi;
  bool ownDaughters;

  PNode() : mother(0), daughters(0), pdgId(0), status(0), mass(0.), pt(0.), eta(0.), phi(0.), ownDaughters(false) {}
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
    return *this;
  }
  std::string print(bool _showMass = false)
  {
    using namespace std;

    string spaceHolder;
    for(unsigned i(0); i != 31; ++i) spaceHolder += " ";
    if(_showMass)
      for(unsigned i(0); i != 9; ++i) spaceHolder += " ";

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
      ss << line << daughters[iD]->print(_showMass);
    }

    return ss.str();
  }
};

#endif
