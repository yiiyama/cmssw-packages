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
  ~PNode() {}

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
  std::string print(std::vector<bool>& _isLastChildAt, bool _showMass = false)
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
  bool hasMother() { return mother != 0; }
};
