#ifndef ToolsetGenTreeViewerPNode_h
#define ToolsetGenTreeViewerPNode_h

#include <vector>
#include <bitset>
#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cmath>

struct PNode {
  PNode* mother;
  std::vector<PNode*> daughters;
  int pdgId;
  int status;
  std::bitset<15> statusBits;
  double mass;
  double pt;
  double eta;
  double phi;
  double px;
  double py;
  double pz;
  double energy;
  bool ownDaughters;

  std::string info;

  static double matchEta;
  static double matchPhi;
  static double matchDR;

  enum MomentumDispMode {
    kShowAllP,
    kShowFinalP,
    kNoP,
    nMomentumDispModes
  };

  enum MassDispMode {
    kShowAllM,
    kShowHardScatM,
    kNoM,
    nMassDispModes
  };

  enum StatusBits {
    kIsPrompt = 0,
    kIsDecayedLeptonHadron,
    kIsTauDecayProduct,
    kIsPromptTauDecayProduct,
    kIsDirectTauDecayProduct,
    kIsDirectPromptTauDecayProduct,
    kIsDirectHadronDecayProduct,
    kIsHardProcess,
    kFromHardProcess,
    kIsHardProcessTauDecayProduct,
    kIsDirectHardProcessTauDecayProduct,
    kFromHardProcessBeforeFSR,
    kIsFirstCopy,
    kIsLastCopy,
    kIsLastCopyBeforeFSR,
    nStatusBits
  };

  PNode() : mother(0), daughters(0), pdgId(0), status(0), statusBits(),
            mass(0.), pt(0.), eta(0.), phi(0.),
            px(0.), py(0.), pz(0.), energy(0.),
            ownDaughters(false), info("") {}
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
    statusBits = _rhs.statusBits;
    mass = _rhs.mass;
    pt = _rhs.pt;
    eta = _rhs.eta;
    phi = _rhs.phi;
    px = _rhs.px;
    py = _rhs.py;
    pz = _rhs.pz;
    energy = _rhs.energy;
    info = "";
    return *this;
  }
  void generateInfo(MomentumDispMode _pMode = kShowAllP, MassDispMode _mMode = kShowAllM, bool _usePtEtaPhi = true)
  {
    if(info != "")
      return;

    using namespace std;

    stringstream ss;
    ss << setw(9) << pdgId;

    if(_mMode == kShowAllM || (_mMode == kShowHardScatM && (statusBits[kIsHardProcess] || statusBits[kFromHardProcess])))
      ss << " [" << setw(6) << fixed << setprecision(2) << mass << "]";

    if(_pMode == kShowAllP || (_pMode == kShowFinalP && status == 1)){
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

    if(status == 1){
      if(matchDR > 0.){
        double dEta(eta - matchEta);
        double dPhi(phi - matchPhi);
        if(dPhi > 3.141593) dPhi -= 6.283185307;
        if(dPhi < -3.141593) dPhi += 6.283185307;
        if(std::sqrt(dEta * dEta + dPhi * dPhi) < matchDR) ss << " *";
      }
    }
    else{
      ss << " " << status << " ";
      bool first(true);
      for(unsigned iB(0); iB != nStatusBits; ++iB){
        if(!statusBits[iB])
          continue;
        if(first){
          ss << "{";
          first = false;
        }
        else
          ss << ",";
        ss << iB;
      }
      if(!first)
        ss << "}";
    }

    info = ss.str();

    for(auto* daughter : daughters)
      daughter->generateInfo(_pMode, _mMode, _usePtEtaPhi);
  }
  std::string print(std::vector<std::string>* spacers = 0, bool firstDaughter = true, bool lastDaughter = true)
  {
    using namespace std;

    stringstream ss;

    if(spacers && !firstDaughter){
      for(string& s : *spacers)
        ss << s;
    }

    ss << "+ " << info << " ";

    string spacer(info.size() + 1, ' ');
    if(lastDaughter)
      spacer = "  " + spacer;
    else
      spacer = "| " + spacer;

    if(spacers)
      spacers->push_back(spacer);
    else
      spacers = new vector<string>(1, spacer);

    if(daughters.size() != 0){
      for(unsigned iD(0); iD != daughters.size(); ++iD)
        ss << daughters[iD]->print(spacers, iD == 0, iD == daughters.size() - 1);
    }
    else{
      ss << endl;
      if(lastDaughter && spacers){
        for(string& s : *spacers)
          ss << s;
        ss << endl;
      }
    }

    spacers->pop_back();
    if(spacers->size() == 0)
      delete spacers;

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
  void cleanDaughters(bool _doDelete = false)
  {
    int motherPdg(std::abs(pdgId));

    for(unsigned iD(0); iD < daughters.size(); ++iD){
      PNode* dnode(daughters[iD]);
      dnode->cleanDaughters(_doDelete);

      unsigned nGD(dnode->daughters.size());
      bool intermediateTerminal(nGD == 0 && dnode->status != 1);
      bool noDecay(nGD == 1 && dnode->pdgId == dnode->daughters.front()->pdgId);
      int pdg(std::abs(dnode->pdgId));
      bool hadronicIntermediate(dnode->status != 1 &&
                                ((pdg / 100) % 10 != 0 || pdg == 21 || (pdg > 80 && pdg < 101)));
      bool firstHeavyHadron((motherPdg / 1000) % 10 < 4 && (motherPdg / 100) % 10 < 4 &&
                            ((pdg / 1000) % 10 >= 4 || (pdg / 100) % 10 >= 4));
      bool lightDecayingToLight(false);
      if(pdg < 4){
        unsigned iGD(0);
        for(; iGD < nGD; ++iGD)
          if(std::abs(dnode->daughters[iGD]->pdgId) > 3) break;
        lightDecayingToLight = iGD == nGD;
      }

      if(intermediateTerminal || noDecay || (hadronicIntermediate && !firstHeavyHadron) || lightDecayingToLight){
        for(unsigned iGD(0); iGD < nGD; ++iGD)
          dnode->daughters[iGD]->mother = this;
        daughters.erase(daughters.begin() + iD);
        daughters.insert(daughters.begin() + iD, dnode->daughters.begin(), dnode->daughters.end());
        dnode->daughters.clear();
        if(_doDelete)
          delete dnode;
        --iD;
      }
    }
  }
};

#endif
