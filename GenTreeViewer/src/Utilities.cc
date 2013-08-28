#include "../interface/Utilities.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <algorithm>

PNode*
setDaughters(reco::GenParticle const* _gen, std::map<reco::GenParticle const*, PNode*>& _nodeMap, double _minPt)
{
  if(_gen->status() == 1 && _gen->pt() < _minPt) return 0;
  unsigned nD(_gen->numberOfDaughters());
  if(nD == 0 && _gen->status() != 1) return 0;

  PNode* node(new PNode);
  node->pdgId = _gen->pdgId();
  node->status = _gen->status();
  node->mass = _gen->mass();
  node->pt = _gen->pt();
  node->eta = _gen->eta();
  node->phi = _gen->phi();
  node->ownDaughters = true;
  _nodeMap[_gen] = node;

  for(unsigned iD(0); iD < nD; ++iD){
    reco::GenParticle const* daughter(static_cast<reco::GenParticle const*>(_gen->daughter(iD)));

    if(_nodeMap.find(daughter) != _nodeMap.end()){
      PNode* dnode(_nodeMap[daughter]);
      PNode* mother(dnode->mother);
      if(mother){
        if(mother == node) continue;

        int mpdg(std::abs(mother->pdgId));
        bool mhad((mpdg / 100) % 10 != 0 || mpdg == 21 || (mpdg > 80 && mpdg < 101));
        int pdg(std::abs(node->pdgId));
        bool had((pdg / 100) % 10 != 0 || pdg == 21 || (pdg > 80 && pdg < 101));
        bool takeAway(false);
        if(had && mhad)
          takeAway = node->pt > mother->pt;
        else if(!had && mhad)
          takeAway = true;

        if(takeAway){
          dnode->mother = node;
          node->daughters.push_back(dnode);
          std::vector<PNode*>::iterator dItr(std::find(mother->daughters.begin(), mother->daughters.end(), dnode));
          mother->daughters.erase(dItr);
        }
      }
      else{
        node->daughters.push_back(dnode);
        dnode->mother = node;
      }
    }
    else{
      PNode* dnode(setDaughters(daughter, _nodeMap, _minPt));
      if(dnode){
        dnode->mother = node;
        node->daughters.push_back(dnode);
      }
    }
  }

  return node;
}

void
cleanDaughters(PNode* node)
{
  int motherPdg(std::abs(node->pdgId));
  std::vector<PNode*>& daughters(node->daughters);

  for(unsigned iD(0); iD < daughters.size(); ++iD){
    PNode* dnode(daughters[iD]);
    cleanDaughters(dnode);

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
        dnode->daughters[iGD]->mother = node;
      daughters.erase(daughters.begin() + iD);
      daughters.insert(daughters.begin() + iD, dnode->daughters.begin(), dnode->daughters.end());
      dnode->daughters.clear();
      delete dnode;
      --iD;
    }
  }
}
