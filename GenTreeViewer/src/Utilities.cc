#include "../interface/Utilities.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include <algorithm>
#include <stdexcept>

PNode*
setDaughters(reco::GenParticle const* _gen, std::map<reco::GenParticle const*, PNode*>& _nodeMap, double _minPt)
{
  if(_gen->status() == 1 && _gen->pt() < _minPt) return 0;
  unsigned nD(_gen->numberOfDaughters());
  if(nD == 0){
    if(_gen->status() != 1) return 0;
  }
  else{
    if(_gen->status() == 1) throw std::logic_error("Status 1 particle with daughters");
  }

  PNode* node(new PNode);
  node->pdgId = _gen->pdgId();
  node->status = _gen->status();
  node->statusBits = _gen->statusFlags().flags_;
  node->mass = _gen->mass();
  node->pt = _gen->pt();
  node->eta = _gen->eta();
  node->phi = _gen->phi();
  node->px = _gen->px();
  node->py = _gen->py();
  node->pz = _gen->pz();
  node->energy = _gen->energy();
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

PNode*
setDaughters(HepMC::GenParticle const* _gen, std::map<HepMC::GenParticle const*, PNode*>& _nodeMap, double _minPt)
{
  if(_gen->status() == 1 && _gen->momentum().perp() < _minPt) return 0;
  HepMC::GenVertex const* vtx(_gen->end_vertex());
  unsigned nD(vtx ? vtx->particles_out_size() : 0);
  if(nD == 0){
    if(_gen->status() != 1) return 0;
  }
  else{
    if(_gen->status() == 1) throw std::logic_error("Status 1 particle with daughters");
  }

  PNode* node(new PNode);
  node->pdgId = _gen->pdg_id();
  node->status = _gen->status();
  node->mass = _gen->generated_mass();
  node->pt = _gen->momentum().perp();
  node->eta = _gen->momentum().eta();
  node->phi = _gen->momentum().phi();
  node->px = _gen->momentum().px();
  node->py = _gen->momentum().py();
  node->pz = _gen->momentum().pz();
  node->energy = _gen->momentum().e();
  node->ownDaughters = true;
  _nodeMap[_gen] = node;

  if(!vtx) return node;

  for(HepMC::GenVertex::particles_out_const_iterator oItr(vtx->particles_out_const_begin()); oItr != vtx->particles_out_const_end(); ++oItr){
    HepMC::GenParticle const* daughter(*oItr);

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
