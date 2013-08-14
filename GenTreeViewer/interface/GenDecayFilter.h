// -*- C++ -*-
//
// Package:    GenDecayFilter
// Class:      GenDecayFilter
// 
/*
 *\class GenDecayFilter GenDecayFilter.cc Toolset/GenTreeViewer/src/GenDecayFilter.cc
 Description: EDFilter based on decay chain information. Filter configuration is done through a single string

 Implementation:
 [Notes on implementation]
*/
//
// Original Author:  Yutaro Iiyama,512 1-005,+41227670489,
//         Created:  Sun Aug 11 14:14:08 CEST 2013
// $Id$
//
//

#ifndef ToolsetGenTreeViewerGenDecayFilter_h
#define ToolsetGenTreeViewerGenDecayFilter_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "PNode.h"

#include "TString.h"

#include <vector>

namespace edm {
  class ParameterSet;
  class ConfigurationDescriptions;
  class Event;
  class EventSetup;
}

struct DecayNode {
  DecayNode() :
    pdgId(0),
    chargeSensitive(false),
    passIfMatch(true)
  {}
  DecayNode(int _pdgId, bool _chargeSensitive, bool _passIfMatch) :
    pdgId(_pdgId),
    chargeSensitive(_chargeSensitive),
    passIfMatch(_passIfMatch)
  {}

  int pdgId;
  bool chargeSensitive;
  bool passIfMatch;
};

class GenFilter {
 public:
  GenFilter(TString const&);
  ~GenFilter();

  TString toString() const;
  bool pass(std::vector<PNode*> const&) const;

  enum Types {
    kAtomic,
    kNOT,
    kOR,
    kAND,
    nTypes
  };

 private:
  bool chainMatch_(unsigned, std::vector<PNode*> const&) const;
  bool oneToOneMatch_(DecayNode const&, PNode&) const;
  std::vector<PNode*> anyDecayMatch_(DecayNode const&, PNode&, std::vector<std::vector<PNode*> >* = 0) const;
  bool vetoMatch_(DecayNode const&, std::vector<PNode*> const&) const;

  std::vector<DecayNode> decayChain_;
  std::vector<GenFilter*> subfilters_;
  Types type_;
};

class GenDecayFilter : public edm::EDFilter {
public:
  explicit GenDecayFilter(edm::ParameterSet const&);
  ~GenDecayFilter() {}

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  virtual bool filter(edm::Event&, edm::EventSetup const&);
      
  edm::InputTag genParticlesTag_;
  GenFilter filter_;
  bool veto_;
};

#endif
