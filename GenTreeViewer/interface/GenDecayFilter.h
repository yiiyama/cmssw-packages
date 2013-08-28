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

#include "GenFilter.h"

namespace edm {
  class ParameterSet;
  class ConfigurationDescriptions;
  class Event;
  class EventSetup;
}

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
