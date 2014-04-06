// -*- C++ -*-
//
// Package:    GenTreeViewer
// Class:      GenTreeViewer
// 
/**\class GenTreeViewer GenTreeViewer.h DataInspection/GenTreeViewer/src/GenTreeViewer.h

 Description: [one line class summary]

*/
//
// Original Author:  Yutaro Iiyama,512 1-005,+41227670489,
//         Created:  Sun Jul  8 21:52:04 CEST 2012
// $Id: GenTreeViewer.cc,v 1.2 2012/10/08 12:31:54 yiiyama Exp $
//
//

#ifndef ToolsetGenTreeViewerGenTreeViewer_h
#define ToolsetGenTreeViewerGenTreeViewer_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Utilities/interface/InputTag.h"

namespace edm {
  class ParameterSet;
  class ConfigurationDescriptions;
  class Event;
  class EventSetup;
}

class GenTreeViewer : public edm::EDAnalyzer {
 public:
  explicit GenTreeViewer(edm::ParameterSet const&);
  ~GenTreeViewer() {}

  static void fillDescriptions(edm::ConfigurationDescriptions&);

 private:
  virtual void analyze(edm::Event const&, edm::EventSetup const&);

  edm::InputTag genParticlesTag_;
  bool showMomentum_;
  bool showMass_;
  bool usePtEtaPhi_;
  int cleaningMode_;
  float minPt_;
};

#endif
