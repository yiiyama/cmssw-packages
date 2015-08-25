#ifndef ToolsetGenTreeViewerUtilities_h
#define ToolsetGenTreeViewerUtilities_h

#include "PNode.h"

#include <map>

namespace reco {
  class GenParticle;
}
namespace HepMC {
  class GenParticle;
}

PNode* setDaughters(reco::GenParticle const*, std::map<reco::GenParticle const*, PNode*>&, double);
PNode* setDaughters(HepMC::GenParticle const*, std::map<HepMC::GenParticle const*, PNode*>&, double);

#endif
