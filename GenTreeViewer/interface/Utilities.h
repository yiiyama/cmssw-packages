#ifndef ToolsetGenTreeViewerUtilities_h
#define ToolsetGenTreeViewerUtilities_h

#include "PNode.h"

#include <map>

namespace reco {
  class GenParticle;
}

PNode* setDaughters(reco::GenParticle const*, std::map<reco::GenParticle const*, PNode*>&, double);
void cleanDaughters(PNode*);

#endif
