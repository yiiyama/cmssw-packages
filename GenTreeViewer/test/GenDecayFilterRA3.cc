#include "SusyEvent.h"

#include "../interface/GenFilter.h"

#include "TString.h"

#include <vector>

class GenDecayFilterRA3 {
public:
  GenDecayFilterRA3(TString const& _filterExpr) : filter_(_filterExpr) {}
  ~GenDecayFilterRA3() {}

  bool pass(susy::Event const&);

private:
  GenFilter filter_;
};

bool
GenDecayFilterRA3::pass(susy::Event const& _event)
{
  susy::ParticleCollection const& particles(_event.genParticles);

  std::vector<PNode> allNodes(particles.size());
  std::vector<PNode*> rootNodes;

  for(unsigned iP(0); iP != particles.size(); ++iP){
    susy::Particle const& gen(particles[iP]);

    PNode& node(allNodes[iP]);
    node.pdgId = gen.pdgId;
    node.status = gen.status;
  }

  for(unsigned iP(0); iP != particles.size(); ++iP){
    susy::Particle const& gen(particles[iP]);

    PNode& node(allNodes[iP]);

    if(gen.motherIndex == -1) rootNodes.push_back(&node);
    else{
      node.mother = &allNodes[gen.motherIndex];
      allNodes[gen.motherIndex].daughters.push_back(&node);
    }
  }

  return filter_.pass(rootNodes);
}
