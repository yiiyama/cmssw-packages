#include "../interface/GenDecayFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "../interface/GenFilter.h"
#include "../interface/Utilities.h"

GenDecayFilter::GenDecayFilter(const edm::ParameterSet& _ps) :
  genParticlesTag_(_ps.getParameter<edm::InputTag>("genParticlesTag")),
  filter_(_ps.getParameter<std::string>("filterExpression")),
  veto_(_ps.getParameter<bool>("veto"))
{
  std::cout << "GenDecayFilter: ";
  if(veto_) std::cout << "!(";
  std::cout << filter_.toString();
  if(veto_) std::cout << ")";
  std::cout << std::endl;
}

bool
GenDecayFilter::filter(edm::Event& _event, const edm::EventSetup&)
{
  if(_event.isRealData()) return true;

  edm::Handle<reco::GenParticleCollection> gpHndl;
  if(!_event.getByLabel(genParticlesTag_, gpHndl)) return false;

  std::map<reco::GenParticle const*, PNode*> nodeMap;
  std::vector<PNode*> rootNodes;

  for(reco::GenParticleCollection::const_iterator genItr(gpHndl->begin()); genItr != gpHndl->end(); ++genItr){
    reco::GenParticle const& gen(*genItr);

    if(gen.numberOfMothers() == 0){
      PNode* node(setDaughters(&gen, nodeMap, 0.));
      if(node) rootNodes.push_back(node);
    }
  }

  for(unsigned iN(0); iN != rootNodes.size(); iN++)
    cleanDaughters(rootNodes[iN]);

  bool pass(filter_.pass(rootNodes));
  if(veto_) pass = !pass;

  for(unsigned iN(0); iN != rootNodes.size(); iN++)
    delete rootNodes[iN];

  return pass;
}

void
GenDecayFilter::fillDescriptions(edm::ConfigurationDescriptions& _descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("genParticlesTag", edm::InputTag("genParticles"));
  desc.add<std::string>("filterExpression", "");
  desc.add<bool>("veto", false);

  _descriptions.add("genDecayFilter", desc);
}

DEFINE_FWK_MODULE(GenDecayFilter);
