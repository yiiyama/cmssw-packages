#include "../interface/GenDecayFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "../interface/GenFilter.h"
#include "../interface/Utilities.h"

#include "TString.h"

GenDecayFilter::GenDecayFilter(const edm::ParameterSet& _ps) :
  sourceTag_(_ps.getParameter<edm::InputTag>("sourceTag")),
  useGenParticles_(_ps.getUntrackedParameter<bool>("useGenParticles")),
  filter_(0),
  veto_(_ps.getParameter<bool>("veto"))
{
  TString expr(_ps.getParameter<std::string>("filterExpression"));
  filter_ = GenFilter::parseExpression(expr);
  std::cout << "GenDecayFilter: ";
  if(veto_) std::cout << "!(";
  std::cout << filter_->toString();
  if(veto_) std::cout << ")";
  std::cout << std::endl;
}

GenDecayFilter::~GenDecayFilter()
{
  delete filter_;
}

bool
GenDecayFilter::filter(edm::Event& _event, const edm::EventSetup&)
{
  if(_event.isRealData()) return true;

  filter_->reset();

  std::vector<PNode*> rootNodes;

  if(useGenParticles_){
    edm::Handle<reco::GenParticleCollection> gpHndl;
    if(!_event.getByLabel(sourceTag_, gpHndl))
      throw cms::Exception("ProductNotFound") << sourceTag_;

    std::map<reco::GenParticle const*, PNode*> nodeMap;

    for(reco::GenParticleCollection::const_iterator genItr(gpHndl->begin()); genItr != gpHndl->end(); ++genItr){
      reco::GenParticle const& gen(*genItr);

      if(gen.numberOfMothers() == 0){
        PNode* node(setDaughters(&gen, nodeMap, 0.));
        if(node) rootNodes.push_back(node);
      }
    }
  }
  else{
    edm::Handle<edm::HepMCProduct> mcHndl;
    if(!_event.getByLabel(sourceTag_, mcHndl))
      throw cms::Exception("ProductNotFound") << sourceTag_;

    HepMC::GenEvent const* event(mcHndl->GetEvent());

    std::map<HepMC::GenParticle const*, PNode*> nodeMap;

    for(HepMC::GenEvent::particle_const_iterator pItr(event->particles_begin()); pItr != event->particles_end(); ++pItr){
      HepMC::GenParticle const& gen(**pItr);
      HepMC::GenVertex const* vtx(gen.production_vertex());

      if(!vtx || vtx->particles_in_size() == 0){
        PNode* node(setDaughters(&gen, nodeMap, 0.));
        if(node) rootNodes.push_back(node);
      }
    }
  }

  for(unsigned iN(0); iN != rootNodes.size(); iN++)
    rootNodes[iN]->cleanDaughters();

  bool pass(filter_->pass(rootNodes));
  if(veto_) pass = !pass;

  for(unsigned iN(0); iN != rootNodes.size(); iN++)
    delete rootNodes[iN];

  return pass;
}

void
GenDecayFilter::fillDescriptions(edm::ConfigurationDescriptions& _descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("sourceTag", edm::InputTag("genParticles"));
  desc.addUntracked<bool>("useGenParticles", true);
  desc.add<std::string>("filterExpression", "");
  desc.add<bool>("veto", false);

  _descriptions.add("genDecayFilter", desc);
}

DEFINE_FWK_MODULE(GenDecayFilter);
