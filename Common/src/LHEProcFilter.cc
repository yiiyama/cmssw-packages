#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include <set>
#include <vector>

class LHEProcFilter : public edm::EDFilter {
public:
  LHEProcFilter(edm::ParameterSet const&);
  ~LHEProcFilter() {}

private:
  bool filter(edm::Event&, edm::EventSetup const&) override;
  
  edm::InputTag lheEventsTag_;
  std::set<int> allowedProcs_;
};

LHEProcFilter::LHEProcFilter(edm::ParameterSet const& _ps) :
  lheEventsTag_(_ps.getUntrackedParameter<edm::InputTag>("lheEventsTag")),
  allowedProcs_()
{
  std::vector<int> procs(_ps.getUntrackedParameter<std::vector<int> >("allowedProcs"));
  for(unsigned iP(0); iP != procs.size(); ++iP)
    allowedProcs_.insert(procs[iP]);
}

bool
LHEProcFilter::filter(edm::Event& _event, edm::EventSetup const&)
{
  edm::Handle<LHEEventProduct> lheEvtH;
  if(!_event.getByLabel(lheEventsTag_, lheEvtH)) return false;

  return allowedProcs_.find(lheEvtH->hepeup().IDPRUP) != allowedProcs_.end();
}

DEFINE_FWK_MODULE(LHEProcFilter);
