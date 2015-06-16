#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include <iostream>

class EcalTPFilter : public edm::EDFilter {
public:
  explicit EcalTPFilter(edm::ParameterSet const&);
  ~EcalTPFilter() {}

private:
  bool filter(edm::Event&, edm::EventSetup const&);

  edm::EDGetTokenT<EcalTrigPrimDigiCollection> tpCollectionToken_;
  unsigned minN_;
  double minEt_;
};

EcalTPFilter::EcalTPFilter(edm::ParameterSet const& _ps) :
  tpCollectionToken_(consumes<EcalTrigPrimDigiCollection>(_ps.getParameter<edm::InputTag>("source"))),
  minN_(_ps.getParameter<int>("minN")),
  minEt_(_ps.getParameter<double>("minEt"))
{
}

bool
EcalTPFilter::filter(edm::Event& _event, edm::EventSetup const&)
{
  if (minN_ == 0)
    return true;
  
  edm::Handle<EcalTrigPrimDigiCollection> handle;
  if(!_event.getByToken(tpCollectionToken_, handle))
    return false;

  auto& tps(*handle);

  std::cout << "filter" << std::endl;

  unsigned nPass(0);
  for (auto&& tp : tps) {
    if (tp.compressedEt() > minEt_)
      ++nPass;
    if (nPass == minN_)
      return true;
  }

  std::cout << "fail" << std::endl;
    
  return false;
}

DEFINE_FWK_MODULE(EcalTPFilter);
