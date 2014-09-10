#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

class LHEPtFilter : public edm::EDFilter {
public:
  explicit LHEPtFilter(edm::ParameterSet const&);
  ~LHEPtFilter() {}

private:
  bool filter(edm::Event&, edm::EventSetup const&);

  edm::InputTag lheEventProductTag_;
  unsigned pdgId_;
  double ptLow_;
  double ptHigh_;
};

LHEPtFilter::LHEPtFilter(edm::ParameterSet const& _ps) :
  lheEventProductTag_(_ps.getParameter<edm::InputTag>("lheEventProductTag")),
  pdgId_(_ps.getParameter<unsigned>("pdgId")),
  ptLow_(_ps.getParameter<double>("ptLow")),
  ptHigh_(_ps.getParameter<double>("ptHigh"))
{
}

bool
LHEPtFilter::filter(edm::Event& _event, edm::EventSetup const&)
{
  edm::Handle<LHEEventProduct> prodHndl;
  if(!_event.getByLabel(lheEventProductTag_, prodHndl)) return false;

  lhef::HEPEUP const& lhe(prodHndl->hepeup());
  int iP(0);
  for(; iP != lhe.NUP; ++iP){
    if(std::abs(lhe.IDUP[iP]) != pdgId_) continue;
    double pt(std::sqrt(lhe.PUP[iP][0] * lhe.PUP[iP][0] + lhe.PUP[iP][1] * lhe.PUP[iP][1]));
    if(pt > ptLow_ && pt < ptHigh_) break;
  }
  
  return iP != lhe.NUP;
}

DEFINE_FWK_MODULE(LHEPtFilter);
