#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

class PhotonPtFilter : public edm::EDFilter {
public:
  explicit PhotonPtFilter(edm::ParameterSet const&);
  ~PhotonPtFilter() {}

private:
  bool filter(edm::Event&, edm::EventSetup const&);

  edm::InputTag photonCollectionTag_;
  double threshold_;
};

PhotonPtFilter::PhotonPtFilter(edm::ParameterSet const& _ps) :
  photonCollectionTag_(_ps.getParameter<edm::InputTag>("photonCollectionTag")),
  threshold_(_ps.getParameter<double>("threshold"))
{
}

bool
PhotonPtFilter::filter(edm::Event& _event, edm::EventSetup const&)
{
  edm::Handle<reco::PhotonCollection> phHndl;
  if(!_event.getByLabel(photonCollectionTag_, phHndl)) return false;

  unsigned iP(0);
  for(; iP != phHndl->size(); ++iP)
    if(phHndl->at(iP).pt() > threshold_) break;
  
  return iP != phHndl->size();
}

DEFINE_FWK_MODULE(PhotonPtFilter);
