#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include <iostream>

class DijetFilter : public edm::stream::EDFilter<> {
public:
  explicit DijetFilter(edm::ParameterSet const&);
  ~DijetFilter() {}

private:
  bool filter(edm::Event&, edm::EventSetup const&);

  edm::EDGetTokenT<reco::PFJetCollection> jetCollectionToken_;
  double minPt_;
  double maxEta_;
  double minMjj_;
};

DijetFilter::DijetFilter(edm::ParameterSet const& _ps) :
  jetCollectionToken_(consumes<reco::PFJetCollection>(_ps.getParameter<edm::InputTag>("source"))),
  minPt_(_ps.getParameter<double>("minPt")),
  maxEta_(_ps.getParameter<double>("maxEta")),
  minMjj_(_ps.getParameter<double>("minMjj"))
{
}

bool
DijetFilter::filter(edm::Event& _event, edm::EventSetup const&)
{
  edm::Handle<reco::PFJetCollection> handle;
  if(!_event.getByToken(jetCollectionToken_, handle))
    return false;

  auto& jets(*handle);

  double maxMjj(0.);
  for (unsigned iJ1(0); iJ1 != jets.size(); ++iJ1) {
    auto& jet1(jets.at(iJ1));
    if (jet1.pt() < minPt_ || std::abs(jet1.eta()) > maxEta_)
      continue;

    for (unsigned iJ2(iJ1 + 1); iJ2 != jets.size(); ++iJ2) {
      auto& jet2(jets.at(iJ2));
      if (jet2.pt() < minPt_ || std::abs(jet2.eta()) > maxEta_)
        continue;
      
      double mjj((jet1.p4() + jet2.p4()).mass());
      if (mjj > maxMjj)
        maxMjj = mjj;
    }
  }

  return maxMjj > minMjj_;
}

DEFINE_FWK_MODULE(DijetFilter);
