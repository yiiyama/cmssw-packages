#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "TH1D.h"

class EcalSimHitFilter : public edm::EDFilter {
public:
  explicit EcalSimHitFilter(edm::ParameterSet const&);
  ~EcalSimHitFilter() {}

private:
  void beginJob() override;
  bool filter(edm::Event&, edm::EventSetup const&) override;

  TH1* hNHitEB_ = 0;
  TH1* hNHitEE_ = 0;
  TH1* hEEB_ = 0;
  TH1* hEEE_ = 0;

  edm::EDGetTokenT<edm::PCaloHitContainer> ebContainerToken_;
  edm::EDGetTokenT<edm::PCaloHitContainer> eeContainerToken_;
  unsigned minN_ = 0;
  unsigned minNEB_ = 0;
  unsigned minNEE_ = 0;
  double thresholdEB_ = 0.;
  double thresholdEE_ = 0.;
  double minE_ = 0.;
  double minEEB_ = 0.;
  double minEEE_ = 0.;
};

EcalSimHitFilter::EcalSimHitFilter(edm::ParameterSet const& _cfg) :
  ebContainerToken_(consumes<edm::PCaloHitContainer>(_cfg.getParameter<edm::InputTag>("ebHits"))),
  eeContainerToken_(consumes<edm::PCaloHitContainer>(_cfg.getParameter<edm::InputTag>("eeHits"))),
  minN_(_cfg.getParameter<int>("minN")),
  minNEB_(_cfg.getParameter<int>("minNEB")),
  minNEE_(_cfg.getParameter<int>("minNEE")),
  thresholdEB_(_cfg.getParameter<double>("thresholdEB")),
  thresholdEE_(_cfg.getParameter<double>("thresholdEE")),
  minE_(_cfg.getParameter<double>("minE")),
  minEEB_(_cfg.getParameter<double>("minEEB")),
  minEEE_(_cfg.getParameter<double>("minEEE"))
{
}

void
EcalSimHitFilter::beginJob()
{
  edm::Service<TFileService> fs;
  hNHitEB_ = fs->make<TH1D>("hNHitEB", "Number of EB hits", 100, 0., 1000.);
  hNHitEE_ = fs->make<TH1D>("hNHitEE", "Number of EE hits", 100, 0., 1000.);
  hEEB_ = fs->make<TH1D>("hEnergyEB", "EB", 100, 0., 100.);
  hEEE_ = fs->make<TH1D>("hEnergyEE", "EE", 100, 0., 100.);
}

bool
EcalSimHitFilter::filter(edm::Event& _event, edm::EventSetup const&)
{
  edm::Handle<edm::PCaloHitContainer> ebHandle;
  if (!_event.getByToken(ebContainerToken_, ebHandle))
    throw cms::Exception("ProductNotFound") << "EB Hits";

  edm::Handle<edm::PCaloHitContainer> eeHandle;
  if (!_event.getByToken(eeContainerToken_, eeHandle))
    throw cms::Exception("ProductNotFound") << "EE Hits";

  unsigned nEB(0);
  double eEB(0.);

  unsigned nEE(0);
  double eEE(0.);

  auto& ebHits(*ebHandle);
  for (auto&& hit : ebHits) {
    double energy(hit.energy());
    if (energy > thresholdEB_)
      ++nEB;

    eEB += energy;
  }
 
  auto& eeHits(*eeHandle);
  for (auto&& hit : eeHits) {
    double energy(hit.energy());
    if (energy > thresholdEE_)
      ++nEE;

    eEE += energy;
  }

  hEEB_->Fill(eEB);
  hEEE_->Fill(eEE);
  hNHitEB_->Fill(nEB);
  hNHitEE_->Fill(nEE);

  if (nEB < minNEB_)
    return false;
  if (eEB < minEEB_)
    return false;
  if (nEE < minNEE_)
    return false;
  if (eEE < minEEE_)
    return false;

  if (nEB + nEE < minN_)
    return false;
  if (eEB + eEE < minE_)
    return false;

  return true;
}

DEFINE_FWK_MODULE(EcalSimHitFilter);
