#include "../interface/SimpleEventDisplay.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TString.h"
#include "TMarker.h"
#include "TEllipse.h"
#include "TMath.h"
#include "TLegendEntry.h"
#include "TList.h"
#include "TVector2.h"

#include <iostream>
#include <map>
#include <utility>
#include <cmath>

double const sizeNorm(2.);
double const etaMax(4.);

SimpleEventDisplay::SimpleEventDisplay() :
  ptThreshold_(1.),
  canvas_("evdisp", "Simple Event Display", 800, 600),
  etaPhiField_("evdisp", ";#eta;#phi", 100, -etaMax, etaMax, 100, -TMath::Pi(), TMath::Pi()),
  eventInfo_(0.82, 0.8, 0.98, 0.95, "brNDC"),
  notes_(0.82, 0.7, 0.98, 0.75, "brNDC"),
  sizeLegend_(0.82, 0.4, 0.98, 0.55),
  objLegend_(0.82, 0.05, 0.98, 0.35)
{
  canvas_.SetGrid(1, 1);
  canvas_.SetMargin(0.05, 0.2, 0.05, 0.05);

  etaPhiField_.SetStats(false);
  etaPhiField_.SetTitleSize(0.15);
  etaPhiField_.GetXaxis()->SetNdivisions(20, 1, 0);
  etaPhiField_.GetYaxis()->SetNdivisions(20, 1, 0);
  etaPhiField_.GetXaxis()->SetTitleOffset(-0.3);
  etaPhiField_.GetYaxis()->SetTitleOffset(-0.3);
  etaPhiField_.GetXaxis()->SetTitleSize(0.05);
  etaPhiField_.GetYaxis()->SetTitleSize(0.05);

  eventInfo_.SetBorderSize(1);

  notes_.SetBorderSize(1);
  notes_.SetTextAlign(12); // left, center
  notes_.SetFillStyle(0);

  formLegend_();
}

SimpleEventDisplay::~SimpleEventDisplay()
{
  formLegend_(true);
}

bool
SimpleEventDisplay::showEvent(unsigned _eventNumber, reco::GenParticleCollection const& _particles)
{
  canvas_.Clear();

  etaPhiField_.SetTitle(TString::Format("Event %d", _eventNumber));
  etaPhiField_.Draw();

  objLegend_.Draw();
  sizeLegend_.Draw();

  for(reco::GenParticleCollection::const_iterator pItr(_particles.begin()); pItr != _particles.end(); ++pItr){
    reco::GenParticle const& part(*pItr);
    if(part.status() != 1 || part.pt() < ptThreshold_) continue;

    TMarker marker;

    switch(std::abs(part.pdgId())){
    case 22:
      marker.SetMarkerColor(kGreen);
      marker.SetMarkerStyle(30);
      break;
    case 11:
      marker.SetMarkerColor(kBlue);
      marker.SetMarkerStyle(3);
      break;
    case 13:
      marker.SetMarkerColor(kRed);
      marker.SetMarkerStyle(2);
      break;
    case 12:
      marker.SetMarkerColor(kCyan);
      marker.SetMarkerStyle(3);
      break;
    case 14:
      marker.SetMarkerColor(kCyan);
      marker.SetMarkerStyle(2);
      break;
    case 16:
      marker.SetMarkerColor(kCyan);
      marker.SetMarkerStyle(26);
      break;
    default:
      marker.SetMarkerColor(kBlack);
      marker.SetMarkerStyle(20);
      break;
    }
    marker.SetMarkerSize(sizeNorm * TMath::Log10(part.pt()));

    marker.DrawMarker(part.eta(), part.phi());
  }

  return true;
}

void
SimpleEventDisplay::formLegend_(bool _clearOnly/* = false*/)
{
  TList* colorList(objLegend_.GetListOfPrimitives());
  for(int iL(0); iL != colorList->GetEntries(); ++iL){
    TLegendEntry* entry(static_cast<TLegendEntry*>(colorList->At(iL)));
    delete entry->GetObject();
  }

  TList* sizeList(sizeLegend_.GetListOfPrimitives());
  for(int iL(0); iL != sizeList->GetEntries(); ++iL){
    TLegendEntry* entry(static_cast<TLegendEntry*>(sizeList->At(iL)));
    delete entry->GetObject();
  }

  objLegend_.Clear();
  sizeLegend_.Clear();

  if(_clearOnly) return;

  TMarker* tenGeV(new TMarker(0., 0., 20));
  TMarker* hundredGeV(new TMarker(0., 0., 20));
  tenGeV->SetMarkerSize(sizeNorm);
  hundredGeV->SetMarkerSize(sizeNorm * 2.);

  sizeLegend_.SetBorderSize(0);
  sizeLegend_.SetFillStyle(0);
  sizeLegend_.SetTextAlign(32); // right, center
  sizeLegend_.AddEntry(tenGeV, "10 GeV", "P");
  sizeLegend_.AddEntry(hundredGeV, "100 GeV", "P");

  TMarker* photon(new TMarker(0., 0., 30));
  photon->SetMarkerColor(kGreen);
  objLegend_.AddEntry(photon, "#gamma", "P");

  TMarker* electron(new TMarker(0., 0., 3));
  electron->SetMarkerColor(kBlue);
  objLegend_.AddEntry(electron, "e", "P");

  TMarker* muon(new TMarker(0., 0., 2));
  muon->SetMarkerColor(kRed);
  objLegend_.AddEntry(muon, "#mu", "P");

  TMarker* eNeutrino(new TMarker(0., 0., 3));
  eNeutrino->SetMarkerColor(kCyan);
  objLegend_.AddEntry(eNeutrino, "#nu_{e}", "P");

  TMarker* muNeutrino(new TMarker(0., 0., 2));
  muNeutrino->SetMarkerColor(kCyan);
  objLegend_.AddEntry(muNeutrino, "#nu_{#mu}", "P");

  TMarker* tauNeutrino(new TMarker(0., 0., 26));
  tauNeutrino->SetMarkerColor(kCyan);
  objLegend_.AddEntry(tauNeutrino, "#nu_{#tau}", "P");

  TMarker* strayHadron(new TMarker(0., 0., 20));
  strayHadron->SetMarkerColor(kBlack);
  objLegend_.AddEntry(strayHadron, "h", "P");
}
