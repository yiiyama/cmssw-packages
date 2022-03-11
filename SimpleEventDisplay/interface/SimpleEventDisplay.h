#ifndef SimpleEventDisplay_h
#define SimpleEventDisplay_h

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TCanvas.h"
#include "TH2F.h"
#include "TPaveText.h"
#include "TLegend.h"

class SimpleEventDisplay {
public:
  SimpleEventDisplay();
  ~SimpleEventDisplay();

  void setPtThreshold(double _pt) { ptThreshold_ = _pt; }

  bool showEvent(unsigned, reco::GenParticleCollection const&);

  void print(char const* _fileName) { canvas_.Print(_fileName); }

private:
  void formLegend_(bool = false);

  double ptThreshold_;

  TCanvas canvas_;
  TH2F etaPhiField_;
  TPaveText eventInfo_;
  TPaveText notes_;

  TLegend sizeLegend_;
  TLegend objLegend_;
};

#endif
