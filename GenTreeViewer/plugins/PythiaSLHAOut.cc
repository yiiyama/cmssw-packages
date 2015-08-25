#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

extern "C" {
  void fioopn_(int*, char const*, int);
  void fiocls_(int*);
}

class PythiaSLHAOut : public edm::EDAnalyzer {
public:
  explicit PythiaSLHAOut(edm::ParameterSet const&);
  ~PythiaSLHAOut();

private:
  void analyze(edm::Event const&, edm::EventSetup const&) {}

  int lun_;
};

PythiaSLHAOut::PythiaSLHAOut(edm::ParameterSet const& _ps) :
  lun_(_ps.getUntrackedParameter<int>("LUN"))
{
  std::string fileName(_ps.getUntrackedParameter<std::string>("fileName"));
  fioopn_(&lun_, fileName.c_str(), fileName.size());
}

PythiaSLHAOut::~PythiaSLHAOut()
{
  fiocls_(&lun_);
}

DEFINE_FWK_MODULE(PythiaSLHAOut);
