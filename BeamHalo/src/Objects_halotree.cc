#include "Toolset/BeamHalo/interface/Objects_halotree.h"
#include "TTree.h"

halotree::Met::Met(Met const& _src) :
  name_(_src.name_),
  met(_src.met),
  phi(_src.phi),
  sumEt(_src.sumEt)
{
}

void
halotree::Met::setStatus(TTree& _tree, Bool_t _status, flatutils::BranchList const& _branches/* = {"*"}*/, Bool_t _whitelist/* = kTRUE*/)
{
  flatutils::setStatus(_tree, name_, "met", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, name_, "phi", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, name_, "sumEt", _status, _branches, _whitelist);
}

void
halotree::Met::setAddress(TTree& _tree, flatutils::BranchList const& _branches/* = {"*"}*/, Bool_t _whitelist/* = kTRUE*/)
{
  flatutils::setStatusAndAddress(_tree, name_, "met", &met, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, name_, "phi", &phi, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, name_, "sumEt", &sumEt, _branches, _whitelist);
}

void
halotree::Met::book(TTree& _tree, flatutils::BranchList const& _branches/* = {"*"}*/, Bool_t _whitelist)
{
  flatutils::book(_tree, name_, "met", "", 'F', &met, _branches, _whitelist);
  flatutils::book(_tree, name_, "phi", "", 'F', &phi, _branches, _whitelist);
  flatutils::book(_tree, name_, "sumEt", "", 'F', &sumEt, _branches, _whitelist);
}

halotree::Met&
halotree::Met::operator=(Met const& _rhs)
{
  met = _rhs.met;
  phi = _rhs.phi;
  sumEt = _rhs.sumEt;
  return *this;
}

void
halotree::Photon::array_data::setStatus(TTree& _tree, TString const& _name, Bool_t _status, flatutils::BranchList const& _branches/* = {"*"}*/, Bool_t _whitelist/* = kTRUE*/)
{
  flatutils::setStatus(_tree, _name, "pt", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, _name, "eta", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, _name, "phi", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, _name, "chIso", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, _name, "nhIso", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, _name, "phIso", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, _name, "sieie", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, _name, "hOverE", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, _name, "seedTime", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, _name, "pixelVeto", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, _name, "loose", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, _name, "medium", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, _name, "tight", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, _name, "isEB", _status, _branches, _whitelist);
}

void
halotree::Photon::array_data::setAddress(TTree& _tree, TString const& _name, flatutils::BranchList const& _branches/* = {"*"}*/, Bool_t _whitelist/* = kTRUE*/)
{
  flatutils::setStatusAndAddress(_tree, _name, "pt", pt, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, _name, "eta", eta, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, _name, "phi", phi, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, _name, "chIso", chIso, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, _name, "nhIso", nhIso, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, _name, "phIso", phIso, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, _name, "sieie", sieie, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, _name, "hOverE", hOverE, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, _name, "seedTime", seedTime, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, _name, "pixelVeto", pixelVeto, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, _name, "loose", loose, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, _name, "medium", medium, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, _name, "tight", tight, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, _name, "isEB", isEB, _branches, _whitelist);
}

void
halotree::Photon::array_data::book(TTree& _tree, TString const& _name, flatutils::BranchList const& _branches/* = {"*"}*/, Bool_t _whitelist/* = kTRUE*/)
{
  flatutils::book(_tree, _name, "pt", _name + ".size", 'F', pt, _branches, _whitelist);
  flatutils::book(_tree, _name, "eta", _name + ".size", 'F', eta, _branches, _whitelist);
  flatutils::book(_tree, _name, "phi", _name + ".size", 'F', phi, _branches, _whitelist);
  flatutils::book(_tree, _name, "chIso", _name + ".size", 'F', chIso, _branches, _whitelist);
  flatutils::book(_tree, _name, "nhIso", _name + ".size", 'F', nhIso, _branches, _whitelist);
  flatutils::book(_tree, _name, "phIso", _name + ".size", 'F', phIso, _branches, _whitelist);
  flatutils::book(_tree, _name, "sieie", _name + ".size", 'F', sieie, _branches, _whitelist);
  flatutils::book(_tree, _name, "hOverE", _name + ".size", 'F', hOverE, _branches, _whitelist);
  flatutils::book(_tree, _name, "seedTime", _name + ".size", 'F', seedTime, _branches, _whitelist);
  flatutils::book(_tree, _name, "pixelVeto", _name + ".size", 'O', pixelVeto, _branches, _whitelist);
  flatutils::book(_tree, _name, "loose", _name + ".size", 'O', loose, _branches, _whitelist);
  flatutils::book(_tree, _name, "medium", _name + ".size", 'O', medium, _branches, _whitelist);
  flatutils::book(_tree, _name, "tight", _name + ".size", 'O', tight, _branches, _whitelist);
  flatutils::book(_tree, _name, "isEB", _name + ".size", 'O', isEB, _branches, _whitelist);
}

halotree::Photon::Photon(array_data& _data, UInt_t _idx) :
  pt(_data.pt[_idx]),
  eta(_data.eta[_idx]),
  phi(_data.phi[_idx]),
  chIso(_data.chIso[_idx]),
  nhIso(_data.nhIso[_idx]),
  phIso(_data.phIso[_idx]),
  sieie(_data.sieie[_idx]),
  hOverE(_data.hOverE[_idx]),
  seedTime(_data.seedTime[_idx]),
  pixelVeto(_data.pixelVeto[_idx]),
  loose(_data.loose[_idx]),
  medium(_data.medium[_idx]),
  tight(_data.tight[_idx]),
  isEB(_data.isEB[_idx])
{
}

halotree::Photon::Photon(Photon const& _src) :
  pt(_src.pt),
  eta(_src.eta),
  phi(_src.phi),
  chIso(_src.chIso),
  nhIso(_src.nhIso),
  phIso(_src.phIso),
  sieie(_src.sieie),
  hOverE(_src.hOverE),
  seedTime(_src.seedTime),
  pixelVeto(_src.pixelVeto),
  loose(_src.loose),
  medium(_src.medium),
  tight(_src.tight),
  isEB(_src.isEB)
{
}

halotree::Photon&
halotree::Photon::operator=(Photon const& _rhs)
{
  pt = _rhs.pt;
  eta = _rhs.eta;
  phi = _rhs.phi;
  chIso = _rhs.chIso;
  nhIso = _rhs.nhIso;
  phIso = _rhs.phIso;
  sieie = _rhs.sieie;
  hOverE = _rhs.hOverE;
  seedTime = _rhs.seedTime;
  pixelVeto = _rhs.pixelVeto;
  loose = _rhs.loose;
  medium = _rhs.medium;
  tight = _rhs.tight;
  isEB = _rhs.isEB;
  return *this;
}

void
halotree::Muon::array_data::setStatus(TTree& _tree, TString const& _name, Bool_t _status, flatutils::BranchList const& _branches/* = {"*"}*/, Bool_t _whitelist/* = kTRUE*/)
{
  flatutils::setStatus(_tree, _name, "pt", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, _name, "eta", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, _name, "phi", _status, _branches, _whitelist);
}

void
halotree::Muon::array_data::setAddress(TTree& _tree, TString const& _name, flatutils::BranchList const& _branches/* = {"*"}*/, Bool_t _whitelist/* = kTRUE*/)
{
  flatutils::setStatusAndAddress(_tree, _name, "pt", pt, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, _name, "eta", eta, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, _name, "phi", phi, _branches, _whitelist);
}

void
halotree::Muon::array_data::book(TTree& _tree, TString const& _name, flatutils::BranchList const& _branches/* = {"*"}*/, Bool_t _whitelist/* = kTRUE*/)
{
  flatutils::book(_tree, _name, "pt", _name + ".size", 'F', pt, _branches, _whitelist);
  flatutils::book(_tree, _name, "eta", _name + ".size", 'F', eta, _branches, _whitelist);
  flatutils::book(_tree, _name, "phi", _name + ".size", 'F', phi, _branches, _whitelist);
}

halotree::Muon::Muon(array_data& _data, UInt_t _idx) :
  pt(_data.pt[_idx]),
  eta(_data.eta[_idx]),
  phi(_data.phi[_idx])
{
}

halotree::Muon::Muon(Muon const& _src) :
  pt(_src.pt),
  eta(_src.eta),
  phi(_src.phi)
{
}

halotree::Muon&
halotree::Muon::operator=(Muon const& _rhs)
{
  pt = _rhs.pt;
  eta = _rhs.eta;
  phi = _rhs.phi;
  return *this;
}

halotree::MIP::MIP(MIP const& _src) :
  name_(_src.name_),
  energy(_src.energy),
  phi(_src.phi),
  r(_src.r),
  nclus(_src.nclus)
{
}

void
halotree::MIP::setStatus(TTree& _tree, Bool_t _status, flatutils::BranchList const& _branches/* = {"*"}*/, Bool_t _whitelist/* = kTRUE*/)
{
  flatutils::setStatus(_tree, name_, "energy", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, name_, "phi", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, name_, "r", _status, _branches, _whitelist);
  flatutils::setStatus(_tree, name_, "nclus", _status, _branches, _whitelist);
}

void
halotree::MIP::setAddress(TTree& _tree, flatutils::BranchList const& _branches/* = {"*"}*/, Bool_t _whitelist/* = kTRUE*/)
{
  flatutils::setStatusAndAddress(_tree, name_, "energy", &energy, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, name_, "phi", &phi, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, name_, "r", &r, _branches, _whitelist);
  flatutils::setStatusAndAddress(_tree, name_, "nclus", &nclus, _branches, _whitelist);
}

void
halotree::MIP::book(TTree& _tree, flatutils::BranchList const& _branches/* = {"*"}*/, Bool_t _whitelist)
{
  flatutils::book(_tree, name_, "energy", "", 'F', &energy, _branches, _whitelist);
  flatutils::book(_tree, name_, "phi", "", 'F', &phi, _branches, _whitelist);
  flatutils::book(_tree, name_, "r", "", 'F', &r, _branches, _whitelist);
  flatutils::book(_tree, name_, "nclus", "", 's', &nclus, _branches, _whitelist);
}

halotree::MIP&
halotree::MIP::operator=(MIP const& _rhs)
{
  energy = _rhs.energy;
  phi = _rhs.phi;
  r = _rhs.r;
  nclus = _rhs.nclus;
  return *this;
}

