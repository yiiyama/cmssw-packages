#include "Toolset/BeamHalo/interface/TreeEntries_halotree.h"
#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"

void
halotree::Event::setStatus(TTree& _tree, Bool_t _status, flatutils::BranchList const& _branches/* = {"*"}*/, Bool_t _whitelist/* = kTRUE*/)
{

  photons.setStatus(_tree, _status, flatutils::subBranchList(_branches, "photons"), _whitelist);
  met.setStatus(_tree, _status, flatutils::subBranchList(_branches, "met"), _whitelist);
  muons.setStatus(_tree, _status, flatutils::subBranchList(_branches, "muons"), _whitelist);
  mip.setStatus(_tree, _status, flatutils::subBranchList(_branches, "mip"), _whitelist);
}

void
halotree::Event::setAddress(TTree& _tree, flatutils::BranchList const& _branches/* = {"*"}*/, Bool_t _whitelist/* = kTRUE*/)
{

  photons.setAddress(_tree, flatutils::subBranchList(_branches, "photons"), _whitelist);
  met.setAddress(_tree, flatutils::subBranchList(_branches, "met"), _whitelist);
  muons.setAddress(_tree, flatutils::subBranchList(_branches, "muons"), _whitelist);
  mip.setAddress(_tree, flatutils::subBranchList(_branches, "mip"), _whitelist);
}

void
halotree::Event::book(TTree& _tree, flatutils::BranchList const& _branches/* = {"*"}*/, Bool_t _whitelist/* = kTRUE*/)
{

  photons.book(_tree, flatutils::subBranchList(_branches, "photons"), _whitelist);
  met.book(_tree, flatutils::subBranchList(_branches, "met"), _whitelist);
  muons.book(_tree, flatutils::subBranchList(_branches, "muons"), _whitelist);
  mip.book(_tree, flatutils::subBranchList(_branches, "mip"), _whitelist);
}

