#ifndef TreeEntries_halotree_h
#define TreeEntries_halotree_h
#include "MitFlat/DataFormats/interface/Collection.h"
#include "Toolset/BeamHalo/interface/Objects_halotree.h"

namespace halotree {

  typedef flatutils::Collection<Photon, flatutils::BaseCollection<kFALSE>> PhotonCollection;
  typedef flatutils::Collection<Muon, flatutils::BaseCollection<kFALSE>> MuonCollection;

  class Event {
  public:
    PhotonCollection photons = PhotonCollection("photons");
    Met met = Met("met");
    MuonCollection muons = MuonCollection("muons");
    MIP mip = MIP("mip");

    void setStatus(TTree&, Bool_t, flatutils::BranchList const& = {"*"}, Bool_t whitelist = kTRUE);
    void setAddress(TTree&, flatutils::BranchList const& = {"*"}, Bool_t whitelist = kTRUE);
    void book(TTree&, flatutils::BranchList const& = {"*"}, Bool_t whitelist = kTRUE);
  };

}

#endif
