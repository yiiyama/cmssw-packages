#ifndef Objects_halotree_h
#define Objects_halotree_h
#include "MitFlat/DataFormats/interface/Utils.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include "TVector2.h"
#include "TString.h"
#include "Rtypes.h"
class TTree;

namespace halotree {

  typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> LorentzVectorM;

  class Met {
  public:
    Met(TString const& name) : name_(name) {}
    Met(Met const&);
    virtual ~Met() {}
    Met& operator=(Met const&);

    void setName(TString const& name) { name_ = name; }
    virtual void setStatus(TTree&, Bool_t, flatutils::BranchList const& = {"*"}, Bool_t whitelist = kTRUE);
    virtual void setAddress(TTree&, flatutils::BranchList const& = {"*"}, Bool_t whitelist = kTRUE);
    virtual void book(TTree&, flatutils::BranchList const& = {"*"}, Bool_t whitelist = kTRUE);

    TVector2 v() const { TVector2 vec; vec.SetMagPhi(met, phi); return vec; }

  protected:
    TString name_;

  public:
    Float_t met{};
    Float_t phi{};
    Float_t sumEt{};
  };

  class Photon {
  public:
    constexpr static double const chIsoCuts[2][3]{{2.67, 1.79, 1.66}, {1.79, 1.09, 1.04}};
    constexpr static double const nhIsoCuts[2][3]{{7.23, 0.16, 0.14}, {8.89, 4.31, 3.89}};
    constexpr static double const phIsoCuts[2][3]{{2.11, 1.9, 1.4}, {3.09, 1.9, 1.4}};
    constexpr static double const sieieCuts[2][3]{{0.0107, 0.01, 0.01}, {0.0272, 0.0267, 0.0265}};
    constexpr static double const hOverECuts[2][3]{{0.028, 0.012, 0.01}, {0.093, 0.023, 0.015}};

    struct array_data {
      static UInt_t const NMAX{32};

      Float_t pt[NMAX]{};
      Float_t eta[NMAX]{};
      Float_t phi[NMAX]{};
      Float_t chIso[NMAX]{};
      Float_t nhIso[NMAX]{};
      Float_t phIso[NMAX]{};
      Float_t sieie[NMAX]{};
      Float_t hOverE[NMAX]{};
      Float_t seedTime[NMAX]{};
      Bool_t pixelVeto[NMAX]{};
      Bool_t loose[NMAX]{};
      Bool_t medium[NMAX]{};
      Bool_t tight[NMAX]{};
      Bool_t isEB[NMAX]{};

      void setStatus(TTree&, TString const&, Bool_t, flatutils::BranchList const& = {"*"}, Bool_t whitelist = kTRUE);
      void setAddress(TTree&, TString const&, flatutils::BranchList const& = {"*"}, Bool_t whitelist = kTRUE);
      void book(TTree&, TString const&, flatutils::BranchList const& = {"*"}, Bool_t whitelist = kTRUE);
    };

    Photon(array_data&, UInt_t idx);
    Photon(Photon const&);
    virtual ~Photon() {}
    Photon& operator=(Photon const&);

    virtual LorentzVectorM p4() const { return LorentzVectorM(pt, eta, phi, 0.); }
    bool passCHIso(UInt_t wp) const { return chIso < (isEB ? chIsoCuts[0][wp] : chIsoCuts[1][wp]); }
    bool passNHIso(UInt_t wp) const { return nhIso < (isEB ? nhIsoCuts[0][wp] : nhIsoCuts[1][wp]); }
    bool passPhIso(UInt_t wp) const { return phIso < (isEB ? phIsoCuts[0][wp] : phIsoCuts[1][wp]); }
    bool passSieie(UInt_t wp) const { return sieie < (isEB ? sieieCuts[0][wp] : sieieCuts[1][wp]); }
    bool passHOverE(UInt_t wp) const { return hOverE < (isEB ? hOverECuts[0][wp] : hOverECuts[1][wp]); }

  public:
    Float_t& pt;
    Float_t& eta;
    Float_t& phi;
    Float_t& chIso;
    Float_t& nhIso;
    Float_t& phIso;
    Float_t& sieie;
    Float_t& hOverE;
    Float_t& seedTime;
    Bool_t& pixelVeto;
    Bool_t& loose;
    Bool_t& medium;
    Bool_t& tight;
    Bool_t& isEB;
  };

  class Muon {
  public:
    struct array_data {
      static UInt_t const NMAX{32};

      Float_t pt[NMAX]{};
      Float_t eta[NMAX]{};
      Float_t phi[NMAX]{};

      void setStatus(TTree&, TString const&, Bool_t, flatutils::BranchList const& = {"*"}, Bool_t whitelist = kTRUE);
      void setAddress(TTree&, TString const&, flatutils::BranchList const& = {"*"}, Bool_t whitelist = kTRUE);
      void book(TTree&, TString const&, flatutils::BranchList const& = {"*"}, Bool_t whitelist = kTRUE);
    };

    Muon(array_data&, UInt_t idx);
    Muon(Muon const&);
    virtual ~Muon() {}
    Muon& operator=(Muon const&);

  public:
    Float_t& pt;
    Float_t& eta;
    Float_t& phi;
  };

  class MIP {
  public:
    MIP(TString const& name) : name_(name) {}
    MIP(MIP const&);
    virtual ~MIP() {}
    MIP& operator=(MIP const&);

    void setName(TString const& name) { name_ = name; }
    virtual void setStatus(TTree&, Bool_t, flatutils::BranchList const& = {"*"}, Bool_t whitelist = kTRUE);
    virtual void setAddress(TTree&, flatutils::BranchList const& = {"*"}, Bool_t whitelist = kTRUE);
    virtual void book(TTree&, flatutils::BranchList const& = {"*"}, Bool_t whitelist = kTRUE);

  protected:
    TString name_;

  public:
    Float_t energy{};
    Float_t phi{};
    Float_t r{};
    UShort_t nclus{};
  };
}

#endif
