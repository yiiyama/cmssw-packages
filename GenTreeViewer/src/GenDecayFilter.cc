#include "../interface/GenDecayFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "../interface/Utilities.h"

#include "TPRegexp.h"
#include "TObjArray.h"

GenFilter::GenFilter(TString const& _expr) :
  decayChain_(0),
  subfilters_(0),
  type_(nTypes)
{
  // Example expression: (neutralino -> Z -> e OR chargino -> W -> e) AND neutralino -> gamma
  // (1000022>23>!15>11 || 1000024>24>!15>11) && 1000022>*>22
  TPRegexp particlePat("([!]?)([+-]?)([0-9]+|b|j|[*])");
  TPRegexp chainPat("(?:" + particlePat.GetPattern() + "+[>]?)+");
  TPRegexp wrappingPat("([!]?)(\\(.*\\))");

  TString expr(_expr.Strip(TString::kBoth));

  if(!expr.Contains("&&") && !expr.Contains("||")){
    if(expr.Contains(" ")) throw cms::Exception(("incorrect syntax in " + expr).Data());

    if(wrappingPat.MatchB(expr)){
      if(expr.Index(wrappingPat) != 0) throw cms::Exception(("incorrect syntax in " + expr).Data());

      TObjArray* matches(wrappingPat.MatchS(expr));
      TString notExpr(matches->At(1)->GetName());
      TString subExpr(matches->At(2)->GetName());
      delete matches;

      expr = subExpr(1, subExpr.Length() - 2);

      if(notExpr == "!"){
        type_ = kNOT;
        subfilters_.push_back(new GenFilter(expr));
        return;
      }
    }

    type_ = kAtomic;

    if(!chainPat.MatchB(expr)) throw cms::Exception(("incorrect syntax in " + expr).Data());

    decayChain_.push_back(DecayNode(0, false, true));

    TString chainStr(expr(chainPat));
    TObjArray* particles(chainStr.Tokenize(">"));

    for(int iP(0); iP < particles->GetEntries(); ++iP){
      TString particle(particles->At(iP)->GetName());

      TObjArray* matches(particlePat.MatchS(particle));

      bool passIfMatch(TString(matches->At(1)->GetName()).Length() == 0);
      TString sign(matches->At(2)->GetName());
      TString idStr(matches->At(3)->GetName());

      delete matches;

      int pdgId(0);
      if(idStr == "j")
	pdgId = 100;
      else if(idStr == "b")
	pdgId = 500;
      else if(idStr == "*")
	pdgId = 0;
      else
	pdgId = (sign + idStr).Atoi();

      if(pdgId == 0 && (iP == 0 || iP == particles->GetEntries() - 1))
	throw cms::Exception("wildcard as the end node");
      if((pdgId == 0 || !passIfMatch) && (decayChain_.back().pdgId == 0 || !decayChain_.back().passIfMatch))
	throw cms::Exception("consequtive wildcarding");

      decayChain_.push_back(DecayNode(pdgId, sign.Length() != 0, passIfMatch));
    }

    delete particles;
    return;
  }

  while(expr.Length() > 0){
    if(expr.Index(chainPat) == 0){
      TString chainStr(expr(chainPat));
      subfilters_.push_back(new GenFilter(chainStr));
      expr = expr(chainStr.Length(), expr.Length());
      expr = expr.Strip(TString::kLeading);
    }
    else if(expr.Index(wrappingPat) == 0){
      TObjArray* matches(wrappingPat.MatchS(expr));

      TString notExpr(matches->At(1)->GetName());
      TString subExpr(matches->At(2)->GetName());

      delete matches;

      int nestDepth(1);
      int iS(1);
      while(iS < subExpr.Length() && nestDepth > 0){
	int nextOpen(subExpr.Index("(", iS));
	int nextClose(subExpr.Index(")", iS));
	if(nextOpen == -1 || nextClose == -1){
	  iS = subExpr.Length();
	  break;
	}

	if(nextOpen < nextClose){
	  iS = nextOpen + 1;
	  ++nestDepth;
	}
	else{
	  iS = nextClose + 1;
	  --nestDepth;
	}
      }
      if(nestDepth != 0) throw cms::Exception(("nested expression incorrect in " + subExpr).Data());

      if(expr.Length() == iS){ // when the whole expression is in the brackets for no reason
        expr = expr(1, iS - 1);
        expr.Strip(TString::kLeading);
        continue;
      }

      subExpr = subExpr(1, iS - 1);

      if(notExpr == "!") subfilters_.push_back(new GenFilter("!(" + subExpr + ")"));
      else subfilters_.push_back(new GenFilter(subExpr));

      expr = expr(notExpr.Length() + iS, expr.Length());
    }
    else if(expr.Index("&&") == 0){
      if(type_ == kOR) throw cms::Exception("mixed and/or");
      type_ = kAND;
      expr = expr(2, expr.Length());
    }
    else if(expr.Index("||") == 0){
      if(type_ == kAND) throw cms::Exception("mixed and/or");
      type_ = kOR;
      expr = expr(2, expr.Length());
    }

    expr = expr.Strip(TString::kLeading);
  }
}

GenFilter::~GenFilter()
{
  for(unsigned iF(0); iF < subfilters_.size(); ++iF)
    delete subfilters_[iF];
}

TString
GenFilter::toString() const
{
  TString result("");

  switch(type_){
  case kAtomic:
    for(unsigned iD(0); iD < decayChain_.size(); ++iD){
      int pdgId(decayChain_[iD].pdgId);
      TString pdgStr;
      switch(pdgId){
      case 100:
        pdgStr = "j";
	break;
      case 500:
        pdgStr = "b";
	break;
      case 0:
	pdgStr = "*";
	break;
      default:
        if(decayChain_[iD].chargeSensitive)
          pdgStr = TString::Format("%+d", pdgId);
        else
          pdgStr = TString::Format("%d", pdgId);
	break;
      }

      if(!decayChain_[iD].passIfMatch) result += "!";
      result += pdgStr;

      if(iD != decayChain_.size() - 1)
        result += ">";
    }
    break;
  case kNOT:
    result = "!(" + subfilters_[0]->toString() + ")";
    break;
  case kOR:
  case kAND:
    for(unsigned iF(0); iF < subfilters_.size(); ++iF){
      if(subfilters_[iF]->type_ != kAtomic) result += "(";
      result += subfilters_[iF]->toString();
      if(subfilters_[iF]->type_ != kAtomic) result += ")";

      if(type_ == kAND && iF < subfilters_.size() - 1)
        result += " && ";
      else if(type_ == kOR && iF < subfilters_.size() - 1)
        result += " || ";
    }
    break;
  default:
    break;
  }

  return result;
}

bool
GenFilter::pass(std::vector<PNode*> const& _rootNodes) const
{
  switch(type_){
  case kAtomic:
    if(decayChain_.size() == 0) return true;
    return chainMatch_(0, _rootNodes);
  case kNOT:
    return !subfilters_[0]->pass(_rootNodes);
  case kOR:
  case kAND:
    for(unsigned iF(0); iF < subfilters_.size(); ++iF){
      if(subfilters_[iF]->pass(_rootNodes)){
        if(type_ == kOR) return true;
      }
      else{
        if(type_ == kAND) return false;
      }
    }

    if(type_ == kAND) return true;
    else return false;
  default:
    return false;
  }
}

bool
GenFilter::chainMatch_(unsigned _decayStage, std::vector<PNode*> const& _list) const
{
  enum MatchType {
    kOneToOne,
    kAny,
    kVeto,
    nMatchType
  };

  if(_decayStage == decayChain_.size()) return true;

  if(_list.size() == 0) return false;

  MatchType type;
  if(decayChain_[_decayStage].pdgId == 0) type = kAny;
  else if(!decayChain_[_decayStage].passIfMatch) type = kVeto;
  else type = kOneToOne;

  if(type == kAny || type == kVeto) ++_decayStage;

  DecayNode const& dnode(decayChain_[_decayStage]);

  for(unsigned iP(0); iP != _list.size(); ++iP){
    PNode& pnode(*_list[iP]);

    switch(type){
    case kOneToOne:
      if(oneToOneMatch_(dnode, pnode) && chainMatch_(_decayStage + 1, pnode.daughters)) return true;
      break;
    case kAny:
      {
        std::vector<PNode*> allMatched(anyDecayMatch_(dnode, pnode));
        for(unsigned iM(0); iM != allMatched.size(); ++iM)
          if(chainMatch_(_decayStage + 1, allMatched[iM]->daughters)) return true;
      }
      break;
    case kVeto:
      {
        std::vector<std::vector<PNode*> > chain;
        std::vector<PNode*> allMatched(anyDecayMatch_(dnode, pnode, &chain));
        for(unsigned iM(0); iM != allMatched.size(); ++iM){
          if(vetoMatch_(decayChain_[_decayStage - 1], chain[iM])) continue;
          if(chainMatch_(_decayStage + 1, allMatched[iM]->daughters)) return true;
        }
      }
      break;
    default:
      break;
    }
  }

  return false;
}

bool
GenFilter::oneToOneMatch_(DecayNode const& _dnode, PNode& _pnode) const
{
  int id(_pnode.pdgId);
  int absId(std::abs(id));

  if(_dnode.pdgId == 500 && (absId == 5 || (absId / 100) % 10 == 5 || (absId / 1000) % 10 == 5)) return true;
  if(_dnode.pdgId == 100 && (absId / 100) % 10 != 0) return true;
  if(_dnode.chargeSensitive && _dnode.pdgId == id) return true;
  if(!_dnode.chargeSensitive && _dnode.pdgId == absId) return true;

  return false;
}

std::vector<PNode*>
GenFilter::anyDecayMatch_(DecayNode const& _dnode, PNode& _pnode, std::vector<std::vector<PNode*> >* _chain/* = 0*/) const
{
  std::vector<PNode*> allMatches;

  if(_chain) _chain->clear();

  if(oneToOneMatch_(_dnode, _pnode)){
    allMatches.push_back(&_pnode);
    if(_chain) _chain->push_back(std::vector<PNode*>(0));
  }
  
  for(unsigned iN(0); iN != _pnode.daughters.size(); ++iN){
    PNode& daughter(*_pnode.daughters[iN]);

    std::vector<std::vector<PNode*> > daughterChain;
    std::vector<PNode*> daughterMatches(anyDecayMatch_(_dnode, daughter, _chain ? &daughterChain : 0));

    for(unsigned iM(0); iM != daughterMatches.size(); ++iM){
      allMatches.push_back(daughterMatches[iM]);
      if(_chain){
        daughterChain[iM].insert(daughterChain[iM].begin(), &_pnode);
        _chain->push_back(daughterChain[iM]);
      }
    }
  }

  return allMatches;
}

bool
GenFilter::vetoMatch_(DecayNode const& _dnode, std::vector<PNode*> const& _chain) const
{
  for(unsigned iC(0); iC != _chain.size(); ++iC)
    if(oneToOneMatch_(_dnode, *_chain[iC])) return true;
  return false;
}

GenDecayFilter::GenDecayFilter(const edm::ParameterSet& _ps) :
  genParticlesTag_(_ps.getParameter<edm::InputTag>("genParticlesTag")),
  filter_(_ps.getParameter<std::string>("filterExpression")),
  veto_(_ps.getParameter<bool>("veto"))
{
  std::cout << "GenDecayFilter: ";
  if(veto_) std::cout << "!(";
  std::cout << filter_.toString();
  if(veto_) std::cout << ")";
  std::cout << std::endl;
}

bool
GenDecayFilter::filter(edm::Event& _event, const edm::EventSetup&)
{
  if(_event.isRealData()) return true;

  edm::Handle<reco::GenParticleCollection> gpHndl;
  if(!_event.getByLabel(genParticlesTag_, gpHndl)) return false;

  std::map<reco::GenParticle const*, PNode*> nodeMap;
  std::vector<PNode*> rootNodes;

  for(reco::GenParticleCollection::const_iterator genItr(gpHndl->begin()); genItr != gpHndl->end(); ++genItr){
    reco::GenParticle const& gen(*genItr);

    if(gen.numberOfMothers() == 0){
      PNode* node(setDaughters(&gen, nodeMap, 0.));
      if(node) rootNodes.push_back(node);
    }
  }

  for(unsigned iN(0); iN != rootNodes.size(); iN++)
    cleanDaughters(rootNodes[iN]);

  bool pass(filter_.pass(rootNodes));
  if(veto_) pass = !pass;

  for(unsigned iN(0); iN != rootNodes.size(); iN++)
    delete rootNodes[iN];

  return pass;
}

void
GenDecayFilter::fillDescriptions(edm::ConfigurationDescriptions& _descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("genParticlesTag", edm::InputTag("genParticles"));
  desc.add<std::string>("filterExpression", "");
  desc.add<bool>("veto", false);

  _descriptions.add("genDecayFilter", desc);
}

DEFINE_FWK_MODULE(GenDecayFilter);
