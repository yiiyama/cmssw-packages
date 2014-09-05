#ifndef ToolsetGenTreeViewerGenFilter_h
#define ToolsetGenTreeViewerGenFilter_h

#include "PNode.h"

#include "TString.h"
#include "TPRegexp.h"
#include "TObjArray.h"

#include <vector>
#include <stdexcept>
#include <iostream>

struct DecayNode {
  DecayNode() :
    pdgId(0),
    chargeSensitive(false),
    passIfMatch(true)
  {}
  DecayNode(int _pdgId, bool _chargeSensitive, bool _passIfMatch) :
    pdgId(_pdgId),
    chargeSensitive(_chargeSensitive),
    passIfMatch(_passIfMatch)
  {}

  int pdgId;
  bool chargeSensitive;
  bool passIfMatch;
};

typedef std::vector<DecayNode const*> DecayChain;

class GenFilter {
 public:
  static GenFilter* parseExpression(TString const&, unsigned* = 0);

  virtual TString toString() const = 0;
  virtual bool pass(std::vector<PNode*> const&) const = 0;
};

class DecayChainFilter : public GenFilter {
 public:
  DecayChainFilter(DecayChain const&);
  ~DecayChainFilter();

  static DecayChainFilter* parseExpression(TString const&, unsigned* = 0);

  TString toString() const;
  bool pass(std::vector<PNode*> const&) const;

 protected:
  bool chainMatch_(unsigned, std::vector<PNode*> const&) const;
  bool oneToOneMatch_(DecayNode const&, PNode&) const;
  std::vector<PNode*> anyDecayMatch_(DecayNode const&, PNode&, std::vector<std::vector<PNode*> >* = 0) const;
  bool vetoMatch_(DecayNode const&, std::vector<PNode*> const&) const;

  DecayChain chain_;
};

class NestedGenFilter : public GenFilter {
 public:
  NestedGenFilter(GenFilter const*, bool);
  ~NestedGenFilter();

  static NestedGenFilter* parseExpression(TString const&, unsigned* = 0);

  TString toString() const;
  bool pass(std::vector<PNode*> const&) const;

  bool invert() const { return invert_; }
  void setInvert(bool _inv) { invert_ = _inv; }

 protected:
  GenFilter const* content_;
  bool invert_;
};

class ChainedGenFilter : public GenFilter {
 public:
  enum Operator {
    kOR,
    kAND,
    nOperators
  };

  ChainedGenFilter(Operator);
  ~ChainedGenFilter();

  TString toString() const;
  bool pass(std::vector<PNode*> const&) const;

  void addTerm(GenFilter const* _term) { terms_.push_back(_term); }

  Operator getOperator() const { return operator_; }

 protected:
  std::vector<GenFilter const*> terms_;
  Operator operator_;
};

/*static*/
GenFilter*
GenFilter::parseExpression(TString const& _expr, unsigned* _nParsed/* = 0*/)
{
  // Example expression: (neutralino -> Z -> e OR chargino -> W -> e) AND neutralino -> gamma
  // (1000022>23>!15>11 || 1000024>24>!15>11) && 1000022>*>22

  TString expr(_expr.Strip(TString::kLeading));

  if(expr.Length() == 0)
    throw std::runtime_error(("Empty expression block in " + _expr).Data());

  unsigned nParsed(0);
  
  GenFilter* filter(0);
  if(expr(0) == '(' || expr(0, 2) == "!(")
    filter = NestedGenFilter::parseExpression(expr, &nParsed);
  else
    filter = DecayChainFilter::parseExpression(expr, &nParsed);

  expr = expr.Remove(0, nParsed).Strip(TString::kLeading);

  if(expr.Length() == 0 || expr(0) == ')'){
    if(_nParsed) *_nParsed = _expr.Length() - expr.Length();
    return filter;
  }

  // next expression must be an operator

  ChainedGenFilter* chainFilter(0);

  if(expr(0, 2) == "||")
    chainFilter = new ChainedGenFilter(ChainedGenFilter::kOR);
  else if(expr(0, 2) == "&&")
    chainFilter = new ChainedGenFilter(ChainedGenFilter::kAND);
  else
    throw std::runtime_error(("Syntax error in " + _expr).Data());

  chainFilter->addTerm(filter);

  expr = expr.Remove(0, 2).Strip(TString::kLeading);

  bool expectTerm(true);
  while(expr.Length() != 0){
    ChainedGenFilter::Operator op(ChainedGenFilter::nOperators);
    if(expr(0, 2) == "||")
      op = ChainedGenFilter::kOR;
    else if(expr(0, 2) == "&&")
      op = ChainedGenFilter::kAND;

    if(op == ChainedGenFilter::nOperators){
      if(expr(0) == ')') break;

      if(!expectTerm)
        throw std::runtime_error(("Syntax error in " + _expr).Data());

      if(expr(0) == '(' || expr(0, 2) == "!(")
        chainFilter->addTerm(NestedGenFilter::parseExpression(expr, &nParsed));
      else
        chainFilter->addTerm(DecayChainFilter::parseExpression(expr, &nParsed));

      expr = expr.Remove(0, nParsed).Strip(TString::kLeading);

      expectTerm = false;
    }
    else{
      if(expectTerm)
        throw std::runtime_error(("Syntax error in " + _expr).Data());
      if(op != chainFilter->getOperator())
        throw std::runtime_error(("Mixed and/or in " + _expr).Data());

      expr = expr.Remove(0, 2).Strip(TString::kLeading);

      expectTerm = true;
    }
  }
  
  if(_nParsed) *_nParsed = _expr.Length() - expr.Length();

  return chainFilter;
}

/*static*/
DecayChainFilter*
DecayChainFilter::parseExpression(TString const& _expr, unsigned* _nParsed/* = 0*/)
{
  TPRegexp particlePat("([!]?)([+-]?)([0-9]+|b|j|[*])");
  TPRegexp chainPat("^[ ]*(?:" + particlePat.GetPattern() + "+[>]?)+");

  if(!chainPat.MatchB(_expr))
    throw std::runtime_error(("Invalid decay chain in " + _expr).Data());

  DecayChain chain;

  chain.push_back(new DecayNode(0, false, true));

  TString chainStr(_expr(chainPat));
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
      throw std::runtime_error(("Wildcard as the end node in " + _expr).Data());
    if((pdgId == 0 || !passIfMatch) && (chain.back()->pdgId == 0 || !chain.back()->passIfMatch))
      throw std::runtime_error(("Consequtive wildcarding " + _expr).Data());

    chain.push_back(new DecayNode(pdgId, sign.Length() != 0, passIfMatch));
  }

  delete particles;

  if(_nParsed) *_nParsed = chainStr.Length();

  return new DecayChainFilter(chain);
}

/*static*/
NestedGenFilter*
NestedGenFilter::parseExpression(TString const& _expr, unsigned* _nParsed/* = 0*/)
{
  TString expr(_expr.Strip(TString::kLeading));

  bool invert(expr(0) == '!');
  expr = expr.Remove(0, expr.Index("(") + 1).Strip(TString::kLeading);

  unsigned nParsed(0);
  GenFilter* content(GenFilter::parseExpression(expr, &nParsed));

  expr = expr.Remove(0, nParsed).Strip(TString::kLeading);

  if(expr(0) != ')')
    throw std::runtime_error(("Syntax error in " + _expr).Data());

  expr = expr.Remove(0, 1).Strip(TString::kLeading);

  if(_nParsed) *_nParsed = _expr.Length() - expr.Length();

  NestedGenFilter* nested(dynamic_cast<NestedGenFilter*>(content));
  if(nested){
    nested->setInvert(nested->invert() != invert);
    return nested;
  }

  return new NestedGenFilter(content, invert);
}


DecayChainFilter::DecayChainFilter(DecayChain const& _chain) :
  chain_(_chain)
{
}

DecayChainFilter::~DecayChainFilter()
{
  for(unsigned iD(0); iD != chain_.size(); ++iD)
    delete chain_[iD];
}

TString
DecayChainFilter::toString() const
{
  TString result("");

  for(unsigned iD(0); iD < chain_.size(); ++iD){
    int pdgId(chain_[iD]->pdgId);
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
      if(chain_[iD]->chargeSensitive)
        pdgStr = TString::Format("%+d", pdgId);
      else
        pdgStr = TString::Format("%d", pdgId);
      break;
    }

    if(!chain_[iD]->passIfMatch) result += "!";
    result += pdgStr;

    if(iD != chain_.size() - 1)
      result += ">";
  }

  return result;
}

bool
DecayChainFilter::pass(std::vector<PNode*> const& _rootNodes) const
{
  return chainMatch_(0, _rootNodes);
}

bool
DecayChainFilter::chainMatch_(unsigned _decayStage, std::vector<PNode*> const& _list) const
{
  enum MatchType {
    kOneToOne,
    kAny,
    kVeto,
    nMatchType
  };

  if(_decayStage == chain_.size()) return true;

  if(_list.size() == 0) return false;

  MatchType type;
  if(chain_[_decayStage]->pdgId == 0) type = kAny;
  else if(!chain_[_decayStage]->passIfMatch) type = kVeto;
  else type = kOneToOne;

  if(type == kAny || type == kVeto) ++_decayStage;

  DecayNode const& dnode(*chain_[_decayStage]);

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
          if(vetoMatch_(*chain_[_decayStage - 1], chain[iM])) continue;
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
DecayChainFilter::oneToOneMatch_(DecayNode const& _dnode, PNode& _pnode) const
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
DecayChainFilter::anyDecayMatch_(DecayNode const& _dnode, PNode& _pnode, std::vector<std::vector<PNode*> >* _chain/* = 0*/) const
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
DecayChainFilter::vetoMatch_(DecayNode const& _dnode, std::vector<PNode*> const& _chain) const
{
  for(unsigned iC(0); iC != _chain.size(); ++iC)
    if(oneToOneMatch_(_dnode, *_chain[iC])) return true;
  return false;
}


NestedGenFilter::NestedGenFilter(GenFilter const* _content, bool _invert) :
  content_(_content),
  invert_(_invert)
{
}

NestedGenFilter::~NestedGenFilter()
{
}

TString
NestedGenFilter::toString() const
{
  TString result("");

  if(invert_) result += "!";
  result += "(";
  result += content_->toString();
  result += ")";

  return result;
}

bool
NestedGenFilter::pass(std::vector<PNode*> const& _rootNodes) const
{
  bool result(content_->pass(_rootNodes));
  return invert_ ? !result : result;
}


ChainedGenFilter::ChainedGenFilter(Operator _operator) :
  terms_(0),
  operator_(_operator)
{
}

ChainedGenFilter::~ChainedGenFilter()
{
  for(unsigned iT(0); iT != terms_.size(); ++iT)
    delete terms_[iT];
}

TString
ChainedGenFilter::toString() const
{
  TString result("");

  for(unsigned iT(0); iT != terms_.size(); ++iT){
    result += terms_[iT]->toString();
    if(iT != terms_.size() - 1){
      if(operator_ == kOR) result += " || ";
      else result += " && ";
    }
  }

  return result;
}

bool
ChainedGenFilter::pass(std::vector<PNode*> const& _rootNodes) const
{
  unsigned iT(0);
  for(; iT != terms_.size(); ++iT){
    if(operator_ == kOR && terms_[iT]->pass(_rootNodes)) break;
    else if(operator_ == kAND && !terms_[iT]->pass(_rootNodes)) break;
  }
  return (operator_ == kOR && iT != terms_.size()) || (operator_ == kAND && iT == terms_.size());
}

#endif
