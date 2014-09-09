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
  virtual ~GenFilter() {}
  static GenFilter* parseExpression(TString const&, unsigned* = 0);

  virtual TString toString() const = 0;
  virtual bool pass(std::vector<PNode*> const&) const = 0;
  virtual void reset() const = 0;

 protected:
  GenFilter() {}
};

class DecayChainFilter : public GenFilter {
 public:
  DecayChainFilter(DecayChain const&);
  ~DecayChainFilter();

  static DecayChainFilter* parseExpression(TString const&, unsigned* = 0);

  TString toString() const;
  bool pass(std::vector<PNode*> const&) const;
  void reset() const;

 protected:
  bool chainMatch_(unsigned, std::vector<PNode*> const&, std::vector<PNode*>* = 0) const;
  // simply match one DNode with one PNode
  bool oneToOneMatch_(DecayNode const&, PNode&) const;
  // trace the PNode sequence all the way down and collect matching PNodes, return sequences leading to the PNode
  std::vector<std::vector<PNode*> > anyDecayMatch_(DecayNode const&, PNode&) const;
  // returns true if DNode match exists in the sequence of PNode
  bool isInSequence_(DecayNode const&, std::vector<PNode*> const&) const;

  DecayChain chain_;
  mutable std::vector<std::vector<PNode*> > matchedSequences_;
};

class NestedGenFilter : public GenFilter {
 public:
  NestedGenFilter(GenFilter const*, int);
  ~NestedGenFilter();

  static NestedGenFilter* parseExpression(TString const&, unsigned* = 0);

  TString toString() const;
  bool pass(std::vector<PNode*> const&) const;
  void reset() const;

  bool repeat() const { return repeat_; }
  void setRepeat(int _repeat) { repeat_ = _repeat; }

 protected:
  GenFilter const* content_;
  int repeat_;
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
  void reset() const;

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

  // try nested pattern first
  GenFilter* filter(NestedGenFilter::parseExpression(expr, &nParsed));
  if(!filter)
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

      GenFilter* term(NestedGenFilter::parseExpression(expr, &nParsed));
      if(!term)
        term = DecayChainFilter::parseExpression(expr, &nParsed);

      chainFilter->addTerm(term);

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

  TPRegexp nestPattern("^([0-9]+[ ]*[*]|!|)[ ]*\\(");
  if(!nestPattern.MatchB(expr)) return 0;

  TObjArray* matches(nestPattern.MatchS(expr));
  TString modifier(matches->At(1)->GetName());

  int repeat(1);
  if(modifier(0) == '!') repeat = -1;
  else if(modifier.Length() != 0) repeat = TString(modifier(0, modifier.Index("*"))).Atoi();

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
    nested->setRepeat(nested->repeat() * repeat);
    return nested;
  }

  return new NestedGenFilter(content, repeat);
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

void
DecayChainFilter::reset() const
{
  matchedSequences_.clear();
}

bool
DecayChainFilter::chainMatch_(unsigned _decayStage, std::vector<PNode*> const& _list, std::vector<PNode*>* _currentSequence/* = 0*/) const
{
  enum MatchType {
    kOneToOne,
    kAny,
    kVeto,
    nMatchType
  };

  if(_decayStage == chain_.size()){
    if(!_currentSequence) return true;

    // parsed down to the end of the chain; does the sequence match any of the already-matched sequences?
    unsigned iS(0);
    for(; iS != matchedSequences_.size(); ++iS){
      if(_currentSequence->size() != matchedSequences_[iS].size()) continue;
      unsigned iP(0);
      for(; iP != _currentSequence->size(); ++iP)
        if((*_currentSequence)[iP] != matchedSequences_[iS][iP]) break;
      if(iP == _currentSequence->size()) break;
    }
    if(iS == matchedSequences_.size()){
      matchedSequences_.push_back(*_currentSequence);
      return true;
    }

    return false;
  }

  if(_list.size() == 0) return false;

  MatchType type;
  if(chain_[_decayStage]->pdgId == 0) type = kAny;
  else if(!chain_[_decayStage]->passIfMatch) type = kVeto;
  else type = kOneToOne;

  // match the next particle in the chain
  if(type == kAny || type == kVeto) ++_decayStage;

  DecayNode const& dnode(*chain_[_decayStage]);

  for(unsigned iP(0); iP != _list.size(); ++iP){
    PNode& pnode(*_list[iP]);

    if(type == kOneToOne){
      if(oneToOneMatch_(dnode, pnode)){
        std::vector<PNode*> sequence;
        if(_currentSequence) sequence.assign(_currentSequence->begin(), _currentSequence->end());
        sequence.push_back(&pnode);
        if(chainMatch_(_decayStage + 1, pnode.daughters, &sequence)) return true;
      }
    }
    else{
      // collect sequences that lead to descendants that match the dnode (which is one down from the level passed as argument to this function)
      std::vector<std::vector<PNode*> > matchedSequences(anyDecayMatch_(dnode, pnode));
      // check whether subsequent decays of any of the descendants matches the decay chain
      for(unsigned iM(0); iM != matchedSequences.size(); ++iM){
        if(type == kVeto){
          // check for vetoed particle in the matched sequences
          if(isInSequence_(*chain_[_decayStage - 1], matchedSequences[iM])) continue;
        }

        std::vector<PNode*> sequence;
        if(_currentSequence) sequence.assign(_currentSequence->begin(), _currentSequence->end());
        sequence.insert(sequence.end(), matchedSequences[iM].begin(), matchedSequences[iM].end());
        if(chainMatch_(_decayStage + 1, matchedSequences[iM].back()->daughters, &sequence)) return true;
      }
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

std::vector<std::vector<PNode*> >
DecayChainFilter::anyDecayMatch_(DecayNode const& _dnode, PNode& _pnode) const
{
  std::vector<std::vector<PNode*> > matchedSequences;

  // first check for a direct match
  if(oneToOneMatch_(_dnode, _pnode)) matchedSequences.push_back(std::vector<PNode*>(1, &_pnode));
  
  for(unsigned iN(0); iN != _pnode.daughters.size(); ++iN){
    PNode& daughter(*_pnode.daughters[iN]);

    std::vector<std::vector<PNode*> > daughterMatches(anyDecayMatch_(_dnode, daughter));

    for(unsigned iM(0); iM != daughterMatches.size(); ++iM)
      daughterMatches[iM].insert(daughterMatches[iM].begin(), &_pnode);

    matchedSequences.insert(matchedSequences.end(), daughterMatches.begin(), daughterMatches.end());
  }

  return matchedSequences;
}

bool
DecayChainFilter::isInSequence_(DecayNode const& _dnode, std::vector<PNode*> const& _sequence) const
{
  unsigned iC(0);
  for(; iC != _sequence.size(); ++iC)
    if(oneToOneMatch_(_dnode, *_sequence[iC])) break;
  return iC != _sequence.size();
}


NestedGenFilter::NestedGenFilter(GenFilter const* _content, int _repeat) :
  content_(_content),
  repeat_(_repeat)
{
}

NestedGenFilter::~NestedGenFilter()
{
}

TString
NestedGenFilter::toString() const
{
  TString result("");

  if(repeat_ < 0) result += "!";
  else if(repeat_ > 1) result += TString::Format("%d * ", repeat_);

  result += "(";
  result += content_->toString();
  result += ")";

  return result;
}

bool
NestedGenFilter::pass(std::vector<PNode*> const& _rootNodes) const
{
  bool result(content_->pass(_rootNodes));
  if(result && repeat_ > 1){
    int iR(1);
    for(; iR != repeat_; ++iR)
      if(!content_->pass(_rootNodes)) break;
    result = (iR == repeat_);
  }

  return repeat_ > 0 ? result : !result;
}

void
NestedGenFilter::reset() const
{
  content_->reset();
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

void
ChainedGenFilter::reset() const
{
  for(unsigned iT(0); iT != terms_.size(); ++iT)
    terms_[iT]->reset();
}

#endif
