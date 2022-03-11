#include "TFile.h"
#include "TSystem.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "../interface/SimpleEventDisplay.h"

#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <sstream>

int
main(int argc, char** argv)
{
  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();

  TFile* source(TFile::Open(argv[1]));
  if(!source){
    std::cerr << "Usage: gendisp source.root [lumi1,]event1 [[lumi2,]event2 ..]" << std::endl;
    return 1;
  }

  std::vector<std::pair<unsigned, unsigned> > eventNumbers;
  for(int iArg(2); iArg != argc; ++iArg){
    std::string arg(argv[iArg]);
    std::pair<unsigned, unsigned> eventNumber(0, 0);
    size_t colon(arg.find(":"));
    if(colon != std::string::npos){
      eventNumber.first = std::atoi(arg.substr(0, colon).c_str());
      eventNumber.second = std::atoi(arg.substr(colon + 1).c_str());
      eventNumbers.push_back(eventNumber);
    }
    else{
      eventNumber.second = std::atoi(arg.c_str());
      eventNumbers.push_back(eventNumber);
    }
  }

  SimpleEventDisplay display;
  display.setPtThreshold(0.);

  fwlite::Event event(source);
  fwlite::Handle<reco::GenParticleCollection> genParticles;

  for(event.toBegin(); !event.atEnd() && eventNumbers.size() != 0; ++event){
    unsigned eventNumber(event.id().event());
    std::vector<std::pair<unsigned, unsigned> >::iterator eItr(eventNumbers.begin());
    for(; eItr != eventNumbers.end(); ++eItr){
      if(eItr->second == eventNumber){
        if(eItr->first == 0) break;
        else if(eItr->first == event.id().luminosityBlock()) break;
      }
    }
    if(eItr == eventNumbers.end()) continue;

    eventNumbers.erase(eItr);

    genParticles.getByLabel(event, "genParticles");

    display.showEvent(event.id().event(), *genParticles.product());

    std::stringstream outputName;
    outputName << event.id().event() << ".pdf";

    display.print(outputName.str().c_str());
  }

  delete source;

  return 0;
}
