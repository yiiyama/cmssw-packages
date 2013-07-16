{
  gROOT->LoadMacro("GenTreeViewerRA3.cc+");

  susy::Event* event = new susy::Event;

  TChain* susyTree = new TChain("susyTree");
  //  susyTree->Add("/store/RA3Ntuples/SusyNtuples/cms533v1/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/PU_S10_START53_V7A-v2/susyEvents_*.root");
  //  susyTree->Add("/store/RA3Ntuples/SusyNtuples/cms533v1/Summer12_DR53X-PrivateSkim/QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/LoosePhoton15/susyEvents_*.root");
  //  susyTree->Add("/store/RA3Ntuples/SusyNtuples/cms538pre1/Summer12_DR53X/GVJets_Incl_8TeV-madgraph/PU_S10_START53_V7C-v1/susyEvents_1_1_YUf.root");
  //  susyTree->Add("/store/RA3Ntuples/SusyNtuples/cms538v0p1/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/PU_S10_START53_V7A-v2/susyEvents_1_1_*.root");
  //  susyTree->Add("root://eoscms//eos/cms/store/user/yiiyama/RA3Ntuples/TChiwg/TChiwg_800_700/TChiwg_800_700_21_1_2HH.root");
  //  susyTree->Add("/data/disk0/RA3Ntuples/SusyNtuples/cms538v0p1/Summer12_DR53X/DYToEE_M_20_TuneZ2star_8TeV_pythia6_v2/PU_S10_START53_V7A-v1/susyEvents_*.root");
  susyTree->Add("/store/RA3Ntuples/SusyNtuples/cms538v0p1/Summer12_DR53X/DYToEE_M_20_TuneZ2star_8TeV_pythia6_v2/PU_S10_START53_V7A-v1//susyEvents_560_1_c19.root");
  susyTree->Add("/store/RA3Ntuples/SusyNtuples/cms538v0p1/Summer12_DR53X/DYToEE_M_20_TuneZ2star_8TeV_pythia6_v2/PU_S10_START53_V7A-v1//susyEvents_633_1_XyZ.root");
  susyTree->Add("/store/RA3Ntuples/SusyNtuples/cms538v0p1/Summer12_DR53X/DYToEE_M_20_TuneZ2star_8TeV_pythia6_v2/PU_S10_START53_V7A-v1//susyEvents_767_1_B1O.root");
  susyTree->Add("/store/RA3Ntuples/SusyNtuples/cms538v0p1/Summer12_DR53X/DYToEE_M_20_TuneZ2star_8TeV_pythia6_v2/PU_S10_START53_V7A-v1//susyEvents_163_1_zXZ.root");

  susyTree->SetBranchStatus("*", 0);
  susyTree->SetBranchStatus("genParticles*", 1);
  susyTree->SetBranchStatus("eventNumber", 1);
  //  susyTree->SetBranchStatus("luminosityBlockNumber", 1);

  event->setInput(*susyTree);

  int iFound = 0;

  int max = 100;

  int iTree = -1;
  long iEvent = 0;
  long iDisplay = 0;
  while(event->getEntry(iEvent++) != 0){
    if(susyTree->GetTreeNumber() != iTree){
      iTree = susyTree->GetTreeNumber();
      cout << "Tree " << iTree << ": " << susyTree->GetCurrentFile()->GetName() << endl;
    }

    switch(event->eventNumber){
    case 16287557:
    case 16715561:
    case 5069132:
      cout << susyTree->GetCurrentFile()->GetName() << endl;
      cout << event->eventNumber << endl;
      viewGenTreeRA3(*event, true);
      ++iDisplay;
      break;
    default:
      break;
    }

    if(iDisplay == 3) break;

    // Look for 5 tau decays
//     unsigned nP = event->genParticles.size();
//     bool hasE(false);
//     bool hasMu(false);
//     bool hasTau(false);
//     for(unsigned iP = 0; iP < nP; ++iP){
//       if(TMath::Abs(event->genParticles.at(iP).pdgId) == 11) hasE = true;
//       if(TMath::Abs(event->genParticles.at(iP).pdgId) == 13) hasMu = true;
//       if(TMath::Abs(event->genParticles.at(iP).pdgId) == 15) hasTau = true;
//     }

//     if((hasE || hasMu) && hasTau){
//       viewGenTreeRA3(*event);
//       if(++iFound == 5) break;
//     }

// select final states with photon
//     vector<susy::Particle>* particles(&event->genParticles);
//     unsigned iP = 0;
//     for(; iP != particles->size(); ++iP){
//       susy::Particle* particle = &particles->at(iP);
//       if(particle->momentum.Pt() < 10.) continue;
//       if(particle->pdgId == 22 && particle->status == 1) break;
//     }
//     if(iP == particles->size()) continue;

// select final states with electron
//     vector<susy::Particle>* particles(&event->genParticles);
//     unsigned iP = 0;
//     for(; iP != particles->size(); ++iP){
//       susy::Particle* particle = &particles->at(iP);
//       if(particle->momentum.Pt() < 10.) continue;
//       if(TMath::Abs(particle->pdgId) == 11 && particle->status == 1) break;
//     }
//     if(iP == particles->size()) continue;

    // Just dump events
//     viewGenTreeRA3(*event, true);
//     ++iDisplay;

//     if(iDisplay == max) break;
  }

}

