{
  gSystem->Load("libSusyEvent.so");
  gSystem->AddIncludePath("-I" + TString(gSystem->Getenv("CMSSW_BASE")) + "/src");

  gROOT->LoadMacro("GenTreeViewerRA3.cc+");

  susy::Event* event = new susy::Event;

  //  TFile* source = TFile::Open("/store/RA3Ntuples/PrivateMC/FastSim/525p1/Spectra_gsq_W/SusyNtuple/cms533v1_v1/tree_1000_1020_375.root");
  //  TFile* source = TFile::Open("/store/RA3Ntuples/PrivateMC/FastSim/533p3_full/naturalHiggsinoNLSP/SusyNtuple/cms533v1_v1/tree_naturalHiggsinoNLSPout_mst_100_M3_5025_mu_175.root");
  //  TFile* source = TFile::Open("http://dcmu00/store/RA3Ntuples/SusyNtuples/cms533v1/Summer12_DR53X/WGToLNuG_TuneZ2star_8TeV-madgraph-tauola/PU_S10_START53_V7A-v1/susyEvents_1_1_aUD.root");
  //  TFile* source = TFile::Open("/store/RA3Ntuples/SusyNtuples/cms533v1/Summer12_DR53X/ZGToLLG_8TeV-madgraph/PU_S10_START53_V7A-v1/susyEvents_1_1_bwB.root");
  //  TFile* source = TFile::Open("/store/RA3Ntuples/SusyNtuples/cms533v1/Summer12_DR53X-PrivateSkim/ZGToLLG_NoFSR/PU_S10_START53_V7A-v1/susyEvents_disk0_0_2.root");
  //  TFile* source = TFile::Open("/store/RA3Ntuples/SusyNtuples/cms533v1/Summer12_DR53X/ZGToLLG_8TeV-madgraph/PU_S10_START53_V7A-v1/susyEvents_100_1_thO.root");
  //  TFile* source = TFile::Open("/data/yiiyama/signalSkim/skim_disk0_0_20_WGToLNuG_Signal.root");
  //  TFile* source = TFile::Open("/store/RA3Ntuples/SusyNtuples/cms533v1/Summer12-PrivateSkim/SMS-T5wg_Mgluino-400to2000_Mchargino-100to2000_8TeV-Pythia6Z/Sorted/susyEvents_1000_1000_disk0_218_220.root");
  //  TFile* source = TFile::Open("/store/RA3Ntuples/SusyNtuples/cms533v1/Summer12_DR53X/QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/PU_S10_START53_V7A-v1/susyEvents_1_1_Yh1.root");
  //  TFile* source = TFile::Open("http://dcmu00/store/RA3Ntuples/SusyNtuples/cms533v1/Summer12_DR53X/DYToEE_M_20_TuneZ2star_8TeV_pythia6_v2/PU_S10_START53_V7A-v1/susyEvents_106_1_2os.root");
  //  TFile* source = TFile::Open("/store/RA3Ntuples/SusyNtuples/cms533v1/Summer12_DR53X-PrivateSkim/GVJets_Incl_8TeV-madgraph/GenLepton/susyEvents_1_1_ADg.root");
  //  TFile* source = TFile::Open("/store/RA3Ntuples/SusyNtuples/cms533v1/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/PU_S10_START53_V7A-v2/susyEvents_1_1_dMf.root");
  //  TFile* source = TFile::Open("/store/RA3Ntuples/SusyNtuples/cms533v1/Summer12_DR53X-PrivateSkim/TT_CT10_TuneZ2star_8TeV-powheg-tauola/GenLepton/susyEvents_disk0_0_4.root");

  //  TTree* susyTree = (TTree*)source->Get("susyTree");
  TChain* susyTree = new TChain("susyTree");
  susyTree->Add("/store/RA3Ntuples/SusyNtuples/cms533v1/Summer12_DR53X/DYToEE_M_20_TuneZ2star_8TeV_pythia6_v2/PU_S10_START53_V7A-v1/susyEvents*.root");

  susyTree->SetBranchStatus("*", 0);
  susyTree->SetBranchStatus("genParticles*", 1);
  susyTree->SetBranchStatus("eventNumber", 1);
  susyTree->SetBranchStatus("luminosityBlockNumber", 1);

  susyTree->SetBranchAddress("susyEvent", &event);

  int iFound = 0;

  int iTree = -1;
  long iEvent = 0;
  long iDisplay = 0;
  while(susyTree->GetEntry(iEvent++) != 0){
    if(susyTree->GetTreeNumber() != iTree){
      iTree = susyTree->GetTreeNumber();
      cout << "Tree " << iTree << ": " << susyTree->GetCurrentFile()->GetName() << endl;
    }

    switch(event->eventNumber){
//     case  596563:
//     case  596613:
//     case  596683:
//     case  729305:
//     case  729549:
//     case  730327:
//     case 5689968:
//     case 5690495:
//     case 5700493:
//     case 5700493:
//     case 5711359:
//     case 5711429:
//     case 5711957:
//     case 5726619:
//     case 5737564:
//     case 5738026:
//     case 5755061:
//     case 5755061:
//     case 5755460:
//     case 5760589:
//     case 5761047:
//     case 5761218:
//     case 5816172:
//     case  597020:
//     case  597218:
    case   730327:
    case  5711957:
    case  5726619:
    case  5816172:
    case   662130:
    case   662339:
    case  3411907:
    case  5693073:
    case  5695295:
    case  5695518:
    case  5695608:
    case  5705941:
    case  6204390:
    case   620373:
    case   621110:
    case  7869374:
    case  7990112:
    case  7990665:
    case 12516056:
    case   636003:
    case   798604:
    case  8705976:
    case   838709:
    case  1552102:
    case  1571968:
      cout << susyTree->GetCurrentFile()->GetName() << endl;
      cout << event->eventNumber << endl;
      viewGenTreeRA3(*event, 2., true);
      ++iDisplay;
      break;
    default:
      break;
    }

    if(iDisplay == 10) break;

    // search for hard FSR in DYToEE
//     switch(event->eventNumber){
//     case  5692103:
//     case 19796384:
//     case  4369043:
//     case  8887813:
//     case 11801987:
//     case  5857306:
//     case  2976895:
//     case  2985458:
//     case  3797650:
//     case  3798158:
//     case  5164093:
//     case  1261097:
//     case  1261877:
//     case  4456123:
//     case  1033091:
//     case  7585757:
//     case  1281969:
//     case  3116714:
//     case  3349775:
//     case  3854058:
//     case  2482344:
//     case 12389630:
//     case  1287527:
//     case  2897859:
//     case 12548191:
//     case  12501397:
//     case  18301718:
//     case  19391341:
//     case  14639217:
//     case   1229683:
//     case   1230331:
//     case   4492423:
//     case   4600551:
//     case  10643065:
//     case   6119542:
//     case  14149816:
//     case  14925889:
//     case     59022:
//     case  19908073:
//     case   4152700:
//     case     97300:
//     case    299217:
//     case   2833696:
//     case   3358941:
//     case   6371866:
//     case   6734755:
//     case   1459949:
//     case  16173358:
//     case  17796065:
//     case   4999468:
//     case   1513567:
//     case   1642374:
//     case   1902129:
//     case   2944142:
//     case   4491937:
//     case   1648999:
//     case   1905067:
//     case  10928516:
//     case   5917318:
//     case   5917885:
//     case   6138396:
//     case   9114712:
//     case   2209529:
//     case  12387261:
//     case   5573693:
//     case   1934715:
//     case   8140168:
//     case  12523232:
//     case   7503045:
//     case   1180334:
//     case   7779086:
//     case   7908328:
//     case   4313902:
//     case   7894155:
//     case  17931165:
//     case   2766544:
//     case   4547447:
//     case   5924156:
//       cout << susyTree->GetCurrentFile()->GetName() << endl;
//       cout << event->luminosityBlockNumber << endl;
//       viewGenTreeRA3(*event, 2., true);
//       ++iDisplay;
//       break;
//     default:
//       break;
//     }

    // search for hard FSR in DYToMuMu
//     switch(event->eventNumber){
//     case   845304:
//     case  7799215:
//     case   840376:
//     case  5091754:
//     case   869564:
//     case   872190:
//     case   873916:
//     case   909189:
//     case   919307:
//     case   948410:
//     case   970843:
//     case  1496493:
//     case  1057347:
//     case  2729925:
//     case  2751390:
//     case  1463193:
//     case  2382543:
//     case  2422854:
//     case 11293415:
//     case  5121340:
//     case 11650513:
//     case  2735500:
//     case  2737926:
//     case  5050390:
//     case  5096268:
//     case  1270072:
//     case  1300295:
//     case  6206333:
//     case  6206619:
//     case  6249250:
//     case   929323:
//     case   931649:
//     case   992629:
//     case  1414198:
//     case   931398:
//     case  1206922:
//     case  2347455:
//     case   934650:
//     case  1445193:
//     case  1445451:
//     case  6975304:
//     case 16816397:
//     case   954067:
//     case  1436918:
//     case  1511010:
//     case  2442808:
//     case    91065:
//     case 10862601:
//     case 10981439:
//     case 14801100:
//     case  8233706:
//     case  1076252:
//     case  3550109:
//     case  1180261:
//     case  1197418:
//     case  1133179:
//     case  3670227:
//     case  7548776:
//     case  1085120:
//     case  1796822:
//     case  1805235:
//     case  1109505:
//     case  1134629:
//     case  3542455:
//     case  3959555:
//     case  9597897:
//     case  1138280:
//     case  1108350:
//     case  1108786:
//     case  1679526:
//     case  3636179:
//     case  1110532:
//     case  8224192:
//     case 13246839:
//     case  1668657:
//     case  3564500:
//     case  3619117:
//     case  5000350:
//     case  3536680:
//     case  6153291:
//       cout << susyTree->GetCurrentFile()->GetName() << endl;
//       cout << event->luminosityBlockNumber << endl;
//       viewGenTreeRA3(*event, 2., true);
//       ++iDisplay;
//       break;
//     default:
//       break;
//     }

//     if(iDisplay == 3) break;

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
//       viewGenTreeRA3(*event, 2.);
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
//      viewGenTreeRA3(*event, 2., true);
//      if(++iDisplay == 10) break;
  }

}

