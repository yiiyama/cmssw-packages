void
exec(TString const& path, long max = 20)
{
  gROOT->LoadMacro("/afs/cern.ch/user/y/yiiyama/cmssw/Ntuplizer538/src/Toolset/GenTreeViewer/test/GenTreeViewerRA3.cc+");
  gROOT->LoadMacro("/afs/cern.ch/user/y/yiiyama/cmssw/Ntuplizer538/src/Toolset/GenTreeViewer/test/GenDecayFilterRA3.cc+");

  TChain input("susyTree");
  input.Add(path);

  input.SetBranchStatus("*", 0);
  input.SetBranchStatus("genParticles*", 1);

  susy::Event* event = new susy::Event;
  event->setInput(input);

  GenDecayFilterRA3 filter("15>24>11");

  int iTree = -1;
  long iEntry = 0;
  long iDisplay = 0;
  while(event->getEntry(iEntry++) != 0){
    if(input.GetTreeNumber() != iTree){
      iTree = input.GetTreeNumber();
      cout << "Tree " << iTree << ": " << input.GetCurrentFile()->GetName() << endl;
    }

    if(!filter.pass(*event)) continue;

    viewGenTreeRA3(*event);
    ++iDisplay;

    if(iDisplay == max) break;
  }

  delete event;
}

