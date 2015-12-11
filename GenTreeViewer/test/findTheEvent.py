import os
import sys
import ROOT
ROOT.gROOT.SetBatch(True)

ROOT.gSystem.Load('libMitAnaDataTree.so')
ROOT.gROOT.LoadMacro(os.environ['CMSSW_BASE'] + '/src/Toolset/GenTreeViewer/test/GenTreeViewerBambu.cc+')

files = []
with open('/home/cmsprod/catalog/t2mit/filefi/042/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM/Files') as fileList:
    for line in fileList:
        words = line.split()
        if words[0] == '0000':
            files.append(words[1])

sourceDir = '/mnt/hadoop/cms/store/user/paus/filefi/042/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM/'

for path in files:
    source = ROOT.TFile.Open(sourceDir + path)
    tree = source.Get('Events')

    hdrBranch = tree.GetBranch('EventHeader')
    header = ROOT.mithep.EventHeader()
    tree.SetBranchAddress('EventHeader', header)

    branch = tree.GetBranch('MCParticles')
    particles = ROOT.mithep.MCParticleArr()
    tree.SetBranchAddress('MCParticles', particles)
    
    iEntry = 0
    while hdrBranch.GetEntry(iEntry) > 0:
        if header.EvtNum() == 32459237:
            branch.GetEntry(iEntry)
            ROOT.GenTreeViewerBambu(particles)
            break

        iEntry += 1

    else:
        continue

    break
