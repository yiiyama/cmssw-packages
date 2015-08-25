import os
import sys
import ROOT
ROOT.gROOT.SetBatch(True)

inputPath = sys.argv[1]
try:
    nevents = int(sys.argv[2])
except:
    nevents = 1

ROOT.gSystem.Load('libMitAnaDataTree.so')
ROOT.gROOT.LoadMacro(os.environ['CMSSW_BASE'] + '/src/Toolset/GenTreeViewer/test/GenTreeViewerBambu.cc+')

source = ROOT.TFile.Open(inputPath)

tree = source.Get('Events')

branch = tree.GetBranch('MCParticles')

particles = ROOT.mithep.MCParticleArr()
tree.SetBranchAddress('MCParticles', particles)

for iE in range(nevents):
    branch.GetEntry(iE)
    ROOT.GenTreeViewerBambu(particles)
