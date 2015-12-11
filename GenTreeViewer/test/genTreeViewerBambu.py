#!/usr/bin/env python

import os
import sys
from argparse import ArgumentParser
import ROOT
ROOT.gROOT.SetBatch(True)

argParser = ArgumentParser(description = 'Dump gen-particle tree.')
argParser.add_argument('sourcePath', metavar = 'FILE', nargs = '+', help = 'Path to the input Bambu file.')
argParser.add_argument('--event', '-e', metavar = 'SPEC', dest = 'events', nargs = '+', help = 'Run:Event')
argParser.add_argument('--nentries', '-n', metavar = 'N', dest = 'nentries', type = int, default = 1, help = 'Number of entries to process.')

argParser.add_argument('--min-pt', '-t', metavar = 'MIN', dest = 'minPt', type = float, default = 0., help = 'pT threshold for final-state particles.')
argParser.add_argument('--show-momenta', '-p', metavar = 'MODE', dest = 'showP', default = 'Final', help = 'Print momenta of the particles. Options are All, Final, and None.')
argParser.add_argument('--show-mass', '-m', metavar = 'MODE', dest = 'showM', default = 'HardScat', help = 'Print mass of the particles. Options are All, HardScat, and None.')
argParser.add_argument('--no-cleaning', '-C', action = 'store_true', dest = 'noCleaning', help = 'Do not apply cleaning.')

args = argParser.parse_args()
sys.argv = []

ROOT.gSystem.Load('libMitAnaDataTree.so')
ROOT.gROOT.LoadMacro(os.environ['CMSSW_BASE'] + '/src/Toolset/GenTreeViewer/test/GenTreeViewerBambu.cc+')

if args.showP == 'All':
    showP = ROOT.PNode.kShowAllP
elif args.showP == 'Final':
    showP = ROOT.PNode.kShowFinalP
elif args.showP == 'None':
    showP = ROOT.PNode.kNoP
else:
    print 'Undefined show-momenta option:', args.showP
    sys.exit(1)

if args.showM == 'All':
    showM = ROOT.PNode.kShowAllM
elif args.showM == 'HardScat':
    showM = ROOT.PNode.kShowHardScatM
elif args.showM == 'None':
    showM = ROOT.PNode.kNoM
else:
    print 'Undefined show-mass option:', args.showM
    sys.exit(1)


tree = ROOT.TChain('Events')
for path in args.sourcePath:
    tree.Add(path)

if args.events is not None:
    eventsToProcess = map(lambda x: tuple(map(int, x.split(':'))), args.events)

    hdBranch = tree.GetBranch('EventHeader')
    header = ROOT.mithep.EventHeader()
    tree.SetBranchAddress('EventHeader', header)

mcBranch = tree.GetBranch('MCParticles')
particles = ROOT.mithep.MCParticleArr()
tree.SetBranchAddress('MCParticles', particles)

iEntry = -1
iDisplayed = 0
while True:
    iEntry += 1

    if args.events is not None:
        if hdBranch.GetEntry(iEntry) <= 0:
            print 'End of input'
            break

        if (header.RunNum(), header.EvtNum()) not in eventsToProcess:
            continue

    mcBranch.GetEntry(iEntry)
    ROOT.GenTreeViewerBambu(particles, args.minPt, showP, showM, not args.noCleaning)

    iDisplayed += 1

    if iDisplayed == args.nentries:
        break
