#!/usr/bin/env python
import ROOT
import os
import numpy as np
from array import array

import argparse
parser = argparse.ArgumentParser(description = "A simple ttree plotter")
parser.add_argument("-t", "--ttree", dest="ttree", default="Ana/passedEvents", help="TTree Name")
parser.add_argument("-i", "--inputfiles", dest="inputfiles", default="Sync_1031_2018_ttH_v2.root", help="List of input root files")
parser.add_argument("-o", "--outputfile", dest="outputfile", default="plots.root", help="Output file containing plots")
parser.add_argument("-s","--substring", dest="substring",default="",help='only submit datasets with this string in the name')
args = parser.parse_args()

print("get file " + str(args.inputfiles))
#input Ntuple
chain = ROOT.TChain(args.ttree)
f = open(args.inputfiles)
for line in f:
    dataset = line.strip('\n')
    print("dataset after strip " + str(dataset))
    chain.Add(dataset)
print 'Total number of events: ' + str(chain.GetEntries())

#variables
lep1_pt = array('f',[0.])
'''
lep_Hindex = array('l'[0])
lep_id = array('l',[0])
lep_pt = array('f',[0.])
lep_phi = array('f',[0.])
lep_eta = array('f',[0.])
lep_mass = array('f',[0.])
lepFSR_pt = array('f',[0.])
lepFSR_phi = array('f',[0.])
lepFSR_eta = array('f',[0.])
lepFSR_mass = array('f',[0.])
Z_pt = array('f',[0.])
Z_eta = array('f',[0.0])
Z_phi = array('f',[0.])
Z_mass = array('f',[0.])
Z_noFSR_pt = array('f',[0.])
Z_noFSR_eta = array('f',[0.0])
Z_noFSR_phi = array('f',[0.])
Z_noFSR_mass = array('f',[0.])
H_pt = array('f',[0.])
H_eta = array('f',[0.])
H_phi = array('f',[0.])
H_mass = array('f',[0.])
passedFullSelection = True
passedTrig = True
passedZ1LSelection = True
passedZ4lSelection = True
'''

#Output file and any Branch we want
file_out = ROOT.TFile(args.outputfile, 'recreate')
passedEvents = ROOT.TTree("passedEvents","passedEvents")
passedEvents.Branch("lep_pt1",lep1_pt,"lep1_pt/F")
'''
passedEvents.Branch("lep_Hindex",lep_Hindex,"lep_Hindex/L")
passedEvents.Branch("lep_id",lep_id,"lep_id/L")
passedEvents.Branch("lep_pt",lep_pt,"lep_pt/F")
passedEvents.Branch("lep_eta",lep_eta,"lep_eta/F")
passedEvents.Branch("lep_phi",lep_phi,"lep_phi/F")
passedEvents.Branch("lep_mass",lep_mass,"lep_mass/F")
passedEvents.Branch("lepFSR_pt",lepFSR_pt,"lepFSR_pt/F")
passedEvents.Branch("lepFSR_eta",lepFSR_eta,"lepFSR_eta/F")
passedEvents.Branch("lepFSR_phi",lepFSR_phi,"lepFSR_phi/F")
passedEvents.Branch("lepFSR_mass",lepFSR_mass,"lepFSR_mass/F")
passedEvents.Branch("Z_pt",Z_pt,"Z_pt/F")
passedEvents.Branch("Z_eta",Z_eta,"Z_eta/F")
passedEvents.Branch("Z_phi",Z_phi,"Z_phi/F")
passedEvents.Branch("Z_mass",Z_mass,"Z_mass/F")
passedEvents.Branch("Z_noFSR_pt",Z_noFSR_pt,"Z_noFSR_pt/F")
passedEvents.Branch("Z_noFSR_eta",Z_noFSR_eta,"Z_noFSR_eta/F")
passedEvents.Branch("Z_noFSR_phi",Z_noFSR_phi,"Z_noFSR_phi")
passedEvents.Branch("Z_noFSR_mass",Z_noFSR_mass,"Z_noFSR_mass")
passedEvents.Branch("H_pt",H_pt,"H_pt/F")
passedEvents.Branch("H_eta",H_eta,"H_eta/F")
passedEvents.Branch("H_phi",H_phi,"H_phi/F")
passedEvents.Branch("H_mass",H_mass,"H_mass")
passedEvents.Branch("passedTrig",passedTrig,"passedTrig/O")
passedEvents.Branch("passedZ1LSelection",passedZ1LSelection,"passedZ1LSelection/O")
passedEvents.Branch("passedZ4lSelection",passedZ4lSelection,"passedZ4lSelection/O")
passedEvents.Branch("passedFullSelection",passedFullSelection,"passedFullSelection/O")
'''

#Loop over all the events in the input ntuple
for ievent,event in enumerate(chain):
    Nlep = event.lep_pt.size()
    for i in range(Nlep):
    #fill tree
       lep1_pt[0] = event.lep_pt[event.lep_Hindex[0]]
#      lep_Hindex[i] = event.lep_Hindex[i]
       passedEvents.Fill()

file_out.Write()
file_out.Close()
