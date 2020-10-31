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
    print("found dataset " + str(dataset))
    chain.Add(dataset)
print 'Total number of events: ' + str(chain.GetEntries())

#variables
lep1_pt = array('f',[0.])
lep1_eta = array('f',[0.])
lep1_phi = array('f',[0.])
lep1_mass = array('f',[0.])
lep2_pt = array('f',[0.])
lep2_eta = array('f',[0.])
lep2_phi = array('f',[0.])
lep2_mass = array('f',[0.])
lep3_pt = array('f',[0.])
lep3_eta = array('f',[0.])
lep3_phi = array('f',[0.])
lep3_mass = array('f',[0.])
lep4_pt = array('f',[0.])
lep4_eta = array('f',[0.])
lep4_phi = array('f',[0.])
lep4_mass = array('f',[0.])
lep1FSR_pt = array('f',[0.])
lep1FSR_eta = array('f',[0.])
lep1FSR_phi = array('f',[0.])
lep1FSR_mass = array('f',[0.])
lep2FSR_pt = array('f',[0.])
lep2FSR_eta = array('f',[0.])
lep2FSR_phi = array('f',[0.])
lep2FSR_mass = array('f',[0.])
lep3FSR_pt = array('f',[0.])
lep3FSR_eta = array('f',[0.])
lep3FSR_phi = array('f',[0.])
lep3FSR_mass = array('f',[0.])
lep4FSR_pt = array('f',[0.])
lep4FSR_eta = array('f',[0.])
lep4FSR_phi = array('f',[0.])
lep4FSR_mass = array('f',[0.])
h_pt = array('f',[0.])
h_eta = array('f',[0.])
h_phi = array('f',[0.])
h_mass = array('f',[0.])

#Output file and any Branch we want
file_out = ROOT.TFile(args.outputfile, 'recreate')
passedEvents = ROOT.TTree("passedEvents","passedEvents")
passedEvents.Branch("lep1_pt",lep1_pt,"lep1_pt/F")
passedEvents.Branch("lep1_eta",lep1_eta,"lep1_eta/F")
passedEvents.Branch("lep1_phi",lep1_phi,"lep1_phi/F")
passedEvents.Branch("lep1_mass",lep1_mass,"lep1_mass/F")
passedEvents.Branch("lep2_pt",lep2_pt,"lep2_pt/F")
passedEvents.Branch("lep2_eta",lep2_eta,"lep2_eta/F")
passedEvents.Branch("lep2_phi",lep2_phi,"lep2_phi/F")
passedEvents.Branch("lep2_mass",lep2_mass,"lep2_mass/F")
passedEvents.Branch("lep3_pt",lep3_pt,"lep3_pt/F")
passedEvents.Branch("lep3_eta",lep3_eta,"lep3_eta/F")
passedEvents.Branch("lep3_phi",lep3_phi,"lep3_phi/F")
passedEvents.Branch("lep3_mass",lep3_mass,"lep3_mass/F")
passedEvents.Branch("lep4_pt",lep4_pt,"lep4_pt/F")
passedEvents.Branch("lep4_eta",lep4_eta,"lep4_eta/F")
passedEvents.Branch("lep4_phi",lep4_phi,"lep4_phi/F")
passedEvents.Branch("lep4_mass",lep4_mass,"lep4_mass/F")
passedEvents.Branch("lep1FSR_pt",lep1FSR_pt,"lep1FSR_pt/F")
passedEvents.Branch("lep1FSR_eta",lep1FSR_eta,"lep1FSR_eta/F")
passedEvents.Branch("lep1FSR_phi",lep1FSR_phi,"lep1FSR_phi/F")
passedEvents.Branch("lep1FSR_mass",lep1FSR_mass,"lep1FSR_mass/F")
passedEvents.Branch("lep2FSR_pt",lep2FSR_pt,"lep2FSR_pt/F")
passedEvents.Branch("lep2FSR_eta",lep2FSR_eta,"lep2FSR_eta/F")
passedEvents.Branch("lep2FSR_phi",lep2FSR_phi,"lep2_phi/F")
passedEvents.Branch("lep2FSR_mass",lep2FSR_mass,"lep2FSR_mass/F")
passedEvents.Branch("lep3FSR_pt",lep3FSR_pt,"lep3FSR_pt/F")
passedEvents.Branch("lep3FSR_eta",lep3FSR_eta,"lep3FSR_eta/F")
passedEvents.Branch("lep3FSR_phi",lep3FSR_phi,"lep3FSR_phi/F")
passedEvents.Branch("lep3FSR_mass",lep3FSR_mass,"lep3FSR_mass/F")
passedEvents.Branch("lep4FSR_pt",lep4FSR_pt,"lep4FSR_pt/F")
passedEvents.Branch("lep4FSR_eta",lep4FSR_eta,"lep4FSR_eta/F")
passedEvents.Branch("lep4FSR_phi",lep4FSR_phi,"lep4FSR_phi/F")
passedEvents.Branch("lep4FSR_mass",lep4FSR_mass,"lep4FSR_mass/F")
passedEvents.Branch("h_pt",h_pt,"h_pt/F")
passedEvents.Branch("h_eta",h_eta,"h_eta/F")
passedEvents.Branch("h_phi",h_phi,"h_phi/F")
passedEvents.Branch("h_mass",h_mass,"h_mass/F")

#Loop over all the events in the input ntuple
for ievent,event in enumerate(chain):
    if(not event.passedTrig): continue
    if(not event.passedFullSelection): continue
    Nlep = event.lep_pt.size()
    for i in range(Nlep):

    #fill tree

       lep1_pt[0] = event.lep_pt[event.lep_Hindex[0]]
       lep1_eta[0] = event.lep_eta[event.lep_Hindex[0]]
       lep1_phi[0] = event.lep_phi[event.lep_Hindex[0]]
       lep1_mass[0] = event.lep_mass[event.lep_Hindex[0]]


       lep2_pt[0] = event.lep_pt[event.lep_Hindex[1]]
       lep2_eta[0] = event.lep_eta[event.lep_Hindex[1]]
       lep2_phi[0] = event.lep_phi[event.lep_Hindex[1]]
       lep2_mass[0] = event.lep_mass[event.lep_Hindex[1]]

       lep3_pt[0] = event.lep_pt[event.lep_Hindex[2]]
       lep3_eta[0] = event.lep_eta[event.lep_Hindex[2]]
       lep3_phi[0] = event.lep_phi[event.lep_Hindex[2]]
       lep3_mass[0] = event.lep_mass[event.lep_Hindex[2]]

       lep4_pt[0] = event.lep_pt[event.lep_Hindex[3]]
       lep4_eta[0] = event.lep_eta[event.lep_Hindex[3]]
       lep4_phi[0] = event.lep_phi[event.lep_Hindex[3]]
       lep4_mass[0] = event.lep_mass[event.lep_Hindex[3]]

       lep1FSR_pt[0] = event.lepFSR_pt[event.lep_Hindex[0]]
       lep1FSR_eta[0] = event.lepFSR_eta[event.lep_Hindex[0]]
       lep1FSR_phi[0] = event.lepFSR_phi[event.lep_Hindex[0]]
       lep1FSR_mass[0] = event.lepFSR_mass[event.lep_Hindex[0]]

       lep2FSR_pt[0] = event.lepFSR_pt[event.lep_Hindex[1]]
       lep2FSR_eta[0] = event.lepFSR_eta[event.lep_Hindex[1]]
       lep2FSR_phi[0] = event.lepFSR_phi[event.lep_Hindex[1]]
       lep2FSR_mass[0] = event.lepFSR_mass[event.lep_Hindex[1]]

       lep3FSR_pt[0] = event.lepFSR_pt[event.lep_Hindex[2]]
       lep3FSR_eta[0] = event.lepFSR_eta[event.lep_Hindex[2]]
       lep3FSR_phi[0] = event.lepFSR_phi[event.lep_Hindex[2]]
       lep3FSR_mass[0] = event.lepFSR_mass[event.lep_Hindex[2]]

       lep4FSR_pt[0] = event.lepFSR_pt[event.lep_Hindex[3]]
       lep4FSR_eta[0] = event.lepFSR_eta[event.lep_Hindex[3]]
       lep4FSR_phi[0] = event.lepFSR_phi[event.lep_Hindex[3]]
       lep4FSR_mass[0] = event.lepFSR_mass[event.lep_Hindex[3]]

    Hhiggs = event.H_pt.size()
    for i in range(Hhiggs):
        h_pt[0] = event.H_pt[i]
        h_eta[0] =event.H_eta[i]
        h_phi[0] = event.H_phi[i]
        h_mass[0] = event.H_mass[i]
#      lep_Hindex[i] = event.lep_Hindex[i]
    passedEvents.Fill()

file_out.Write()
file_out.Close()
