#!/usr/bin/env python
import ROOT
from array import array

# input file and get tree
ggTozz = ROOT.TChain(passedEvents)
ggTozz.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC/ GluGluToContinToZZ*.root")

qqTozz = ROOT.TChain(passedEvents)
qqTozz.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC/ZZTo4L*.root")

DataSim = ROOT.TChain(passedEvents)
DataSim.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC/GluGluHToZZ*.root")
DataSim.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC/VBF_HToZZTo4L*.root")

Signal = ROOT.TFile('/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/2016_allsignal.root')
t = Signal.Get('passedEvents')

# book histogram and canvas
c = ROOT.TCanvas()

gg = ROOT.TH1D("gg->zz","gg->zz",50,70,170)
gg.SetFillColor(KBlue)

qq = ROOT.TH1D("qq->zz","qq->zz",50,70,170)
qq.SetFillColor(7)

Sim = ROOT.TH1D("sim","sim",50,70,170)
Sim.SetFillColor(KRed)

Data = ROOT.TH1D("Data","Data",50,70,170)
Data.GetYaxis().SetTitle("Events / 2 GeV")
Data.SetMarkerStyle(20);
Data.SetMarkerColor(kBlack);
Data.SetMarkerSize(1.2);
Data.SetLineColor(kBlack);
Data.SetLineWidth(1);
Data.SetStats(kFALSE);

leg = ROOT.TLegend(0.7, 0.7, 0.85, 0.85)
leg.AddEntry(Data,"Data","PE1")
leg.AddEntry(gg,"gg->zz","f")
leg.AddEntry(qq,"qq->zz","f")
leg.AddEntry(Sim,"H(125)","f")
leg.SetTextSize(0.038)
leg.SetFillColor(10)
leg.SetLineColor(10)

#Loop over all the events and fill histogram
for ievent,event in enumerate(ggTozz):
    gg.Fill(event.H_FSR)

for ievent,event in enumerate(qqTozz):
    qq.Fill(event.H_FSR)

for ievent,event in enumerate(DataSim):
    Sim.Fill(event.H_FSR)

for ievent,event in enumerate(t):
    Data.Fill(event.H_FSR)

#set histoand drew
qq.Draw()
gg.Draw("SAME")
Sim.Draw("SAME")
leg.Draw("SAME")
Data.Draw("PE1")
c.SaveAs("reult.png")
