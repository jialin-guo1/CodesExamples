#!/usr/bin/env python
import ROOT
from array import array

ROOT.gStyle.SetOptStat(False)

# input file and get tree
ggTozz = ROOT.TChain("passedEvents")
ggTozz.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC/GluGluToContinToZZTo2e2mu_FL01.root")

qqTozz = ROOT.TChain("passedEvents")
qqTozz.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC/ZZTo4L_FL01.root")

DataSim = ROOT.TChain("passedEvents")
DataSim.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC/GluGluHToZZTo4L_CN01.root")
#DataSim.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC/VBF_HToZZTo4L*.root")

Signal = ROOT.TFile('/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/2016_allsignal.root')
t = Signal.Get('passedEvents')

# book histogram and canvas
c = ROOT.TCanvas()

gg = ROOT.TH1D("gg->zz","Backgrund(2016)",50,70,170)
gg.SetFillColor(ROOT.kBlue)
gg.GetYaxis().SetTitle("Events / 2 GeV")

qq = ROOT.TH1D("qq->zz","Backgrund(2016)",50,70,170)
qq.SetFillColor(7)
qq.GetYaxis().SetTitle("Events / 2 GeV")

Sim = ROOT.TH1D("sim","Backgrund(2016)",50,70,170)
Sim.SetFillColor(ROOT.kRed)
Sim.GetYaxis().SetTitle("Events / 2 GeV")

Data = ROOT.TH1D("Data","Data",50,70,170)
Data.GetYaxis().SetTitle("Events / 2 GeV")
Data.SetMarkerStyle(20)
Data.SetMarkerColor(ROOT.kBlack)
Data.SetMarkerSize(1.2)
Data.SetLineColor(ROOT.kBlack)
Data.SetLineWidth(1)
#Data.SetStats(ROOT.kFALSE)

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

#normal
gg.Scale((0.359*1000)/gg.Integral())
qq.Scale((0.359*1000)/qq.Integral())
Sim.Scale((0.359*1000)/Sim.Integral())

#set histo and drew
Data.Draw("PE1")
hstack = ROOT.THStack("hstack","2016reuslt")
hstack.Add(gg)
hstack.Add(qq)
hstack.Add(Sim)
hstack.Draw("histo")
leg.Draw("same")
Data.Draw("samePE1")
c.SaveAs("result.png")
