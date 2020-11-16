#!/usr/bin/env python
import ROOT
from array import array

ROOT.gStyle.SetOptStat(False)

# input file and get tree
ggTozz = ROOT.TChain("passedEvents")
ggTozz.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC3/GluGluToContinToZZ*.root")

qqTozz = ROOT.TChain("passedEvents")
qqTozz.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/ZZTo4L*.root")

DataSim_gg = ROOT.TChain("passedEvents")
DataSim_gg.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC3/GluGluHToZZTo4L*.root")

DataSim_qq = ROOT.TChain("passedEvents")
DataSim_qq.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC3/VBF_HToZZTo4L*.root")

DataSim_WplusH = ROOT.TChain("passedEvents")
DataSim_WplusH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC3/WplusH*.root")

DataSim_WminH = ROOT.TChain("passedEvents")
DataSim_WminH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC3/WminusH*.root")

DataSim_ZH = ROOT.TChain("passedEvents")
DataSim_ZH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC3/ZH_HToZZ*.root")

DataSim_ttH = ROOT.TChain("passedEvents")
DataSim_ttH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC3/ZH_HToZZ*.root")

Signal = ROOT.TFile('/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/2016_allsignal_new.root')
t = Signal.Get('passedEvents')

# book histogram and canvas
c = ROOT.TCanvas()

gg = ROOT.TH1D("gg->zz","Backgrund(2016)",50,70,170)
gg.SetFillColor(ROOT.kBlue)
gg.GetYaxis().SetTitle("Events / 2 GeV")

qq = ROOT.TH1D("qq->zz","Backgrund(2016)",50,70,170)
qq.SetFillColor(7)
qq.GetYaxis().SetTitle("Events / 2 GeV")

Sim_gg = ROOT.TH1D("Sim_gg","Backgrund(2016)",50,70,170)
Sim_gg.SetFillColor(ROOT.kRed)
Sim_gg.GetYaxis().SetTitle("Events / 2 GeV")

Sim_qq = ROOT.TH1D("Sim_qq","Backgrund(2016)",50,70,170)
Sim_qq.SetFillColor(ROOT.kRed)
Sim_qq.GetYaxis().SetTitle("Events / 2 GeV")

Sim_WplusH = ROOT.TH1D("Sim_WplusH","Backgrund(2016)",50,70,170)
Sim_WplusH.SetFillColor(ROOT.kRed)
Sim_WplusH.GetYaxis().SetTitle("Events / 2 GeV")

Sim_WminH = ROOT.TH1D("Sim_WminH","Backgrund(2016)",50,70,170)
Sim_WminH.SetFillColor(ROOT.kRed)
Sim_WminH.GetYaxis().SetTitle("Events / 2 GeV")

Sim_ttH = ROOT.TH1D("Sim_ttH","Backgrund(2016)",50,70,170)
Sim_ttH.SetFillColor(ROOT.kRed)
Sim_ttH.GetYaxis().SetTitle("Events / 2 GeV")

Sim_ZH = ROOT.TH1D("Sim_ZH","Backgrund(2016)",50,70,170)
Sim_ZH.SetFillColor(ROOT.kRed)
Sim_ZH.GetYaxis().SetTitle("Events / 2 GeV")


Data = ROOT.TH1D("Data","Data(2016)",50,70,170)
Data.GetYaxis().SetTitle("Events / 2 GeV")
Data.SetMarkerStyle(20)
Data.SetMarkerColor(ROOT.kBlack)
Data.SetMarkerSize(1.2)
Data.SetLineColor(ROOT.kBlack)
Data.SetLineWidth(1)
#Data.SetStats(ROOT.kFALSE)

#Loop over all the events and fill histogram
for ievent,event in enumerate(ggTozz):
    gg.Fill(event.H_FSR,35.9*1000*0.01434*event.weight/event.cross*event.k_gg)

for ievent,event in enumerate(qqTozz):
    qq.Fill(event.H_FSR,35.9*1000*1.256*event.weight/event.cross*event.k_qq_qcd_pt*event.k_qq_ewk)

for ievent,event in enumerate(DataSim_gg):
    Sim_gg.Fill(event.H_FSR,35.9*12.18*event.weight/event.cross)

for ievent,event in enumerate(DataSim_qq):
    Sim_qq.Fill(event.H_FSR,35.9*1.044*event.weight/event.cross)

for ievent,event in enumerate(DataSim_WplusH):
    Sim_WplusH.Fill(event.H_FSR,35.9*0.232*event.weight/event.cross)

for ievent,event in enumerate(DataSim_WminH):
    Sim_WminH.Fill(event.H_FSR,35.9*0.147*event.weight/event.cross)

for ievent,event in enumerate(DataSim_ZH):
    Sim_ZH.Fill(event.H_FSR,35.9*0.668*event.weight/event.cross)

for ievent,event in enumerate(DataSim_ttH):
    Sim_ttH.Fill(event.H_FSR,35.9*0.393*event.weight/event.cross)

for ievent,event in enumerate(t):
    Data.Fill(event.H_FSR)


print "number of ggHToZZ = " + str(Sim_gg.Integral())
print "number of VBFHToZZ = " + str(Sim_qq.Integral())
print "number of WplusH = " + str(Sim_WplusH.Integral())
print "number of WminH = " + str(Sim_WminH.Integral())
print "number of ZH = " + str(Sim_ZH.Integral())
print "number of ttH = " + str(Sim_ttH.Integral())
print "number of data = " + str(Data.Integral())

#normal
#gg.Scale(35.9*1000*0.00637)
#qq.Scale(35.9*1000*2.468)
#Sim_gg.Scale(35.9*12.18)
#Sim_qq.Scale(35.9*1.044)

Sim = ROOT.TH1D("Sim","Backgrund(2016)",50,70,170)
Sim.SetFillColor(ROOT.kRed)
Sim.Sumw2()
Sim.Add(Sim_gg,Sim_qq)
Sim.Add(Sim,Sim_ZH)
Sim.Add(Sim,Sim_WplusH)
Sim.Add(Sim,Sim_WminH)
Sim.Add(Sim,Sim_ttH)

print "number of all MC = " + str(Sim.Integral())


# set leg
leg = ROOT.TLegend(0.7, 0.7, 0.85, 0.85)
leg.AddEntry(Data,"Data","PE1")
leg.AddEntry(gg,"gg->zz","f")
leg.AddEntry(qq,"qq->zz","f")
leg.AddEntry(Sim_gg,"H(125)","f")
leg.SetTextSize(0.038)
leg.SetFillColor(10)
leg.SetLineColor(10)


#set histo and drew
Data.Draw("E1")
#hstack = ROOT.THStack("hstack","2016reuslt")
#hstack.Add(gg)
#hstack.Add(qq)
#hstack.Add(Sim)
#hstack.Draw("same histo")
Sim.Draw("same histo")
#qq.Draw("same histo")
gg.Draw("same histo")
leg.Draw()
Data.Draw("same E1")
c.SaveAs("All.png")

qq.Draw("histo")
Data.Draw("same E1")
c.SaveAs("qqZZ.png")
