#!/usr/bin/env python
import ROOT
from array import array

ROOT.gStyle.SetOptStat(False)

# input file and get tree
ggTozz2e2u = ROOT.TChain("passedEvents")
ggTozz2e2u.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC1/GluGluToContinToZZTo2e2mu_FL01.root")
ggTozz2e2u.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2017/MC/GluGluToContinToZZTo2e2mu_CN01.root")
ggTozz2e2u.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2018/MC/GluGluToContinToZZTo2e2mu_FL01.root")

ggTozz2u2t = ROOT.TChain("passedEvents")
ggTozz2u2t.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC1/GluGluToContinToZZTo2mu2tau_FL02.root")
ggTozz2u2t.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2017/MC/GluGluToContinToZZTo2mu2tau_CN01.root")
ggTozz2u2t.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2018/MC/GluGluToContinToZZTo2mu2tau_CN01.root")

ggTozz2e2t = ROOT.TChain("passedEvents")
ggTozz2e2t.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC1/GluGluToContinToZZTo2e2tau_FL01.root")
ggTozz2e2t.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2017/MC/GluGluToContinToZZTo2e2tau_CN01.root")
ggTozz2e2t.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2018/MC/GluGluToContinToZZTo2e2tau_CN01.root")

ggTozz4e = ROOT.TChain("passedEvents")
ggTozz4e.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC1/GluGluToContinToZZTo4e_FL01.root")
ggTozz4e.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2017/MC/GluGluToContinToZZTo4e_CN01.root")
ggTozz4e.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2018/MC/GluGluToContinToZZTo4e_CN01.root")

ggTozz4u = ROOT.TChain("passedEvents")
ggTozz4u.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC1/GluGluToContinToZZTo4mu_FL02.root")
ggTozz4u.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2017/MC/GluGluToContinToZZTo4mu_CN01.root")
ggTozz4u.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2018/MC/GluGluToContinToZZTo4mu_CN01.root")

ggTozz4t = ROOT.TChain("passedEvents")
ggTozz4t.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC1/GluGluToContinToZZTo4tau_FL01.root")
ggTozz4t.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2017/MC/GluGluToContinToZZTo4tau_CN01.root")
ggTozz4t.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2018/MC/GluGluToContinToZZTo4tau_CN01.root")

qqTozz = ROOT.TChain("passedEvents")
#qqTozz.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/ZZTo4L_FL04.root")
#qqTozz.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/ZZTo4L_FL02.root")
qqTozz.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC1/ZZTo4L_FLALL.root")
qqTozz.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2017/MC/ZZTo4L_FL01.root")
qqTozz.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2018/MC/ZZTo4L_CN01.root")

#qqTozz0 = ROOT.TChain("passedEvents")
#qqTozz0.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/ZZTo4L_FL02.root")

DataSim_gg = ROOT.TChain("passedEvents")
DataSim_gg.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC1/GluGluHToZZTo4L*.root")
DataSim_gg.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2017/MC/GluGluHToZZTo4L_CN01.root")
DataSim_gg.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2018/MC/GluGluHToZZTo4L_CN01.root")

DataSim_qq = ROOT.TChain("passedEvents")
DataSim_qq.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC1/VBF_HToZZTo4L*.root")
DataSim_qq.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2017/MC/VBF_HToZZTo4L_CN01.root")
#DataSim_qq.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2018/MC/VBF_HToZZTo4L_CN01.root")

DataSim_WplusH = ROOT.TChain("passedEvents")
DataSim_WplusH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC1/WplusH*.root")
DataSim_WplusH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2017/MC/WplusH_HToZZTo4L_CN01.root")
DataSim_WplusH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2018/MC/WplusH_HToZZTo4L_CN01.root")

DataSim_WminH = ROOT.TChain("passedEvents")
DataSim_WminH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC1/WminusH*.root")
DataSim_WminH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2017/MC/WminusH_HToZZTo4L_CN01.root")
DataSim_WminH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2018/MC/WminusH_HToZZTo4L_CN01.root")


DataSim_ZH = ROOT.TChain("passedEvents")
DataSim_ZH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC1/ZH_HToZZ*.root")
DataSim_ZH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2017/MC/ZH_HToZZ_4L_CN01.root")
DataSim_ZH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2018/MC/ZH_HToZZ_4L_CN01.root")

DataSim_ttH = ROOT.TChain("passedEvents")
DataSim_ttH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/MC1/ZH_HToZZ*.root")
DataSim_ttH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2017/MC/ttH_HToZZ_CN01.root")
DataSim_ttH.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2018/MC/ttH_HToZZ_CN01.root")


#Signal = ROOT.TFile('/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/2016_allsignal_new.root')
#t = Signal.Get('passedEvents')
t = ROOT.TChain("passedEvents")
t.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/2016_noDuplicates_new.root")
t.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2017/2017_noDuplicates_new.root")
t.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2018/2018_noDuplicates_new.root")

# book histogram and canvas
c = ROOT.TCanvas()
gg2e2u = ROOT.TH1D("gg->ggTozz2e2u","Backgrund(2016)",50,70,170)
gg2e2u.SetFillColor(ROOT.kBlue)
gg2e2u.GetYaxis().SetTitle("Events / 2 GeV")

gg2u2t = ROOT.TH1D("gg->ggTozz2u2t","Backgrund(2016)",50,70,170)
gg2u2t.SetFillColor(ROOT.kBlue)
gg2u2t.GetYaxis().SetTitle("Events / 2 GeV")

gg2e2t = ROOT.TH1D("gg->ggTozz2e2t","Backgrund(2016)",50,70,170)
gg2e2t.SetFillColor(ROOT.kBlue)
gg2e2t.GetYaxis().SetTitle("Events / 2 GeV")

gg4e = ROOT.TH1D("gg->ggTozz4e","Backgrund(2016)",50,70,170)
gg4e.SetFillColor(ROOT.kBlue)
gg4e.GetYaxis().SetTitle("Events / 2 GeV")

gg4u = ROOT.TH1D("gg->ggTozz4u","Backgrund(2016)",50,70,170)
gg4u.SetFillColor(ROOT.kBlue)
gg4u.GetYaxis().SetTitle("Events / 2 GeV")

gg4t = ROOT.TH1D("gg->ggTozz4t","Backgrund(2016)",50,70,170)
gg4t.SetFillColor(ROOT.kBlue)
gg4t.GetYaxis().SetTitle("Events / 2 GeV")

qq = ROOT.TH1D("qq->zz","Backgrund(2016)",50,70,170)
qq.SetFillColor(7)
qq.GetYaxis().SetTitle("Events / 2 GeV")

qq0 = ROOT.TH1D("qq->zz0","Backgrund(2016)",50,70,170)
qq0.SetFillColor(7)
qq0.GetYaxis().SetTitle("Events / 2 GeV")

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


#Data = ROOT.TH1D("Data","Data(2016)",50,70,170)
Data = ROOT.TH1D("Data","result(2017)",50,70,170)
Data.GetYaxis().SetTitle("Events / 2 GeV")
Data.SetMarkerStyle(20)
Data.SetMarkerColor(ROOT.kBlack)
Data.SetMarkerSize(1.2)
Data.SetLineColor(ROOT.kBlack)
Data.SetLineWidth(1)
#Data.SetStats(ROOT.kFALSE)


#Loop over all the events and fill histogram
for ievent,event in enumerate(ggTozz2e2u):
    gg2e2u.Fill(event.H_FSR,35.9*1000*0.0039*event.weight*event.k_gg/event.cross)

for ievent,event in enumerate(ggTozz2e2t):
    gg2e2t.Fill(event.H_FSR,35.9*1000*0.00319*event.weight*event.k_gg/event.cross)

for ievent,event in enumerate(ggTozz2u2t):
    gg2u2t.Fill(event.H_FSR,35.9*1000*0.00319*event.weight*event.k_gg/event.cross)

for ievent,event in enumerate(ggTozz4e):
    gg4e.Fill(event.H_FSR,35.9*1000*0.00159*event.weight*event.k_gg/event.cross)

for ievent,event in enumerate(ggTozz4u):
    gg4u.Fill(event.H_FSR,35.9*1000*0.00159*event.weight*event.k_gg/event.cross)

for ievent,event in enumerate(ggTozz4t):
    gg4t.Fill(event.H_FSR,35.9*1000*0.00159*event.weight*event.k_gg/event.cross)

for ievent,event in enumerate(qqTozz):
    qq.Fill(event.H_FSR,35.9*1000*1.256*event.weight*event.k_qq_qcd_M*event.k_qq_ewk/event.cross)

#for ievent,event in enumerate(qqTozz0):
#    qq0.Fill(event.H_FSR,35.9*1000*event.weight*event.k_qq_qcd_M*event.k_qq_ewk)

for ievent,event in enumerate(DataSim_gg):
    Sim_gg.Fill(event.H_FSR,35.9*12.31*event.weight/event.cross)

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

Data.Draw("E1")
Sim.Draw("same hito")
c.SaveAs("Sum.png")

ggSum = ROOT.TH1D("qqSum","Backgrund(2016)",50,70,170)
ggSum.SetFillColor(ROOT.kBlue)
ggSum.Sumw2()
ggSum.Add(gg2e2u,gg2u2t)
print "sum 1 = "+ str(ggSum.Integral())
ggSum.Add(ggSum,gg2e2t)
print "sum 2 = "+ str(ggSum.Integral())
ggSum.Add(ggSum,gg4e)
print "sum 3 = "+ str(ggSum.Integral())
ggSum.Add(ggSum,gg4u)
print "sum 4 = "+ str(ggSum.Integral())
ggSum.Add(ggSum,gg4t)
print "sum 5 = "+ str(ggSum.Integral())


#qqSum = ROOT.TH1D("qqSum","Backgrund(2016)",50,70,170)
#qqSum.SetFillColor(7)
#qqSum.Sumw2()
#qqSum.Add(qq,qq0)

print "number of all MC = " + str(Sim.Integral())


# set leg
leg = ROOT.TLegend(0.7, 0.7, 0.85, 0.85)
leg.AddEntry(Data,"Data","PE1")
leg.AddEntry(ggSum,"gg->zz","f")
leg.AddEntry(qq,"qq->zz","f")
leg.AddEntry(Sim_gg,"H(125)","f")
leg.SetTextSize(0.038)
leg.SetFillColor(10)
leg.SetLineColor(10)


#set histo and drew
Data.Draw("E1")
hstack = ROOT.THStack("hstack","2016reuslt")
hstack.Add(ggSum)
hstack.Add(qq)
hstack.Add(Sim)
#hstack.Add(Sim_gg)
#hstack.Add(Sim_qq)
#hstack.Add(Sim_ZH)
#hstack.Add(Sim_ttH)
#hstack.Add(Sim_WminH)
#hstack.Add(Sim_WplusH)
hstack.Draw("same histo")
#Sim.Draw("same histo")
#qq.Draw("same histo")
#ggSum.Draw("same histo")
leg.Draw()
Data.Draw("same E1")
c.SaveAs("result2017.png")

#qq.Draw("histo")
#Data.Draw("same E1")
#c.SaveAs("qqZZ.png")
