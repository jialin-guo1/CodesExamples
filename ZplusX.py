import ROOT

t = ROOT.TChain("passedEvents")
t.Add("/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/2016_noDuplicates_new.root")

#Z+X
c = ROOT.TCanvas()

ZX1 = ROOT.TH1D("ZX1","ZX1",10,0,1)
ZX1.GetYaxis().SetTitle("Events / 2 GeV")
ZX1.SetLineColor(ROOT.kRed)

ZX2 = ROOT.TH1D("ZX2","ZX2",10,0,1)
ZX2.GetYaxis().SetTitle("Events / 2 GeV")
ZX2.SetLineColor(ROOT.kBlue)

ZX3 = ROOT.TH1D("ZX3","ZX3",10,0,1)
ZX3.GetYaxis().SetTitle("Events / 2 GeV")
ZX3.SetLineColor(ROOT.kBlack)

ZX4 = ROOT.TH1D("ZX4","ZX4",10,0,1)
ZX4.GetYaxis().SetTitle("Events / 2 GeV")
ZX4.SetLineColor(7)

for ievent,event in enumerate(t):
    ZX1.Fill(event.lep_RelIsoNoFSR1)
    ZX2.Fill(event.lep_RelIsoNoFSR2)
    ZX3.Fill(event.lep_RelIsoNoFSR3)
    ZX4.Fill(event.lep_RelIsoNoFSR4)

ZX1.Draw()
ZX2.Draw("same")
ZX3.Draw("same")
ZX4.Draw("same")
c.SaveAs("ZX.png")
