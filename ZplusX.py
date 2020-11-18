import ROOT

Signal = ROOT.TFile('/afs/cern.ch/work/g/guoj/XToZZ_FullRunII/Data2016/2016_noDuplicates.root')
t = Signal.Get('passedEvents')

#Z+X
c = ROOT.TCanvas()

#ZX1 = ROOT.TH1D("ZX1","ZX1",10,0.3,1)
#ZX1.GetYaxis().SetTitle("Events / 2 GeV")
#ZX1.SetLineColor(ROOT.kRed)

#ZX2 = ROOT.TH1D("ZX2","ZX2",10,0.3,1)
#ZX2.GetYaxis().SetTitle("Events / 2 GeV")
#ZX2.SetLineColor(ROOT.kBlue)

#ZX3 = ROOT.TH1D("ZX3","ZX3",10,0.3,1)
#ZX3.GetYaxis().SetTitle("Events / 2 GeV")
#ZX3.SetLineColor(ROOT.kBlack)

#ZX4 = ROOT.TH1D("ZX4","ZX4",10,0.3,1)
#ZX4.GetYaxis().SetTitle("Events / 2 GeV")
#ZX4.SetLineColor(7)

#nZX = ROOT.TH1D("nZX","nZX",10,0,3)

dataR1 =ROOT.TH1D("dataR1","dataR1",50,70,170)
dataR1.SetLineColor(ROOT.kRed)

dataR2 =ROOT.TH1D("dataR2","dataR2",50,70,170)
dataR2.SetLineColor(ROOT.kBlue)

for ievent,event in enumerate(t):
    nlep = event.lep_pt.size()
    if(not event.passedZXCRSelection): continue
    if(event.nZXCRFailedLeptons==1):
        for i in range(nlep):
            l1FSR = ROOT.TLorentzVector()
            l2FSR = ROOT.TLorentzVector()
            l3FSR = ROOT.TLorentzVector()
            l4FSR = ROOT.TLorentzVector()
            l1FSR.SetPtEtaPhiM(event.lepFSR_pt[event.lep_Hindex[0]],event.lepFSR_eta[event.lep_Hindex[0]],event.lepFSR_phi[event.lep_Hindex[0]],event.lepFSR_mass[event.lep_Hindex[0]])
            l2FSR.SetPtEtaPhiM(event.lepFSR_pt[event.lep_Hindex[1]],event.lepFSR_eta[event.lep_Hindex[1]],event.lepFSR_phi[event.lep_Hindex[1]],event.lepFSR_mass[event.lep_Hindex[1]])
            l3FSR.SetPtEtaPhiM(event.lepFSR_pt[event.lep_Hindex[2]],event.lepFSR_eta[event.lep_Hindex[2]],event.lepFSR_phi[event.lep_Hindex[2]],event.lepFSR_mass[event.lep_Hindex[2]])
            l4FSR.SetPtEtaPhiM(event.lepFSR_pt[event.lep_Hindex[3]],event.lepFSR_eta[event.lep_Hindex[3]],event.lepFSR_phi[event.lep_Hindex[3]],event.lepFSR_mass[event.lep_Hindex[3]])
            H4massFSR = ROOT.TLorentzVector()
            H4massFSR = l1FSR+l2FSR+l3FSR+l4FSR
            dataR1.Fill(H4massFSR.M())
    if(event.nZXCRFailedLeptons==2):
        for i in range(nlep):
            l1FSR = ROOT.TLorentzVector()
            l2FSR = ROOT.TLorentzVector()
            l3FSR = ROOT.TLorentzVector()
            l4FSR = ROOT.TLorentzVector()
            l1FSR.SetPtEtaPhiM(event.lepFSR_pt[event.lep_Hindex[0]],event.lepFSR_eta[event.lep_Hindex[0]],event.lepFSR_phi[event.lep_Hindex[0]],event.lepFSR_mass[event.lep_Hindex[0]])
            l2FSR.SetPtEtaPhiM(event.lepFSR_pt[event.lep_Hindex[1]],event.lepFSR_eta[event.lep_Hindex[1]],event.lepFSR_phi[event.lep_Hindex[1]],event.lepFSR_mass[event.lep_Hindex[1]])
            l3FSR.SetPtEtaPhiM(event.lepFSR_pt[event.lep_Hindex[2]],event.lepFSR_eta[event.lep_Hindex[2]],event.lepFSR_phi[event.lep_Hindex[2]],event.lepFSR_mass[event.lep_Hindex[2]])
            l4FSR.SetPtEtaPhiM(event.lepFSR_pt[event.lep_Hindex[3]],event.lepFSR_eta[event.lep_Hindex[3]],event.lepFSR_phi[event.lep_Hindex[3]],event.lepFSR_mass[event.lep_Hindex[3]])
            H4massFSR = ROOT.TLorentzVector()
            H4massFSR = l1FSR+l2FSR+l3FSR+l4FSR
            dataR2.Fill(H4massFSR.M()


dataR1.Draw()
dataR2.Draw("same")
c.SaveAs("dataCR.png")




#    nZX.Fill(event.nZXCRFailedLeptons)
#    ZX1.Fill(event.lep_RelIsoNoFSR1)
#    ZX2.Fill(event.lep_RelIsoNoFSR2)
#    ZX3.Fill(event.lep_RelIsoNoFSR3)
#    ZX4.Fill(event.lep_RelIsoNoFSR4)

#ZX1.Draw()
#ZX2.Draw("same")
#ZX3.Draw("same")
#ZX4.Draw("same")
#c.SaveAs("ZX01.png")

#nZX.Draw()
#c.SaveAs("passedZXCR.png")
