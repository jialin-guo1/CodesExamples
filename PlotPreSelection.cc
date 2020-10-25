#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TH1D.h"
#include "TMath.h"
#include "TAxis.h"
#include <iostream>
using namespace std;

//declar some variables
bool passedZ1LSelection;
bool passedFullSelection;
bool passedZ4lSelection;
int final4mu=1;
int final4e=2;
int final2mu2e=3;
int final2e2mu=4;
Int_t finalState;
Int_t lep_Hindex[4]; //position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z2 sub
std::vector<int> *lep_id;
std::vector<float> *lep_pt,*lep_eta, *lep_phi, *lep_mass, *lepFSR_pt, *lepFSR_eta, *lepFSR_phi;
std::vector<float> *H_pt, *H_eta, *H_phi ,*H_Mass;
TH1D *Muon_pt = new TH1D("Muon_pt","Muon_pt",50,0,100);
TH1D *Muon_phi = new TH1D("Muon_phi","Muon_phi",50,-TMath::Pi(),TMath::Pi());
TH1D *Muon_eta = new TH1D("Muon_eta","Muon_eta",50,-3,3);
TH1D *Muon_mass = new TH1D("Muon_mass","Muon_mass",50,0,1);
TH1D *Muon_id = new TH1D("Muon_id","Muon_id",50,-15,15);
TH1D *Electron_pt = new TH1D("Electron_pt","Electron_pt",50,0,100);
TH1D *Electron_eta = new TH1D("Electron_eta","Electron_eta",50,-3,3);
TH1D *Electron_phi = new TH1D("Electron_phi","Electron_phi",50,-TMath::Pi(),TMath::Pi());
TH1D *H_mass = new TH1D("H_mass","H_mass",50,70,170);
TH1D *H_mass_noZinfo = new TH1D("H_mass_noZinfo","H_mass_noZinfo",50,70,170);
TH1D *Fourlep_mass = new TH1D("4l_mass","4l_mass",1000,70,3000);
TH1D *H_MassMC = new TH1D("H_MassMC","H_MassMC",50,70,170);


//function to set TTress address
void SetAddressHisto(TTree* t);

void PlotPreSelection(){
  TCanvas* c = new TCanvas();
//chain all the files for signal
  TChain* chain = new TChain("Ana/passedEvents");
  chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/SingleMuon/crab_SingleMuon_Run2018A-17Sep2018-v2/201018_142705/0000/SingleMuon_Run2018A-17Sep2018-v2_19*.root");
  //chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/SingleMuon/crab_SingleMuon_Run2018A-17Sep2018-v2/201018_142705/0001/*.root");
  //chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/SingleMuon/crab_SingleMuon_Run2018A-17Sep2018-v2/201018_142705/0002/*.root");
  //chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/DoubleMuon/crab_DoubleMuon_Run2018A-17Sep2018-v2/201018_142029/0000/*.root");
  //chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/EGamma/crab_EGamma_Run2018A-17Sep2018-v2/201018_142452/0000/*.root");
  //chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/EGamma/crab_EGamma_Run2018A-17Sep2018-v2/201018_142452/0001/*.root");
  //chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/EGamma/crab_EGamma_Run2018A-17Sep2018-v2/201018_142452/0002/*.root");
  //chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/EGamma/crab_EGamma_Run2018A-17Sep2018-v2/201018_142452/0003/*.root");
  //chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/MuonEG/crab_MuonEG_Run2018A-17Sep2018-v1/201018_142240/0000/*.root");
//chain all the files for MC
  TChain* chainMC = new TChain("Ana/passedEvents");
  chainMC->Add("/pnfs/ihep.ac.cn/data/cms/store/user/chenguan/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/crab_2018VBF_125/200510_155950/0000/2018VBF_125_13.root");

//loop over the MC TFiles that have been added to this chain
  SetAddressHisto(chainMC);
  for(int i=0;i<chainMC->GetEntries(); i++){
    chainMC->GetEntry(i);
    if(passedFullSelection==1&&(finalState==final4mu||finalState==final4e||finalState==final2mu2e||final2e2mu)){
      for(int j=0; j<H_pt->size(); j++){
        TLorentzVector H_massMC;
        H_massMC.SetPtEtaPhiM((*H_pt)[j],(*H_eta)[j],(*H_phi)[j],(*H_Mass)[j]);
        H_MassMC->Fill(H_massMC.M());
      }
    }
  }


// loop over the signal TFiles that have been added to this chain
  SetAddressHisto(chain);
  for(int i=0;i<chain->GetEntries(); i++){
      chain->GetEntry(i);
      for(int j=0; j<lep_id->size(); j++){
         if(lep_id->at(j)==13 || lep_id->at(j)==-13){
        // cout<<"find a pair muon"<<i<<endl;
            Muon_pt->Fill(lep_pt->at(j));
            Muon_eta->Fill((*lep_eta)[j]);
            Muon_phi->Fill((*lep_phi)[j]);
            Muon_id->Fill((*lep_id)[j]);
         }
         if(lep_id->at(j)==11 || lep_id->at(j)==-11){
           Electron_pt->Fill(lep_pt->at(j));
           Electron_eta->Fill(lep_eta->at(j));
           Electron_phi->Fill(lep_phi->at(j));
         }
      }

      // h mass spectrum
      if(lep_id->size()<4) continue;
      TLorentzVector lep1_noZ, lep2_noZ, lep3_noZ, lep4_noZ;
      lep1_noZ.SetPtEtaPhiM((*lep_pt)[0],(*lep_eta)[0],(*lep_phi)[0],(*lep_mass)[0]);
      lep2_noZ.SetPtEtaPhiM((*lep_pt)[1],(*lep_eta)[1],(*lep_phi)[1],(*lep_mass)[1]);
      lep3_noZ.SetPtEtaPhiM((*lep_pt)[2],(*lep_eta)[2],(*lep_phi)[2],(*lep_mass)[2]);
      lep4_noZ.SetPtEtaPhiM((*lep_pt)[3],(*lep_eta)[3],(*lep_phi)[3],(*lep_mass)[3]);
      TLorentzVector Higgs_noZ = lep1_noZ+lep2_noZ+lep3_noZ+lep4_noZ;
      H_mass_noZinfo->Fill(Higgs_noZ.M());
      Fourlep_mass->Fill(Higgs_noZ.M());

      if(passedFullSelection==1&&(finalState==final4mu||finalState==final4e||finalState==final2mu2e||final2e2mu)){
      TLorentzVector lep1, lep2, lep3, lep4;
      lep1.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[0]],(*lepFSR_eta)[lep_Hindex[0]],(*lepFSR_phi)[lep_Hindex[0]],(*lep_mass)[lep_Hindex[0]]);
      lep2.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[1]],(*lepFSR_eta)[lep_Hindex[1]],(*lepFSR_phi)[lep_Hindex[1]],(*lep_mass)[lep_Hindex[1]]);
      lep3.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[2]],(*lepFSR_eta)[lep_Hindex[2]],(*lepFSR_phi)[lep_Hindex[2]],(*lep_mass)[lep_Hindex[2]]);
      lep4.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[3]],(*lepFSR_eta)[lep_Hindex[3]],(*lepFSR_phi)[lep_Hindex[3]],(*lep_mass)[lep_Hindex[3]]);
      TLorentzVector Higgs=lep1+lep2+lep3+lep4;
      H_mass->Fill(Higgs.M());
    }
  }
      Muon_pt->Draw();  c->SaveAs("Mu_pt.png");
      Muon_eta->Draw(); c->SaveAs("Mu_eta.png");
      Muon_phi->Draw(); c->SaveAs("Mu_phi.png");
      Muon_id->Draw();  c->SaveAs("Mu_id.png");
      Electron_pt->Draw(); c->SaveAs("e_pt.png");
      Electron_eta->Draw(); c->SaveAs("e_eta.png");
      Electron_phi->Draw(); c->SaveAs("e_phi.png");
      H_mass_noZinfo->Draw(); c->SaveAs("H_mass_noZinfo.png");
      Fourlep_mass->Draw(); c->SaveAs("4l_mass.png");
      TCanvas* c1 = new TCanvas();
      c1->cd();
      TLegend *leg=new TLegend(0.62, 0.70, 0.82, 0.88);
      leg->AddEntry(H_mass,"Data","PE1");
      leg->AddEntry(H_MassMC,"MC","f");
      leg->SetFillColor(kRed);
      leg->Draw();
      H_MassMC->SetFillColor(kRed);
      H_MassMC->GetYaxis()->SetTitle("Events / 2 GeV");
      H_mass->GetYaxis()->SetTitle("Events / 2 GeV");
      H_mass->SetMarkerStyle(20);
      H_mass->SetMarkerColor(kBlack);
      H_mass->SetMarkerSize(1.2);
      H_mass->SetLineColor(kBlack);
      H_mass->SetLineWidth(1);
      H_mass->Draw("P E1");
      THStack hs("hs","Hmass");
      hs.Add(H_MassMC);
      hs.Add(H_mass);
      hs.Draw("SAME");
      H_mass->Draw("SAME P E1");
      c1->SaveAs("H_mass.png");
}


//define function
void SetAddressHisto(TTree* t){
  t->SetBranchAddress("lep_Hindex",&lep_Hindex);
  t->SetBranchAddress("lep_id",&lep_id);
  t->SetBranchAddress("lep_pt",&lep_pt);
  t->SetBranchAddress("lep_eta",&lep_eta);
  t->SetBranchAddress("lep_phi",&lep_phi);
  t->SetBranchAddress("lep_mass",&lep_mass);
  t->SetBranchAddress("passedZ1LSelection",&passedZ1LSelection);
  t->SetBranchAddress("passedFullSelection",&passedFullSelection);
  t->SetBranchAddress("passedZ4lSelection",&passedZ4lSelection);
  t->SetBranchAddress("finalState",&finalState);
  t->SetBranchAddress("lepFSR_pt",&lepFSR_pt);
  t->SetBranchAddress("lepFSR_eta",&lepFSR_eta);
  t->SetBranchAddress("lepFSR_phi",&lepFSR_phi);
  t->SetBranchAddress("H_pt",&H_pt);
  t->SetBranchAddress("H_eta",&H_eta);
  t->SetBranchAddress("H_phi",&H_phi);
  t->SetBranchAddress("H_mass",&H_Mass);
}
