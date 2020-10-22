#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TH1D.h"
#include "TMath.h"
#include <iostream>
using namespace std;

//declar some variables
bool passedZ1LSelection;
bool passedFullSelection;
int finalid1=1;
int finalid2=2;
Int_t finalState;
std::vector<int> *lep_id;
std::vector<float> *lep_pt,*lep_eta, *lep_phi, *lep_mass;
TH1D *Muon_pt = new TH1D("Muon_pt","Muon_pt",50,0,100);
TH1D *Muon_phi = new TH1D("Muon_phi","Muon_phi",50,-TMath::Pi(),TMath::Pi());
TH1D *Muon_eta = new TH1D("Muon_eta","Muon_eta",50,-3,3);
TH1D *Muon_mass = new TH1D("Muon_mass","Muon_mass",50,0,1);
TH1D *Muon_id = new TH1D("Muon_id","Muon_id",50,-15,15);
TH1D *Electron_pt = new TH1D("Electron_pt","Electron_pt",50,0,100);
TH1D *Electron_eta = new TH1D("Electron_eta","Electron_eta",50,-3,3);
TH1D *Electron_phi = new TH1D("Electron_phi","Electron_phi",50,-TMath::Pi(),TMath::Pi());
TH1D *H_mass = new TH1D("H_mass","H_mass",50,70,170);
TH1D *Fourlep_mass = new TH1D("4l_mass","4l_mass",1000,70,3000);


//function to set TTress address
void SetAddressHisto(TTree* t); 

void PlotPreSelection(){
  TCanvas* c = new TCanvas();
//chain all the files
  TChain* chain = new TChain("Ana/passedEvents");
  chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/SingleMuon/crab_SingleMuon_Run2018A-17Sep2018-v2/201018_142705/0000/*.root");
  SetAddressHisto(chain);
  cout<<chain->GetEntries()<<endl;
// loop over the TFiles that have been added to this chain
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
      if(lep_id->size()<4) continue;
      if(passedZ1LSelection==!0) continue;
      TLorentzVector lep1, lep2, lep3, lep4;
      lep1.SetPtEtaPhiM((*lep_pt)[0],(*lep_eta)[0],(*lep_phi)[0],(*lep_mass)[0]);
      lep2.SetPtEtaPhiM((*lep_pt)[1],(*lep_eta)[1],(*lep_phi)[1],(*lep_mass)[1]);
      lep3.SetPtEtaPhiM((*lep_pt)[2],(*lep_eta)[2],(*lep_phi)[2],(*lep_mass)[2]);
      lep4.SetPtEtaPhiM((*lep_pt)[3],(*lep_eta)[3],(*lep_phi)[3],(*lep_mass)[3]);
      TLorentzVector Higgs=lep1+lep2+lep3+lep4;
      H_mass->Fill(Higgs.M());
      Fourlep_mass->Fill(Higgs.M()); 
  }
      Muon_pt->Draw();  c->SaveAs("Mu_pt.png");
      Muon_eta->Draw(); c->SaveAs("Mu_eta.png");
      Muon_phi->Draw(); c->SaveAs("Mu_phi.png");
      Muon_id->Draw();  c->SaveAs("Mu_id.png");
      Electron_pt->Draw(); c->SaveAs("e_pt.png");
      Electron_eta->Draw(); c->SaveAs("e_eta.png");
      Electron_phi->Draw(); c->SaveAs("e_phi.png");
      H_mass->Draw(); c->SaveAs("H_mass.png");
      Fourlep_mass->Draw(); c->SaveAs("4l_mass.png");
}


//define function
void SetAddressHisto(TTree* t){
  t->SetBranchAddress("lep_id",&lep_id);
  t->SetBranchAddress("lep_pt",&lep_pt);
  t->SetBranchAddress("lep_eta",&lep_eta);
  t->SetBranchAddress("lep_phi",&lep_phi);
  t->SetBranchAddress("lep_mass",&lep_mass);
  t->SetBranchAddress("passedZ1LSelection",&passedZ1LSelection);
  t->SetBranchAddress("passedFullSelection",&passedFullSelection);
  t->SetBranchAddress("finalState",&finalState);
}
