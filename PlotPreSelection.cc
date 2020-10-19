#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TH1D.h"
#include "TMath.h"

void SetAddressHisto(TTree* t); //function to set TTress address

void PlotPreSelection(){
  TCanvas* c= new TCanvas("c","c",1000,1000);
  c->cd();

//chain all the files
  TChain* chain = new TChain("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/SingleMuon/crab_SingleMuon_Run2018A-17Sep2018-v2/201018_142705/0000
");
  chain.Add("*.root");

// loop over the TFiles that have been added to this chain
  TObjArray *fileElements=chain->GetListOfFiles();
  TIter next(fileElements);
  TChainElement *chEl=0;
  while (( chEl=(TChainElement*)next() )) {
     TFile f(chEl->GetTitle());
     TTree* t = (TTree*)f->Get("passedEvents");
     SetAddress(t);
     int entries=t->GetEntries();
     for(int i=0; i<entries; i++){
       if(lep_id[i]==13 || lep_id[i]==-13){
         Muon_pt->Fill(lep_pt[i]);
         Muon_eta->Fill(lep_eta[i]);
         Muon_phi->Fill(lep_phi[i]);
       }
     }
  }
  Muon_pt->Draw();
  Muon_eta->Draw();
  Muon_phi->Draw();
  c->SaveAs("PlotPreSelection.png");
}

//define function SetAddress()
void SetAddressHisto(TTree* t){
  t->SetBranchAddress("lep_id",&lep_id);
  t->SetBranchAddress("lep_pt",&lep_pt);
  t->SetBranchAddress("lep_eta",&lep_eta);
  t->SetBranchAddress("lep_phi",&lep_phi);
  t->SetBranchAddress("lep_mass",&lep_mass);

  TH1D *Muon_pt("Muon_pt","Muon_pt",50,15,100);
  TH1D *Muon_phi("Muon_phi","Muon_phi",50,-TMath::Pi(),TMath::Pi());
  TH1D *Muon_eta("Muon_eta","Muon_eta",50,-3,3);
}
