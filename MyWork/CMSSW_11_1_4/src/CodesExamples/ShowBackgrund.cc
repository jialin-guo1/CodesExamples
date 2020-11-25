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

bool passedZ1LSelection;
bool passedFullSelection;
bool passedZ4lSelection;
int final4mu=1;
int final4e=2;
int final2mu2e=3;
int final2e2mu=4;
Int_t finalState;
Int_t lep_Hindex[4];
std::vector<int> *lep_id;
std::vector<float> *lep_pt,*lep_eta, *lep_phi, *lep_mass, *lepFSR_pt, *lepFSR_eta, *lepFSR_phi;
std::vector<float> *H_pt, *H_eta, *H_phi ,*H_Mass;

TH1D* ZX = new TH1d("ZX","ZX",50,70,170);


//function to set TTress address
void SetAddressHisto(TTree* t);

void ShowBackgrund(){
  TCanvas* c = new TCanvas();
  TChain* chain = new TChain("Ana/passedEvents");
  chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_MC/GluGluHToZZTo4L_M120_13TeV_powheg2_JHUgenV6_pythia8/crab_GluGluHToZZTo4L_M120_13TeV_powheg2_JHUgenV6_pythia8_RunIISummer16MiniAODv2/201020_054859/0000/GluGluHToZZTo4L_M120_13TeV_powheg2_JHUgenV6_pythia8_RunIISummer16MiniAODv2_10.root")
  SetAddressHisto(chain);
  for(int i=0; i<chain->GetEntries(); i++){
    chain->GetEntry(i);
    if(lep_id.size()>3) continue;
    if(!passedZ1LSelection) continue;
    for(int j=0; j<lep_id.size()){
      TLorentzVector l1, l2 , X;
      l1.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[0]],(*lepFSR_eta)[lep_Hindex[0]],(*lepFSR_phi)[lep_Hindex[0]],(*lep_mass)[lep_Hindex[0]]);
      l2.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[1]],(*lepFSR_eta)[lep_Hindex[1]],(*lepFSR_phi)[lep_Hindex[1]],(*lep_mass)[lep_Hindex[1]]);
      X.SetPtEtaPhiM((*lepFSR_pt)[lep_Hindex[2]],(*lepFSR_eta)[lep_Hindex[2]],(*lepFSR_phi)[lep_Hindex[2]],(*lep_mass)[lep_Hindex[2]]);

    }

  }
}



void SetAddressHisto(TTree* t){
  t->SetBranchAddress("passedZ1LSelection",&passedZ1LSelection);
  t->SetBranchAddress("lep_Hindex",&lep_Hindex);
  t->SetBranchAddress("lep_id",&lep_id);
  t->SetBranchAddress("lep_pt",&lep_pt);
  t->SetBranchAddress("lep_eta",&lep_eta);
  t->SetBranchAddress("lep_phi",&lep_phi);
  t->SetBranchAddress("lep_mass"&lep_mass);
  t->SetBranchAddress("finalState",&finalState);
  t->SetBranchAddress("lepFSR_pt",&lepFSR_pt);
  t->SetBranchAddress("lepFSR_eta",&lepFSR_eta);
  t->SetBranchAddress("lepFSR_phi",&lepFSR_phi);

}
