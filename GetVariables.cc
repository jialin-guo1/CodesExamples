#include <iostream>
#include <string>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
using namespace std;

void GetVariables(){
  bool passedTring;
  bool passedZ1LSelection;
  bool passedZ4lSelection;
  bool passedFullSelection;
  std::vector<int> *lep_id;
  std::vector<float> *lep_pt;
  std::vector<float> *lep_eta;
  std::vector<float> *lep_phi;
  std::vector<float> *lep_mass;
  std::vector<float> *lepFSR_pt;
  std::vector<float> *lepFSR_eta;
  std::vector<float> *lepFSR_phi;
  std::vector<float> *lepFSR_mass;
  Int_t finalState;
  Int_t lep_Hindex;
  vector<float>   *Z_pt;
  vector<float>   *Z_eta;
  vector<float>   *Z_phi;
  vector<float>   *Z_mass;
  vector<float>   *Z_noFSR_pt;
  vector<float>   *Z_noFSR_eta;
  vector<float>   *Z_noFSR_phi;
  vector<float>   *Z_noFSR_mass;
  std::vector<float> *H_pt;
  std::vector<float> *H_eta;
  std::vector<float> *H_phi;
  std::vector<float> *H_mass;

  TString outputFile = "Fourlep_2018_01.root";
  TFile ofile (outputFile,"recreate");
  TTree* newtree = new TTree("newtree","newtree");
  newtree->Branch("passedZ1LSelection", &passedZ1LSelection);
  newtree->Branch("passedFullSelection", &passedFullSelection);
  newtree->Branch("finalState",&finalState);
  newtree->Branch("passedZ4lSelection", &passedZ4lSelection);
  newtree->Branch("lep_id", &lep_id);
  newtree->Branch("lep_pt", &lep_pt);
  newtree->Branch("lep_eta", &lep_eta);
  newtree->Branch("lep_phi", &lep_phi);
  newtree->Branch("lep_mass", &lep_mass);
  newtree->Branch("lepFSR_pt", &lepFSR_pt);
  newtree->Branch("lepFSR_eta", &lepFSR_eta);
  newtree->Branch("lepFSR_phi", &lepFSR_phi);
  newtree->Branch("lepFSR_mass", &lepFSR_mass);
  newtree->Branch("Z_pt",&Z_pt);
  newtree->Branch("Z_eta",&Z_eta);
  newtree->Branch("Z_phi",&Z_phi);
  newtree->Branch("Z_mass",&Z_mass);
  newtree->Branch("Z_noFSR_pt",&Z_noFSR_pt);
  newtree->Branch("Z_noFSR_eta",&Z_noFSR_eta);
  newtree->Branch("Z_noFSR_phi",&Z_noFSR_phi);
  newtree->Branch("Z_noFSR_mass",&Z_noFSR_mass);
  newtree->Branch("H_pt",&H_pt);
  newtree->Branch("H_eta",&H_eta);
  newtree->Branch("H_phi",&H_phi);
  newtree->Branch("H_mass",&H_mass);

  TChain* chain = new TChain("Ana/passedEvents");
  chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/SingleMuon/crab_SingleMuon_Run2018A-17Sep2018-v2/201018_142705/0000/*.root");
  chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/SingleMuon/crab_SingleMuon_Run2018A-17Sep2018-v2/201018_142705/0001/*.root");
  chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/SingleMuon/crab_SingleMuon_Run2018A-17Sep2018-v2/201018_142705/0002/*.root");
  chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/DoubleMuon/crab_DoubleMuon_Run2018A-17Sep2018-v2/201018_142029/0000/*.root");
  chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/EGamma/crab_EGamma_Run2018A-17Sep2018-v2/201018_142452/0000/*.root");
  chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/EGamma/crab_EGamma_Run2018A-17Sep2018-v2/201018_142452/0001/*.root");
  chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/EGamma/crab_EGamma_Run2018A-17Sep2018-v2/201018_142452/0002/*.root");
  chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/EGamma/crab_EGamma_Run2018A-17Sep2018-v2/201018_142452/0003/*.root");
  chain->Add("/pnfs/ihep.ac.cn/data/cms/store/user/guoj/2018data/UFHZZAnalysisRun2/myTask_Data/MuonEG/crab_MuonEG_Run2018A-17Sep2018-v1/201018_142240/0000/*.root");

  long64_t nentries = chain->GetEntries();
  lep_id = 0;
  lep_pt = 0;
  lep_eta = 0;
  lep_phi = 0;
  lep_mass = 0;
  lepFSR_pt = 0;
  lepFSR_eta = 0;
  lepFSR_phi = 0;
  lepFSR_mass = 0;
  Z_pt = 0;
  Z_eta = 0;
  Z_phi = 0;
  Z_mass = 0;
  Z_noFSR_pt = 0;
  Z_noFSR_eta = 0;
  Z_noFSR_phi = 0;
  Z_noFSR_mass = 0;
  H_pt = 0;
  H_eta = 0;
  H_phi = 0;
  H_mass = 0;

  chain->SetBranchAddress("passedZ1LSelection",&passedZ1LSelection);
  chain->SetBranchAddress("passedZ4lSelection",&passedZ4lSelection);
  chain->SetBranchAddress("passedFullSelection",&passedFullSelection);
  chain->SetBranchAddress("passedTrig", &passedTrig);
  chain->SetBranchAddress("finalState",&finalState);
  chain->SetBranchAddress("lep_Hindex",&lep_Hindex);
  chain->SetBranchAddress("lep_id", &lep_id);
  chain->SetBranchAddress("lep_pt", &lep_pt);
  chain->SetBranchAddress("lep_eta", &lep_eta);
  chain->SetBranchAddress("lep_phi", &lep_phi);
  chain->SetBranchAddress("lep_mass", &lep_mass);
  chain->SetBranchAddress("lepFSR_pt", &lepFSR_pt);
  chain->SetBranchAddress("lepFSR_eta", &lepFSR_eta);
  chain->SetBranchAddress("lepFSR_phi", &lepFSR_phi);
  chain->SetBranchAddress("lepFSR_mass", &lepFSR_mass);
  chain->SetBranchAddress("Z_pt", &Z_pt);
  chain->SetBranchAddress("Z_eta", &Z_eta);
  chain->SetBranchAddress("Z_phi", &Z_phi);
  chain->SetBranchAddress("Z_mass", &Z_mass);
  chain->SetBranchAddress("Z_noFSR_pt", &Z_noFSR_pt);
  chain->SetBranchAddress("Z_noFSR_eta", &Z_noFSR_eta);
  chain->SetBranchAddress("Z_noFSR_phi", &Z_noFSR_phi);
  chain->SetBranchAddress("Z_noFSR_mass", &Z_noFSR_mass);
  chain->SetBranchAddress("H_pt",&H_pt);
  chain->SetBranchAddress("H_eta",&H_eta);
  chain->SetBranchAddress("H_phi",&H_phi);
  chain->SetBranchAddress("H_mass",&H_mass);

  for(long64_t i == 0; i<nentries; i++){
    chain->GetEntry(i);
    if (!passedTrig)    continue;
    newtree->Fill();
  }

  ofile.cd();
  newtree->Write();
  ofile.Close();

}
