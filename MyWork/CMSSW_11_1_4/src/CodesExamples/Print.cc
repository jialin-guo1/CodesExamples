#include <vector>
#include "TH1D.h"
#include <iostream>
#include "TLorentzVector.h"
void Print() {
    using namespace std;

    TFile *f1 = TFile::Open("/afs/cern.ch/user/l/lipe/public/WWPTMASS/GF_HH_SM_slc6_LHEBqrk.root");
    TFile *f2 = TFile::Open("/afs/cern.ch/user/l/lipe/public/WWPTMASS/GF_HH_SM_slc6_LHEBqrk_ZZ.root");

    TTree *t2 = (TTree*)f2->Get("otree");
    TTree *t1 = (TTree*)f1->Get("otree");
// Define some variables
    double gen_leading_WBoson_Pt;
    double gen_Subleading_WBoson_Pt;
    double gen_leading_WBoson_M;
    double gen_Subleading_WBoson_M;

    double gen_leading_ZBoson_Pt;
    double gen_Subleading_ZBoson_Pt;
    double gen_leading_ZBoson_M;
    double gen_Subleading_ZBoson_M;

    int Pass_select_num = 0;
    int numOfLargerMass = 0;
    int Pass_select_numZ = 0;
    int numOfLargerMassZ = 0;

// Get variables 
    t1->SetBranchAddress("gen_leading_WBoson_Pt",&gen_leading_WBoson_Pt);
    t1->SetBranchAddress("gen_Subleading_WBoson_Pt",&gen_Subleading_WBoson_Pt);
    t1->SetBranchAddress("gen_leading_WBoson_M",&gen_leading_WBoson_M);
    t1->SetBranchAddress("gen_Subleading_WBoson_M",&gen_Subleading_WBoson_M);

    t2->SetBranchAddress("gen_leading_ZBoson_Pt",&gen_leading_ZBoson_Pt);
    t2->SetBranchAddress("gen_Subleading_ZBoson_Pt",&gen_Subleading_ZBoson_Pt);
    t2->SetBranchAddress("gen_leading_ZBoson_M",&gen_leading_ZBoson_M);
    t2->SetBranchAddress("gen_Subleading_ZBoson_M",&gen_Subleading_ZBoson_M);


    int entries1=t1->GetEntries();
    int entries2=t2->GetEntries();
    cout<<"Totle number of events (W) = "<<entries1<<endl;
    cout<<"Totle number of events (Z) = "<<entries2<<endl;

//create histogram
    TH1D *leadZ_Mass = new TH1D("leadZ_Mass","leadZ_Mass",50,60,110);
    TH1D *subleadingZ_Mass = new TH1D("subleadingZ_Mass","subleadingZ_Mass",50,10,60);
    TH1D *DiZ_Mass = new TH1D("DiZ_Mass","DiZ_Mass",50,80,150);

//loop and print informations

    for(int i=0; i<entries1; i++){
        t1->GetEntry(i);
        if(gen_leading_WBoson_M>gen_Subleading_WBoson_M){
            numOfLargerMass ++;
            if(gen_leading_WBoson_Pt>gen_Subleading_WBoson_Pt){
                Pass_select_num ++;
            }
        }
    }

    for(int j=0; j<entries2; j++){
        t2->GetEntry(j);
        if(gen_leading_ZBoson_M>gen_Subleading_ZBoson_M){
            numOfLargerMassZ ++;
            if(gen_leading_ZBoson_Pt>gen_Subleading_ZBoson_Pt){
                Pass_select_numZ ++;
    
            }
        }
        leadZ_Mass->Fill(gen_leading_ZBoson_M);
        subleadingZ_Mass->Fill(gen_Subleading_ZBoson_M);
        DiZ_Mass->Fill(gen_leading_ZBoson_M+gen_Subleading_ZBoson_M);
    }


//Print the result
    cout<<"Number of event passing selection (W) = "<<Pass_select_num<<endl;
    cout<<"Pass rate (W) = "<<(double) Pass_select_num/entries1 <<endl;
    cout<<"number of leadingMass larger than subleadingMss (W) = "<< numOfLargerMass<<endl;
    cout<<"Number of event passing selection (Z) = "<<Pass_select_numZ<<endl;
    cout<<"number of leadingMass larger than subleadingMss (Z) = "<< numOfLargerMassZ<<endl;
    cout<<"Pass rate (Z) = "<<(double) Pass_select_numZ/entries2 <<endl;
//
// Draw
    TCanvas *c = new TCanvas("c","c",1000,800);
    c->Divide(3,1);
    c->cd(1);
    leadZ_Mass->Draw();
    c->cd(2);
    subleadingZ_Mass->Draw();
    c->cd(3);
    DiZ_Mass->Draw();
}
