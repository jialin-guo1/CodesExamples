// Package:    UFHZZ4LAna
// Class:      UFHZZ4LAna

// system include files
#include <memory>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iomanip>
#include <set>

#define PI 3.14159
// user include files
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TSpline.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
#include "TClonesArray.h"
#include "TCanvas.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

//HTXS
#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
//#include "SimDataFormats/HZZFiducial/interface/HZZFiducialVolume.h"

// PAT
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// Reco
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

// KD's
#include "ZZMatrixElement/MELA/interface/Mela.h"

//Helper
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LHelper.h"
//Muons
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LMuonAna.h"
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LMuonTree.h"
//Electrons
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LElectronTree.h"
//Photons
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LPhotonTree.h"
//Jets
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LJetTree.h"
//Final Leps
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LFinalLepTree.h"
//Sip
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LSipAna.h"
//PU
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LPileUp.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//GEN
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LGENAna.h"
//VBF Jets
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LJets.h"

// Jet energy correction
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"


// Kinematic Fit
#include "KinZfitter/KinZfitter/interface/KinZfitter.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"

// EWK corrections
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/EwkCorrections.h"

// JEC related
#include "PhysicsTools/PatAlgos/plugins/PATJetUpdater.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"

//JER related
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

//BTag Calibration

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

//Muon MVA
//#include "MuonMVAReader/Reader/interface/MuonGBRForestReader.hpp"

//
// class declaration
//
using namespace EwkCorrections;

class UFHZZ4LAna : public edm::EDAnalyzer {
public:
    explicit UFHZZ4LAna(const edm::ParameterSet&);
    ~UFHZZ4LAna();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    static bool sortByPt( const reco::GenParticle &p1, const reco::GenParticle &p2 ){ return (p1.pt() > p2.pt()); };

private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    virtual void beginRun(edm::Run const&, const edm::EventSetup& iSetup);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const& lumiSeg,edm::EventSetup const& eSetup);

    void findHiggsCandidate(std::vector< pat::Muon > &selectedMuons, std::vector< pat::Electron > &selectedElectrons, const edm::Event& iEvent);
    void findZ1LCandidate(const edm::Event& iEvent);

    //Helper Class
    HZZ4LHelper helper;
    //GEN
    HZZ4LGENAna genAna;
    //VBF
    HZZ4LJets jetHelper;
    //PU Reweighting
    edm::LumiReWeighting *lumiWeight;
    HZZ4LPileUp pileUp;
    //JES Uncertainties
    std::unique_ptr<JetCorrectionUncertainty> jecunc;
    // kfactors
    TSpline3 *kFactor_ggzz;
    std::vector<std::vector<float> > tableEwk;
    // data/MC scale factors
    TH2F *hElecScaleFac;
    TH2F *hElecScaleFac_Cracks;
    TH2F *hElecScaleFacGsf;
    TH2F *hElecScaleFacGsfLowET;
    TH2F *hMuScaleFac;
    TH2F *hMuScaleFacUnc;
    TH1D *h_pileup;
    TH1D *h_pileupUp;
    TH1D *h_pileupDn;
    std::vector<TH1F*> h_medians;
    TH2F *hbTagEffi;
    TH2F *hcTagEffi;
    TH2F *hudsgTagEffi;

    BTagCalibrationReader* reader;

    //Saved Events Trees
    TTree *passedEventsTree_All;

    void bookPassedEventTree(TString treeName, TTree *tree);
    void setTreeVariables( const edm::Event&, const edm::EventSetup&,
                           std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons,
                           std::vector<pat::Muon> recoMuons, std::vector<pat::Electron> recoElectrons,
                           std::vector<pat::Jet> goodJets, std::vector<float> goodJetQGTagger,
                           std::vector<float> goodJetaxis2, std::vector<float> goodJetptD, std::vector<int> goodJetmult,
                           std::vector<pat::Jet> selectedMergedJets);
    void setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles,
                         edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles,
                         edm::Handle<edm::View<reco::GenJet> > genJets);
    bool mZ1_mZ2(unsigned int& L1, unsigned int& L2, unsigned int& L3, unsigned int& L4, bool makeCuts);

    // -------------------------
    // RECO level information
    // -------------------------

    // Event Variables
    ULong64_t Run, Event, LumiSect;
    int nVtx, nInt;
    int finalState;
    std::string triggersPassed;
    bool passedTrig, passedFullSelection, passedZ4lSelection, passedQCDcut;
    bool passedZ1LSelection, passedZ4lZ1LSelection, passedZ4lZXCRSelection, passedZXCRSelection;
    int nZXCRFailedLeptons;


    // Event Weights
    float genWeight, pileupWeight, pileupWeightUp, pileupWeightDn, dataMCWeight, eventWeight, prefiringWeight;
    float k_ggZZ, k_qqZZ_qcd_dPhi, k_qqZZ_qcd_M, k_qqZZ_qcd_Pt, k_qqZZ_ewk;
    // pdf weights
    vector<float> qcdWeights;
    vector<float> nnloWeights;
    vector<float> pdfWeights;
    int posNNPDF;
    float pdfRMSup, pdfRMSdown, pdfENVup, pdfENVdown;
    // lepton variables
    vector<double> lep_pt; vector<double> lep_pterr; vector<double> lep_pterrold;
    vector<double> lep_p; vector<double> lep_ecalEnergy; vector<int> lep_isEB; vector<int> lep_isEE;
    vector<double> lep_eta; vector<double> lep_phi; vector<double> lep_mass;
    vector<double> lepFSR_pt; vector<double> lepFSR_eta; vector<double> lepFSR_phi; vector<double> lepFSR_mass;
    int lep_Hindex[4];//position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z2 sub
    float pTL1, pTL2, pTL3, pTL4;
    float etaL1, etaL2, etaL3, etaL4;
    int idL1, idL2, idL3, idL4;
    float mL1, mL2, mL3, mL4;
    float pTErrL1, pTErrL2, pTErrL3, pTErrL4;
    float phiL1, phiL2, phiL3, phiL4;
    float pTL1FSR, pTL2FSR, pTL3FSR, pTL4FSR;
    vector<float> lep_dataMC; vector<float> lep_dataMCErr;
    vector<int> lep_genindex; //position of lepton in GENlep_p4 (if gen matched, -1 if not gen matched)
    vector<int> lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId; // gen matching even if not in GENlep_p4
    vector<int> lep_id;
    vector<float> lep_mva; vector<int> lep_ecalDriven;
    vector<int> lep_tightId; vector<int> lep_tightIdSUS; vector<int> lep_tightIdHiPt; //vector<int> lep_tightId_old;
    vector<float> lep_Sip; vector<float> lep_IP; vector<float> lep_isoNH; vector<float> lep_isoCH; vector<float> lep_isoPhot;
    vector<float> lep_isoPU; vector<float> lep_isoPUcorr;
    vector<float> lep_RelIso; vector<float> lep_RelIsoNoFSR; vector<float> lep_MiniIso;
    vector<float> lep_ptRatio; vector<float> lep_ptRel;
    vector<int> lep_missingHits;
    vector<string> lep_filtersMatched; // for each lepton, all filters it is matched to
    int nisoleptons;
    double muRho, elRho, rhoSUS;

    // tau variables
    vector<int> tau_id;
    vector<double> tau_pt, tau_eta, tau_phi, tau_mass;

    // photon variables
    vector<double> pho_pt, pho_eta, pho_phi, photonCutBasedIDLoose;

    // Higgs candidate variables
    vector<double> H_pt; vector<double> H_eta; vector<double> H_phi; vector<double> H_mass;
    vector<double> H_noFSR_pt; vector<double> H_noFSR_eta; vector<double> H_noFSR_phi; vector<double> H_noFSR_mass;
    float mass4l, mass4l_noFSR, mass4e, mass4mu, mass2e2mu, pT4l, eta4l, phi4l, rapidity4l;
    float cosTheta1, cosTheta2, cosThetaStar, Phi, Phi1;
    float mass3l;

    // kin fit
    float mass4lREFIT, massZ1REFIT, massZ2REFIT, mass4lErr, mass4lErrREFIT;

    // Z candidate variables
    vector<double> Z_pt; vector<double> Z_eta; vector<double> Z_phi; vector<double> Z_mass;
    vector<double> Z_noFSR_pt; vector<double> Z_noFSR_eta; vector<double> Z_noFSR_phi; vector<double> Z_noFSR_mass;
    int Z_Hindex[2]; // position of Z1 and Z2 in Z_p4
    float massZ1, massZ1_Z1L, massZ2, pTZ1, pTZ2;

    // MET
    float met; float met_phi;
    float met_jesup, met_phi_jesup, met_jesdn, met_phi_jesdn;
    float met_uncenup, met_phi_uncenup, met_uncendn, met_phi_uncendn;

    // Jets
    vector<int>    jet_iscleanH4l;
    int jet1index, jet2index;
    vector<double> jet_pt; vector<double> jet_eta; vector<double> jet_phi; vector<double> jet_mass; vector<double> jet_pt_raw;
    vector<float>  jet_pumva, jet_csvv2, jet_csvv2_; vector<int> jet_isbtag;
    vector<int>    jet_hadronFlavour, jet_partonFlavour;
    vector<float>  jet_QGTagger, jet_QGTagger_jesup, jet_QGTagger_jesdn;
    vector<float>  jet_axis2, jet_ptD; vector<int> jet_mult;
    vector<float>  jet_relpterr; vector<float>  jet_phierr;
    vector<float>  jet_bTagEffi;
    vector<float>  jet_cTagEffi;
    vector<float>  jet_udsgTagEffi;
    vector<int>    jet_jesup_iscleanH4l;
    vector<double> jet_jesup_pt; vector<double> jet_jesup_eta;
    vector<double> jet_jesup_phi; vector<double> jet_jesup_mass;
    vector<int>    jet_jesdn_iscleanH4l;
    vector<double> jet_jesdn_pt; vector<double> jet_jesdn_eta;
    vector<double> jet_jesdn_phi; vector<double> jet_jesdn_mass;
    vector<int>    jet_jerup_iscleanH4l;
    vector<double> jet_jerup_pt; vector<double> jet_jerup_eta;
    vector<double> jet_jerup_phi; vector<double> jet_jerup_mass;
    vector<int>    jet_jerdn_iscleanH4l;
    vector<double> jet_jerdn_pt; vector<double> jet_jerdn_eta;
    vector<double> jet_jerdn_phi; vector<double> jet_jerdn_mass;
    int njets_pt30_eta4p7; int njets_pt30_eta4p7_jesup; int njets_pt30_eta4p7_jesdn;
    int njets_pt30_eta4p7_jerup; int njets_pt30_eta4p7_jerdn;
    int njets_pt30_eta2p5; int njets_pt30_eta2p5_jesup; int njets_pt30_eta2p5_jesdn;
    int njets_pt30_eta2p5_jerup; int njets_pt30_eta2p5_jerdn;
    int nbjets_pt30_eta4p7; int nvjets_pt40_eta2p4;
    float pt_leadingjet_pt30_eta4p7;
    float pt_leadingjet_pt30_eta4p7_jesup; float pt_leadingjet_pt30_eta4p7_jesdn;
    float pt_leadingjet_pt30_eta4p7_jerup; float pt_leadingjet_pt30_eta4p7_jerdn;
    float pt_leadingjet_pt30_eta2p5;
    float pt_leadingjet_pt30_eta2p5_jesup; float pt_leadingjet_pt30_eta2p5_jesdn;
    float pt_leadingjet_pt30_eta2p5_jerup; float pt_leadingjet_pt30_eta2p5_jerdn;
    float absrapidity_leadingjet_pt30_eta4p7;
    float absrapidity_leadingjet_pt30_eta4p7_jesup; float absrapidity_leadingjet_pt30_eta4p7_jesdn;
    float absrapidity_leadingjet_pt30_eta4p7_jerup; float absrapidity_leadingjet_pt30_eta4p7_jerdn;
    float absdeltarapidity_hleadingjet_pt30_eta4p7;
    float absdeltarapidity_hleadingjet_pt30_eta4p7_jesup; float absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn;
    float absdeltarapidity_hleadingjet_pt30_eta4p7_jerup; float absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn;
    float DijetMass, DijetDEta, DijetFisher;

    // merged jets
    vector<int>   mergedjet_iscleanH4l;
    vector<float> mergedjet_pt; vector<float> mergedjet_eta; vector<float> mergedjet_phi; vector<float> mergedjet_mass;

    vector<float> mergedjet_tau1; vector<float> mergedjet_tau2;
    vector<float> mergedjet_btag;

    vector<float> mergedjet_L1;
    vector<float> mergedjet_prunedmass; vector<float> mergedjet_softdropmass;

    vector<int> mergedjet_nsubjet;
    vector<vector<float> > mergedjet_subjet_pt; vector<vector<float> > mergedjet_subjet_eta;
    vector<vector<float> > mergedjet_subjet_phi; vector<vector<float> > mergedjet_subjet_mass;
    vector<vector<float> > mergedjet_subjet_btag;
    vector<vector<int> > mergedjet_subjet_partonFlavour, mergedjet_subjet_hadronFlavour;

    // FSR Photons
    int nFSRPhotons;
    vector<int> fsrPhotons_lepindex;
    vector<double> fsrPhotons_pt; vector<double> fsrPhotons_pterr;
    vector<double> fsrPhotons_eta; vector<double> fsrPhotons_phi;
    vector<double> fsrPhotons_mass;
    vector<float> fsrPhotons_dR; vector<float> fsrPhotons_iso;
    vector<float> allfsrPhotons_dR; vector<float> allfsrPhotons_pt; vector<float> allfsrPhotons_iso;

    // Z4l? FIXME
    float theta12, theta13, theta14;
    float minM3l, Z4lmaxP, minDeltR, m3l_soft;
    float minMass2Lep, maxMass2Lep;
    float thetaPhoton, thetaPhotonZ;

    // Event Category
    int EventCat;

    // -------------------------
    // GEN level information
    // -------------------------

    //Event variables
    int GENfinalState;
    bool passedFiducialSelection;

    // lepton variables
    vector<double> GENlep_pt; vector<double> GENlep_eta; vector<double> GENlep_phi; vector<double> GENlep_mass;
    vector<int> GENlep_id; vector<int> GENlep_status;
    vector<int> GENlep_MomId; vector<int> GENlep_MomMomId;
    int GENlep_Hindex[4];//position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z3 sub
    vector<float> GENlep_isoCH; vector<float> GENlep_isoNH; vector<float> GENlep_isoPhot; vector<float> GENlep_RelIso;

    // Higgs candidate variables (calculated using selected gen leptons)
    vector<double> GENH_pt; vector<double> GENH_eta; vector<double> GENH_phi; vector<double> GENH_mass;
    float GENmass4l, GENmass4e, GENmass4mu, GENmass2e2mu, GENpT4l, GENeta4l, GENrapidity4l;
    float GENMH; //mass directly from gen particle with id==25
    float GENcosTheta1, GENcosTheta2, GENcosThetaStar, GENPhi, GENPhi1;

    // Z candidate variables
    vector<double> GENZ_pt; vector<double> GENZ_eta; vector<double> GENZ_phi; vector<double> GENZ_mass;
    vector<int> GENZ_DaughtersId; vector<int> GENZ_MomId;
    float  GENmassZ1, GENmassZ2, GENpTZ1, GENpTZ2, GENdPhiZZ, GENmassZZ, GENpTZZ;

    // Higgs variables directly from GEN particle
    float GENHmass;

    // Jets
    vector<double> GENjet_pt; vector<double> GENjet_eta; vector<double> GENjet_phi; vector<double> GENjet_mass;
    int GENnjets_pt30_eta4p7; float GENpt_leadingjet_pt30_eta4p7;
    int GENnjets_pt30_eta2p5; float GENpt_leadingjet_pt30_eta2p5;
    float GENabsrapidity_leadingjet_pt30_eta4p7; float GENabsdeltarapidity_hleadingjet_pt30_eta4p7;
    int lheNb, lheNj, nGenStatus2bHad;

    // STXS info
    int stage0cat;
    int stage1cat;
    int stage1p1cat;
    int stage1p2cat;
    // Fiducial Rivet
    int passedFiducialRivet;
    float GENpT4lRivet;
    int GENnjets_pt30_eta4p7Rivet;
    float GENpt_leadingjet_pt30_eta4p7Rivet;

    //KinZfitter
    KinZfitter *kinZfitter;

    // MEM
    Mela* mela;

    float me_0plus_JHU, me_qqZZ_MCFM, p0plus_m4l, bkg_m4l;
    float D_bkg_kin, D_bkg, D_g4, D_g1g4;

    float p0minus_VAJHU, Dgg10_VAMCFM, pg1g4_VAJHU;

    // old but working
    float phjj_VAJHU, pvbf_VAJHU;
    float pwh_hadronic_VAJHU, pzh_hadronic_VAJHU;
    float pAux_vbf_VAJHU, phj_VAJHU;

    float p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal;
    float p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal;
    float p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal;
    float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal;
    float p_HadWH_SIG_ghw1_1_JHUGen_JECNominal;
    float p_HadZH_SIG_ghz1_1_JHUGen_JECNominal;
    float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
    float p_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
    float pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
    float p_JQCD_SIG_ghv1_1_JHUGen_JECNominal;
    float p_JQCD_SIG_ghg2_1_JHUGen_JECNominal;

    float p_JJVBF_BKG_MCFM_JECNominal;
    float p_HadZH_BKG_MCFM_JECNominal;
    float p_HadWH_BKG_MCFM_JECNominal;
    float p_JJQCD_BKG_MCFM_JECNominal;

    float p_HadZH_mavjj_JECNominal;
    float p_HadZH_mavjj_true_JECNominal;
    float p_HadWH_mavjj_JECNominal;
    float p_HadWH_mavjj_true_JECNominal;

    float pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal;
    float pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal;
    float pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal;
    float pConst_JJVBF_BKG_MCFM_JECNominal;
    float pConst_HadZH_BKG_MCFM_JECNominal;
    float pConst_HadWH_BKG_MCFM_JECNominal;
    float pConst_JJQCD_BKG_MCFM_JECNominal;
    float pConst_JJVBF_SIG_ghv1_1_JHUGen_JECNominal;
    float pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal;
    float pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal;
    float pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;

    float D_VBF, D_VBF1j, D_HadWH, D_HadZH;
    float D_VBF_QG, D_VBF1j_QG, D_HadWH_QG, D_HadZH_QG;
    float D_bkg_VBFdec, D_bkg_VHdec;

    // a vector<float> for each vector<double>
    vector<float> lep_pt_float, lep_pterr_float, lep_pterrold_float;
    vector<float> lep_p_float, lep_ecalEnergy_float;
    vector<float> lep_eta_float, lep_phi_float, lep_mass_float;
    vector<float> lepFSR_pt_float, lepFSR_eta_float;
    vector<float> lepFSR_phi_float, lepFSR_mass_float;
    vector<float> tau_pt_float, tau_eta_float, tau_phi_float, tau_mass_float;
    vector<float> pho_pt_float, pho_eta_float, pho_phi_float, photonCutBasedIDLoose_float;
    vector<float> H_pt_float, H_eta_float, H_phi_float, H_mass_float;
    vector<float> H_noFSR_pt_float, H_noFSR_eta_float;
    vector<float> H_noFSR_phi_float, H_noFSR_mass_float;
    vector<float> Z_pt_float, Z_eta_float, Z_phi_float, Z_mass_float;
    vector<float> Z_noFSR_pt_float, Z_noFSR_eta_float;
    vector<float> Z_noFSR_phi_float, Z_noFSR_mass_float;
    vector<float> jet_pt_float, jet_eta_float, jet_phi_float, jet_mass_float, jet_pt_raw_float;
    vector<float> jet_jesup_pt_float, jet_jesup_eta_float;
    vector<float> jet_jesup_phi_float, jet_jesup_mass_float;
    vector<float> jet_jesdn_pt_float, jet_jesdn_eta_float;
    vector<float> jet_jesdn_phi_float, jet_jesdn_mass_float;
    vector<float> jet_jerup_pt_float, jet_jerup_eta_float;
    vector<float> jet_jerup_phi_float, jet_jerup_mass_float;
    vector<float> jet_jerdn_pt_float, jet_jerdn_eta_float;
    vector<float> jet_jerdn_phi_float, jet_jerdn_mass_float;
    vector<float> fsrPhotons_pt_float, fsrPhotons_pterr_float;
    vector<float> fsrPhotons_eta_float, fsrPhotons_phi_float, fsrPhotons_mass_float;
    vector<float> GENlep_pt_float, GENlep_eta_float;
    vector<float> GENlep_phi_float, GENlep_mass_float;
    vector<float> GENH_pt_float, GENH_eta_float;
    vector<float> GENH_phi_float, GENH_mass_float;
    vector<float> GENZ_pt_float, GENZ_eta_float;
    vector<float> GENZ_phi_float, GENZ_mass_float;
    vector<float> GENjet_pt_float, GENjet_eta_float;
    vector<float> GENjet_phi_float, GENjet_mass_float;

    // Global Variables but not stored in the tree
    vector<double> lep_ptreco;
    vector<int> lep_ptid; vector<int> lep_ptindex;
    vector<pat::Muon> recoMuons; vector<pat::Electron> recoElectrons; vector<pat::Electron> recoElectronsUnS;
    vector<pat::Tau> recoTaus; vector<pat::Photon> recoPhotons;
    vector<pat::PFParticle> fsrPhotons;
    TLorentzVector HVec, HVecNoFSR, Z1Vec, Z2Vec;
    TLorentzVector GENZ1Vec, GENZ2Vec;
    bool RecoFourMuEvent, RecoFourEEvent, RecoTwoETwoMuEvent, RecoTwoMuTwoEEvent;
    bool foundHiggsCandidate; bool foundZ1LCandidate; bool firstEntry;
    float jet1pt, jet2pt;

    // hist container
    std::map<std::string,TH1F*> histContainer_;

    //Input edm
    edm::EDGetTokenT<edm::View<pat::Electron> > elecSrc_;
    edm::EDGetTokenT<edm::View<pat::Electron> > elecUnSSrc_;
    edm::EDGetTokenT<edm::View<pat::Muon> > muonSrc_;
    edm::EDGetTokenT<edm::View<pat::Tau> > tauSrc_;
    edm::EDGetTokenT<edm::View<pat::Photon> > photonSrc_;
    edm::EDGetTokenT<edm::View<pat::Jet> > jetSrc_;
    edm::EDGetTokenT<edm::ValueMap<float> > qgTagSrc_;
    edm::EDGetTokenT<edm::ValueMap<float> > axis2Src_;
    edm::EDGetTokenT<edm::ValueMap<int> > multSrc_;
    edm::EDGetTokenT<edm::ValueMap<float> > ptDSrc_;
    edm::EDGetTokenT<edm::View<pat::Jet> > mergedjetSrc_;
    edm::EDGetTokenT<edm::View<pat::MET> > metSrc_;
    //edm::InputTag triggerSrc_;
    edm::EDGetTokenT<edm::TriggerResults> triggerSrc_;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
    edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
    edm::EDGetTokenT<std::vector<reco::Conversion> > conversionSrc_;
    edm::EDGetTokenT<double> muRhoSrc_;
    edm::EDGetTokenT<double> elRhoSrc_;
    edm::EDGetTokenT<double> rhoSrcSUS_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSrc_;
    edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandsSrc_;
    edm::EDGetTokenT<edm::View<pat::PFParticle> > fsrPhotonsSrc_;
    edm::EDGetTokenT<reco::GenParticleCollection> prunedgenParticlesSrc_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedgenParticlesSrc_;
    edm::EDGetTokenT<edm::View<reco::GenJet> > genJetsSrc_;
    edm::EDGetTokenT<GenEventInfoProduct> generatorSrc_;
    edm::EDGetTokenT<LHEEventProduct> lheInfoSrc_;
    edm::EDGetTokenT<LHERunInfoProduct> lheRunInfoToken_;
    edm::EDGetTokenT<HTXS::HiggsClassification> htxsSrc_;
    //edm::EDGetTokenT<HZZFid::FiducialSummary> fidRivetSrc_;
    edm::EDGetTokenT< double > prefweight_token_;

    // Configuration
    const float Zmass;
    float mZ1Low, mZ2Low, mZ1High, mZ2High, m4lLowCut;
    float jetpt_cut, jeteta_cut;
    std::string elecID;
    bool isMC, isSignal;
    float mH;
    float crossSection;
    bool weightEvents;
    float isoCutEl, isoCutMu;
    double isoConeSizeEl, isoConeSizeMu;
    float sip3dCut, leadingPtCut, subleadingPtCut;
    float genIsoCutEl, genIsoCutMu;
    double genIsoConeSizeEl, genIsoConeSizeMu;
    float _elecPtCut, _muPtCut, _tauPtCut, _phoPtCut;
    float BTagCut;
    bool reweightForPU;
    std::string PUVersion;
    bool doFsrRecovery,bestCandMela, doMela, GENbestM4l;
    bool doPUJetID;
    int jetIDLevel;
    bool doJER;
    bool doJEC;
    bool doRefit;
    bool doTriggerMatching;
    bool checkOnlySingle;
    std::vector<std::string> triggerList;
    int skimLooseLeptons, skimTightLeptons;
    bool verbose;

    int year;///use to choose Muon BDT

    // register to the TFileService
    edm::Service<TFileService> fs;

    // Counters
    float nEventsTotal;
    float sumWeightsTotal;
    float sumWeightsTotalPU;

    // JER
    JME::JetResolution resolution_pt, resolution_phi;
    JME::JetResolutionScaleFactor resolution_sf;

    string EleBDT_name_161718;
    string heepID_name_161718;

};

UFHZZ4LAna::UFHZZ4LAna(const edm::ParameterSet& iConfig):
    histContainer_(),
    elecSrc_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc"))),
    elecUnSSrc_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("electronUnSSrc"))),
    muonSrc_(consumes<edm::View<pat::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc"))),
    tauSrc_(consumes<edm::View<pat::Tau> >(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc"))),
    photonSrc_(consumes<edm::View<pat::Photon> >(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc"))),
    jetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc"))),
    qgTagSrc_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"))),
    axis2Src_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "axis2"))),
    multSrc_(consumes<edm::ValueMap<int>>(edm::InputTag("QGTagger", "mult"))),
    ptDSrc_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "ptD"))),
     mergedjetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("mergedjetSrc"))),
    metSrc_(consumes<edm::View<pat::MET> >(iConfig.getUntrackedParameter<edm::InputTag>("metSrc"))),
    triggerSrc_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerSrc"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
    vertexSrc_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc"))),
    beamSpotSrc_(consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotSrc"))),
    conversionSrc_(consumes<std::vector<reco::Conversion> >(iConfig.getUntrackedParameter<edm::InputTag>("conversionSrc"))),
    muRhoSrc_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("muRhoSrc"))),
    elRhoSrc_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("elRhoSrc"))),
    rhoSrcSUS_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("rhoSrcSUS"))),
    pileupSrc_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getUntrackedParameter<edm::InputTag>("pileupSrc"))),
    pfCandsSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfCandsSrc"))),
    fsrPhotonsSrc_(consumes<edm::View<pat::PFParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("fsrPhotonsSrc"))),
    prunedgenParticlesSrc_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("prunedgenParticlesSrc"))),
    packedgenParticlesSrc_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("packedgenParticlesSrc"))),
    genJetsSrc_(consumes<edm::View<reco::GenJet> >(iConfig.getUntrackedParameter<edm::InputTag>("genJetsSrc"))),
    generatorSrc_(consumes<GenEventInfoProduct>(iConfig.getUntrackedParameter<edm::InputTag>("generatorSrc"))),
    lheInfoSrc_(consumes<LHEEventProduct>(iConfig.getUntrackedParameter<edm::InputTag>("lheInfoSrc"))),
    lheRunInfoToken_(consumes<LHERunInfoProduct,edm::InRun>(edm::InputTag("externalLHEProducer",""))),
    htxsSrc_(consumes<HTXS::HiggsClassification>(edm::InputTag("rivetProducerHTXS","HiggsClassification"))),
    prefweight_token_(consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"))),//2017 2017 need this. not for 2018
    //fidRivetSrc_(consumes<HZZFid::FiducialSummary>(edm::InputTag("rivetProducerHZZFid","FiducialSummary"))),
    Zmass(91.1876),
    mZ1Low(iConfig.getUntrackedParameter<double>("mZ1Low",40.0)),
    mZ2Low(iConfig.getUntrackedParameter<double>("mZ2Low",12.0)), // was 12
    mZ1High(iConfig.getUntrackedParameter<double>("mZ1High",120.0)),
    mZ2High(iConfig.getUntrackedParameter<double>("mZ2High",120.0)),
    m4lLowCut(iConfig.getUntrackedParameter<double>("m4lLowCut",70.0)),
    jetpt_cut(iConfig.getUntrackedParameter<double>("jetpt_cut",10.0)),
    jeteta_cut(iConfig.getUntrackedParameter<double>("eta_cut",4.7)),
    elecID(iConfig.getUntrackedParameter<std::string>("elecID","NonTrig")),
    isMC(iConfig.getUntrackedParameter<bool>("isMC",true)),
    isSignal(iConfig.getUntrackedParameter<bool>("isSignal",false)),
    mH(iConfig.getUntrackedParameter<double>("mH",0.0)),
    crossSection(iConfig.getUntrackedParameter<double>("CrossSection",1.0)),
    weightEvents(iConfig.getUntrackedParameter<bool>("weightEvents",false)),
    isoCutEl(iConfig.getUntrackedParameter<double>("isoCutEl",9999.0)),
    isoCutMu(iConfig.getUntrackedParameter<double>("isoCutMu",0.35)),/////ios is applied to new Muon BDT //previous 0.35///Qianying
    //isoCutMu(iConfig.getUntrackedParameter<double>("isoCutMu",9999.0)),/////ios is applied to new Muon BDT //previous 0.35///Qianying
    isoConeSizeEl(iConfig.getUntrackedParameter<double>("isoConeSizeEl",0.3)),
    isoConeSizeMu(iConfig.getUntrackedParameter<double>("isoConeSizeMu",0.3)),
    sip3dCut(iConfig.getUntrackedParameter<double>("sip3dCut",4)),
    leadingPtCut(iConfig.getUntrackedParameter<double>("leadingPtCut",20.0)),
    subleadingPtCut(iConfig.getUntrackedParameter<double>("subleadingPtCut",10.0)),
    genIsoCutEl(iConfig.getUntrackedParameter<double>("genIsoCutEl",0.35)),
    genIsoCutMu(iConfig.getUntrackedParameter<double>("genIsoCutMu",0.35)),
    genIsoConeSizeEl(iConfig.getUntrackedParameter<double>("genIsoConeSizeEl",0.3)),
    genIsoConeSizeMu(iConfig.getUntrackedParameter<double>("genIsoConeSizeMu",0.3)),
    _elecPtCut(iConfig.getUntrackedParameter<double>("_elecPtCut",7.0)),
    _muPtCut(iConfig.getUntrackedParameter<double>("_muPtCut",5.0)),
    _tauPtCut(iConfig.getUntrackedParameter<double>("_tauPtCut",20.0)),
    _phoPtCut(iConfig.getUntrackedParameter<double>("_phoPtCut",10.0)),
    //BTagCut(iConfig.getUntrackedParameter<double>("BTagCut",0.4184)),/////2016: 0.6321; 2017: 0.4941; 2018: 0.4184
    reweightForPU(iConfig.getUntrackedParameter<bool>("reweightForPU",true)),
    PUVersion(iConfig.getUntrackedParameter<std::string>("PUVersion","Summer16_80X")),
    doFsrRecovery(iConfig.getUntrackedParameter<bool>("doFsrRecovery",true)),
    bestCandMela(iConfig.getUntrackedParameter<bool>("bestCandMela",true)),
    doMela(iConfig.getUntrackedParameter<bool>("doMela",true)),
    GENbestM4l(iConfig.getUntrackedParameter<bool>("GENbestM4l",false)),
    doPUJetID(iConfig.getUntrackedParameter<bool>("doPUJetID",true)),
    jetIDLevel(iConfig.getUntrackedParameter<int>("jetIDLevel",2)),
    doJER(iConfig.getUntrackedParameter<bool>("doJER",true)),
    doJEC(iConfig.getUntrackedParameter<bool>("doJEC",true)),
    doRefit(iConfig.getUntrackedParameter<bool>("doRefit",true)),
    doTriggerMatching(iConfig.getUntrackedParameter<bool>("doTriggerMatching",!isMC)),
    checkOnlySingle(iConfig.getUntrackedParameter<bool>("checkOnlySingle",false)),
    triggerList(iConfig.getUntrackedParameter<std::vector<std::string>>("triggerList")),
    skimLooseLeptons(iConfig.getUntrackedParameter<int>("skimLooseLeptons",2)),
    skimTightLeptons(iConfig.getUntrackedParameter<int>("skimTightLeptons",2)),
    verbose(iConfig.getUntrackedParameter<bool>("verbose",false)),
    year(iConfig.getUntrackedParameter<int>("year",2018))////for year put 2016,2017, or 2018 to select correct training
{

        if(!isMC){reweightForPU = false;}

        nEventsTotal=0.0;
        sumWeightsTotal=0.0;
        sumWeightsTotalPU=0.0;
        histContainer_["NEVENTS"]=fs->make<TH1F>("nEvents","nEvents in Sample",2,0,2);
        histContainer_["SUMWEIGHTS"]=fs->make<TH1F>("sumWeights","sum Weights of Sample",2,0,2);
        histContainer_["SUMWEIGHTSPU"]=fs->make<TH1F>("sumWeightsPU","sum Weights and PU of Sample",2,0,2);
        histContainer_["NVTX"]=fs->make<TH1F>("nVtx","Number of Vertices",36,-0.5,35.5);
        histContainer_["NVTX_RW"]=fs->make<TH1F>("nVtx_ReWeighted","Number of Vertices",36,-0.5,35.5);
        histContainer_["NINTERACT"]=fs->make<TH1F>("nInteractions","Number of True Interactions",61,-0.5,60.5);
        histContainer_["NINTERACT_RW"]=fs->make<TH1F>("nInteraction_ReWeighted","Number of True Interactions",61,-0.5,60.5);

        passedEventsTree_All = new TTree("passedEvents","passedEvents");

        edm::FileInPath kfacfileInPath("UFHZZAnalysisRun2/UFHZZ4LAna/data/Kfactor_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
        TFile *fKFactor = TFile::Open(kfacfileInPath.fullPath().c_str());
        kFactor_ggzz = (TSpline3*) fKFactor->Get("sp_Kfactor");
        fKFactor->Close();
        delete fKFactor;

        tableEwk = readFile_and_loadEwkTable("ZZBG");

        kinZfitter = new KinZfitter(!isMC);

       if(doMela){
        mela = new Mela(13.0, 125.0, TVar::SILENT);
        mela->setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    	}

        //string elec_scalefac_Cracks_name_161718[3] = {"egammaEffi.txt_EGM2D_cracks.root", "egammaEffi.txt_EGM2D_Moriond2018v1_gap.root", "egammaEffi.txt_EGM2D_Moriond2019_v1_gap.root"};
        string elec_scalefac_Cracks_name_161718[3] = {"ElectronSF_Legacy_2016_Gap.root", "ElectronSF_Legacy_2017_Gap.root", "ElectronSF_Legacy_2018_Gap.root"};
        edm::FileInPath elec_scalefacFileInPathCracks(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+elec_scalefac_Cracks_name_161718[year-2016]).c_str());
        TFile *fElecScalFacCracks = TFile::Open(elec_scalefacFileInPathCracks.fullPath().c_str());
        hElecScaleFac_Cracks = (TH2F*)fElecScalFacCracks->Get("EGamma_SF2D");

        //string elec_scalefac_name_161718[3] = {"egammaEffi.txt_EGM2D.root", "egammaEffi.txt_EGM2D_Moriond2018v1.root", "egammaEffi.txt_EGM2D_Moriond2019_v1.root"};
        string elec_scalefac_name_161718[3] = {"ElectronSF_Legacy_2016_NoGap.root", "ElectronSF_Legacy_2017_NoGap.root", "ElectronSF_Legacy_2018_NoGap.root"};
        edm::FileInPath elec_scalefacFileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+elec_scalefac_name_161718[year-2016]).c_str());
        TFile *fElecScalFac = TFile::Open(elec_scalefacFileInPath.fullPath().c_str());
        hElecScaleFac = (TH2F*)fElecScalFac->Get("EGamma_SF2D");

        //string elec_Gsfscalefac_name_161718[3] = {"egammaEffi.txt_EGM2D_GSF.root", "egammaEffi.txt_EGM2D_Moriond2018v1_runBCDEF_passingRECO.root", "Ele_Reco_2018.root"};//was previous;
        string elec_Gsfscalefac_name_161718[3] = {"Ele_Reco_2016.root", "Ele_Reco_2017.root", "Ele_Reco_2018.root"};
        edm::FileInPath elec_GsfscalefacFileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+elec_Gsfscalefac_name_161718[year-2016]).c_str());
        TFile *fElecScalFacGsf = TFile::Open(elec_GsfscalefacFileInPath.fullPath().c_str());
        hElecScaleFacGsf = (TH2F*)fElecScalFacGsf->Get("EGamma_SF2D");

        //string elec_GsfLowETscalefac_name_161718[3]= {"", "egammaEffi.txt_EGM2D_Moriond2018v1_runBCDEF_passingRECO_lowEt.root", "Ele_Reco_LowEt_2018.root"};//was previous
        string elec_GsfLowETscalefac_name_161718[3]= {"Ele_Reco_LowEt_2016.root", "Ele_Reco_LowEt_2017.root", "Ele_Reco_LowEt_2018.root"};
        edm::FileInPath elec_GsfLowETscalefacFileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+elec_GsfLowETscalefac_name_161718[year-2016]).c_str());
        TFile *fElecScalFacGsfLowET = TFile::Open(elec_GsfLowETscalefacFileInPath.fullPath().c_str());
        hElecScaleFacGsfLowET = (TH2F*)fElecScalFacGsfLowET->Get("EGamma_SF2D");

        //string mu_scalefac_name_161718[3] = {"final_HZZ_Moriond17Preliminary_v4.root", "ScaleFactors_mu_Moriond2018_final.root", "final_HZZ_muon_SF_2018RunA2D_ER_2702.root"};//was previous;
        string mu_scalefac_name_161718[3] = {"final_HZZ_SF_2016_legacy_mupogsysts.root", "final_HZZ_SF_2017_rereco_mupogsysts_3010.root", "final_HZZ_SF_2018_rereco_mupogsysts_3010.root"};
        edm::FileInPath mu_scalefacFileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+mu_scalefac_name_161718[year-2016]).c_str());
        TFile *fMuScalFac = TFile::Open(mu_scalefacFileInPath.fullPath().c_str());
        hMuScaleFac = (TH2F*)fMuScalFac->Get("FINAL");
        hMuScaleFacUnc = (TH2F*)fMuScalFac->Get("ERROR");

        //string pileup_name_161718[3] = {"puWeightsMoriond17_v2.root", "puWeightsMoriond18.root", "pu_weights_2018.root"};///was previous
        string pileup_name_161718[3] = {"pu_weights_2016.root", "pu_weights_2017.root", "pu_weights_2018.root"};
        edm::FileInPath pileup_FileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+pileup_name_161718[year-2016]).c_str());
        TFile *f_pileup = TFile::Open(pileup_FileInPath.fullPath().c_str());
        h_pileup = (TH1D*)f_pileup->Get("weights");
        h_pileupUp = (TH1D*)f_pileup->Get("weights_varUp");
        h_pileupDn = (TH1D*)f_pileup->Get("weights_varDn");


        string bTagEffi_name_161718[3] = {"bTagEfficiencies_2016.root", "bTagEfficiencies_2017.root", "bTagEfficiencies_2018.root"};
        edm::FileInPath BTagEffiInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+bTagEffi_name_161718[year-2016]).c_str());
        TFile *fbTagEffi = TFile::Open(BTagEffiInPath.fullPath().c_str());
        hbTagEffi = (TH2F*)fbTagEffi->Get("eff_b_M_ALL");
        hcTagEffi = (TH2F*)fbTagEffi->Get("eff_c_M_ALL");
        hudsgTagEffi = (TH2F*)fbTagEffi->Get("eff_udsg_M_ALL");


        //BTag calibration
        string csv_name_161718[3] = {"DeepCSV_2016LegacySF_V1.csv", "DeepCSV_94XSF_V4_B_F.csv", "DeepCSV_102XSF_V1.csv"};
        edm::FileInPath btagfileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+csv_name_161718[year-2016]).c_str());

        BTagCalibration calib("DeepCSV", btagfileInPath.fullPath().c_str());
        reader = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,  // operating point
                                           "central",             // central sys type
                                           {"up", "down"});      // other sys types


        reader->load(calib,                // calibration instance
                    BTagEntry::FLAV_B,    // btag flavour
                    "comb");               // measurement type

        if(year==2018)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Autumn18IdIsoValues"; BTagCut=0.4184; heepID_name_161718 = "heepElectronID-HEEPV70";}
        if(year==2017)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Fall17IsoV2Values"; BTagCut=0.4941; heepID_name_161718 = "heepElectronID-HEEPV70";}
        if(year==2016)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Summer16IdIsoValues"; BTagCut=0.6321; heepID_name_161718 = "heepElectronID-HEEPV70";}

}

UFHZZ4LAna::~UFHZZ4LAna()
{
    //destructor --- don't do anything here
}

// ------------ method called for each event  ------------
void UFHZZ4LAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace pat;
  using namespace trigger;
  using namespace EwkCorrections;

  nEventsTotal += 1.0;

  Run=iEvent.id().run();
  Event = iEvent.id().event();
  LumiSect = iEvent.id().luminosityBlock();

  if(verbose){cout<<"Run: " << Run << ",Event: " << Event << ",LumiSect: "<<LumiSect<<endl;}

  // ======= Get Collections ======= //
  if (verbose) {cout<<"getting collections"<<endl;}

  // trigger collection
  edm::Handle<edm::TriggerResults> trigger;
  iEvent.getByToken(triggerSrc_,trigger);
  const edm::TriggerNames trigNames = iEvent.triggerNames(*trigger);

  // trigger Objects
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  // vertex collection
  edm::Handle<reco::VertexCollection> vertex;
  iEvent.getByToken(vertexSrc_,vertex);

  // electron collection
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByToken(elecSrc_,electrons);
  if (verbose) cout<<electrons->size()<<" total electrons in the collection"<<endl;

  // electron before scale/smearing corrections
  edm::Handle<edm::View<pat::Electron> > electronsUnS;
  iEvent.getByToken(elecUnSSrc_,electronsUnS);

  // muon collection
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muonSrc_,muons);
  if (verbose) cout<<muons->size()<<" total muons in the collection"<<endl;

  // tau collection
  edm::Handle<edm::View<pat::Tau> > taus;
  iEvent.getByToken(tauSrc_,taus);
  if (verbose) cout<<taus->size()<<" total taus in the collection"<<endl;

  // photon collection
  edm::Handle<edm::View<pat::Photon> > photons;
  iEvent.getByToken(photonSrc_,photons);
  if (verbose) cout<<photons->size()<<" total photons in the collection"<<endl;

  // met collection
  edm::Handle<edm::View<pat::MET> > mets;
  iEvent.getByToken(metSrc_,mets);

  // Rho Correction
  edm::Handle<double> eventRhoMu;
  iEvent.getByToken(muRhoSrc_,eventRhoMu);
  muRho = *eventRhoMu;

  edm::Handle<double> eventRhoE;
  iEvent.getByToken(elRhoSrc_,eventRhoE);
  elRho = *eventRhoE;

  edm::Handle<double> eventRhoSUS;
  iEvent.getByToken(rhoSrcSUS_,eventRhoSUS);
  rhoSUS = *eventRhoSUS;

  // Conversions
  edm::Handle< std::vector<reco::Conversion> > theConversions;
  iEvent.getByToken(conversionSrc_, theConversions);

  // Beam Spot
  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken(beamSpotSrc_,beamSpot);
  const reco::BeamSpot BS = *beamSpot;

  // Particle Flow Cands
  edm::Handle<pat::PackedCandidateCollection> pfCands;
  iEvent.getByToken(pfCandsSrc_,pfCands);

  // FSR Photons
  edm::Handle<edm::View<pat::PFParticle> > photonsForFsr;
  iEvent.getByToken(fsrPhotonsSrc_,photonsForFsr);

  // Jets
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetSrc_,jets);

  if (!jecunc) {
      edm::ESHandle<JetCorrectorParametersCollection> jetCorrParameterSet;
      iSetup.get<JetCorrectionsRecord>().get("AK4PFchs", jetCorrParameterSet);
      const JetCorrectorParameters& jetCorrParameters = (*jetCorrParameterSet)["Uncertainty"];
      jecunc.reset(new JetCorrectionUncertainty(jetCorrParameters));
  }

  resolution_pt = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
  resolution_phi = JME::JetResolution::get(iSetup, "AK4PFchs_phi");
  resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");

  edm::Handle<edm::ValueMap<float> > qgHandle;
  iEvent.getByToken(qgTagSrc_, qgHandle);

  edm::Handle<edm::ValueMap<float> > axis2Handle;
  iEvent.getByToken(axis2Src_, axis2Handle);

  edm::Handle<edm::ValueMap<int> > multHandle;
  iEvent.getByToken(multSrc_, multHandle);

  edm::Handle<edm::ValueMap<float> > ptDHandle;
  iEvent.getByToken(ptDSrc_, ptDHandle);

  edm::Handle<edm::View<pat::Jet> > mergedjets;
  iEvent.getByToken(mergedjetSrc_,mergedjets;

  // GEN collections
  edm::Handle<reco::GenParticleCollection> prunedgenParticles;
  iEvent.getByToken(prunedgenParticlesSrc_, prunedgenParticles);

  edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles;
  iEvent.getByToken(packedgenParticlesSrc_, packedgenParticles);

  edm::Handle<edm::View<reco::GenJet> > genJets;
  iEvent.getByToken(genJetsSrc_, genJets);

  edm::Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByToken(generatorSrc_,genEventInfo);

  edm::Handle<LHEEventProduct> lheInfo;
  iEvent.getByToken(lheInfoSrc_, lheInfo);

  if (isMC) {
      edm::Handle<HTXS::HiggsClassification> htxs;
      iEvent.getByToken(htxsSrc_,htxs);
      stage0cat = htxs->stage0_cat;
      stage1cat = htxs->stage1_cat_pTjet30GeV;
      //stage1p1cat = htxs->stage1p1_cat;
      stage1p1cat = htxs->stage1_1_cat_pTjet30GeV;
      stage1p2cat = htxs->stage1_2_cat_pTjet30GeV;
      if (verbose) cout<<"stage1cat "<<stage1cat<<endl;
  }

  if (isMC) {
      edm::Handle< double > theprefweight;
      iEvent.getByToken(prefweight_token_, theprefweight ) ;
      if (year == 2016 || year == 2017)
        prefiringWeight =(*theprefweight);
      else if (year == 2018)
        prefiringWeight =1.0;
    }
    else
      prefiringWeight =1.0;

      // ============ Initialize Variables ============= //

      // Event Variables
      if (verbose) {cout<<"clear variables"<<endl;}
      nVtx = -1.0; nInt = -1.0;
      finalState = -1;
      triggersPassed="";
      passedTrig=false; passedFullSelection=false; passedZ4lSelection=false; passedQCDcut=false;
      passedZ1LSelection=false; passedZ4lZ1LSelection=false; passedZ4lZXCRSelection=false; passedZXCRSelection=false;
      nZXCRFailedLeptons=0;

      // Event Weights
      genWeight=1.0; pileupWeight=1.0; pileupWeightUp=1.0; pileupWeightDn=1.0; dataMCWeight=1.0; eventWeight=1.0;
      k_ggZZ=1.0; k_qqZZ_qcd_dPhi = 1.0; k_qqZZ_qcd_M = 1.0; k_qqZZ_qcd_Pt = 1.0; k_qqZZ_ewk = 1.0;

      qcdWeights.clear(); nnloWeights.clear(); pdfWeights.clear();
      pdfRMSup=1.0; pdfRMSdown=1.0; pdfENVup=1.0; pdfENVdown=1.0;

      //lepton variables
      lep_pt.clear(); lep_pterr.clear(); lep_pterrold.clear();
      lep_p.clear(); lep_ecalEnergy.clear(); lep_isEB.clear(); lep_isEE.clear();
      lep_eta.clear(); lep_phi.clear(); lep_mass.clear();
      lepFSR_pt.clear(); lepFSR_eta.clear(); lepFSR_phi.clear(); lepFSR_mass.clear();
      for (int i=0; i<4; ++i) {lep_Hindex[i]=-1;}
      pTL1=-1.0; pTL2=-1.0; pTL3=-1.0; pTL4=-1.0;
      etaL1=9999.0; etaL2=9999.0; etaL3=9999.0; etaL4=9999.0;
      idL1=9999; idL2=9999; idL3=9999; idL4=9999;
      mL1=-1.0; mL2=-1.0; mL3=-1.0; mL4=-1.0;
      pTErrL1=-1.0; pTErrL2=-1.0; pTErrL3=-1.0; pTErrL4=-1.0;
      phiL1=9999.0; phiL2=9999.0; phiL3=9999.0; phiL4=9999.0;
      pTL1FSR=-1.0; pTL2FSR=-1.0; pTL3FSR=-1.0; pTL4FSR=-1.0;
      lep_genindex.clear(); lep_id.clear(); lep_dataMC.clear(); lep_dataMCErr.clear();
      lep_matchedR03_PdgId.clear(); lep_matchedR03_MomId.clear(); lep_matchedR03_MomMomId.clear();
      lep_mva.clear(); lep_ecalDriven.clear();
      lep_tightId.clear(); lep_tightIdSUS.clear(); lep_tightIdHiPt.clear(); //lep_tightId_old.clear();
      lep_Sip.clear(); lep_IP.clear();
      lep_isoNH.clear(); lep_isoCH.clear(); lep_isoPhot.clear(); lep_isoPU.clear(); lep_isoPUcorr.clear();
      lep_RelIso.clear(); lep_RelIsoNoFSR.clear(); lep_MiniIso.clear();
      lep_ptRatio.clear(); lep_ptRel.clear();
      lep_missingHits.clear();
      lep_filtersMatched.clear();
      nisoleptons=0;

      //  tau variables
      tau_id.clear(); tau_pt.clear(); tau_eta.clear(); tau_phi.clear(); tau_mass.clear();

      // photon variables
      pho_pt.clear(); pho_eta.clear(); pho_phi.clear(); photonCutBasedIDLoose.clear();

      // Higgs candidate variables
      H_pt.clear(); H_eta.clear(); H_phi.clear(); H_mass.clear();
      H_noFSR_pt.clear(); H_noFSR_eta.clear(); H_noFSR_phi.clear(); H_noFSR_mass.clear();
      mass4l=-1.0; mass4l_noFSR=-1.0; mass4e=-1.0; mass4mu=-1.0; mass2e2mu=-1.0; pT4l=-1.0; eta4l=9999.0; phi4l=9999.0; rapidity4l=9999.0;
      cosTheta1=9999.0; cosTheta2=9999.0; cosThetaStar=9999.0; Phi=9999.0; Phi1=9999.0;
      mass3l=-1.0;

      // kin fitter
      mass4lREFIT = -999.0; massZ1REFIT = -999.0; massZ2REFIT = -999.0; mass4lErr = -999.0; mass4lErrREFIT = -999.0;

      // Z candidate variables
      Z_pt.clear(); Z_eta.clear(); Z_phi.clear(); Z_mass.clear();
      Z_noFSR_pt.clear(); Z_noFSR_eta.clear(); Z_noFSR_phi.clear(); Z_noFSR_mass.clear();
      for (int i=0; i<2; ++i) {Z_Hindex[i]=-1;}
      massZ1=-1.0; massZ1_Z1L=-1.0; massZ2=-1.0; pTZ1=-1.0; pTZ2=-1.0;

      // MET
      met=-1.0; met_phi=9999.0;
      met_jesup=-1.0; met_phi_jesup=9999.0; met_jesdn=-1.0; met_phi_jesdn=9999.0;
      met_uncenup=-1.0; met_phi_uncenup=9999.0; met_uncendn=-1.0; met_phi_uncendn=9999.0;

      // Jets
      jet_pt.clear(); jet_eta.clear(); jet_phi.clear(); jet_mass.clear(); jet_pt_raw.clear();
      jet_jesup_pt.clear(); jet_jesup_eta.clear(); jet_jesup_phi.clear(); jet_jesup_mass.clear();
      jet_jesdn_pt.clear(); jet_jesdn_eta.clear(); jet_jesdn_phi.clear(); jet_jesdn_mass.clear();
      jet_jerup_pt.clear(); jet_jerup_eta.clear(); jet_jerup_phi.clear(); jet_jerup_mass.clear();
      jet_jerdn_pt.clear(); jet_jerdn_eta.clear(); jet_jerdn_phi.clear(); jet_jerdn_mass.clear();
      jet_pumva.clear(); jet_csvv2.clear(); jet_isbtag.clear();
      jet_csvv2_.clear();
      jet_hadronFlavour.clear(); jet_partonFlavour.clear();
      jet_QGTagger.clear(); jet_QGTagger_jesup.clear(); jet_QGTagger_jesdn.clear();
      jet_relpterr.clear(); jet_phierr.clear();
      jet_bTagEffi.clear();
      jet_cTagEffi.clear();
      jet_udsgTagEffi.clear();
      jet_axis2.clear(); jet_ptD.clear(); jet_mult.clear();

      jet_iscleanH4l.clear();
      jet1index=-1; jet2index=-1;
      jet_jesup_iscleanH4l.clear(); jet_jesdn_iscleanH4l.clear();
      jet_jerup_iscleanH4l.clear(); jet_jerdn_iscleanH4l.clear();

      njets_pt30_eta4p7=0;
      njets_pt30_eta4p7_jesup=0; njets_pt30_eta4p7_jesdn=0;
      njets_pt30_eta4p7_jerup=0; njets_pt30_eta4p7_jerdn=0;

      njets_pt30_eta2p5=0;
      njets_pt30_eta2p5_jesup=0; njets_pt30_eta2p5_jesdn=0;
      njets_pt30_eta2p5_jerup=0; njets_pt30_eta2p5_jerdn=0;

      nbjets_pt30_eta4p7=0; nvjets_pt40_eta2p4=0;

      pt_leadingjet_pt30_eta4p7=-1.0;
      pt_leadingjet_pt30_eta4p7_jesup=-1.0; pt_leadingjet_pt30_eta4p7_jesdn=-1.0;
      pt_leadingjet_pt30_eta4p7_jerup=-1.0; pt_leadingjet_pt30_eta4p7_jerdn=-1.0;

      pt_leadingjet_pt30_eta2p5=-1.0;
      pt_leadingjet_pt30_eta2p5_jesup=-1.0; pt_leadingjet_pt30_eta2p5_jesdn=-1.0;
      pt_leadingjet_pt30_eta2p5_jerup=-1.0; pt_leadingjet_pt30_eta2p5_jerdn=-1.0;

      absrapidity_leadingjet_pt30_eta4p7=-1.0;
      absrapidity_leadingjet_pt30_eta4p7_jesup=-1.0; absrapidity_leadingjet_pt30_eta4p7_jesdn=-1.0;
      absrapidity_leadingjet_pt30_eta4p7_jerup=-1.0; absrapidity_leadingjet_pt30_eta4p7_jerdn=-1.0;

      absdeltarapidity_hleadingjet_pt30_eta4p7=-1.0;
      absdeltarapidity_hleadingjet_pt30_eta4p7_jesup=-1.0; absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn=-1.0;
      absdeltarapidity_hleadingjet_pt30_eta4p7_jerup=-1.0; absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn=-1.0;

      DijetMass=-1.0; DijetDEta=9999.0; DijetFisher=9999.0;

      mergedjet_iscleanH4l.clear();
      mergedjet_pt.clear(); mergedjet_eta.clear(); mergedjet_phi.clear(); mergedjet_mass.clear();
      mergedjet_L1.clear();
      mergedjet_softdropmass.clear(); mergedjet_prunedmass.clear();
      mergedjet_tau1.clear(); mergedjet_tau2.clear();
      mergedjet_btag.clear();

      mergedjet_nsubjet.clear();
      mergedjet_subjet_pt.clear(); mergedjet_subjet_eta.clear();
      mergedjet_subjet_phi.clear(); mergedjet_subjet_mass.clear();
      mergedjet_subjet_btag.clear();
      mergedjet_subjet_partonFlavour.clear(); mergedjet_subjet_hadronFlavour.clear();

      // FSR Photons
      nFSRPhotons=0;
      fsrPhotons_lepindex.clear(); fsrPhotons_pt.clear(); fsrPhotons_pterr.clear();
      fsrPhotons_eta.clear(); fsrPhotons_phi.clear();
      fsrPhotons_dR.clear(); fsrPhotons_iso.clear();
      allfsrPhotons_dR.clear(); allfsrPhotons_pt.clear(); allfsrPhotons_iso.clear();

      // Z4l? FIXME
      theta12=9999.0; theta13=9999.0; theta14=9999.0;
      minM3l=-1.0; Z4lmaxP=-1.0; minDeltR=9999.0; m3l_soft=-1.0;
      minMass2Lep=-1.0; maxMass2Lep=-1.0;
      thetaPhoton=9999.0; thetaPhotonZ=9999.0;

      // -------------------------
      // GEN level information
      // -------------------------

      //Event variables
      GENfinalState=-1;
      passedFiducialSelection=false;

      // lepton variables
      GENlep_pt.clear(); GENlep_eta.clear(); GENlep_phi.clear(); GENlep_mass.clear();
      GENlep_id.clear(); GENlep_status.clear(); GENlep_MomId.clear(); GENlep_MomMomId.clear();
      for (int i=0; i<4; ++i) {GENlep_Hindex[i]=-1;};//position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z3 sub
      GENlep_isoCH.clear(); GENlep_isoNH.clear(); GENlep_isoPhot.clear(); GENlep_RelIso.clear();

      // Higgs candidate variables (calculated using selected gen leptons)
      GENH_pt.clear(); GENH_eta.clear(); GENH_phi.clear(); GENH_mass.clear();
      GENmass4l=-1.0; GENmassZ1=-1.0; GENmassZ2=-1.0; GENpT4l=-1.0; GENeta4l=9999.0; GENrapidity4l=9999.0; GENMH=-1.0;
      GENcosTheta1=9999.0; GENcosTheta2=9999.0; GENcosThetaStar=9999.0; GENPhi=9999.0; GENPhi1=9999.0;

      // Z candidate variables
      GENZ_DaughtersId.clear(); GENZ_MomId.clear();
      GENZ_pt.clear(); GENZ_eta.clear(); GENZ_phi.clear(); GENZ_mass.clear();
      GENmassZ1=-1.0; GENmassZ2=-1.0; GENpTZ1=-1.0; GENpTZ2=-1.0, GENdPhiZZ=9999.0, GENmassZZ=-1.0, GENpTZZ=-1.0;

      // Higgs variables directly from GEN particle
      GENHmass=-1.0;

      // Jets
      GENjet_pt.clear(); GENjet_eta.clear(); GENjet_phi.clear(); GENjet_mass.clear();
      GENnjets_pt30_eta4p7=0;
      GENnjets_pt30_eta2p5=0;
      GENpt_leadingjet_pt30_eta4p7=-1.0; GENabsrapidity_leadingjet_pt30_eta4p7=-1.0; GENabsdeltarapidity_hleadingjet_pt30_eta4p7=-1.0;
      GENpt_leadingjet_pt30_eta2p5=-1.0;
      lheNb=0; lheNj=0; nGenStatus2bHad=0;

      //ME
      me_0plus_JHU=999.0; me_qqZZ_MCFM=999.0; p0plus_m4l=999.0; bkg_m4l=999.0; D_bkg_kin=999.0; D_bkg=999.0;
      p0minus_VAJHU=999.0; pg1g4_VAJHU=999.0; Dgg10_VAMCFM=999.0, D_g4=999.0; D_g1g4=999.0;

      // OLD but working
      phjj_VAJHU=999.0; pvbf_VAJHU=999.0; pAux_vbf_VAJHU=999.0;
      pwh_hadronic_VAJHU=999.0; pwh_hadronic_VAJHU=999.0;
      pzh_hadronic_VAJHU=999.0; pzh_hadronic_VAJHU=999.0;

      p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal=999.0;
      p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal=999.0;
      p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal=999.0;
      p_JJVBF_BKG_MCFM_JECNominal=999.0;
      p_HadZH_BKG_MCFM_JECNominal=999.0;
      p_HadWH_BKG_MCFM_JECNominal=999.0;
      p_JJQCD_BKG_MCFM_JECNominal=999.0;
      p_HadZH_mavjj_JECNominal=999.0;
      p_HadZH_mavjj_true_JECNominal=999.0;
      p_HadWH_mavjj_JECNominal=999.0;
      p_HadWH_mavjj_true_JECNominal=999.0;
      pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal=999.0;
      pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal=999.0;
      pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal=999.0;
      pConst_JJVBF_BKG_MCFM_JECNominal=999.0;
      pConst_HadZH_BKG_MCFM_JECNominal=999.0;
      pConst_HadWH_BKG_MCFM_JECNominal=999.0;
      pConst_JJQCD_BKG_MCFM_JECNominal=999.0;

      p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal=999.0;
      pConst_JJVBF_SIG_ghv1_1_JHUGen_JECNominal=999.0;
      p_HadWH_SIG_ghw1_1_JHUGen_JECNominal=999.0;
      pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal=999.0;
      p_HadZH_SIG_ghz1_1_JHUGen_JECNominal=999.0;
      pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal=999.0;

      D_HadWH=999.0; D_HadZH=999.0;
      D_VBF=999.0; D_VBF1j=999.0; D_HadWH=999.0; D_HadZH=999.0;
      D_VBF_QG=999.0; D_VBF1j_QG=999.0; D_HadWH_QG=999.0; D_HadZH_QG=999.0;
      D_bkg_VBFdec=999.0; D_bkg_VHdec=999.0;

      if (verbose) {cout<<"clear other variables"<<endl; }
      // Resolution
      //massErrorUCSD=-1.0; massErrorUCSDCorr=-1.0; massErrorUF=-1.0; massErrorUFCorr=-1.0; massErrorUFADCorr=-1.0;

      // Event
      EventCat=-1;

      // Global variables not stored in tree
      lep_ptreco.clear(); lep_ptid.clear(); lep_ptindex.clear();
      recoMuons.clear(); recoElectrons.clear(); fsrPhotons.clear(); recoElectronsUnS.clear();
      HVec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
      HVecNoFSR.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
      Z1Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
      Z2Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
      GENZ1Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
      GENZ2Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
      RecoFourMuEvent = false; RecoFourEEvent = false;
      RecoTwoETwoMuEvent = false; RecoTwoMuTwoEEvent = false;
      foundHiggsCandidate = false; foundZ1LCandidate = false;
      jet1pt=-1.0; jet2pt=-1.0;

      // Float vectors
      lep_pt_float.clear(); lep_pterr_float.clear(); lep_pterrold_float.clear();
      lep_p_float.clear(); lep_ecalEnergy_float.clear();
      lep_eta_float.clear(); lep_phi_float.clear(); lep_mass_float.clear();
      lepFSR_pt_float.clear(); lepFSR_eta_float.clear(); lepFSR_phi_float.clear(); lepFSR_mass_float.clear();
      tau_pt_float.clear(); tau_eta_float.clear(); tau_phi_float.clear(); tau_mass_float.clear();
      pho_pt_float.clear(); pho_eta_float.clear(); pho_phi_float.clear(); photonCutBasedIDLoose_float.clear();
      H_pt_float.clear(); H_eta_float.clear(); H_phi_float.clear(); H_mass_float.clear();
      H_noFSR_pt_float.clear(); H_noFSR_eta_float.clear(); H_noFSR_phi_float.clear(); H_noFSR_mass_float.clear();
      Z_pt_float.clear(); Z_eta_float.clear(); Z_phi_float.clear(); Z_mass_float.clear();
      Z_noFSR_pt_float.clear(); Z_noFSR_eta_float.clear(); Z_noFSR_phi_float.clear(); Z_noFSR_mass_float.clear();
      jet_pt_float.clear(); jet_eta_float.clear(); jet_phi_float.clear(); jet_mass_float.clear(); jet_pt_raw_float.clear();
      jet_jesup_pt_float.clear(); jet_jesup_eta_float.clear(); jet_jesup_phi_float.clear(); jet_jesup_mass_float.clear();
      jet_jesdn_pt_float.clear(); jet_jesdn_eta_float.clear(); jet_jesdn_phi_float.clear(); jet_jesdn_mass_float.clear();
      jet_jerup_pt_float.clear(); jet_jerup_eta_float.clear(); jet_jerup_phi_float.clear(); jet_jerup_mass_float.clear();
      jet_jerdn_pt_float.clear(); jet_jerdn_eta_float.clear(); jet_jerdn_phi_float.clear();  jet_jerdn_mass_float.clear();
      fsrPhotons_pt_float.clear(); fsrPhotons_pterr_float.clear(); fsrPhotons_eta_float.clear(); fsrPhotons_phi_float.clear(); fsrPhotons_mass_float.clear();

      // ====================== Do Analysis ======================== //

      std::map<int, TLorentzVector> fsrmap;
      vector<reco::Candidate*> selectedLeptons;
      std::map<unsigned int, TLorentzVector> selectedFsrMap;

      fsrmap.clear(); selectedFsrMap.clear(); selectedLeptons.clear();

      if (verbose) cout<<"start pileup reweighting"<<endl;
      // PU information
      if(isMC && reweightForPU) {
        edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
        iEvent.getByToken(pileupSrc_, PupInfo);

        if (verbose) cout<<"got pileup info"<<endl;
        std::vector<PileupSummaryInfo>::const_iterator PVI;
        int npv = -1;
        for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
          int BX = PVI->getBunchCrossing();
          if(BX == 0) { npv = PVI->getTrueNumInteractions(); continue;}
        }
        if (verbose) cout<<"N true interations = "<<npv<<endl;
        nInt = npv;
        //pileupWeight = pileUp.getPUWeight(npv,PUVersion);
        pileupWeight = pileUp.getPUWeight(h_pileup,npv);
        pileupWeightUp = pileUp.getPUWeight(h_pileupUp,npv);
        pileupWeightDn = pileUp.getPUWeight(h_pileupDn,npv);
        if (verbose) cout<<"pileup weight = "<<pileupWeight<<", filling histograms"<<endl;
        histContainer_["NINTERACT"]->Fill(npv);
        histContainer_["NINTERACT_RW"]->Fill(npv,pileupWeight);
      } else { pileupWeight = 1.0;}

      if (verbose) {cout<<"finished pileup reweighting"<<endl; }

      if(isMC){
        float tmpWeight = genEventInfo->weight();
        genWeight = (tmpWeight > 0 ? 1.0 : -1.0);
        if (verbose) {cout<<"tmpWeight: "<<tmpWeight<<"; genWeight: "<<genWeight<<endl;}
        double rms = 0.0;

        //std::cout<<"tmpWeight: "<<tmpWeight<<std::endl;

        if(lheInfo.isValid()){
          for(unsigned int i = 0; i < lheInfo->weights().size(); i++){

            tmpWeight = genEventInfo->weight();
            tmpWeight *= lheInfo->weights()[i].wgt/lheInfo->originalXWGTUP();
            pdfWeights.push_back(tmpWeight);

            if (i<=8 or int(i)>=posNNPDF){
              tmpWeight = genEventInfo->weight();
              tmpWeight *= lheInfo->weights()[i].wgt/lheInfo->originalXWGTUP();
              if (int(i)<posNNPDF) {qcdWeights.push_back(tmpWeight);}
            }
            else{
              tmpWeight = lheInfo->weights()[i].wgt;
              tmpWeight /= lheInfo->originalXWGTUP();
              //if (i==9) genWeight = tmpWeight;
              if (int(i)<posNNPDF) {nnloWeights.push_back(tmpWeight);}
            }
            // NNPDF30 variations
            if (int(i)>=posNNPDF && int(i)<=(posNNPDF+100)){
              rms += tmpWeight*tmpWeight;
              if (tmpWeight>pdfENVup) pdfENVup=tmpWeight;
              if (tmpWeight<pdfENVdown) pdfENVdown=tmpWeight;
            }
          }
          pdfRMSup=sqrt(rms/100.0); pdfRMSdown=1.0/pdfRMSup;
          if (verbose) cout<<"pdfRMSup "<<pdfRMSup<<" pdfRMSdown "<<pdfRMSdown<<endl;

          const lhef::HEPEUP& lheEvent = lheInfo->hepeup();
          std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
          for ( size_t idxParticle = 0; idxParticle < lheParticles.size(); ++idxParticle ){
            int id = std::abs(lheEvent.IDUP[idxParticle]);
            int status = lheEvent.ISTUP[idxParticle];
            if ( status == 1 && id==5 ) {
              lheNb += 1;
            }
            if ( status == 1 && ((id >= 1 && id <= 6) || id == 21) ) {
              lheNj += 1;
            }
          }
        }

        if (verbose) cout<<"setting gen variables"<<endl;
        setGENVariables(prunedgenParticles,packedgenParticles,genJets);
        if (verbose) { cout<<"finshed setting gen variables"<<endl;  }

        if (int(GENZ_pt.size()) == 2){
          GENZ1Vec.SetPtEtaPhiM(GENZ_pt[0], GENZ_eta[0], GENZ_phi[0], GENZ_mass[0]);
          GENZ2Vec.SetPtEtaPhiM(GENZ_pt[1], GENZ_eta[1], GENZ_phi[1], GENZ_mass[1]);
          k_qqZZ_ewk = getEwkCorrections(prunedgenParticles, tableEwk, genEventInfo, GENZ1Vec, GENZ2Vec);
          k_ggZZ = kFactor_ggzz->Eval((GENZ1Vec+GENZ2Vec).M());
          if (verbose) cout<<"ZZmass: "<< (GENZ1Vec+GENZ2Vec).M() <<"k factor qqZZ ewk: "<<k_qqZZ_ewk<<" ggZZ qcd"<<k_ggZZ << endl;
        }
      }

      sumWeightsTotal += genWeight;
      sumWeightsTotalPU += pileupWeight*genWeight;

      eventWeight = pileupWeight*genWeight;
      unsigned int _tSize = trigger->size();
      // create a string with all passing trigger names
      for (unsigned int i=0; i<_tSize; ++i){
        std::string triggerName = trigNames.triggerName(i);
        if (strstr(triggerName.c_str(),"_step")) continue;
        if (strstr(triggerName.c_str(),"MC_")) continue;
        if (strstr(triggerName.c_str(),"AlCa_")) continue;
        if (strstr(triggerName.c_str(),"DST_")) continue;
        if (strstr(triggerName.c_str(),"HLT_HI")) continue;
        if (strstr(triggerName.c_str(),"HLT_Physics")) continue;
        if (strstr(triggerName.c_str(),"HLT_Random")) continue;
        if (strstr(triggerName.c_str(),"HLT_ZeroBias")) continue;
        if (strstr(triggerName.c_str(),"HLT_IsoTrack")) continue;
        if (strstr(triggerName.c_str(),"Hcal")) continue;
        if (strstr(triggerName.c_str(),"Ecal")) continue;
        if (trigger->accept(i)) triggersPassed += triggerName;
      }
      if (firstEntry) cout<<"triggersPassed: "<<triggersPassed<<endl;
      firstEntry = false;
      // check if any of the triggers in the user list have passed
      bool passedSingleEl=false;
      bool passedSingleMu=false;
      bool passedAnyOther=false;
      for (unsigned int i=0; i<triggerList.size(); ++i){
        if (strstr(triggersPassed.c_str(),triggerList.at(i).c_str())){
          passedTrig=true;
          if (!isMC){
            if (strstr(triggerList.at(i).c_str(),"_WP")) passedSingleEl=true;
            if (strstr(triggerList.at(i).c_str(),"HLT_Iso")) passedSingleMu=true;
            if (strstr(triggerList.at(i).c_str(),"CaloIdL")) passedAnyOther=true;
            if (strstr(triggerList.at(i).c_str(),"TrkIsoVVL")) passedAnyOther=true;
            if (strstr(triggerList.at(i).c_str(),"Triple")) passedAnyOther=true;
          }
        }
      }

      bool passedOnlySingle=((passedSingleEl && !passedAnyOther) || (passedSingleMu && !passedSingleEl && !passedAnyOther));
      bool trigConditionData = true;

      if (verbose) cout<<"checking PV"<<endl;
      const reco::Vertex *PV = 0;
      int theVertex = -1;
      for (unsigned int i=0; i<vertex->size(); i++){
        PV = &(vertex->at(i));
        if (verbose) std::cout<<"isFake: "<<PV->isFake()<<" chi2 "<<PV->chi2()<<" ndof "<<PV->ndof()<<" rho "<<PV->position().Rho()<<" Z "<<PV->position().Z()<<endl;
        if (PV->isFake()) continue;
        if (PV->ndof()<=4 || PV->position().Rho()>2.0 || fabs(PV->position().Z())>24.0) continue;
        theVertex=(int)i; break;
      }

      if (verbose) std::cout<<"vtx: "<<theVertex<<" trigConditionData "<<trigConditionData<<" passedTrig "<<passedTrig<<std::endl;
      if(theVertex >= 0 && (isMC || (!isMC && trigConditionData)) ){

        if (verbose) cout<<"good PV "<<theVertex<<endl;
        //N Vertex
        if (verbose) {cout<<"fill nvtx histogram"<<endl;}
        nVtx = vertex->size();
        histContainer_["NVTX"]->Fill(nVtx);
        histContainer_["NVTX_RW"]->Fill(nVtx,pileupWeight);

        //MET
        if (verbose) {cout<<"get met value"<<endl;}
        if (!mets->empty()){
          met = (*mets)[0].et();
          met_phi = (*mets)[0].phi();
          met_jesup = (*mets)[0].shiftedPt(pat::MET::JetEnUp);
          met_phi_jesup = (*mets)[0].shiftedPhi(pat::MET::JetEnUp);
          met_jesdn = (*mets)[0].shiftedPt(pat::MET::JetEnDown);
          met_phi_jesdn = (*mets)[0].shiftedPhi(pat::MET::JetEnDown);
          met_uncenup = (*mets)[0].shiftedPt(pat::MET::UnclusteredEnUp);
          met_phi_uncenup = (*mets)[0].shiftedPhi(pat::MET::UnclusteredEnUp);
          met_uncendn = (*mets)[0].shiftedPt(pat::MET::UnclusteredEnDown);
          met_phi_uncendn = (*mets)[0].shiftedPhi(pat::MET::UnclusteredEnDown);
        }

        if (verbose) cout<<"start lepton analysis"<<endl;
        vector<pat::Electron> AllElectrons; vector<pat::Muon> AllMuons;
        vector<pat::Electron> AllElectronsUnS;////uncorrected electron
        vector<pat::Tau> AllTaus; vector<pat::Photon> AllPhotons;
        AllElectrons = helper.goodLooseElectrons2012(electrons,_elecPtCut); // _elecPtCut > 7Gev
        AllElectronsUnS = helper.goodLooseElectrons2012(electrons,electronsUnS,_elecPtCut);
        AllMuons = helper.goodLooseMuons2012(muons,_muPtCut); //_muPtCut>5GeV
        AllTaus = helper.goodLooseTaus2015(taus,_tauPtCut);
        AllPhotons = helper.goodLoosePhotons2015(photons,_phoPtCut);//_phoPtCut>10Gev

        helper.cleanOverlappingLeptons(AllMuons,AllElectrons,PV);
        helper.cleanOverlappingLeptons(AllMuons,AllElectronsUnS,PV);
        recoMuons = helper.goodMuons2015_noIso_noPf(AllMuons,_muPtCut,PV,sip3dCut); //sip3dCut = 4
        recoElectrons = helper.goodElectrons2015_noIso_noBdt(AllElectrons,_elecPtCut,elecID,PV,iEvent,sip3dCut);
        recoElectronsUnS = helper.goodElectrons2015_noIso_noBdt(AllElectronsUnS,_elecPtCut,elecID,PV,iEvent,sip3dCut);
        helper.cleanOverlappingTaus(recoMuons,recoElectrons,AllTaus,isoCutMu,isoCutEl,muRho,elRho);
        recoTaus = helper.goodTaus2015(AllTaus,_tauPtCut); //_tauPtCut>20Gev
        recoPhotons = helper.goodPhotons2015(AllPhotons,_phoPtCut,year);

        if (verbose) cout<<AllMuons.size()<<" loose muons "<<AllElectrons.size()<<" loose electrons"<<endl;

        //sort electrons and muons by pt
        if (verbose) cout<<recoMuons.size()<<" good muons and "<<recoElectrons.size()<<" good electrons to be sorted"<<endl;
        if (verbose) cout<<"start pt-sorting leptons"<<endl;
        if (verbose) cout<<"adding muons to sorted list"<<endl;

        if( (recoMuons.size() + recoElectrons.size()) >= (uint)skimLooseLeptons ) { //skimLooseLeptons=2

          if (verbose) cout<<"found two leptons"<<endl;
          for(unsigned int i = 0; i < recoMuons.size(); i++){
            if (lep_ptreco.size()==0 || recoMuons[i].pt()<lep_ptreco[lep_ptreco.size()-1]){
              lep_ptreco.push_back(recoMuons[i].pt());
              lep_ptid.push_back(recoMuons[i].pdgId());
              lep_ptindex.push_back(i);
              continue;
            }
            for (unsigned int j=0; j<lep_ptreco.size(); j++){
              if (recoMuons[i].pt()>lep_ptreco[j]){
                lep_ptreco.insert(lep_ptreco.begin()+j,recoMuons[i].pt());
                lep_ptid.insert(lep_ptid.begin()+j,recoMuons[i].pdgId());
                lep_ptindex.insert(lep_ptindex.begin()+j,i);
                break;
              }
            }
          }

          if (verbose) cout<<"adding electrons to sorted list"<<endl;
          for(unsigned int i = 0; i < recoElectrons.size(); i++){
            if (lep_ptreco.size()==0 || recoElectrons[i].pt()<lep_ptreco[lep_ptreco.size()-1]) {
              lep_ptreco.push_back(recoElectrons[i].pt());
              lep_ptid.push_back(recoElectrons[i].pdgId());
              lep_ptindex.push_back(i);
              continue;
            }
            for (unsigned int j=0; j<lep_ptreco.size(); j++){
              if (recoElectrons[i].pt()>lep_ptreco[j]) {
                lep_ptreco.insert(lep_ptreco.begin()+j,recoElectrons[i].pt());
                lep_ptid.insert(lep_ptid.begin()+j,recoElectrons[i].pdgId());
                lep_ptindex.insert(lep_ptindex.begin()+j,i);
                break;
              }
            }
          }

          for(unsigned int i = 0; i < lep_ptreco.size(); i++){
            if (verbose) cout<<"sorted lepton "<<i<<" pt "<<lep_ptreco[i]<<" id "<<lep_ptid[i]<<" index "<<lep_ptindex[i]<<endl;

            if (abs(lep_ptid[i])==11){
              lep_id.push_back(recoElectrons[lep_ptindex[i]].pdgId());
              lep_pt.push_back(recoElectrons[lep_ptindex[i]].pt());
              lep_pterrold.push_back(recoElectrons[lep_ptindex[i]].p4Error(reco::GsfElectron::P4_COMBINATION));
              lep_isEB.push_back(recoElectrons[lep_ptindex[i]].isEB());
              lep_isEE.push_back(recoElectrons[lep_ptindex[i]].isEE());
              lep_p.push_back(recoElectrons[lep_ptindex[i]].p());
              lep_ecalEnergy.push_back(recoElectrons[lep_ptindex[i]].correctedEcalEnergy());
              double perr = 0.0;
              if (recoElectrons[lep_ptindex[i]].ecalDriven()){
                perr = recoElectrons[lep_ptindex[i]].p4Error(reco::GsfElectron::P4_COMBINATION);
              }
              else{
                double ecalEnergy = recoElectrons[lep_ptindex[i]].correctedEcalEnergy();
                double err2 = 0.0;
                if (recoElectrons[lep_ptindex[i]].isEB()){
                  err2 += (5.24e-02*5.24e-02)/ecalEnergy;
                  err2 += (2.01e-01*2.01e-01)/(ecalEnergy*ecalEnergy);
                  err2 += 1.00e-02*1.00e-02;
                } else if (recoElectrons[lep_ptindex[i]].isEE()){
                  err2 += (1.46e-01*1.46e-01)/ecalEnergy;
                  err2 += (9.21e-01*9.21e-01)/(ecalEnergy*ecalEnergy);
                  err2 += 1.94e-03*1.94e-03;
                }
                perr = ecalEnergy * sqrt(err2);
              }
              double pterr = perr*recoElectrons[lep_ptindex[i]].pt()/recoElectrons[lep_ptindex[i]].p();
              lep_pterr.push_back(pterr);
              lep_eta.push_back(recoElectrons[lep_ptindex[i]].eta());
              lep_phi.push_back(recoElectrons[lep_ptindex[i]].phi());
              lep_mass.push_back(recoElectrons[lep_ptindex[i]].mass());
              lepFSR_pt.push_back(recoElectrons[lep_ptindex[i]].pt());
              lepFSR_eta.push_back(recoElectrons[lep_ptindex[i]].eta());
              lepFSR_phi.push_back(recoElectrons[lep_ptindex[i]].phi());
              lepFSR_mass.push_back(recoElectrons[lep_ptindex[i]].mass());
              if (isoConeSizeEl==0.4){
                lep_RelIso.push_back(helper.pfIso(recoElectrons[lep_ptindex[i]],elRho));
                lep_RelIsoNoFSR.push_back(helper.pfIso(recoElectrons[lep_ptindex[i]],elRho));
                lep_isoCH.push_back(recoElectrons[lep_ptindex[i]].chargedHadronIso());
                lep_isoNH.push_back(recoElectrons[lep_ptindex[i]].neutralHadronIso());
                lep_isoPhot.push_back(recoElectrons[lep_ptindex[i]].photonIso());
                lep_isoPU.push_back(recoElectrons[lep_ptindex[i]].puChargedHadronIso());
                lep_isoPUcorr.push_back(helper.getPUIso(recoElectrons[lep_ptindex[i]],elRho));
              } else if(isoConeSizeEl==0.3){
                lep_RelIso.push_back(helper.pfIso03(recoElectrons[lep_ptindex[i]],elRho));
                lep_RelIsoNoFSR.push_back(helper.pfIso03(recoElectrons[lep_ptindex[i]],elRho));
                lep_isoCH.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumChargedHadronPt);
                lep_isoNH.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumNeutralHadronEt);
                lep_isoPhot.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumPhotonEt);
                lep_isoPU.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumPUPt);
                lep_isoPUcorr.push_back(helper.getPUIso03(recoElectrons[lep_ptindex[i]],elRho));
              }
              lep_MiniIso.push_back(helper.miniIsolation(pfCands, dynamic_cast<const reco::Candidate *>(&recoElectrons[lep_ptindex[i]]), 0.05, 0.2, 10., rhoSUS, false));
              lep_Sip.push_back(helper.getSIP3D(recoElectrons[lep_ptindex[i]]));
              //lep_mva.push_back(recoElectrons[lep_ptindex[i]].userFloat("ElectronMVAEstimatorRun2Autumn18IdIsoValues"));
              //cout<<EleBDT_name_161718<<endl;
              lep_mva.push_back(recoElectrons[lep_ptindex[i]].userFloat(EleBDT_name_161718.c_str()));
              lep_ecalDriven.push_back(recoElectrons[lep_ptindex[i]].ecalDriven());
              //lep_tightId.push_back(helper.passTight_BDT_Id(recoElectrons[lep_ptindex[i]],recoElectrons[lep_ptindex[i]].userFloat("ElectronMVAEstimatorRun2Autumn18IdIsoValues"), year));
              lep_tightId.push_back(helper.passTight_BDT_Id(recoElectronsUnS[lep_ptindex[i]],year));
              //lep_tightId_old.push_back(helper.passTight_BDT_Id(recoElectronsUnS[lep_ptindex[i]],year));
              //lep_tightId_old.push_back(helper.passTight_BDT_Id(0);
              //cout<<"old "<<recoElectrons[lep_ptindex[i]].userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values") <<" new" <<recoElectrons[lep_ptindex[i]].userFloat("ElectronMVAEstimatorRun2Spring16HZZV1Values")<<endl;
              lep_tightIdSUS.push_back(helper.passTight_Id_SUS(recoElectrons[lep_ptindex[i]],elecID,PV,BS,theConversions, year));
              //lep_tightIdHiPt.push_back(recoElectrons[lep_ptindex[i]].electronID("heepElectronID-HEEPV70"));
              lep_tightIdHiPt.push_back(recoElectrons[lep_ptindex[i]].electronID(heepID_name_161718.c_str()));
              lep_ptRatio.push_back(helper.ptRatio(recoElectrons[lep_ptindex[i]],jets,true)); // no L2L3 yet
              lep_ptRel.push_back(helper.ptRel(recoElectrons[lep_ptindex[i]],jets,true)); // no L2L3 yet
              lep_dataMC.push_back(helper.dataMC(recoElectrons[lep_ptindex[i]],hElecScaleFac,hElecScaleFac_Cracks,hElecScaleFacGsf,hElecScaleFacGsfLowET));
              lep_dataMCErr.push_back(helper.dataMCErr(recoElectrons[lep_ptindex[i]],hElecScaleFac,hElecScaleFac_Cracks));
              lep_genindex.push_back(-1.0);
            }
            if (abs(lep_ptid[i])==13){
              lep_id.push_back(recoMuons[lep_ptindex[i]].pdgId());
              lep_pt.push_back(recoMuons[lep_ptindex[i]].pt());
              lep_pterrold.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->ptError());
              lep_isEB.push_back(0);
              lep_isEE.push_back(0);
              lep_p.push_back(recoMuons[lep_ptindex[i]].p());
              lep_ecalEnergy.push_back(0);
              if (recoMuons[lep_ptindex[i]].hasUserFloat("correctedPtError")){
                lep_pterr.push_back(recoMuons[lep_ptindex[i]].userFloat("correctedPtError"));
              } else{
                lep_pterr.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->ptError());
              }
              lep_eta.push_back(recoMuons[lep_ptindex[i]].eta());
              lep_phi.push_back(recoMuons[lep_ptindex[i]].phi());
              if (recoMuons[lep_ptindex[i]].mass()<0.105) cout<<"muon mass: "<<recoMuons[lep_ptindex[i]].mass()<<endl;
              lep_mass.push_back(recoMuons[lep_ptindex[i]].mass());
              lepFSR_pt.push_back(recoMuons[lep_ptindex[i]].pt());
              lepFSR_eta.push_back(recoMuons[lep_ptindex[i]].eta());
              lepFSR_phi.push_back(recoMuons[lep_ptindex[i]].phi());
              lepFSR_mass.push_back(recoMuons[lep_ptindex[i]].mass());
              if (isoConeSizeMu==0.4){
                lep_RelIso.push_back(helper.pfIso(recoMuons[lep_ptindex[i]],muRho));
                lep_RelIsoNoFSR.push_back(helper.pfIso(recoMuons[lep_ptindex[i]],muRho));
                lep_isoCH.push_back(recoMuons[lep_ptindex[i]].chargedHadronIso());
                lep_isoNH.push_back(recoMuons[lep_ptindex[i]].neutralHadronIso());
                lep_isoPhot.push_back(recoMuons[lep_ptindex[i]].photonIso());
                lep_isoPU.push_back(recoMuons[lep_ptindex[i]].puChargedHadronIso());
                lep_isoPUcorr.push_back(helper.getPUIso(recoMuons[lep_ptindex[i]],muRho));
              } else if(isoConeSizeMu==0.3){
                lep_RelIso.push_back(helper.pfIso03(recoMuons[lep_ptindex[i]],muRho));
                lep_RelIsoNoFSR.push_back(helper.pfIso03(recoMuons[lep_ptindex[i]],muRho));
                lep_isoCH.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumChargedHadronPt);
                lep_isoNH.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumNeutralHadronEt);
                lep_isoPhot.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumPhotonEt);
                lep_isoPU.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumPUPt);
                lep_isoPUcorr.push_back(helper.getPUIso03(recoMuons[lep_ptindex[i]],muRho));
              }
              lep_MiniIso.push_back(helper.miniIsolation(pfCands, dynamic_cast<const reco::Candidate *>(&recoMuons[lep_ptindex[i]]), 0.05, 0.2, 10., rhoSUS, false));
              lep_Sip.push_back(helper.getSIP3D(recoMuons[lep_ptindex[i]]));
              lep_mva.push_back(recoMuons[lep_ptindex[i]].isPFMuon());
              //lep_mva.push_back(helper.get_Muon_MVA_Value(recoMuons[lep_ptindex[i]],vertex,muRho,year,PV));
              lep_ecalDriven.push_back(0);
              lep_tightId.push_back(helper.passTight_Id(recoMuons[lep_ptindex[i]],PV));
              //lep_tightId.push_back(helper.passTight_BDT_Id(recoMuons[lep_ptindex[i]],vertex,muRho,year,PV));
              //lep_tightId_old.push_back(helper.passTight_Id(recoMuons[lep_ptindex[i]],PV));
              //lep_tightId_old.push_back(helper.passTight_Id(0);
              lep_tightIdSUS.push_back(helper.passTight_Id_SUS(recoMuons[lep_ptindex[i]],PV));
              lep_tightIdHiPt.push_back(recoMuons[lep_ptindex[i]].isHighPtMuon(*PV));
              lep_ptRatio.push_back(helper.ptRatio(recoMuons[lep_ptindex[i]],jets,true)); // no L2L3 yet
              lep_ptRel.push_back(helper.ptRel(recoMuons[lep_ptindex[i]],jets,true)); // no L2L3 yet
              lep_dataMC.push_back(helper.dataMC(recoMuons[lep_ptindex[i]],hMuScaleFac));
              lep_dataMCErr.push_back(helper.dataMCErr(recoMuons[lep_ptindex[i]],hMuScaleFacUnc));
              lep_genindex.push_back(-1.0);
            }
            if (verbose) {cout<<" eta: "<<lep_eta[i]<<" phi: "<<lep_phi[i];
            if(abs(lep_ptid[i])==11)  cout<<" eSuperClusterOverP: "<<recoElectrons[lep_ptindex[i]].eSuperClusterOverP()<<" ecalEnergy: "<<recoElectrons[lep_ptindex[i]].ecalEnergy()<<" p: "<<recoElectrons[lep_ptindex[i]].p();
            cout<<" RelIso: "<<lep_RelIso[i]<<" isoCH: "<<lep_isoCH[i]<<" isoNH: "<<lep_isoNH[i]
                <<" isoPhot: "<<lep_isoPhot[i]<<" lep_isoPU: "<<lep_isoPU[i]<<" isoPUcorr: "<<lep_isoPUcorr[i]<<" Sip: "<<lep_Sip[i]
                <<" MiniIso: "<<lep_MiniIso[i]<<" ptRatio: "<<lep_ptRatio[i]<<" ptRel: "<<lep_ptRel[i]<<" lep_mva: "<<lep_mva[i];
                if(abs(lep_ptid[i])==11)    cout<<" SCeta: "<<recoElectrons[lep_ptindex[i]].superCluster()->eta()<<" dxy: "<<recoElectrons[lep_ptindex[i]].gsfTrack()->dxy(PV->position())<<" dz: "<<recoElectrons[lep_ptindex[i]].gsfTrack()->dz(PV->position());
                if(abs(lep_ptid[i])==11)    cout<<" Rho: "<<elRho<<" EleBDT_name: "<<EleBDT_name_161718<<" Uncorrected electron pt: "<<recoElectronsUnS[lep_ptindex[i]].pt();
                if(abs(lep_ptid[i])==13)    cout<<" Rho: "<<muRho;
                cout<<" dataMC: "<<lep_dataMC[i]<<" dataMCErr: "<<lep_dataMCErr[i];
                cout<<" lep_pterr: "<<lep_pterr[i]<<" lep_pterrold: "<<lep_pterrold[i]<<" lep_tightIdHiPt: "<<lep_tightIdHiPt[i]<<endl;
                if((abs(lep_ptid[i])==13)&&lep_pt[i]>200)    cout<<"Muon pt over 200 isTrackerHighPtID? "<<helper.isTrackerHighPt(recoMuons[lep_ptindex[i]],PV)<<endl;}
          }
          if (verbose) cout<<"adding taus to sorted list"<<endl;
          for(int i = 0; i < (int)recoTaus.size(); i++){
            tau_id.push_back(recoTaus[i].pdgId());
            tau_pt.push_back(recoTaus[i].pt());
            tau_eta.push_back(recoTaus[i].eta());
            tau_phi.push_back(recoTaus[i].phi());
            tau_mass.push_back(recoTaus[i].mass());
          }

          if (verbose) cout<<"adding photons to sorted list"<<endl;
          for(int i = 0; i < (int)recoPhotons.size(); i++){
            pho_pt.push_back(recoPhotons[i].pt());
            pho_eta.push_back(recoPhotons[i].eta());
            pho_phi.push_back(recoPhotons[i].phi());
            photonCutBasedIDLoose.push_back(recoPhotons[i].photonID("cutBasedPhotonID-Fall17-94X-V2-loose"));
          }

          if (doTriggerMatching){
            if (verbose) cout<<"start trigger matching"<<endl;
            //trigger matching
            for(unsigned int i = 0; i < lep_pt.size(); i++){

              TLorentzVector reco;
              reco.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);

              double reco_eta = reco.Eta();
              double reco_phi = reco.Phi();

              std::string filtersMatched = "";
              for (pat::TriggerObjectStandAlone obj : *triggerObjects){
                double hlt_eta = obj.eta();
                double hlt_phi = obj.phi();
                double dR =  deltaR(reco_eta,reco_phi,hlt_eta,hlt_phi);
                if (dR<0.5){
                  obj.unpackFilterLabels(iEvent, *trigger);
                  for (unsigned h = 0; h < obj.filterLabels().size(); ++h) filtersMatched += obj.filterLabels()[h];
                }
              }
              if (verbose) cout<<"Trigger matching lep id: "<<lep_id[i]<<" pt: "<<reco.Pt()<<" filters: "<<filtersMatched<<endl;
              lep_filtersMatched.push_back(filtersMatched);
            }
          }
        }
      }





}

// ------------ method called once each job just before starting event loop  ------------
void UFHZZ4LAna::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void UFHZZ4LAna::endJob()
{}

void UFHZZ4LAna::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{}

// ------------ method called when ending the processing of a run  ------------
void UFHZZ4LAna::endRun(const edm::Run& iRun, edm::EventSetup const&)
{}

// ------------ method called when starting to processes a luminosity block  ------------
void UFHZZ4LAna::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a luminosity block  ------------
void UFHZZ4LAna::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg,edm::EventSetup const& eSetup)
{}

  // ============================ UF Functions ============================= //


  //Find Z1,Z2, and Higgs candidate
  //Pass good leptons for analysis as candMuons and candElectrons
  //Pass empty vectors of pat leptons as selectedMuons and selectedElectrons
  // these will be filled in the function and then useable for more analysis.
void UFHZZ4LAna::findHiggsCandidate(std::vector< pat::Muon > &selectedMuons, std::vector< pat::Electron > &selectedElectrons,const edm::Event& iEvent )
{
  using namespace edm;
  using namespace pat;
  using namespace std;

  const double Zmass = 91.1876;

  unsigned int Nlep = lepFSR_pt.size();
  if (verbose) cout<<Nlep<<" leptons in total"<<endl;

  // First, make all Z candidates including any FSR photons
  int n_Zs = 0;
  std::vector<int> Z_lepindex1;
  std::vector<int> Z_lepindex2;
  for(unsigned int i=0; i<Nlep; i++){
    for(unsigned int j=i+1; j<Nlep; j++){

      if((lep_id[i]+lep[j]!=0)) continue;

      TLorentzVector li, lj;
      li.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);
      lj.SetPtEtaPhiM(lep_pt[j],lep_eta[j],lep_phi[j],lep_mass[j]);

      TLorentzVector lifsr, ljfsr;
      lifsr.SetPtEtaPhiM(lepFSR_pt[i],lepFSR_eta[i],lepFSR_phi[i],lepFSR_mass[i]);
      ljfsr.SetPtEtaPhiM(lepFSR_pt[j],lepFSR_eta[j],lepFSR_phi[j],lepFSR_mass[j]);

      TLorentzVector liljfsr = lifsr + ljfsr;

      if (verbose) {
                cout<<"OSSF pair: i="<<i<<" id1="<<lep_id[i]<<" j="<<j<<" id2="<<lep_id[j]<<" pt1: "<<lifsr.Pt()<<" pt2: "<<ljfsr.Pt()<<" M: "<<liljfsr.M()<<endl;
      }

      TLorentzVector Z, Z_noFSR;
      Z = lifsr+ljfsr;
      Z_noFSR = li+lj;

      if (verbose) cout<<"this Z mass: "<<Z.M()<<" mZ2Low: "<<mZ2Low<<endl;

      if(Z.M()>0.0){
        n_Zs++;
        Z_pt.push_back(Z.Pt());
        Z_eta.push_back(Z.Eta());
        Z_phi.push_back(Z.Phi());
        Z_mass.push_back(Z.M());
        Z_noFSR_pt.push_back(Z_noFSR.Pt());
        Z_noFSR_eta.push_back(Z_noFSR.Eta());
        Z_noFSR_phi.push_back(Z_noFSR.Phi());
        Z_noFSR_mass.push_back(Z_noFSR.M());
        Z_lepindex1.push_back(i);
        Z_lepindex2.push_back(j);
        if (verbose) cout<<" add Z_lepindex1: "<<i<<" Z_lepindex2: "<<j<<endl;
      }
    } //lep j
  }//lep i

  if( (recoMuons.size() + recoElectrons.size()) < 4 ) return;
  if (verbose) cout<<"found four leptons"<<endl;

  bool properLep_ID = false; int Nmm = 0; int Nmp = 0; int Nem = 0; int Nep = 0;
  for(unsigned int i=0; i<recoMuons.size(); i++){
    if(recoMuons[i].pdgId()<0) Nmm = Nmm+1;
    if(recoMuons[i].pdgId()>0) Nmp = Nmp+1;
  }
  for(unsigned int i =0; i<recoElectrons.size(); i++){
    if(recoElectrons[i].pdgId()<0) Nem = Nem+1;
    if(recoElectrons[i].pdgId()>0) Nep = Nep+1;
  }
  if(Nmm>=2 && Nmp>=2) properLep_ID = true; //4mu
  if(Nem>=2 && Nep>=2) properLep_ID = true; //4e
  if(Nmm>0 && Nmp>0 && Nem>0 && Nep>0) properLep_ID = true; //2e2mu

  // four proper charge flavor combination
  if(!properLep_ID) return;

  // Consider all ZZ candidates
  double minZ1DeltaM_SR=9999.9; double minZ1DeltaM_CR=99999.9;
  double maxZ2SumPt_SR=0.0; double maxZ2SumPt_CR=0.0;
  double max_D_bkg_kin_SR=0.0; double max_D_bkg_kin_CR=0.0;
  bool foundSRCandidate=false;

  for(int i=0; i<n_Zs; i++){
    for(int j=i+1; j<n_Zs; j++){

      int i1 = Z_lepindex1[i]; int i2 = Z_lepindex2[i];
      int j1 = Z_lepindex1[j]; int j2 = Z_lepindex2[j];

      if(i1==j1 || i1==j2 || i2==j1 || i2==j2) continue;

      TLorentzVector lep_i1, lep_i2, lep_j1, lep_j2;
      lep_i1.SetPtEtaPhiM(lepFSR_pt[i1],lepFSR_eta[i1],lepFSR_phi[i1],lepFSR_mass[i1]);
      lep_i2.SetPtEtaPhiM(lepFSR_pt[i2],lepFSR_eta[i2],lepFSR_phi[i2],lepFSR_mass[i2]);
      lep_j1.SetPtEtaPhiM(lepFSR_pt[j1],lepFSR_eta[j1],lepFSR_phi[j1],lepFSR_mass[j1]);
      lep_j2.SetPtEtaPhiM(lepFSR_pt[j2],lepFSR_eta[j2],lepFSR_phi[j2],lepFSR_mass[j2]);

      TLorentzVector lep_i1_nofsr, lep_i2_nofsr, lep_j1_nofsr, lep_j2_nofsr;
      lep_i1_nofsr.SetPtEtaPhiM(lep_pt[i1],lep_eta[i1],lep_phi[i1],lep_mass[i1]);
      lep_i2_nofsr.SetPtEtaPhiM(lep_pt[i2],lep_eta[i2],lep_phi[i2],lep_mass[i2]);
      lep_j1_nofsr.SetPtEtaPhiM(lep_pt[j1],lep_eta[j1],lep_phi[j1],lep_mass[j1]);
      lep_j2_nofsr.SetPtEtaPhiM(lep_pt[j2],lep_eta[j2],lep_phi[j2],lep_mass[j2]);

      TLorentzVector Zi, Zj;
      Zi.SetPtEtaPhiM(Z_pt[i],Z_eta[i],Z_phi[i],Z_mass[i]);
      Zj.SetPtEtaPhiM(Z_pt[j],Z_eta[j],Z_phi[j],Z_mass[j]);

      if (verbose) {cout<<"ZZ candidate Zi->M() "<<Zi.M()<<" Zj->M() "<<Zj.M()<<endl;}

      TLorentzVector Z1, Z2;
      int Z1index, Z2index;
      int Z1_lepindex[2] = {0,0};
      int Z2_lepindex[2] = {0,0};
      double Z1DeltaM, Z2SumPt;

      if (abs(Zi.M()-Zmass)<abs(Zj.M()-Zmass)){ //Zi mass closer to Zmass
        Z1index = i; Z2index = j;
        Z1 = Zi; Z2 = Zj;
        if(lep_i1.Pt()>lepi2.Pt()){ Z1_lepindex[0] = i1; Z1_lepindex[1] = i2;}
        else { Z1_lepindex[0] = i2;  Z1_lepindex[1] = i1; }
        if(lep_j1.Pt()>lep_j2.Pt()){ Z2_lepindex[0] = j1;  Z2_lepindex[1] = j2;}
        else{ Z2_lepindex[0] = j2;  Z2_lepindex[1] = j1;}
        Z1DeltaM = abs(Zi.M()-Zmass);
        Z2SumPt = lep_j1_nofsr.Pt()+lep_j2_nofsr.Pt();
      }
      else{
        Z1index = j; Z2index = i;
        if (lep_j1.Pt()>lep_j2.Pt()) { Z1_lepindex[0] = j1;  Z1_lepindex[1] = j2; }
        else { Z1_lepindex[0] = j2;  Z1_lepindex[1] = j1; }
        if (lep_i1.Pt()>lep_i2.Pt()) { Z2_lepindex[0] = i1;  Z2_lepindex[1] = i2; }
        else { Z2_lepindex[0] = i2;  Z2_lepindex[1] = i1; }
        Z1DeltaM = abs(Zj.M()-Zmass);
        Z2SumPt = lep_i1_nofsr.Pt()+lep_i2_nofsr.Pt();
      }

      // Check isolation cut (without FSR ) for Z1 leptons
      if (lep_RelIsoNoFSR[Z1_lepindex[0]]>((abs(lep_id[Z1_lepindex[0]])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed
      if (lep_RelIsoNoFSR[Z1_lepindex[1]]>((abs(lep_id[Z1_lepindex[1]])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed
      // Check tight ID cut for Z1 leptons
      if (!(lep_tightId[Z1_lepindex[0]])) continue; // checking tight lepton ID
      if (!(lep_tightId[Z1_lepindex[1]])) continue; // checking tight lepton ID

      // Check Leading and Subleading pt Cut
      vector<double> allPt;
      allPt.push_back(lep_i1_nofsr.Pt()); allPt.push_back(lep_i2_nofsr.Pt());
      allPt.push_back(lep_j1_nofsr.Pt()); allPt.push_back(lep_j2_nofsr.Pt());
      std::sort(allPt.begin(), allPt.end());
      if (verbose) cout<<" leading pt: "<<allPt[3]<<" cut: "<<leadingPtCut<<" subleadingPt: "<<allPt[2]<<" cut: "<<subleadingPtCut<<endl;
      if (allPt[3]<leadingPtCut || allPt[2]<subleadingPtCut ) continue;

      // Check dR(li,lj)>0.02 for any i,j
      vector<double> alldR;
      alldR.push_back(deltaR(lep_i1_nofsr.Eta(),lep_i1_nofsr.Phi(),lep_i2_nofsr.Eta(),lep_i2_nofsr.Phi())); //dR(i1,i2)
      alldR.push_back(deltaR(lep_i1_nofsr.Eta(),lep_i1_nofsr.Phi(),lep_j1_nofsr.Eta(),lep_j1_nofsr.Phi())); //dR(i1,j1)
      alldR.push_back(deltaR(lep_i1_nofsr.Eta(),lep_i1_nofsr.Phi(),lep_j2_nofsr.Eta(),lep_j2_nofsr.Phi()));
      alldR.push_back(deltaR(lep_i2_nofsr.Eta(),lep_i2_nofsr.Phi(),lep_j1_nofsr.Eta(),lep_j1_nofsr.Phi()));
      alldR.push_back(deltaR(lep_i2_nofsr.Eta(),lep_i2_nofsr.Phi(),lep_j2_nofsr.Eta(),lep_j2_nofsr.Phi()));
      alldR.push_back(deltaR(lep_j1_nofsr.Eta(),lep_j1_nofsr.Phi(),lep_j2_nofsr.Eta(),lep_j2_nofsr.Phi()));
      if (verbose) cout<<" minDr: "<<*min_element(alldR.begin(),alldR.end())<<endl;
      if (*min_element(alldR.begin(),alldR.end())<0.02) continue;

      // Check M(l+,l-)>4.0 GeV for any OS pair
      // Do not include FSR photons
      vector<double> allM;
      TLorentzVector i1i2;
      TLorentzVector _4l_temp;
      i1i2 = lep_i1_nofsr + lep_i2_nofsr; allM.push_back(i1i2.M());
      TLorentzVector j1j2;
      j1j2 = lep_j1_nofsr+lep_j2_nofsr; allM.push_back(j1j2.M());
      _4l_temp = Z1 + Z2;

      if (lep_id[i1]*lep_id[j1]<0){ //OS pair for i1j1 and i2j2
        TLorentzVector i1j1;
        i1j1 = (lep_i1_nofsr)+(lep_j1_nofsr); allM.push_back(i1j1.M());
        TLorentzVector i2j2;
        i2j2 = (lep_i2_nofsr)+(lep_j2_nofsr); allM.push_back(i2j2.M());
      } else{ //OS pair for i1j2 and i2j1
        TLorentzVector i1j2;
        i1j2 = (lep_i1_nofsr)+(lep_j2_nofsr); allM.push_back(i1j2.M());
        TLorentzVector i2j1;
        i2j1 = (lep_i2_nofsr)+(lep_j1_nofsr); allM.push_back(i2j1.M());
      }
      if (verbose) cout<<" min m(l+l-): "<<*min_element(allM.begin(),allM.end())<<endl;
      if (*min_element(allM.begin(),allM.end())<4.0) {passedQCDcut=false; continue;}

      // Check the "smart cut": !( |mZa-mZ| < |mZ1-mZ| && mZb<12)
      // only for 4mu or 4e ZZ candidates
      bool passSmartCut=true;
      if ( abs(lep_id[i1])==abs(lep_id[j1])){
        TLorentzVector Za, Zb;
        if (lep_id[i1]==lep_id[j1]){
          Za = (lep_i1)+(lep_j2);
          Zb = (lep_i2)+(lep_j1);
        } else{
          Za = (lep_i1)+(lep_j1);
          Zb = (lep_i2)+(lep_j2);
        }
        if ( abs(Za.M()-Zmass)<abs(Zb.M()-Zmass) ){ //Za mass closer to Zmass
          if (verbose) cout<<"abs(Za.M()-Zmass)-abs(Z1.M()-Zmass): "<<abs(Za.M()-Zmass)-abs(Z1.M()-Zmass)<<" Zb.M(): "<<Zb.M()<<endl;
          if ( abs(Za.M()-Zmass)<abs(Z1.M()-Zmass) && Zb.M()<mZ2Low ) passSmartCut=false;
        }
        else{
          if (verbose) cout<<"abs(Zb.M()-Zmass)-abs(Z1.M()-Zmass): "<<abs(Zb.M()-Zmass)-abs(Z1.M()-Zmass)<<" Za.M(): "<<Za.M()<<endl;
          if ( abs(Zb.M()-Zmass)<abs(Z1.M()-Zmass) && Za.M()<mZ2Low ) passSmartCut=false;
        }
      }
      if(!passSmartCut) continue;

      if (verbose) cout<<" massZ1: "<<Z1.M()<<" massZ2: "<<Z2.M()<<endl;
      if ( (Z1.M() < mZ1Low) || (Z1.M() > mZ1High) || (Z2.M() < mZ2Low) || (Z2.M() > mZ2High) ) continue;
      if (verbose) cout<<" mass4l: "<<_4l_temp.M()<<endl;
      if ( _4l_temp.M() < m4lLowCut ) continue;

      // Signal region if Z2 leptons are both tight ID Iso
      bool signalRegion=true;
      if (lep_RelIsoNoFSR[Z2_lepindex[0]]>((abs(lep_id[Z2_lepindex[0]])==11) ? isoCutEl : isoCutMu)) signalRegion=false; // checking iso with FSR removed
      if (lep_RelIsoNoFSR[Z2_lepindex[1]]>((abs(lep_id[Z2_lepindex[1]])==11) ? isoCutEl : isoCutMu)) signalRegion=false; // checking iso with FSR removed
      if (!(lep_tightId[Z2_lepindex[0]])) signalRegion=false; // checking tight lepton ID
      if (!(lep_tightId[Z2_lepindex[1]])) signalRegion=false; // checking tight lepton ID

      // Check if this candidate has the highest D_bkg_kin
      vector<TLorentzVector> P4s;
      P4s.clear();
      vector<int> tmpIDs;
      tmpIDs.clear();

      if (Z1_lepindex[0] == i1){
        P4s.push_back(lep_i1); P4s.push_back(lep_i2);
        if(Z2_lepindex[0]==j1){
          P4s.push_back(lep_j1); P4s.push_back(lep_j2);
        } else{
          P4s.push_back(lep_j2); P4s.push_back(lep_j1);
        }
      } else if (Z1_lepindex[0]==i2){
        P4s.push_back(lep_i2); P4s.push_back(lep_i1);
        if(Z2_lepindex[0]==j1){
          P4s.push_back(lep_j1); P4s.push_back(lep_j2);
        } else{
          P4s.push_back(lep_j2); P4s.push_back(lep_j1);
        }
      } else if(Z1_lepindex[0]==j1){
        P4s.push_back(lep_j1); P4s.push_back(lep_j2);
        if(Z2_lepindex[0]==i1){
          P4s.push_back(lep_i1); P4s.push_back(lep_i2);
        } else{
          P4s.push_back(lep_i2); P4s.push_back(lep_i1);
        }
      } else if(Z1_lepindex[0]==j2){
        P4s.push_back(lep_j2); P4s.push_back(lep_j1);
        if(Z2_lepindex[0]==i1){
          P4s.push_back(lep_i1); P4s.push_back(lep_i2);
        } else{
          P4s.push_back(lep_i2); P4s.push_back(lep_i1);
        }
      }

      tmpIDs.push_back(lep_id[Z1_lepindex[0]]); tmpIDs.push_back(lep_id[Z1_lepindex[1]]);
      tmpIDs.push_back(lep_id[Z2_lepindex[0]]); tmpIDs.push_back(lep_id[Z2_lepindex[1]]);

      SimpleParticleCollection_t daughters;
      daughters.push_back(SimpleParticle_t(tmpIDs[0],P4s[0]));
      daughters.push_back(SimpleParticle_t(tmpIDs[1],P4s[1]));
      daughters.push_back(SimpleParticle_t(tmpIDs[2],P4s[2]));
      daughters.push_back(SimpleParticle_t(tmpIDs[3],P4s[3]));

      SimpleParticleCollection_t associated;
      float D_bkg_kin_tmp;
      if(doMela){
        mela->setInputEvent(&daughters, &associated, 0, 0);
        mela->setCurrentCandidateFromIndex(0);

        float me_0plus_JHU_tmp, me_qqZZ_MCFM_tmp;
        mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
        mela->computeP(me_0plus_JHU_tmp, true);
        mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
        mela->computeP(me_qqZZ_MCFM_tmp, true);
        D_bkg_kin_tmp = me_0plus_JHU_tmp / (me_0plus_JHU_tmp + me_qqZZ_MCFM_tmp);

        mela->resetInputEvent();
      }

      if (verbose) cout<<"good ZZ candidate, D_bkg_kin: "<<D_bkg_kin_tmp<<" max D_bkg_kin SR: "<<max_D_bkg_kin_SR<<" max D_bkg_kin CR: "<<max_D_bkg_kin_CR<<endl;

      bool same4l=false;
      bool foundZ11=false; bool foundZ12=false; bool foundZ21=false; bool foundZ22=false;
      for(int l = 0; l<4; l++){
        if (lep_Hindex[l]==Z1_lepindex[0]) foundZ11 = true;
        if (lep_Hindex[l]==Z1_lepindex[1]) foundZ12 = true;
        if (lep_Hindex[l]==Z2_lepindex[0]) foundZ21 = true;
        if (lep_Hindex[l]==Z2_lepindex[1]) foundZ22 = true;
      }
      same4l = (foundZ11 && foundZ12 && foundZ21 && foundZ22);

      if(singnalRegion){ //Singnal Region has priorty

        if (!foundSRCandidate) same4l=false;
        if ( (bestCandMela && ((!same4l && D_bkg_kin_tmp>max_D_bkg_kin_SR) || (same4l && Z1DeltaM<=minZ1DeltaM_SR)))|| (!bestCandMela && Z1DeltaM<=minZ1DeltaM_SR) ){
          max_D_bkg_kin_SR = D_bkg_kin_tmp;
          minZ1DeltaM_SR = Z1DeltaM;

          if (!bestCandMela && Z_Hindex[0]==Z1index && Z2SumPt<maxZ2SumPt_SR) continue;

          Z_Hindex[0] = Z1index;
          lep_Hindex[0] = Z1_lepindex[0];
          lep_Hindex[1] = Z1_lepindex[1];

          maxZ2SumPt_SR = Z2SumPt;
          Z_Hindex[1] = Z2index;
          lep_Hindex[2] = Z2_lepindex[0];
          lep_Hindex[3] = Z2_lepindex[1];

          Z1Vec = Z1; Z2Vec = Z2; HVec = Z1+Z2;
          massZ1 = Z1Vec.M(); massZ2 = Z2Vec.M(); mass4l = HVec.M();

          if (verbose) cout<<" new best candidate SR: mass4l: "<<HVec.M()<<endl;
          if (HVec.M()>m4lLowCut){ //m4lLowCut move forward
            foundHiggsCandidate=true;
            foundSRCandidate=true;
          }
        }
      } else if(!foundSRCandidate){// Control regions get second priority

        if ( (bestCandMela && ((!same4l && D_bkg_kin_tmp>max_D_bkg_kin_CR) || (same4l && Z1DeltaM<=minZ1DeltaM_CR)))
                     || (!bestCandMela && Z1DeltaM<=minZ1DeltaM_CR) ){

                       max_D_bkg_kin_CR = D_bkg_kin_tmp;
                       minZ1DeltaM_CR = Z1DeltaM;

                       if (!bestCandMela && Z_Hindex[0]==Z1index && Z2SumPt<maxZ2SumPt_CR) continue;

                       Z_Hindex[0] = Z1index;
                       lep_Hindex[0] = Z1_lepindex[0];
                       lep_Hindex[1] = Z1_lepindex[1];

                       maxZ2SumPt_CR = Z2SumPt;
                       Z_Hindex[1] = Z2index;
                       lep_Hindex[2] = Z2_lepindex[0];
                       lep_Hindex[3] = Z2_lepindex[1];

                       Z1Vec = Z1; Z2Vec = Z2; HVec = Z1+Z2;
                       massZ1 = Z1Vec.M(); massZ2 = Z2Vec.M(); mass4l = HVec.M();
                       if (verbose) cout<<" new best candidate CR: mass4l: "<<HVec.M()<<endl;
                       if (HVec.M()>m4lLowCut) foundHiggsCandidate=true;  //m4lLowCut move forward
                     }
      }
      if (verbose) cout<<"Z_Hindex[0]: "<<Z_Hindex[0]<<" lep_Hindex[0]: "<<lep_Hindex[0]<<" lep_Hindex[1]: "<<lep_Hindex[1]
                             <<"Z_Hindex[1]: "<<Z_Hindex[1]<<" lep_Hindex[2]: "<<lep_Hindex[2]<<" lep_Hindex[3]: "<<lep_Hindex[3]<<endl;
    }//Zj
  }//Zi

  if(foundHiggsCandidate){
    if (verbose) cout<<" lep_Hindex[0]: "<<lep_Hindex[0]<<" lep_Hindex[1]: "<<lep_Hindex[1]<<" lep_Hindex[2]: "<<lep_Hindex[2]<<" lep_Hindex[3]: "<<lep_Hindex[3]<<endl;

    if (verbose) cout<<" lep_id[lep_Hindex[0]]: "<<lep_id[lep_Hindex[0]]<<" lep_id[lep_Hindex[1]]: "<<lep_id[lep_Hindex[1]]
                         <<" lep_id[lep_Hindex[2]]: "<<lep_id[lep_Hindex[2]]<<" lep_id[lep_Hindex[3]]: "<<lep_id[lep_Hindex[3]]<<endl;
    if ( abs(lep_id[lep_Hindex[0]])==13 && abs(lep_id[lep_Hindex[2]])==13 ){
      RecoFourMuEvent = true;
      selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[0]]]);
      selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[1]]]);
      selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[2]]]);
      selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[3]]]);
    }
    else if ( abs(lep_id[lep_Hindex[0]])==11 && abs(lep_id[lep_Hindex[2]])==13 ){
      RecoTwoETwoMuEvent = true;
      selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[0]]]);
      selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[1]]]);
      selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[2]]]);
      selectedMuons.push_back(recoMuons[lep_ptindex[lep_Hindex[3]]]);
    }
    else if (abs(lep_id[lep_Hindex[0]])==11 && abs(lep_id[lep_Hindex[2]])==11){
      RecoFourEEvent = true;
      selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[0]]]);
      selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[1]]]);
      selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[2]]]);
      selectedElectrons.push_back(recoElectrons[lep_ptindex[lep_Hindex[3]]]);
    }
  }
}

//Find Z + 1L candidate for fake rate study
void UFHZZ4LAna::findZ1LCandidate(const edm::Event& iEvent )
{
  using namespace edm;
  using namespace pat;
  using namespace std;

  const double Zmss = 91.1876;
  unsigned int Nlep = lepFSR_pt.size();
  if (verbose) cout<<Nlep<<" leptons in total"<<endl;
  if( Nlep != 3 ) return;

  // First, make all Z candidates including any FSR photons
  int n_Zs=0;
  vector<int> Z_Z1L_lepindex1;
  vector<int> Z_Z1L_lepindex2;

  for(unsigned int i=0; i<Nlep; i++){
    for(unsigned int j=i+1; j<Nlep; j++){

      // same flavor opposite charge
      if(lep_id[i]+lep[j])!= 0); continue;

      TLorentzVector li, lj;
      li.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);
      lj.SetPtEtaPhiM(lep_pt[j],lep_eta[j],lep_phi[j],lep_mass[j]);

      TLorentzVector lifsr, ljfsr;
      lifsr.SetPtEtaPhiM(lepFSR_pt[i],lepFSR_eta[i],lepFSR_phi[i],lepFSR_mass[i]);
      ljfsr.SetPtEtaPhiM(lepFSR_pt[j],lepFSR_eta[j],lepFSR_phi[j],lepFSR_mass[j]);

      TLorentzVector liljfsr = lifsr+ljfsr;

      if (verbose) {
                cout<<"OSSF pair: i="<<i<<" id1="<<lep_id[i]<<" j="<<j<<" id2="<<lep_id[j]<<" pt1: "<<lifsr.Pt()<<" pt2: "<<ljfsr.Pt()<<" M: "<<liljfsr.M()<<endl;
      }

      TLorentzVector Z, Z_noFSR;
      Z = lifsr+ljfsr;
      Z_noFSR = li+lj;
      if (verbose) cout<<"this Z mass: "<<Z.M()<<" mZ2Low: "<<mZ2Low<<endl;

      if (Z.M()>0.0){
        n_Zs++;
        Z_Z1L_lepindex1.push_back(i);
        Z_Z1L_lepindex2.push_back(j);
        if (verbose) cout<<" add Z_lepindex1: "<<i<<" Z_lepindex2: "<<j<<endl;
      }
    }// lep j
  }//lep i

  bool properLep_ID = false; int Nmm = 0; int Nmp = 0; int Nem = 0; int Nep = 0;
  for(unsigned int i =0; i<recoMuons.size(); i++){
    if(recoMuons[i].pdgId()<0) Nmm = Nmm+1;
    if(recoMuons[i].pdgId()>0) Nmp = Nmp+1;
  }
  for(unsigned int i=0; i<recoElectrons.size(); i++){
    if(recoElectrons[i].pdgId()<0) Nem = Nem+1;
    if(recoElectrons[i].pdgId()>0) Nep = Nep+1;
  }

  if(Nmm>=1 && Nmp>=1) properLep_ID = true; //2mu + x
  if(Nem>=1 && Nep>=1) properLep_ID = true; //2e + x

  // proper charge flavor combination for Z + 1L
  if(!properLep_ID) return;

  if (verbose) cout<<"found three leptons"<<endl;

  // Consider all Z candidates
  double minZ1DeltaM=9999.9;
  for (int i=0; i<n_Zs; i++){
    int i1 = Z_Z1L_lepindex1[i]; int i2 = Z_Z1L_lepindex2[i];
    int j1 = 3 - i1 - i2; // index of the third lepton (check if this works)

    TLorentzVector lep_i1, lep_i2, lep_j1;
    lep_i1.SetPtEtaPhiM(lepFSR_pt[i1],lepFSR_eta[i1],lepFSR_phi[i1],lepFSR_mass[i1]);
    lep_i2.SetPtEtaPhiM(lepFSR_pt[i2],lepFSR_eta[i2],lepFSR_phi[i2],lepFSR_mass[i2]);
    lep_j1.SetPtEtaPhiM(lepFSR_pt[j1],lepFSR_eta[j1],lepFSR_phi[j1],lepFSR_mass[j1]);

    TLorentzVector lep_i1_nofsr, lep_i2_nofsr, lep_j1_nofsr;
    lep_i1_nofsr.SetPtEtaPhiM(lep_pt[i1],lep_eta[i1],lep_phi[i1],lep_mass[i1]);
    lep_i2_nofsr.SetPtEtaPhiM(lep_pt[i2],lep_eta[i2],lep_phi[i2],lep_mass[i2]);
    lep_j1_nofsr.SetPtEtaPhiM(lep_pt[j1],lep_eta[j1],lep_phi[j1],lep_mass[j1]);

    TLorentzVector Zi;
    Zi = lep_i1+lep_i2;
    if (verbose) {cout<<"Z candidate Zi->M() "<<Zi.M()<<endl;}

    TLorentzVector Z1 = Zi;
    double Z1DeltaM = abs(Zi.M()-Zmass);
    int Z1_lepindex[2] = {0,0};
    if (lep_i1.Pt()>lep_i2.Pt()) { Z1_lepindex[0] = i1;  Z1_lepindex[1] = i2; }
    else { Z1_lepindex[0] = i2;  Z1_lepindex[1] = i1; }

    // Check Leading and Subleading pt Cut
    vector<double> allPt;
    allPt.push_back(lep_i1.Pt()); allPt.push_back(lep_i2.Pt());
    std::sort(allPt.begin(), allPt.end());
    if (verbose) cout<<" leading pt: "<<allPt[1]<<" cut: "<<leadingPtCut<<" subleadingPt: "<<allPt[0]<<" cut: "<<subleadingPtCut<<endl;
    if (allPt[1]<leadingPtCut || allPt[0]<subleadingPtCut ) continue;

    // Check dR(li,lj)>0.02 for any i,j
    vector<double> alldR;
    alldR.push_back(deltaR(lep_i1.Eta(),lep_i1.Phi(),lep_i2.Eta(),lep_i2.Phi()));
    alldR.push_back(deltaR(lep_i1.Eta(),lep_i1.Phi(),lep_j1.Eta(),lep_j1.Phi()));
    alldR.push_back(deltaR(lep_i2.Eta(),lep_i2.Phi(),lep_j1.Eta(),lep_j1.Phi()));
    if (verbose) cout<<" minDr: "<<*min_element(alldR.begin(),alldR.end())<<endl;
    if (*min_element(alldR.begin(),alldR.end())<0.02) continue;

    // Check M(l+,l-)>4.0 GeV for any OS pair
    // Do not include FSR photons
    vector<double> allM;
    TLorentzVector i1i2;
    i1i2 = (lep_i1_nofsr)+(lep_i2_nofsr); allM.push_back(i1i2.M());
    if (lep_id[i1]*lep_id[j1]<0){
      TLorentzVector i1j1;
      i1j1 = (lep_i1_nofsr)+(lep_j1_nofsr); allM.push_back(i1j1.M());
    } else{
      TLorentzVector i2j1;
      i2j1 = (lep_i2_nofsr)+(lep_j1_nofsr); allM.push_back(i2j1.M());
    }
    if (verbose) cout<<" min m(l+l-): "<<*min_element(allM.begin(),allM.end())<<endl;
    if (*min_element(allM.begin(),allM.end())<4.0) {passedQCDcut=false; continue;}

    // Check isolation cut (without FSR ) for Z1 leptons
    if (lep_RelIsoNoFSR[Z1_lepindex[0]]>((abs(lep_id[Z1_lepindex[0]])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed
    if (lep_RelIsoNoFSR[Z1_lepindex[1]]>((abs(lep_id[Z1_lepindex[1]])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed
    // Check tight ID cut for Z1 leptons
    if (!(lep_tightId[Z1_lepindex[0]])) continue; // checking tight lepton ID
    if (!(lep_tightId[Z1_lepindex[1]])) continue; // checking tight lepton ID

    if ( (Z1.M() < mZ1Low) || (Z1.M() > mZ1High) ) continue;
    if (verbose) cout<<"good Z1L candidate, Z1DeltaM: "<<Z1DeltaM<<" minZ1DeltaM: "<<minZ1DeltaM<<endl;

    // Check if this candidate has the best Z1 and highest scalar sum of Z2 lepton pt
    if ( Z1DeltaM<=minZ1DeltaM ){
      minZ1DeltaM = Z1DeltaM;

      TLorentzVector Z1L;
      Z1L = Z1+lep_j1;

      massZ1_Z1L = Z1.M();
      mass3l = Z1L.M();

      lep_Hindex[0] = Z1_lepindex[0];
      lep_Hindex[1] = Z1_lepindex[1];
      lep_Hindex[2] = j1;

      if (verbose) cout<<" new best Z1L candidate: massZ1: "<<massZ1<<" (mass3l: "<<mass3l<<")"<<endl;
      foundZ1LCandidate=true;
    }
  }
}

void UFHZZ4LAna::bookPassedEventTree(TString treeName, TTree *tree)
{}

void UFHZZ4LAna::setTreeVariables( const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                     std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons,
                                     std::vector<pat::Muon> recoMuons, std::vector<pat::Electron> recoElectrons,
                                     std::vector<pat::Jet> goodJets, std::vector<float> goodJetQGTagger,
                                     std::vector<float> goodJetaxis2, std::vector<float> goodJetptD, std::vector<int> goodJetmult,
                                     std::vector<pat::Jet> selectedMergedJets)
{}

void UFHZZ4LAna::setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles,
                                   edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles,
                                   edm::Handle<edm::View<reco::GenJet> > genJets)
{}

bool UFHZZ4LAna::mZ1_mZ2(unsigned int& L1, unsigned int& L2, unsigned int& L3, unsigned int& L4, bool makeCuts)
{}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void UFHZZ4LAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{}

//define this as a plug-in
DEFINE_FWK_MODULE(UFHZZ4LAna);

//  LocalWords:  ecalDriven
