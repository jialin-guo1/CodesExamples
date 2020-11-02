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
    eleSrc
