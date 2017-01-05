//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec 13 00:42:01 2016 by ROOT version 6.06/06
// from TTree BaconTree/Events
// found on file: DMSpin0_ggPhibb1j_150_1000pb_weighted.root
//////////////////////////////////////////////////////////

#ifndef BaconTree_h
#define BaconTree_h
using namespace std;
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class BaconTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          runNum;
   UInt_t          lumiSec;
   UInt_t          evtNum;
   UInt_t          passJson;
   UInt_t          metfilter;
   UInt_t          triggerBits;
   UInt_t          selectBits;
   Double_t        sf_eleTrig;
   Double_t        sf_metTrig;
   Double_t        sf_phoTrig;
   Double_t        sf_lep;
   Double_t        sf_lepTrack;
   UInt_t          npu;
   UInt_t          npv;
   Float_t         puWeight;
   Float_t         scale1fb;
   Float_t         evtWeight;
   Float_t         kfactor;
   Float_t         kfactorNLO;
   Float_t         PDF;
   Float_t         PDF_UP;
   Float_t         PDF_DO;
   Float_t         RenScale_UP;
   Float_t         RenScale_DO;
   Float_t         FacScale_UP;
   Float_t         FacScale_DO;
   Float_t         pfmet;
   Float_t         pfmetphi;
   Float_t         puppet;
   Float_t         puppetphi;
   Double_t        AK8Puppijet0_pt;
   Double_t        AK8Puppijet0_eta;
   Double_t        AK8Puppijet0_phi;
   Double_t        AK8Puppijet1_pt;
   Double_t        AK8Puppijet1_eta;
   Double_t        AK8Puppijet1_phi;
   Double_t        AK8Puppijet2_pt;
   Double_t        AK8Puppijet2_eta;
   Double_t        AK8Puppijet2_phi;
   Double_t        AK8Puppijet0_mass;
   Double_t        AK8Puppijet0_csv;
   Double_t        AK8Puppijet0_CHF;
   Double_t        AK8Puppijet0_NHF;
   Double_t        AK8Puppijet0_NEMF;
   Double_t        AK8Puppijet0_tau21;
   Double_t        AK8Puppijet0_tau32;
   Double_t        AK8Puppijet0_msd;
   Double_t        AK8Puppijet0_rho;
   Double_t        AK8Puppijet0_minsubcsv;
   Double_t        AK8Puppijet0_maxsubcsv;
   Double_t        AK8Puppijet0_doublecsv;
   Double_t        AK8Puppijet0_doublesub;
   Double_t        AK8Puppijet0_ptraw;
   Double_t        AK8Puppijet0_genpt;
   Double_t        AK8Puppijet0_N2sdb1;
   Double_t        AK8Puppijet0_N2sdb2;
   Double_t        AK8Puppijet0_M2sdb1;
   Double_t        AK8Puppijet0_M2sdb2;
   Double_t        AK8Puppijet0_D2sdb1;
   Double_t        AK8Puppijet0_D2sdb2;
   Double_t        AK8Puppijet0_N2b1;
   Double_t        AK8Puppijet0_N2b2;
   Double_t        AK8Puppijet0_M2b1;
   Double_t        AK8Puppijet0_M2b2;
   Double_t        AK8Puppijet0_D2b1;
   Double_t        AK8Puppijet0_D2b2;
   Double_t        AK8Puppijet1_csv;
   Double_t        AK8Puppijet1_CHF;
   Double_t        AK8Puppijet1_NHF;
   Double_t        AK8Puppijet1_NEMF;
   Double_t        AK8Puppijet1_tau21;
   Double_t        AK8Puppijet1_tau32;
   Double_t        AK8Puppijet1_msd;
   Double_t        AK8Puppijet1_rho;
   Double_t        AK8Puppijet1_minsubcsv;
   Double_t        AK8Puppijet1_maxsubcsv;
   Double_t        AK8Puppijet1_doublecsv;
   Double_t        AK8Puppijet1_doublesub;
   Double_t        AK8Puppijet1_ptraw;
   Double_t        AK8Puppijet1_genpt;
   Double_t        AK8Puppijet1_e2_b1;
   Double_t        AK8Puppijet1_e3_b1;
   Double_t        AK8Puppijet1_e3_v1_b1;
   Double_t        AK8Puppijet1_e3_v2_b1;
   Double_t        AK8Puppijet1_e4_v1_b1;
   Double_t        AK8Puppijet1_e4_v2_b1;
   Double_t        AK8Puppijet1_e2_b2;
   Double_t        AK8Puppijet1_e3_b2;
   Double_t        AK8Puppijet1_e3_v1_b2;
   Double_t        AK8Puppijet1_e3_v2_b2;
   Double_t        AK8Puppijet1_e4_v1_b2;
   Double_t        AK8Puppijet1_e4_v2_b2;
   Double_t        AK8Puppijet1_e2_sdb1;
   Double_t        AK8Puppijet1_e3_sdb1;
   Double_t        AK8Puppijet1_e3_v1_sdb1;
   Double_t        AK8Puppijet1_e3_v2_sdb1;
   Double_t        AK8Puppijet1_e4_v1_sdb1;
   Double_t        AK8Puppijet1_e4_v2_sdb1;
   Double_t        AK8Puppijet1_e2_sdb2;
   Double_t        AK8Puppijet1_e3_sdb2;
   Double_t        AK8Puppijet1_e3_v1_sdb2;
   Double_t        AK8Puppijet1_e3_v2_sdb2;
   Double_t        AK8Puppijet1_e4_v1_sdb2;
   Double_t        AK8Puppijet1_e4_v2_sdb2;
   Double_t        AK8Puppijet1_N2sdb1;
   Double_t        AK8Puppijet1_N2sdb2;
   Double_t        AK8Puppijet1_M2sdb1;
   Double_t        AK8Puppijet1_M2sdb2;
   Double_t        AK8Puppijet1_D2sdb1;
   Double_t        AK8Puppijet1_D2sdb2;
   Double_t        AK8Puppijet1_N2b1;
   Double_t        AK8Puppijet1_N2b2;
   Double_t        AK8Puppijet1_M2b1;
   Double_t        AK8Puppijet1_M2b2;
   Double_t        AK8Puppijet1_D2b1;
   Double_t        AK8Puppijet1_D2b2;
   Double_t        AK8Puppijet1_mass;
   Double_t        AK8Puppijet2_csv;
   Double_t        AK8Puppijet2_CHF;
   Double_t        AK8Puppijet2_NHF;
   Double_t        AK8Puppijet2_NEMF;
   Double_t        AK8Puppijet2_tau21;
   Double_t        AK8Puppijet2_tau32;
   Double_t        AK8Puppijet2_msd;
   Double_t        AK8Puppijet2_rho;
   Double_t        AK8Puppijet2_minsubcsv;
   Double_t        AK8Puppijet2_maxsubcsv;
   Double_t        AK8Puppijet2_doublecsv;
   Double_t        AK8Puppijet2_doublesub;
   Double_t        AK8Puppijet2_ptraw;
   Double_t        AK8Puppijet2_genpt;
   Double_t        AK8Puppijet2_e2_b1;
   Double_t        AK8Puppijet2_e3_b1;
   Double_t        AK8Puppijet2_e3_v1_b1;
   Double_t        AK8Puppijet2_e3_v2_b1;
   Double_t        AK8Puppijet2_e4_v1_b1;
   Double_t        AK8Puppijet2_e4_v2_b1;
   Double_t        AK8Puppijet2_e2_b2;
   Double_t        AK8Puppijet2_e3_b2;
   Double_t        AK8Puppijet2_e3_v1_b2;
   Double_t        AK8Puppijet2_e3_v2_b2;
   Double_t        AK8Puppijet2_e4_v1_b2;
   Double_t        AK8Puppijet2_e4_v2_b2;
   Double_t        AK8Puppijet2_e2_sdb1;
   Double_t        AK8Puppijet2_e3_sdb1;
   Double_t        AK8Puppijet2_e3_v1_sdb1;
   Double_t        AK8Puppijet2_e3_v2_sdb1;
   Double_t        AK8Puppijet2_e4_v1_sdb1;
   Double_t        AK8Puppijet2_e4_v2_sdb1;
   Double_t        AK8Puppijet2_e2_sdb2;
   Double_t        AK8Puppijet2_e3_sdb2;
   Double_t        AK8Puppijet2_e3_v1_sdb2;
   Double_t        AK8Puppijet2_e3_v2_sdb2;
   Double_t        AK8Puppijet2_e4_v1_sdb2;
   Double_t        AK8Puppijet2_e4_v2_sdb2;
   Double_t        AK8Puppijet2_N2sdb1;
   Double_t        AK8Puppijet2_N2sdb2;
   Double_t        AK8Puppijet2_M2sdb1;
   Double_t        AK8Puppijet2_M2sdb2;
   Double_t        AK8Puppijet2_D2sdb1;
   Double_t        AK8Puppijet2_D2sdb2;
   Double_t        AK8Puppijet2_N2b1;
   Double_t        AK8Puppijet2_N2b2;
   Double_t        AK8Puppijet2_M2b1;
   Double_t        AK8Puppijet2_M2b2;
   Double_t        AK8Puppijet2_D2b1;
   Double_t        AK8Puppijet2_D2b2;
   Int_t           nAK8Puppijets;
   Int_t           AK8Puppijet0_isTightVJet;
   Int_t           AK8Puppijet0_isHadronicV;
   Double_t        AK8Puppijet0_vMatching;
   Double_t        AK8Puppijet0_vSize;
   Int_t           AK8Puppijet0_partonFlavor;
   Int_t           AK8Puppijet0_hadronFlavor;
   Int_t           AK8Puppijet0_nCharged;
   Int_t           AK8Puppijet0_nNeutrals;
   Int_t           AK8Puppijet0_nParticles;
   Double_t        AK8Puppijet0_ratioCA15_04;
   Double_t        AK8CHSjet0_doublecsv;
   Double_t        AK8CHSjet0_doublesub;
   Double_t        AK8CHSjet0_pt;
   Double_t        AK8CHSjet0_eta;
   Double_t        AK8CHSjet0_phi;
   Double_t        AK8CHSjet1_doublecsv;
   Double_t        AK8CHSjet1_doublesub;
   Double_t        AK8CHSjet1_pt;
   Double_t        AK8CHSjet1_eta;
   Double_t        AK8CHSjet1_phi;
   Double_t        AK8CHSjet2_doublecsv;
   Double_t        AK8CHSjet2_doublesub;
   Double_t        AK8CHSjet2_pt;
   Double_t        AK8CHSjet2_eta;
   Double_t        AK8CHSjet2_phi;
   Int_t           nCA15Puppijets;
   Double_t        CA15CHSjet0_doublecsv;
   Double_t        CA15CHSjet0_doublesub;
   Double_t        CA15CHSjet0_pt;
   Double_t        CA15CHSjet0_eta;
   Double_t        CA15CHSjet0_phi;
   Double_t        CA15CHSjet1_doublecsv;
   Double_t        CA15CHSjet1_doublesub;
   Double_t        CA15CHSjet1_pt;
   Double_t        CA15CHSjet1_eta;
   Double_t        CA15CHSjet1_phi;
   Double_t        CA15CHSjet2_doublecsv;
   Double_t        CA15CHSjet2_doublesub;
   Double_t        CA15CHSjet2_pt;
   Double_t        CA15CHSjet2_eta;
   Double_t        CA15CHSjet2_phi;
   Int_t           nAK4Puppijets;
   Int_t           nAK4Puppijetsfwd;
   Int_t           nAK4PuppijetsbtagL;
   Int_t           nAK4PuppijetsbtagM;
   Int_t           nAK4PuppijetsbtagT;
   Int_t           nAK4PuppijetsdR08;
   Int_t           nAK4PuppijetsLdR08;
   Int_t           nAK4PuppijetsMdR08;
   Int_t           nAK4PuppijetsTdR08;
   Int_t           nAK4PuppijetsLPt100dR08;
   Int_t           nAK4PuppijetsMPt100dR08;
   Int_t           nAK4PuppijetsTPt100dR08;
   Int_t           nAK4PuppijetsLPt150dR08;
   Int_t           nAK4PuppijetsMPt150dR08;
   Int_t           nAK4PuppijetsTPt150dR08;
   Double_t        AK4Puppijet0_pt;
   Double_t        AK4Puppijet0_eta;
   Double_t        AK4Puppijet0_phi;
   Double_t        AK4Puppijet1_pt;
   Double_t        AK4Puppijet1_eta;
   Double_t        AK4Puppijet1_phi;
   Double_t        AK4Puppijet2_pt;
   Double_t        AK4Puppijet2_eta;
   Double_t        AK4Puppijet2_phi;
   Double_t        AK4Puppijet3_pt;
   Double_t        AK4Puppijet3_eta;
   Double_t        AK4Puppijet3_phi;
   Double_t        AK4Puppijet0_mass;
   Double_t        AK4Puppijet0_csv;
   Double_t        AK4Puppijet0_qgid;
   Double_t        AK4Puppijet0_dR08;
   Double_t        AK4Puppijet0_dPhi08;
   Double_t        AK4Puppijet1_mass;
   Double_t        AK4Puppijet1_csv;
   Double_t        AK4Puppijet1_qgid;
   Double_t        AK4Puppijet1_dR08;
   Double_t        AK4Puppijet1_dPhi08;
   Double_t        AK4Puppijet2_mass;
   Double_t        AK4Puppijet2_csv;
   Double_t        AK4Puppijet2_qgid;
   Double_t        AK4Puppijet2_dR08;
   Double_t        AK4Puppijet2_dPhi08;
   Double_t        AK4Puppijet3_mass;
   Double_t        AK4Puppijet3_csv;
   Double_t        AK4Puppijet3_qgid;
   Double_t        AK4Puppijet3_dR08;
   Double_t        AK4Puppijet3_dPhi08;
   Int_t           nmuLoose;
   Int_t           nmuTight;
   Int_t           ismu0Tight;
   Double_t        vmuoLoose0_pt;
   Double_t        vmuoLoose0_eta;
   Double_t        vmuoLoose0_phi;
   Double_t        vmuoLoose1_pt;
   Double_t        vmuoLoose1_eta;
   Double_t        vmuoLoose1_phi;
   Int_t           neleLoose;
   Int_t           neleTight;
   Int_t           isele0Tight;
   Double_t        veleLoose0_pt;
   Double_t        veleLoose0_eta;
   Double_t        veleLoose0_phi;
   Double_t        veleLoose1_pt;
   Double_t        veleLoose1_eta;
   Double_t        veleLoose1_phi;
   Int_t           ntau;
   Int_t           nphoLoose;
   Int_t           nphoTight;
   Int_t           ispho0Tight;
   Double_t        vpho0_pt;
   Double_t        vpho0_eta;
   Double_t        vpho0_phi;
   Float_t         genVPt;
   Float_t         genVPhi;
   Float_t         genVMass;
   Float_t         genVEta;
   Int_t           genVPdfId;
   Int_t           genEleFromW;
   Int_t           genMuFromW;
   Int_t           genTauFromW;

   // List of branches
   TBranch        *b_fRun;   //!
   TBranch        *b_fLumi;   //!
   TBranch        *b_fEvtV;   //!
   TBranch        *b_fPassJson;   //!
   TBranch        *b_fMetFilters;   //!
   TBranch        *b_fITrigger;   //!
   TBranch        *b_fselectBits;   //!
   TBranch        *b_fsf_eleTrig;   //!
   TBranch        *b_fsf_metTrig;   //!
   TBranch        *b_fsf_phoTrig;   //!
   TBranch        *b_fsf_lep;   //!
   TBranch        *b_fsf_lepTrack;   //!
   TBranch        *b_fNPU;   //!
   TBranch        *b_fNVtx;   //!
   TBranch        *b_fPUWeight;   //!
   TBranch        *b_fScale;   //!
   TBranch        *b_fevtWeight;   //!
   TBranch        *b_fkfactor;   //!
   TBranch        *b_fkFactor_CENT;   //!
   TBranch        *b_fPDF;   //!
   TBranch        *b_fPDF_UP;   //!
   TBranch        *b_fPDF_DO;   //!
   TBranch        *b_fRenScale_UP;   //!
   TBranch        *b_fRenScale_DO;   //!
   TBranch        *b_fFacScale_UP;   //!
   TBranch        *b_fFacScale_DO;   //!
   TBranch        *b_fMet;   //!
   TBranch        *b_fMetPhi;   //!
   TBranch        *b_fPuppEt;   //!
   TBranch        *b_fPuppEtPhi;   //!
   TBranch        *b_AK8Puppijet0_pt;   //!
   TBranch        *b_AK8Puppijet0_eta;   //!
   TBranch        *b_AK8Puppijet0_phi;   //!
   TBranch        *b_AK8Puppijet1_pt;   //!
   TBranch        *b_AK8Puppijet1_eta;   //!
   TBranch        *b_AK8Puppijet1_phi;   //!
   TBranch        *b_AK8Puppijet2_pt;   //!
   TBranch        *b_AK8Puppijet2_eta;   //!
   TBranch        *b_AK8Puppijet2_phi;   //!
   TBranch        *b_AK8Puppijet0_mass;   //!
   TBranch        *b_AK8Puppijet0_csv;   //!
   TBranch        *b_AK8Puppijet0_CHF;   //!
   TBranch        *b_AK8Puppijet0_NHF;   //!
   TBranch        *b_AK8Puppijet0_NEMF;   //!
   TBranch        *b_AK8Puppijet0_tau21;   //!
   TBranch        *b_AK8Puppijet0_tau32;   //!
   TBranch        *b_AK8Puppijet0_msd;   //!
   TBranch        *b_AK8Puppijet0_rho;   //!
   TBranch        *b_AK8Puppijet0_minsubcsv;   //!
   TBranch        *b_AK8Puppijet0_maxsubcsv;   //!
   TBranch        *b_AK8Puppijet0_doublecsv;   //!
   TBranch        *b_AK8Puppijet0_doublesub;   //!
   TBranch        *b_AK8Puppijet0_ptraw;   //!
   TBranch        *b_AK8Puppijet0_genpt;   //!
   TBranch        *b_AK8Puppijet0_N2sdb1;   //!
   TBranch        *b_AK8Puppijet0_N2sdb2;   //!
   TBranch        *b_AK8Puppijet0_M2sdb1;   //!
   TBranch        *b_AK8Puppijet0_M2sdb2;   //!
   TBranch        *b_AK8Puppijet0_D2sdb1;   //!
   TBranch        *b_AK8Puppijet0_D2sdb2;   //!
   TBranch        *b_AK8Puppijet0_N2b1;   //!
   TBranch        *b_AK8Puppijet0_N2b2;   //!
   TBranch        *b_AK8Puppijet0_M2b1;   //!
   TBranch        *b_AK8Puppijet0_M2b2;   //!
   TBranch        *b_AK8Puppijet0_D2b1;   //!
   TBranch        *b_AK8Puppijet0_D2b2;   //!
   TBranch        *b_AK8Puppijet1_csv;   //!
   TBranch        *b_AK8Puppijet1_CHF;   //!
   TBranch        *b_AK8Puppijet1_NHF;   //!
   TBranch        *b_AK8Puppijet1_NEMF;   //!
   TBranch        *b_AK8Puppijet1_tau21;   //!
   TBranch        *b_AK8Puppijet1_tau32;   //!
   TBranch        *b_AK8Puppijet1_msd;   //!
   TBranch        *b_AK8Puppijet1_rho;   //!
   TBranch        *b_AK8Puppijet1_minsubcsv;   //!
   TBranch        *b_AK8Puppijet1_maxsubcsv;   //!
   TBranch        *b_AK8Puppijet1_doublecsv;   //!
   TBranch        *b_AK8Puppijet1_doublesub;   //!
   TBranch        *b_AK8Puppijet1_ptraw;   //!
   TBranch        *b_AK8Puppijet1_genpt;   //!
   TBranch        *b_AK8Puppijet1_e2_b1;   //!
   TBranch        *b_AK8Puppijet1_e3_b1;   //!
   TBranch        *b_AK8Puppijet1_e3_v1_b1;   //!
   TBranch        *b_AK8Puppijet1_e3_v2_b1;   //!
   TBranch        *b_AK8Puppijet1_e4_v1_b1;   //!
   TBranch        *b_AK8Puppijet1_e4_v2_b1;   //!
   TBranch        *b_AK8Puppijet1_e2_b2;   //!
   TBranch        *b_AK8Puppijet1_e3_b2;   //!
   TBranch        *b_AK8Puppijet1_e3_v1_b2;   //!
   TBranch        *b_AK8Puppijet1_e3_v2_b2;   //!
   TBranch        *b_AK8Puppijet1_e4_v1_b2;   //!
   TBranch        *b_AK8Puppijet1_e4_v2_b2;   //!
   TBranch        *b_AK8Puppijet1_e2_sdb1;   //!
   TBranch        *b_AK8Puppijet1_e3_sdb1;   //!
   TBranch        *b_AK8Puppijet1_e3_v1_sdb1;   //!
   TBranch        *b_AK8Puppijet1_e3_v2_sdb1;   //!
   TBranch        *b_AK8Puppijet1_e4_v1_sdb1;   //!
   TBranch        *b_AK8Puppijet1_e4_v2_sdb1;   //!
   TBranch        *b_AK8Puppijet1_e2_sdb2;   //!
   TBranch        *b_AK8Puppijet1_e3_sdb2;   //!
   TBranch        *b_AK8Puppijet1_e3_v1_sdb2;   //!
   TBranch        *b_AK8Puppijet1_e3_v2_sdb2;   //!
   TBranch        *b_AK8Puppijet1_e4_v1_sdb2;   //!
   TBranch        *b_AK8Puppijet1_e4_v2_sdb2;   //!
   TBranch        *b_AK8Puppijet1_N2sdb1;   //!
   TBranch        *b_AK8Puppijet1_N2sdb2;   //!
   TBranch        *b_AK8Puppijet1_M2sdb1;   //!
   TBranch        *b_AK8Puppijet1_M2sdb2;   //!
   TBranch        *b_AK8Puppijet1_D2sdb1;   //!
   TBranch        *b_AK8Puppijet1_D2sdb2;   //!
   TBranch        *b_AK8Puppijet1_N2b1;   //!
   TBranch        *b_AK8Puppijet1_N2b2;   //!
   TBranch        *b_AK8Puppijet1_M2b1;   //!
   TBranch        *b_AK8Puppijet1_M2b2;   //!
   TBranch        *b_AK8Puppijet1_D2b1;   //!
   TBranch        *b_AK8Puppijet1_D2b2;   //!
   TBranch        *b_AK8Puppijet1_mass;   //!
   TBranch        *b_AK8Puppijet2_csv;   //!
   TBranch        *b_AK8Puppijet2_CHF;   //!
   TBranch        *b_AK8Puppijet2_NHF;   //!
   TBranch        *b_AK8Puppijet2_NEMF;   //!
   TBranch        *b_AK8Puppijet2_tau21;   //!
   TBranch        *b_AK8Puppijet2_tau32;   //!
   TBranch        *b_AK8Puppijet2_msd;   //!
   TBranch        *b_AK8Puppijet2_rho;   //!
   TBranch        *b_AK8Puppijet2_minsubcsv;   //!
   TBranch        *b_AK8Puppijet2_maxsubcsv;   //!
   TBranch        *b_AK8Puppijet2_doublecsv;   //!
   TBranch        *b_AK8Puppijet2_doublesub;   //!
   TBranch        *b_AK8Puppijet2_ptraw;   //!
   TBranch        *b_AK8Puppijet2_genpt;   //!
   TBranch        *b_AK8Puppijet2_e2_b1;   //!
   TBranch        *b_AK8Puppijet2_e3_b1;   //!
   TBranch        *b_AK8Puppijet2_e3_v1_b1;   //!
   TBranch        *b_AK8Puppijet2_e3_v2_b1;   //!
   TBranch        *b_AK8Puppijet2_e4_v1_b1;   //!
   TBranch        *b_AK8Puppijet2_e4_v2_b1;   //!
   TBranch        *b_AK8Puppijet2_e2_b2;   //!
   TBranch        *b_AK8Puppijet2_e3_b2;   //!
   TBranch        *b_AK8Puppijet2_e3_v1_b2;   //!
   TBranch        *b_AK8Puppijet2_e3_v2_b2;   //!
   TBranch        *b_AK8Puppijet2_e4_v1_b2;   //!
   TBranch        *b_AK8Puppijet2_e4_v2_b2;   //!
   TBranch        *b_AK8Puppijet2_e2_sdb1;   //!
   TBranch        *b_AK8Puppijet2_e3_sdb1;   //!
   TBranch        *b_AK8Puppijet2_e3_v1_sdb1;   //!
   TBranch        *b_AK8Puppijet2_e3_v2_sdb1;   //!
   TBranch        *b_AK8Puppijet2_e4_v1_sdb1;   //!
   TBranch        *b_AK8Puppijet2_e4_v2_sdb1;   //!
   TBranch        *b_AK8Puppijet2_e2_sdb2;   //!
   TBranch        *b_AK8Puppijet2_e3_sdb2;   //!
   TBranch        *b_AK8Puppijet2_e3_v1_sdb2;   //!
   TBranch        *b_AK8Puppijet2_e3_v2_sdb2;   //!
   TBranch        *b_AK8Puppijet2_e4_v1_sdb2;   //!
   TBranch        *b_AK8Puppijet2_e4_v2_sdb2;   //!
   TBranch        *b_AK8Puppijet2_N2sdb1;   //!
   TBranch        *b_AK8Puppijet2_N2sdb2;   //!
   TBranch        *b_AK8Puppijet2_M2sdb1;   //!
   TBranch        *b_AK8Puppijet2_M2sdb2;   //!
   TBranch        *b_AK8Puppijet2_D2sdb1;   //!
   TBranch        *b_AK8Puppijet2_D2sdb2;   //!
   TBranch        *b_AK8Puppijet2_N2b1;   //!
   TBranch        *b_AK8Puppijet2_N2b2;   //!
   TBranch        *b_AK8Puppijet2_M2b1;   //!
   TBranch        *b_AK8Puppijet2_M2b2;   //!
   TBranch        *b_AK8Puppijet2_D2b1;   //!
   TBranch        *b_AK8Puppijet2_D2b2;   //!
   TBranch        *b_nAK8Puppijets;   //!
   TBranch        *b_AK8Puppijet0_isTightVJet;   //!
   TBranch        *b_AK8Puppijet0_isHadronicV;   //!
   TBranch        *b_AK8Puppijet0_vMatching;   //!
   TBranch        *b_AK8Puppijet0_vSize;   //!
   TBranch        *b_AK8Puppijet0_partonFlavor;   //!
   TBranch        *b_AK8Puppijet0_hadronFlavor;   //!
   TBranch        *b_AK8Puppijet0_nCharged;   //!
   TBranch        *b_AK8Puppijet0_nNeutrals;   //!
   TBranch        *b_AK8Puppijet0_nParticles;   //!
   TBranch        *b_AK8Puppijet0_ratioCA15_04;   //!
   TBranch        *b_AK8CHSjet0_doublecsv;   //!
   TBranch        *b_AK8CHSjet0_doublesub;   //!
   TBranch        *b_AK8CHSjet0_pt;   //!
   TBranch        *b_AK8CHSjet0_eta;   //!
   TBranch        *b_AK8CHSjet0_phi;   //!
   TBranch        *b_AK8CHSjet1_doublecsv;   //!
   TBranch        *b_AK8CHSjet1_doublesub;   //!
   TBranch        *b_AK8CHSjet1_pt;   //!
   TBranch        *b_AK8CHSjet1_eta;   //!
   TBranch        *b_AK8CHSjet1_phi;   //!
   TBranch        *b_AK8CHSjet2_doublecsv;   //!
   TBranch        *b_AK8CHSjet2_doublesub;   //!
   TBranch        *b_AK8CHSjet2_pt;   //!
   TBranch        *b_AK8CHSjet2_eta;   //!
   TBranch        *b_AK8CHSjet2_phi;   //!
   TBranch        *b_nCA15Puppijets;   //!
   TBranch        *b_CA15CHSjet0_doublecsv;   //!
   TBranch        *b_CA15CHSjet0_doublesub;   //!
   TBranch        *b_CA15CHSjet0_pt;   //!
   TBranch        *b_CA15CHSjet0_eta;   //!
   TBranch        *b_CA15CHSjet0_phi;   //!
   TBranch        *b_CA15CHSjet1_doublecsv;   //!
   TBranch        *b_CA15CHSjet1_doublesub;   //!
   TBranch        *b_CA15CHSjet1_pt;   //!
   TBranch        *b_CA15CHSjet1_eta;   //!
   TBranch        *b_CA15CHSjet1_phi;   //!
   TBranch        *b_CA15CHSjet2_doublecsv;   //!
   TBranch        *b_CA15CHSjet2_doublesub;   //!
   TBranch        *b_CA15CHSjet2_pt;   //!
   TBranch        *b_CA15CHSjet2_eta;   //!
   TBranch        *b_CA15CHSjet2_phi;   //!
   TBranch        *b_nAK4Puppijets;   //!
   TBranch        *b_nAK4Puppijetsfwd;   //!
   TBranch        *b_nAK4PuppijetsbtagL;   //!
   TBranch        *b_nAK4PuppijetsbtagM;   //!
   TBranch        *b_nAK4PuppijetsbtagT;   //!
   TBranch        *b_nAK4PuppijetsdR08;   //!
   TBranch        *b_nAK4PuppijetsLdR08;   //!
   TBranch        *b_nAK4PuppijetsMdR08;   //!
   TBranch        *b_nAK4PuppijetsTdR08;   //!
   TBranch        *b_nAK4PuppijetsLPt100dR08;   //!
   TBranch        *b_nAK4PuppijetsMPt100dR08;   //!
   TBranch        *b_nAK4PuppijetsTPt100dR08;   //!
   TBranch        *b_nAK4PuppijetsLPt150dR08;   //!
   TBranch        *b_nAK4PuppijetsMPt150dR08;   //!
   TBranch        *b_nAK4PuppijetsTPt150dR08;   //!
   TBranch        *b_AK4Puppijet0_pt;   //!
   TBranch        *b_AK4Puppijet0_eta;   //!
   TBranch        *b_AK4Puppijet0_phi;   //!
   TBranch        *b_AK4Puppijet1_pt;   //!
   TBranch        *b_AK4Puppijet1_eta;   //!
   TBranch        *b_AK4Puppijet1_phi;   //!
   TBranch        *b_AK4Puppijet2_pt;   //!
   TBranch        *b_AK4Puppijet2_eta;   //!
   TBranch        *b_AK4Puppijet2_phi;   //!
   TBranch        *b_AK4Puppijet3_pt;   //!
   TBranch        *b_AK4Puppijet3_eta;   //!
   TBranch        *b_AK4Puppijet3_phi;   //!
   TBranch        *b_AK4Puppijet0_mass;   //!
   TBranch        *b_AK4Puppijet0_csv;   //!
   TBranch        *b_AK4Puppijet0_qgid;   //!
   TBranch        *b_AK4Puppijet0_dR08;   //!
   TBranch        *b_AK4Puppijet0_dPhi08;   //!
   TBranch        *b_AK4Puppijet1_mass;   //!
   TBranch        *b_AK4Puppijet1_csv;   //!
   TBranch        *b_AK4Puppijet1_qgid;   //!
   TBranch        *b_AK4Puppijet1_dR08;   //!
   TBranch        *b_AK4Puppijet1_dPhi08;   //!
   TBranch        *b_AK4Puppijet2_mass;   //!
   TBranch        *b_AK4Puppijet2_csv;   //!
   TBranch        *b_AK4Puppijet2_qgid;   //!
   TBranch        *b_AK4Puppijet2_dR08;   //!
   TBranch        *b_AK4Puppijet2_dPhi08;   //!
   TBranch        *b_AK4Puppijet3_mass;   //!
   TBranch        *b_AK4Puppijet3_csv;   //!
   TBranch        *b_AK4Puppijet3_qgid;   //!
   TBranch        *b_AK4Puppijet3_dR08;   //!
   TBranch        *b_AK4Puppijet3_dPhi08;   //!
   TBranch        *b_fNMuonsLoose;   //!
   TBranch        *b_fNMuonsTight;   //!
   TBranch        *b_fismu0Tight;   //!
   TBranch        *b_vmuoLoose0_pt;   //!
   TBranch        *b_vmuoLoose0_eta;   //!
   TBranch        *b_vmuoLoose0_phi;   //!
   TBranch        *b_vmuoLoose1_pt;   //!
   TBranch        *b_vmuoLoose1_eta;   //!
   TBranch        *b_vmuoLoose1_phi;   //!
   TBranch        *b_fNElectronsLoose;   //!
   TBranch        *b_fNElectronsTight;   //!
   TBranch        *b_fisele0Tight;   //!
   TBranch        *b_veleLoose0_pt;   //!
   TBranch        *b_veleLoose0_eta;   //!
   TBranch        *b_veleLoose0_phi;   //!
   TBranch        *b_veleLoose1_pt;   //!
   TBranch        *b_veleLoose1_eta;   //!
   TBranch        *b_veleLoose1_phi;   //!
   TBranch        *b_fNTaus;   //!
   TBranch        *b_fNPhotonsLoose;   //!
   TBranch        *b_fNPhotonsTight;   //!
   TBranch        *b_fispho0Tight;   //!
   TBranch        *b_vpho0_pt;   //!
   TBranch        *b_vpho0_eta;   //!
   TBranch        *b_vpho0_phi;   //!
   TBranch        *b_fBosonPt;   //!
   TBranch        *b_fBosonPhi;   //!
   TBranch        *b_fBosonMass;   //!
   TBranch        *b_fBosonEta;   //!
   TBranch        *b_fBosonPdgId;   //!
   TBranch        *b_genEleFromW;   //!
   TBranch        *b_genMuFromW;   //!
   TBranch        *b_genTauFromW;   //!

   BaconTree(TTree *tree=0);
   virtual ~BaconTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef BaconTree_cxx
BaconTree::BaconTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("DMSpin0_ggPhibb1j_150_1000pb_weighted.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("DMSpin0_ggPhibb1j_150_1000pb_weighted.root");
      }
      f->GetObject("BaconTree",tree);

   }
   Init(tree);
}

BaconTree::~BaconTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BaconTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BaconTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void BaconTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNum", &runNum, &b_fRun);
   fChain->SetBranchAddress("lumiSec", &lumiSec, &b_fLumi);
   fChain->SetBranchAddress("evtNum", &evtNum, &b_fEvtV);
   fChain->SetBranchAddress("passJson", &passJson, &b_fPassJson);
   fChain->SetBranchAddress("metfilter", &metfilter, &b_fMetFilters);
   fChain->SetBranchAddress("triggerBits", &triggerBits, &b_fITrigger);
   fChain->SetBranchAddress("selectBits", &selectBits, &b_fselectBits);
   fChain->SetBranchAddress("sf_eleTrig", &sf_eleTrig, &b_fsf_eleTrig);
   fChain->SetBranchAddress("sf_metTrig", &sf_metTrig, &b_fsf_metTrig);
   fChain->SetBranchAddress("sf_phoTrig", &sf_phoTrig, &b_fsf_phoTrig);
   fChain->SetBranchAddress("sf_lep", &sf_lep, &b_fsf_lep);
   fChain->SetBranchAddress("sf_lepTrack", &sf_lepTrack, &b_fsf_lepTrack);
   fChain->SetBranchAddress("npu", &npu, &b_fNPU);
   fChain->SetBranchAddress("npv", &npv, &b_fNVtx);
   fChain->SetBranchAddress("puWeight", &puWeight, &b_fPUWeight);
   fChain->SetBranchAddress("scale1fb", &scale1fb, &b_fScale);
   fChain->SetBranchAddress("evtWeight", &evtWeight, &b_fevtWeight);
   fChain->SetBranchAddress("kfactor", &kfactor, &b_fkfactor);
   fChain->SetBranchAddress("kfactorNLO", &kfactorNLO, &b_fkFactor_CENT);
   fChain->SetBranchAddress("PDF", &PDF, &b_fPDF);
   fChain->SetBranchAddress("PDF_UP", &PDF_UP, &b_fPDF_UP);
   fChain->SetBranchAddress("PDF_DO", &PDF_DO, &b_fPDF_DO);
   fChain->SetBranchAddress("RenScale_UP", &RenScale_UP, &b_fRenScale_UP);
   fChain->SetBranchAddress("RenScale_DO", &RenScale_DO, &b_fRenScale_DO);
   fChain->SetBranchAddress("FacScale_UP", &FacScale_UP, &b_fFacScale_UP);
   fChain->SetBranchAddress("FacScale_DO", &FacScale_DO, &b_fFacScale_DO);
   fChain->SetBranchAddress("pfmet", &pfmet, &b_fMet);
   fChain->SetBranchAddress("pfmetphi", &pfmetphi, &b_fMetPhi);
   fChain->SetBranchAddress("puppet", &puppet, &b_fPuppEt);
   fChain->SetBranchAddress("puppetphi", &puppetphi, &b_fPuppEtPhi);
   fChain->SetBranchAddress("AK8Puppijet0_pt", &AK8Puppijet0_pt, &b_AK8Puppijet0_pt);
   fChain->SetBranchAddress("AK8Puppijet0_eta", &AK8Puppijet0_eta, &b_AK8Puppijet0_eta);
   fChain->SetBranchAddress("AK8Puppijet0_phi", &AK8Puppijet0_phi, &b_AK8Puppijet0_phi);
   fChain->SetBranchAddress("AK8Puppijet1_pt", &AK8Puppijet1_pt, &b_AK8Puppijet1_pt);
   fChain->SetBranchAddress("AK8Puppijet1_eta", &AK8Puppijet1_eta, &b_AK8Puppijet1_eta);
   fChain->SetBranchAddress("AK8Puppijet1_phi", &AK8Puppijet1_phi, &b_AK8Puppijet1_phi);
   fChain->SetBranchAddress("AK8Puppijet2_pt", &AK8Puppijet2_pt, &b_AK8Puppijet2_pt);
   fChain->SetBranchAddress("AK8Puppijet2_eta", &AK8Puppijet2_eta, &b_AK8Puppijet2_eta);
   fChain->SetBranchAddress("AK8Puppijet2_phi", &AK8Puppijet2_phi, &b_AK8Puppijet2_phi);
   fChain->SetBranchAddress("AK8Puppijet0_mass", &AK8Puppijet0_mass, &b_AK8Puppijet0_mass);
   fChain->SetBranchAddress("AK8Puppijet0_csv", &AK8Puppijet0_csv, &b_AK8Puppijet0_csv);
   fChain->SetBranchAddress("AK8Puppijet0_CHF", &AK8Puppijet0_CHF, &b_AK8Puppijet0_CHF);
   fChain->SetBranchAddress("AK8Puppijet0_NHF", &AK8Puppijet0_NHF, &b_AK8Puppijet0_NHF);
   fChain->SetBranchAddress("AK8Puppijet0_NEMF", &AK8Puppijet0_NEMF, &b_AK8Puppijet0_NEMF);
   fChain->SetBranchAddress("AK8Puppijet0_tau21", &AK8Puppijet0_tau21, &b_AK8Puppijet0_tau21);
   fChain->SetBranchAddress("AK8Puppijet0_tau32", &AK8Puppijet0_tau32, &b_AK8Puppijet0_tau32);
   fChain->SetBranchAddress("AK8Puppijet0_msd", &AK8Puppijet0_msd, &b_AK8Puppijet0_msd);
   fChain->SetBranchAddress("AK8Puppijet0_rho", &AK8Puppijet0_rho, &b_AK8Puppijet0_rho);
   fChain->SetBranchAddress("AK8Puppijet0_minsubcsv", &AK8Puppijet0_minsubcsv, &b_AK8Puppijet0_minsubcsv);
   fChain->SetBranchAddress("AK8Puppijet0_maxsubcsv", &AK8Puppijet0_maxsubcsv, &b_AK8Puppijet0_maxsubcsv);
   fChain->SetBranchAddress("AK8Puppijet0_doublecsv", &AK8Puppijet0_doublecsv, &b_AK8Puppijet0_doublecsv);
   fChain->SetBranchAddress("AK8Puppijet0_doublesub", &AK8Puppijet0_doublesub, &b_AK8Puppijet0_doublesub);
   fChain->SetBranchAddress("AK8Puppijet0_ptraw", &AK8Puppijet0_ptraw, &b_AK8Puppijet0_ptraw);
   fChain->SetBranchAddress("AK8Puppijet0_genpt", &AK8Puppijet0_genpt, &b_AK8Puppijet0_genpt);
   fChain->SetBranchAddress("AK8Puppijet0_N2sdb1", &AK8Puppijet0_N2sdb1, &b_AK8Puppijet0_N2sdb1);
   fChain->SetBranchAddress("AK8Puppijet0_N2sdb2", &AK8Puppijet0_N2sdb2, &b_AK8Puppijet0_N2sdb2);
   fChain->SetBranchAddress("AK8Puppijet0_M2sdb1", &AK8Puppijet0_M2sdb1, &b_AK8Puppijet0_M2sdb1);
   fChain->SetBranchAddress("AK8Puppijet0_M2sdb2", &AK8Puppijet0_M2sdb2, &b_AK8Puppijet0_M2sdb2);
   fChain->SetBranchAddress("AK8Puppijet0_D2sdb1", &AK8Puppijet0_D2sdb1, &b_AK8Puppijet0_D2sdb1);
   fChain->SetBranchAddress("AK8Puppijet0_D2sdb2", &AK8Puppijet0_D2sdb2, &b_AK8Puppijet0_D2sdb2);
   fChain->SetBranchAddress("AK8Puppijet0_N2b1", &AK8Puppijet0_N2b1, &b_AK8Puppijet0_N2b1);
   fChain->SetBranchAddress("AK8Puppijet0_N2b2", &AK8Puppijet0_N2b2, &b_AK8Puppijet0_N2b2);
   fChain->SetBranchAddress("AK8Puppijet0_M2b1", &AK8Puppijet0_M2b1, &b_AK8Puppijet0_M2b1);
   fChain->SetBranchAddress("AK8Puppijet0_M2b2", &AK8Puppijet0_M2b2, &b_AK8Puppijet0_M2b2);
   fChain->SetBranchAddress("AK8Puppijet0_D2b1", &AK8Puppijet0_D2b1, &b_AK8Puppijet0_D2b1);
   fChain->SetBranchAddress("AK8Puppijet0_D2b2", &AK8Puppijet0_D2b2, &b_AK8Puppijet0_D2b2);
//    fChain->SetBranchAddress("AK8Puppijet0_mass", &AK8Puppijet0_mass, &b_AK8Puppijet0_mass);
   fChain->SetBranchAddress("AK8Puppijet1_csv", &AK8Puppijet1_csv, &b_AK8Puppijet1_csv);
   fChain->SetBranchAddress("AK8Puppijet1_CHF", &AK8Puppijet1_CHF, &b_AK8Puppijet1_CHF);
   fChain->SetBranchAddress("AK8Puppijet1_NHF", &AK8Puppijet1_NHF, &b_AK8Puppijet1_NHF);
   fChain->SetBranchAddress("AK8Puppijet1_NEMF", &AK8Puppijet1_NEMF, &b_AK8Puppijet1_NEMF);
   fChain->SetBranchAddress("AK8Puppijet1_tau21", &AK8Puppijet1_tau21, &b_AK8Puppijet1_tau21);
   fChain->SetBranchAddress("AK8Puppijet1_tau32", &AK8Puppijet1_tau32, &b_AK8Puppijet1_tau32);
   fChain->SetBranchAddress("AK8Puppijet1_msd", &AK8Puppijet1_msd, &b_AK8Puppijet1_msd);
   fChain->SetBranchAddress("AK8Puppijet1_rho", &AK8Puppijet1_rho, &b_AK8Puppijet1_rho);
   fChain->SetBranchAddress("AK8Puppijet1_minsubcsv", &AK8Puppijet1_minsubcsv, &b_AK8Puppijet1_minsubcsv);
   fChain->SetBranchAddress("AK8Puppijet1_maxsubcsv", &AK8Puppijet1_maxsubcsv, &b_AK8Puppijet1_maxsubcsv);
   fChain->SetBranchAddress("AK8Puppijet1_doublecsv", &AK8Puppijet1_doublecsv, &b_AK8Puppijet1_doublecsv);
   fChain->SetBranchAddress("AK8Puppijet1_doublesub", &AK8Puppijet1_doublesub, &b_AK8Puppijet1_doublesub);
   fChain->SetBranchAddress("AK8Puppijet1_ptraw", &AK8Puppijet1_ptraw, &b_AK8Puppijet1_ptraw);
   fChain->SetBranchAddress("AK8Puppijet1_genpt", &AK8Puppijet1_genpt, &b_AK8Puppijet1_genpt);
   fChain->SetBranchAddress("AK8Puppijet1_e2_b1", &AK8Puppijet1_e2_b1, &b_AK8Puppijet1_e2_b1);
   fChain->SetBranchAddress("AK8Puppijet1_e3_b1", &AK8Puppijet1_e3_b1, &b_AK8Puppijet1_e3_b1);
   fChain->SetBranchAddress("AK8Puppijet1_e3_v1_b1", &AK8Puppijet1_e3_v1_b1, &b_AK8Puppijet1_e3_v1_b1);
   fChain->SetBranchAddress("AK8Puppijet1_e3_v2_b1", &AK8Puppijet1_e3_v2_b1, &b_AK8Puppijet1_e3_v2_b1);
   fChain->SetBranchAddress("AK8Puppijet1_e4_v1_b1", &AK8Puppijet1_e4_v1_b1, &b_AK8Puppijet1_e4_v1_b1);
   fChain->SetBranchAddress("AK8Puppijet1_e4_v2_b1", &AK8Puppijet1_e4_v2_b1, &b_AK8Puppijet1_e4_v2_b1);
   fChain->SetBranchAddress("AK8Puppijet1_e2_b2", &AK8Puppijet1_e2_b2, &b_AK8Puppijet1_e2_b2);
   fChain->SetBranchAddress("AK8Puppijet1_e3_b2", &AK8Puppijet1_e3_b2, &b_AK8Puppijet1_e3_b2);
   fChain->SetBranchAddress("AK8Puppijet1_e3_v1_b2", &AK8Puppijet1_e3_v1_b2, &b_AK8Puppijet1_e3_v1_b2);
   fChain->SetBranchAddress("AK8Puppijet1_e3_v2_b2", &AK8Puppijet1_e3_v2_b2, &b_AK8Puppijet1_e3_v2_b2);
   fChain->SetBranchAddress("AK8Puppijet1_e4_v1_b2", &AK8Puppijet1_e4_v1_b2, &b_AK8Puppijet1_e4_v1_b2);
   fChain->SetBranchAddress("AK8Puppijet1_e4_v2_b2", &AK8Puppijet1_e4_v2_b2, &b_AK8Puppijet1_e4_v2_b2);
   fChain->SetBranchAddress("AK8Puppijet1_e2_sdb1", &AK8Puppijet1_e2_sdb1, &b_AK8Puppijet1_e2_sdb1);
   fChain->SetBranchAddress("AK8Puppijet1_e3_sdb1", &AK8Puppijet1_e3_sdb1, &b_AK8Puppijet1_e3_sdb1);
   fChain->SetBranchAddress("AK8Puppijet1_e3_v1_sdb1", &AK8Puppijet1_e3_v1_sdb1, &b_AK8Puppijet1_e3_v1_sdb1);
   fChain->SetBranchAddress("AK8Puppijet1_e3_v2_sdb1", &AK8Puppijet1_e3_v2_sdb1, &b_AK8Puppijet1_e3_v2_sdb1);
   fChain->SetBranchAddress("AK8Puppijet1_e4_v1_sdb1", &AK8Puppijet1_e4_v1_sdb1, &b_AK8Puppijet1_e4_v1_sdb1);
   fChain->SetBranchAddress("AK8Puppijet1_e4_v2_sdb1", &AK8Puppijet1_e4_v2_sdb1, &b_AK8Puppijet1_e4_v2_sdb1);
   fChain->SetBranchAddress("AK8Puppijet1_e2_sdb2", &AK8Puppijet1_e2_sdb2, &b_AK8Puppijet1_e2_sdb2);
   fChain->SetBranchAddress("AK8Puppijet1_e3_sdb2", &AK8Puppijet1_e3_sdb2, &b_AK8Puppijet1_e3_sdb2);
   fChain->SetBranchAddress("AK8Puppijet1_e3_v1_sdb2", &AK8Puppijet1_e3_v1_sdb2, &b_AK8Puppijet1_e3_v1_sdb2);
   fChain->SetBranchAddress("AK8Puppijet1_e3_v2_sdb2", &AK8Puppijet1_e3_v2_sdb2, &b_AK8Puppijet1_e3_v2_sdb2);
   fChain->SetBranchAddress("AK8Puppijet1_e4_v1_sdb2", &AK8Puppijet1_e4_v1_sdb2, &b_AK8Puppijet1_e4_v1_sdb2);
   fChain->SetBranchAddress("AK8Puppijet1_e4_v2_sdb2", &AK8Puppijet1_e4_v2_sdb2, &b_AK8Puppijet1_e4_v2_sdb2);
   fChain->SetBranchAddress("AK8Puppijet1_N2sdb1", &AK8Puppijet1_N2sdb1, &b_AK8Puppijet1_N2sdb1);
   fChain->SetBranchAddress("AK8Puppijet1_N2sdb2", &AK8Puppijet1_N2sdb2, &b_AK8Puppijet1_N2sdb2);
   fChain->SetBranchAddress("AK8Puppijet1_M2sdb1", &AK8Puppijet1_M2sdb1, &b_AK8Puppijet1_M2sdb1);
   fChain->SetBranchAddress("AK8Puppijet1_M2sdb2", &AK8Puppijet1_M2sdb2, &b_AK8Puppijet1_M2sdb2);
   fChain->SetBranchAddress("AK8Puppijet1_D2sdb1", &AK8Puppijet1_D2sdb1, &b_AK8Puppijet1_D2sdb1);
   fChain->SetBranchAddress("AK8Puppijet1_D2sdb2", &AK8Puppijet1_D2sdb2, &b_AK8Puppijet1_D2sdb2);
   fChain->SetBranchAddress("AK8Puppijet1_N2b1", &AK8Puppijet1_N2b1, &b_AK8Puppijet1_N2b1);
   fChain->SetBranchAddress("AK8Puppijet1_N2b2", &AK8Puppijet1_N2b2, &b_AK8Puppijet1_N2b2);
   fChain->SetBranchAddress("AK8Puppijet1_M2b1", &AK8Puppijet1_M2b1, &b_AK8Puppijet1_M2b1);
   fChain->SetBranchAddress("AK8Puppijet1_M2b2", &AK8Puppijet1_M2b2, &b_AK8Puppijet1_M2b2);
   fChain->SetBranchAddress("AK8Puppijet1_D2b1", &AK8Puppijet1_D2b1, &b_AK8Puppijet1_D2b1);
   fChain->SetBranchAddress("AK8Puppijet1_D2b2", &AK8Puppijet1_D2b2, &b_AK8Puppijet1_D2b2);
   fChain->SetBranchAddress("AK8Puppijet1_mass", &AK8Puppijet1_mass, &b_AK8Puppijet1_mass);
   fChain->SetBranchAddress("AK8Puppijet2_csv", &AK8Puppijet2_csv, &b_AK8Puppijet2_csv);
   fChain->SetBranchAddress("AK8Puppijet2_CHF", &AK8Puppijet2_CHF, &b_AK8Puppijet2_CHF);
   fChain->SetBranchAddress("AK8Puppijet2_NHF", &AK8Puppijet2_NHF, &b_AK8Puppijet2_NHF);
   fChain->SetBranchAddress("AK8Puppijet2_NEMF", &AK8Puppijet2_NEMF, &b_AK8Puppijet2_NEMF);
   fChain->SetBranchAddress("AK8Puppijet2_tau21", &AK8Puppijet2_tau21, &b_AK8Puppijet2_tau21);
   fChain->SetBranchAddress("AK8Puppijet2_tau32", &AK8Puppijet2_tau32, &b_AK8Puppijet2_tau32);
   fChain->SetBranchAddress("AK8Puppijet2_msd", &AK8Puppijet2_msd, &b_AK8Puppijet2_msd);
   fChain->SetBranchAddress("AK8Puppijet2_rho", &AK8Puppijet2_rho, &b_AK8Puppijet2_rho);
   fChain->SetBranchAddress("AK8Puppijet2_minsubcsv", &AK8Puppijet2_minsubcsv, &b_AK8Puppijet2_minsubcsv);
   fChain->SetBranchAddress("AK8Puppijet2_maxsubcsv", &AK8Puppijet2_maxsubcsv, &b_AK8Puppijet2_maxsubcsv);
   fChain->SetBranchAddress("AK8Puppijet2_doublecsv", &AK8Puppijet2_doublecsv, &b_AK8Puppijet2_doublecsv);
   fChain->SetBranchAddress("AK8Puppijet2_doublesub", &AK8Puppijet2_doublesub, &b_AK8Puppijet2_doublesub);
   fChain->SetBranchAddress("AK8Puppijet2_ptraw", &AK8Puppijet2_ptraw, &b_AK8Puppijet2_ptraw);
   fChain->SetBranchAddress("AK8Puppijet2_genpt", &AK8Puppijet2_genpt, &b_AK8Puppijet2_genpt);
   fChain->SetBranchAddress("AK8Puppijet2_e2_b1", &AK8Puppijet2_e2_b1, &b_AK8Puppijet2_e2_b1);
   fChain->SetBranchAddress("AK8Puppijet2_e3_b1", &AK8Puppijet2_e3_b1, &b_AK8Puppijet2_e3_b1);
   fChain->SetBranchAddress("AK8Puppijet2_e3_v1_b1", &AK8Puppijet2_e3_v1_b1, &b_AK8Puppijet2_e3_v1_b1);
   fChain->SetBranchAddress("AK8Puppijet2_e3_v2_b1", &AK8Puppijet2_e3_v2_b1, &b_AK8Puppijet2_e3_v2_b1);
   fChain->SetBranchAddress("AK8Puppijet2_e4_v1_b1", &AK8Puppijet2_e4_v1_b1, &b_AK8Puppijet2_e4_v1_b1);
   fChain->SetBranchAddress("AK8Puppijet2_e4_v2_b1", &AK8Puppijet2_e4_v2_b1, &b_AK8Puppijet2_e4_v2_b1);
   fChain->SetBranchAddress("AK8Puppijet2_e2_b2", &AK8Puppijet2_e2_b2, &b_AK8Puppijet2_e2_b2);
   fChain->SetBranchAddress("AK8Puppijet2_e3_b2", &AK8Puppijet2_e3_b2, &b_AK8Puppijet2_e3_b2);
   fChain->SetBranchAddress("AK8Puppijet2_e3_v1_b2", &AK8Puppijet2_e3_v1_b2, &b_AK8Puppijet2_e3_v1_b2);
   fChain->SetBranchAddress("AK8Puppijet2_e3_v2_b2", &AK8Puppijet2_e3_v2_b2, &b_AK8Puppijet2_e3_v2_b2);
   fChain->SetBranchAddress("AK8Puppijet2_e4_v1_b2", &AK8Puppijet2_e4_v1_b2, &b_AK8Puppijet2_e4_v1_b2);
   fChain->SetBranchAddress("AK8Puppijet2_e4_v2_b2", &AK8Puppijet2_e4_v2_b2, &b_AK8Puppijet2_e4_v2_b2);
   fChain->SetBranchAddress("AK8Puppijet2_e2_sdb1", &AK8Puppijet2_e2_sdb1, &b_AK8Puppijet2_e2_sdb1);
   fChain->SetBranchAddress("AK8Puppijet2_e3_sdb1", &AK8Puppijet2_e3_sdb1, &b_AK8Puppijet2_e3_sdb1);
   fChain->SetBranchAddress("AK8Puppijet2_e3_v1_sdb1", &AK8Puppijet2_e3_v1_sdb1, &b_AK8Puppijet2_e3_v1_sdb1);
   fChain->SetBranchAddress("AK8Puppijet2_e3_v2_sdb1", &AK8Puppijet2_e3_v2_sdb1, &b_AK8Puppijet2_e3_v2_sdb1);
   fChain->SetBranchAddress("AK8Puppijet2_e4_v1_sdb1", &AK8Puppijet2_e4_v1_sdb1, &b_AK8Puppijet2_e4_v1_sdb1);
   fChain->SetBranchAddress("AK8Puppijet2_e4_v2_sdb1", &AK8Puppijet2_e4_v2_sdb1, &b_AK8Puppijet2_e4_v2_sdb1);
   fChain->SetBranchAddress("AK8Puppijet2_e2_sdb2", &AK8Puppijet2_e2_sdb2, &b_AK8Puppijet2_e2_sdb2);
   fChain->SetBranchAddress("AK8Puppijet2_e3_sdb2", &AK8Puppijet2_e3_sdb2, &b_AK8Puppijet2_e3_sdb2);
   fChain->SetBranchAddress("AK8Puppijet2_e3_v1_sdb2", &AK8Puppijet2_e3_v1_sdb2, &b_AK8Puppijet2_e3_v1_sdb2);
   fChain->SetBranchAddress("AK8Puppijet2_e3_v2_sdb2", &AK8Puppijet2_e3_v2_sdb2, &b_AK8Puppijet2_e3_v2_sdb2);
   fChain->SetBranchAddress("AK8Puppijet2_e4_v1_sdb2", &AK8Puppijet2_e4_v1_sdb2, &b_AK8Puppijet2_e4_v1_sdb2);
   fChain->SetBranchAddress("AK8Puppijet2_e4_v2_sdb2", &AK8Puppijet2_e4_v2_sdb2, &b_AK8Puppijet2_e4_v2_sdb2);
   fChain->SetBranchAddress("AK8Puppijet2_N2sdb1", &AK8Puppijet2_N2sdb1, &b_AK8Puppijet2_N2sdb1);
   fChain->SetBranchAddress("AK8Puppijet2_N2sdb2", &AK8Puppijet2_N2sdb2, &b_AK8Puppijet2_N2sdb2);
   fChain->SetBranchAddress("AK8Puppijet2_M2sdb1", &AK8Puppijet2_M2sdb1, &b_AK8Puppijet2_M2sdb1);
   fChain->SetBranchAddress("AK8Puppijet2_M2sdb2", &AK8Puppijet2_M2sdb2, &b_AK8Puppijet2_M2sdb2);
   fChain->SetBranchAddress("AK8Puppijet2_D2sdb1", &AK8Puppijet2_D2sdb1, &b_AK8Puppijet2_D2sdb1);
   fChain->SetBranchAddress("AK8Puppijet2_D2sdb2", &AK8Puppijet2_D2sdb2, &b_AK8Puppijet2_D2sdb2);
   fChain->SetBranchAddress("AK8Puppijet2_N2b1", &AK8Puppijet2_N2b1, &b_AK8Puppijet2_N2b1);
   fChain->SetBranchAddress("AK8Puppijet2_N2b2", &AK8Puppijet2_N2b2, &b_AK8Puppijet2_N2b2);
   fChain->SetBranchAddress("AK8Puppijet2_M2b1", &AK8Puppijet2_M2b1, &b_AK8Puppijet2_M2b1);
   fChain->SetBranchAddress("AK8Puppijet2_M2b2", &AK8Puppijet2_M2b2, &b_AK8Puppijet2_M2b2);
   fChain->SetBranchAddress("AK8Puppijet2_D2b1", &AK8Puppijet2_D2b1, &b_AK8Puppijet2_D2b1);
   fChain->SetBranchAddress("AK8Puppijet2_D2b2", &AK8Puppijet2_D2b2, &b_AK8Puppijet2_D2b2);
   fChain->SetBranchAddress("nAK8Puppijets", &nAK8Puppijets, &b_nAK8Puppijets);
   fChain->SetBranchAddress("AK8Puppijet0_isTightVJet", &AK8Puppijet0_isTightVJet, &b_AK8Puppijet0_isTightVJet);
   fChain->SetBranchAddress("AK8Puppijet0_isHadronicV", &AK8Puppijet0_isHadronicV, &b_AK8Puppijet0_isHadronicV);
   fChain->SetBranchAddress("AK8Puppijet0_vMatching", &AK8Puppijet0_vMatching, &b_AK8Puppijet0_vMatching);
   fChain->SetBranchAddress("AK8Puppijet0_vSize", &AK8Puppijet0_vSize, &b_AK8Puppijet0_vSize);
   fChain->SetBranchAddress("AK8Puppijet0_partonFlavor", &AK8Puppijet0_partonFlavor, &b_AK8Puppijet0_partonFlavor);
   fChain->SetBranchAddress("AK8Puppijet0_hadronFlavor", &AK8Puppijet0_hadronFlavor, &b_AK8Puppijet0_hadronFlavor);
   fChain->SetBranchAddress("AK8Puppijet0_nCharged", &AK8Puppijet0_nCharged, &b_AK8Puppijet0_nCharged);
   fChain->SetBranchAddress("AK8Puppijet0_nNeutrals", &AK8Puppijet0_nNeutrals, &b_AK8Puppijet0_nNeutrals);
   fChain->SetBranchAddress("AK8Puppijet0_nParticles", &AK8Puppijet0_nParticles, &b_AK8Puppijet0_nParticles);
   fChain->SetBranchAddress("AK8Puppijet0_ratioCA15_04", &AK8Puppijet0_ratioCA15_04, &b_AK8Puppijet0_ratioCA15_04);
   fChain->SetBranchAddress("AK8CHSjet0_doublecsv", &AK8CHSjet0_doublecsv, &b_AK8CHSjet0_doublecsv);
   fChain->SetBranchAddress("AK8CHSjet0_doublesub", &AK8CHSjet0_doublesub, &b_AK8CHSjet0_doublesub);
   fChain->SetBranchAddress("AK8CHSjet0_pt", &AK8CHSjet0_pt, &b_AK8CHSjet0_pt);
   fChain->SetBranchAddress("AK8CHSjet0_eta", &AK8CHSjet0_eta, &b_AK8CHSjet0_eta);
   fChain->SetBranchAddress("AK8CHSjet0_phi", &AK8CHSjet0_phi, &b_AK8CHSjet0_phi);
   fChain->SetBranchAddress("AK8CHSjet1_doublecsv", &AK8CHSjet1_doublecsv, &b_AK8CHSjet1_doublecsv);
   fChain->SetBranchAddress("AK8CHSjet1_doublesub", &AK8CHSjet1_doublesub, &b_AK8CHSjet1_doublesub);
   fChain->SetBranchAddress("AK8CHSjet1_pt", &AK8CHSjet1_pt, &b_AK8CHSjet1_pt);
   fChain->SetBranchAddress("AK8CHSjet1_eta", &AK8CHSjet1_eta, &b_AK8CHSjet1_eta);
   fChain->SetBranchAddress("AK8CHSjet1_phi", &AK8CHSjet1_phi, &b_AK8CHSjet1_phi);
   fChain->SetBranchAddress("AK8CHSjet2_doublecsv", &AK8CHSjet2_doublecsv, &b_AK8CHSjet2_doublecsv);
   fChain->SetBranchAddress("AK8CHSjet2_doublesub", &AK8CHSjet2_doublesub, &b_AK8CHSjet2_doublesub);
   fChain->SetBranchAddress("AK8CHSjet2_pt", &AK8CHSjet2_pt, &b_AK8CHSjet2_pt);
   fChain->SetBranchAddress("AK8CHSjet2_eta", &AK8CHSjet2_eta, &b_AK8CHSjet2_eta);
   fChain->SetBranchAddress("AK8CHSjet2_phi", &AK8CHSjet2_phi, &b_AK8CHSjet2_phi);
   fChain->SetBranchAddress("nCA15Puppijets", &nCA15Puppijets, &b_nCA15Puppijets);
   fChain->SetBranchAddress("CA15CHSjet0_doublecsv", &CA15CHSjet0_doublecsv, &b_CA15CHSjet0_doublecsv);
   fChain->SetBranchAddress("CA15CHSjet0_doublesub", &CA15CHSjet0_doublesub, &b_CA15CHSjet0_doublesub);
   fChain->SetBranchAddress("CA15CHSjet0_pt", &CA15CHSjet0_pt, &b_CA15CHSjet0_pt);
   fChain->SetBranchAddress("CA15CHSjet0_eta", &CA15CHSjet0_eta, &b_CA15CHSjet0_eta);
   fChain->SetBranchAddress("CA15CHSjet0_phi", &CA15CHSjet0_phi, &b_CA15CHSjet0_phi);
   fChain->SetBranchAddress("CA15CHSjet1_doublecsv", &CA15CHSjet1_doublecsv, &b_CA15CHSjet1_doublecsv);
   fChain->SetBranchAddress("CA15CHSjet1_doublesub", &CA15CHSjet1_doublesub, &b_CA15CHSjet1_doublesub);
   fChain->SetBranchAddress("CA15CHSjet1_pt", &CA15CHSjet1_pt, &b_CA15CHSjet1_pt);
   fChain->SetBranchAddress("CA15CHSjet1_eta", &CA15CHSjet1_eta, &b_CA15CHSjet1_eta);
   fChain->SetBranchAddress("CA15CHSjet1_phi", &CA15CHSjet1_phi, &b_CA15CHSjet1_phi);
   fChain->SetBranchAddress("CA15CHSjet2_doublecsv", &CA15CHSjet2_doublecsv, &b_CA15CHSjet2_doublecsv);
   fChain->SetBranchAddress("CA15CHSjet2_doublesub", &CA15CHSjet2_doublesub, &b_CA15CHSjet2_doublesub);
   fChain->SetBranchAddress("CA15CHSjet2_pt", &CA15CHSjet2_pt, &b_CA15CHSjet2_pt);
   fChain->SetBranchAddress("CA15CHSjet2_eta", &CA15CHSjet2_eta, &b_CA15CHSjet2_eta);
   fChain->SetBranchAddress("CA15CHSjet2_phi", &CA15CHSjet2_phi, &b_CA15CHSjet2_phi);
   fChain->SetBranchAddress("nAK4Puppijets", &nAK4Puppijets, &b_nAK4Puppijets);
   fChain->SetBranchAddress("nAK4Puppijetsfwd", &nAK4Puppijetsfwd, &b_nAK4Puppijetsfwd);
   fChain->SetBranchAddress("nAK4PuppijetsbtagL", &nAK4PuppijetsbtagL, &b_nAK4PuppijetsbtagL);
   fChain->SetBranchAddress("nAK4PuppijetsbtagM", &nAK4PuppijetsbtagM, &b_nAK4PuppijetsbtagM);
   fChain->SetBranchAddress("nAK4PuppijetsbtagT", &nAK4PuppijetsbtagT, &b_nAK4PuppijetsbtagT);
   fChain->SetBranchAddress("nAK4PuppijetsdR08", &nAK4PuppijetsdR08, &b_nAK4PuppijetsdR08);
   fChain->SetBranchAddress("nAK4PuppijetsLdR08", &nAK4PuppijetsLdR08, &b_nAK4PuppijetsLdR08);
   fChain->SetBranchAddress("nAK4PuppijetsMdR08", &nAK4PuppijetsMdR08, &b_nAK4PuppijetsMdR08);
   fChain->SetBranchAddress("nAK4PuppijetsTdR08", &nAK4PuppijetsTdR08, &b_nAK4PuppijetsTdR08);
   fChain->SetBranchAddress("nAK4PuppijetsLPt100dR08", &nAK4PuppijetsLPt100dR08, &b_nAK4PuppijetsLPt100dR08);
   fChain->SetBranchAddress("nAK4PuppijetsMPt100dR08", &nAK4PuppijetsMPt100dR08, &b_nAK4PuppijetsMPt100dR08);
   fChain->SetBranchAddress("nAK4PuppijetsTPt100dR08", &nAK4PuppijetsTPt100dR08, &b_nAK4PuppijetsTPt100dR08);
   fChain->SetBranchAddress("nAK4PuppijetsLPt150dR08", &nAK4PuppijetsLPt150dR08, &b_nAK4PuppijetsLPt150dR08);
   fChain->SetBranchAddress("nAK4PuppijetsMPt150dR08", &nAK4PuppijetsMPt150dR08, &b_nAK4PuppijetsMPt150dR08);
   fChain->SetBranchAddress("nAK4PuppijetsTPt150dR08", &nAK4PuppijetsTPt150dR08, &b_nAK4PuppijetsTPt150dR08);
   fChain->SetBranchAddress("AK4Puppijet0_pt", &AK4Puppijet0_pt, &b_AK4Puppijet0_pt);
   fChain->SetBranchAddress("AK4Puppijet0_eta", &AK4Puppijet0_eta, &b_AK4Puppijet0_eta);
   fChain->SetBranchAddress("AK4Puppijet0_phi", &AK4Puppijet0_phi, &b_AK4Puppijet0_phi);
   fChain->SetBranchAddress("AK4Puppijet1_pt", &AK4Puppijet1_pt, &b_AK4Puppijet1_pt);
   fChain->SetBranchAddress("AK4Puppijet1_eta", &AK4Puppijet1_eta, &b_AK4Puppijet1_eta);
   fChain->SetBranchAddress("AK4Puppijet1_phi", &AK4Puppijet1_phi, &b_AK4Puppijet1_phi);
   fChain->SetBranchAddress("AK4Puppijet2_pt", &AK4Puppijet2_pt, &b_AK4Puppijet2_pt);
   fChain->SetBranchAddress("AK4Puppijet2_eta", &AK4Puppijet2_eta, &b_AK4Puppijet2_eta);
   fChain->SetBranchAddress("AK4Puppijet2_phi", &AK4Puppijet2_phi, &b_AK4Puppijet2_phi);
   fChain->SetBranchAddress("AK4Puppijet3_pt", &AK4Puppijet3_pt, &b_AK4Puppijet3_pt);
   fChain->SetBranchAddress("AK4Puppijet3_eta", &AK4Puppijet3_eta, &b_AK4Puppijet3_eta);
   fChain->SetBranchAddress("AK4Puppijet3_phi", &AK4Puppijet3_phi, &b_AK4Puppijet3_phi);
   fChain->SetBranchAddress("AK4Puppijet0_mass", &AK4Puppijet0_mass, &b_AK4Puppijet0_mass);
   fChain->SetBranchAddress("AK4Puppijet0_csv", &AK4Puppijet0_csv, &b_AK4Puppijet0_csv);
   fChain->SetBranchAddress("AK4Puppijet0_qgid", &AK4Puppijet0_qgid, &b_AK4Puppijet0_qgid);
   fChain->SetBranchAddress("AK4Puppijet0_dR08", &AK4Puppijet0_dR08, &b_AK4Puppijet0_dR08);
   fChain->SetBranchAddress("AK4Puppijet0_dPhi08", &AK4Puppijet0_dPhi08, &b_AK4Puppijet0_dPhi08);
   fChain->SetBranchAddress("AK4Puppijet1_mass", &AK4Puppijet1_mass, &b_AK4Puppijet1_mass);
   fChain->SetBranchAddress("AK4Puppijet1_csv", &AK4Puppijet1_csv, &b_AK4Puppijet1_csv);
   fChain->SetBranchAddress("AK4Puppijet1_qgid", &AK4Puppijet1_qgid, &b_AK4Puppijet1_qgid);
   fChain->SetBranchAddress("AK4Puppijet1_dR08", &AK4Puppijet1_dR08, &b_AK4Puppijet1_dR08);
   fChain->SetBranchAddress("AK4Puppijet1_dPhi08", &AK4Puppijet1_dPhi08, &b_AK4Puppijet1_dPhi08);
   fChain->SetBranchAddress("AK4Puppijet2_mass", &AK4Puppijet2_mass, &b_AK4Puppijet2_mass);
   fChain->SetBranchAddress("AK4Puppijet2_csv", &AK4Puppijet2_csv, &b_AK4Puppijet2_csv);
   fChain->SetBranchAddress("AK4Puppijet2_qgid", &AK4Puppijet2_qgid, &b_AK4Puppijet2_qgid);
   fChain->SetBranchAddress("AK4Puppijet2_dR08", &AK4Puppijet2_dR08, &b_AK4Puppijet2_dR08);
   fChain->SetBranchAddress("AK4Puppijet2_dPhi08", &AK4Puppijet2_dPhi08, &b_AK4Puppijet2_dPhi08);
   fChain->SetBranchAddress("AK4Puppijet3_mass", &AK4Puppijet3_mass, &b_AK4Puppijet3_mass);
   fChain->SetBranchAddress("AK4Puppijet3_csv", &AK4Puppijet3_csv, &b_AK4Puppijet3_csv);
   fChain->SetBranchAddress("AK4Puppijet3_qgid", &AK4Puppijet3_qgid, &b_AK4Puppijet3_qgid);
   fChain->SetBranchAddress("AK4Puppijet3_dR08", &AK4Puppijet3_dR08, &b_AK4Puppijet3_dR08);
   fChain->SetBranchAddress("AK4Puppijet3_dPhi08", &AK4Puppijet3_dPhi08, &b_AK4Puppijet3_dPhi08);
   fChain->SetBranchAddress("nmuLoose", &nmuLoose, &b_fNMuonsLoose);
   fChain->SetBranchAddress("nmuTight", &nmuTight, &b_fNMuonsTight);
   fChain->SetBranchAddress("ismu0Tight", &ismu0Tight, &b_fismu0Tight);
   fChain->SetBranchAddress("vmuoLoose0_pt", &vmuoLoose0_pt, &b_vmuoLoose0_pt);
   fChain->SetBranchAddress("vmuoLoose0_eta", &vmuoLoose0_eta, &b_vmuoLoose0_eta);
   fChain->SetBranchAddress("vmuoLoose0_phi", &vmuoLoose0_phi, &b_vmuoLoose0_phi);
   fChain->SetBranchAddress("vmuoLoose1_pt", &vmuoLoose1_pt, &b_vmuoLoose1_pt);
   fChain->SetBranchAddress("vmuoLoose1_eta", &vmuoLoose1_eta, &b_vmuoLoose1_eta);
   fChain->SetBranchAddress("vmuoLoose1_phi", &vmuoLoose1_phi, &b_vmuoLoose1_phi);
   fChain->SetBranchAddress("neleLoose", &neleLoose, &b_fNElectronsLoose);
   fChain->SetBranchAddress("neleTight", &neleTight, &b_fNElectronsTight);
   fChain->SetBranchAddress("isele0Tight", &isele0Tight, &b_fisele0Tight);
   fChain->SetBranchAddress("veleLoose0_pt", &veleLoose0_pt, &b_veleLoose0_pt);
   fChain->SetBranchAddress("veleLoose0_eta", &veleLoose0_eta, &b_veleLoose0_eta);
   fChain->SetBranchAddress("veleLoose0_phi", &veleLoose0_phi, &b_veleLoose0_phi);
   fChain->SetBranchAddress("veleLoose1_pt", &veleLoose1_pt, &b_veleLoose1_pt);
   fChain->SetBranchAddress("veleLoose1_eta", &veleLoose1_eta, &b_veleLoose1_eta);
   fChain->SetBranchAddress("veleLoose1_phi", &veleLoose1_phi, &b_veleLoose1_phi);
   fChain->SetBranchAddress("ntau", &ntau, &b_fNTaus);
   fChain->SetBranchAddress("nphoLoose", &nphoLoose, &b_fNPhotonsLoose);
   fChain->SetBranchAddress("nphoTight", &nphoTight, &b_fNPhotonsTight);
   fChain->SetBranchAddress("ispho0Tight", &ispho0Tight, &b_fispho0Tight);
   fChain->SetBranchAddress("vpho0_pt", &vpho0_pt, &b_vpho0_pt);
   fChain->SetBranchAddress("vpho0_eta", &vpho0_eta, &b_vpho0_eta);
   fChain->SetBranchAddress("vpho0_phi", &vpho0_phi, &b_vpho0_phi);
   fChain->SetBranchAddress("genVPt", &genVPt, &b_fBosonPt);
   fChain->SetBranchAddress("genVPhi", &genVPhi, &b_fBosonPhi);
   fChain->SetBranchAddress("genVMass", &genVMass, &b_fBosonMass);
   fChain->SetBranchAddress("genVEta", &genVEta, &b_fBosonEta);
   fChain->SetBranchAddress("genVPdfId", &genVPdfId, &b_fBosonPdgId);
   fChain->SetBranchAddress("genEleFromW", &genEleFromW, &b_genEleFromW);
   fChain->SetBranchAddress("genMuFromW", &genMuFromW, &b_genMuFromW);
   fChain->SetBranchAddress("genTauFromW", &genTauFromW, &b_genTauFromW);
   Notify();
}

Bool_t BaconTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BaconTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BaconTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BaconTree_cxx
