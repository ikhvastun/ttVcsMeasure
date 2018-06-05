#ifndef readTreeSync_H
#define readTreeSync_H

#include <iostream>
#include <fstream>

const int nFlavors = 2;

const int leptonSelectionAnalysis = 3;

const int nSamples = 20;
const int dataSample = 0;

const int runB = 0;
const int runCDE = 1;
const int runF = 2;

TString flavorsString[2] = {"el", "mu"};
TString additionalString[2] = {"_NC", ""};


struct BinLabelOptions{
  int index;
  std::string labelSR;
};

std::vector<BinLabelOptions> theSRLabelOptionsFor2L = {

      
      {1, "2j"},
      {2, "3j1b"},
      {3, "3j>1b"},
      {4, ">3j1b"},
      {5, ">3j>1b"},
      {6, "2j"},
      {7, "3j1b"},
      {8, "3j>1b"},
      {9, ">3j1b"},
      {10, ">3j>1b"},
      {11, "2j"},
      {12, "3j1b"},
      {13, "3j>1b"},
      {14, ">3j1b"},
      {15, ">3j>1b"},
      {16, "2j"},
      {17, "3j1b"},
      {18, "3j>1b"},
      {19, ">3j1b"},
      {20, ">3j>1b"},
      {21, "2j"},
      {22, "3j"},
      {23, ">3j"},

    };
      
      
std::vector<BinLabelOptions> theSRLabelOptionsFor3L = {

      {1, "nB0,nJ2"},
      {2, "nB0,nJ3"},
      {3, "nB0,nJ#geq4"},
      {4, "nB1,nJ2"},
      {5, "nB1,nJ3"},
      {6, "nB1,nJ#geq4"},
      {7, "nB#geq2,nJ2"},
      {8, "nB#geq2,nJ3"},
      {9, "nB#geq2,nJ#geq4"},
      
    };

std::vector<BinLabelOptions> flavourLabelOptionsFor2L = {
      
      
      {1, "#mu^{-}#mu^{-}"},
      {2, "#mu^{-}e^{-}"},
      {3, "e^{-}e^{-}"},

      {4, "#mu^{+}#mu^{+}"},
      {5, "#mu^{+}e^{+}"},
      {6, "e^{+}e^{+}"},
      };

std::vector<BinLabelOptions> flavourLabelOptionsFor3L = {
      
      {1, "#mu#mu#mu"},
      {2, "#mu#mu e"},
      {3, "#mu ee"},
      {4, "eee"},
      
    };

const int nVars  = 11;

TString varN[nVars] = {
    "Raw E_{T}^{miss}", "Type I E_{T}^{miss}", 
    "u_{para}", "u_{perp}", 
    "p_{T}^{leading} [GeV]", "p_{T}^{trailing} [GeV]", 
    "#eta^{leading} [GeV]", "#eta^{trailing} [GeV]", 
    "NPV", "mll",
    "E_{T}^{miss} smeared"
}; 

double varMin[nVars] = {
    0, 0, 
    -200, -200, 
    0, 0, 
    -2.5, -2.5, 
    0, 81,
    0
};
    
double varMax[nVars] = {
    300, 300,
    200, 200, 
    100, 100, 
    2.5, 2.5, 
    90, 101,
    300
};
    
int nBins[nVars] = {
    60, 60,
    80, 80, 
    100, 100, 
    50, 50, 
    90, 40,
    60
};

double weightMC[3] = {4.8, 23.1, 13.5};
// Lepton SF
/*
TFile *file_dataMC_runB = TFile::Open("data/pileUpReweighing/puWeights_2017data_2017MC_4p8fb_nTrueVertices_nTrueVgr10.root","READ"); // PU reweighing
TFile *file_dataMC_runCDE = TFile::Open("data/pileUpReweighing/puWeights_2017data_2017MC_23p1fb_nTrueVertices_nTrueVgr10.root","READ");
TFile *file_dataMC_runF = TFile::Open("data/pileUpReweighing/puWeights_2017data_2017MC_13p5fb_nTrueVertices_nTrueVgr10.root","READ");

TH1D *h_dataMC[3] = {(TH1D*)file_dataMC_runB->Get("puw"), 
                     (TH1D*)file_dataMC_runCDE->Get("puw"),
                     (TH1D*)file_dataMC_runF->Get("puw")
                   }; 
*/

TFile *file_dataMC_2017 = TFile::Open("data/pileUpReweighing/puWeights_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_Fall17.root","READ"); // PU reweighing
TH1D *h_dataMC_2017[5] = {(TH1D*)file_dataMC_2017->Get("puw_Run2017B_central"),
                          (TH1D*)file_dataMC_2017->Get("puw_Run2017C_central"),                         
                          (TH1D*)file_dataMC_2017->Get("puw_Run2017D_central"),                         
                          (TH1D*)file_dataMC_2017->Get("puw_Run2017E_central"),                         
                          (TH1D*)file_dataMC_2017->Get("puw_Run2017F_central")
};

TFile *lepSF_ele_file_B = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D_runB_passingTight94X.root","READ");
TFile *lepSF_ele_file_C = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D_runC_passingTight94X.root","READ");
TFile *lepSF_ele_file_D = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D_runD_passingTight94X.root","READ");
TFile *lepSF_ele_file_E = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D_runE_passingTight94X.root","READ");
TFile *lepSF_ele_file_F = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D_runF_passingTight94X.root","READ");

TFile * recoToLoose_leptonSF_mu1 = TFile::Open("data/leptonSF/TnP_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root","read");
TFile * recoToLoose_leptonSF_mu2 = TFile::Open("data/leptonSF/TnP_NUM_MiniIsoLoose_DENOM_LooseID_VAR_map_pt_eta.root","read");
TFile * recoToLoose_leptonSF_mu3 = TFile::Open("data/leptonSF/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root","read");

TFile *lepSF_mu_trackBF = TFile::Open("data/leptonSF/Tracking_EfficienciesAndSF_BCDEF.root", "READ");
TFile *lepSF_mu_trackGH = TFile::Open("data/leptonSF/Tracking_EfficienciesAndSF_GH.root", "READ");

TFile *electronTrack_B = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D_runB_passingRECO.root","READ");
TFile *electronTrack_C = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D_runC_passingRECO.root","READ");
TFile *electronTrack_D = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D_runD_passingRECO.root","READ");
TFile *electronTrack_E = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D_runE_passingRECO.root","READ");
TFile *electronTrack_F = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D_runF_passingRECO.root","READ");
TFile *electronTrack_lowEt = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt.root","READ");

TFile *lepSF_el_LeptonMVAfile = TFile::Open("data/leptonSF/scaleFactors.root","READ");
//TFile *lepSF_el_LeptonMVAfile = leptonSelectionAnalysis == 2 ? TFile::Open("leptonSF/lepMVAEffSF_el_2lss.root","READ") : TFile::Open("leptonSF/lepMVAEffSF_el_3l.root","READ");
TFile *lepSF_mu_LeptonMVAfile = leptonSelectionAnalysis == 2 ? TFile::Open("data/leptonSF/lepMVAEffSF_mu_2lss.root","READ") : TFile::Open("data/leptonSF/lepMVAEffSF_mu_3l.root","READ");


TGraphAsymmErrors* lepSFMaps1DMuon[2] = {
    (TGraphAsymmErrors*) lepSF_mu_trackBF->Get("ratio_eff_eta3_dr030e030_corr"),
    (TGraphAsymmErrors*) lepSF_mu_trackGH->Get("ratio_eff_eta3_dr030e030_corr"),
  };


TH2F* lepSFMapsElectron[11] = {
    (TH2F*) (electronTrack_B->Get("EGamma_SF2D")),
    (TH2F*) (electronTrack_C->Get("EGamma_SF2D")),
    (TH2F*) (electronTrack_D->Get("EGamma_SF2D")),
    (TH2F*) (electronTrack_E->Get("EGamma_SF2D")),
    (TH2F*) (electronTrack_F->Get("EGamma_SF2D")),
    (TH2F*) (electronTrack_lowEt->Get("EGamma_SF2D")),
    (TH2F*) (lepSF_ele_file_B->Get("EGamma_SF2D")),
    (TH2F*) (lepSF_ele_file_C->Get("EGamma_SF2D")),
    (TH2F*) (lepSF_ele_file_D->Get("EGamma_SF2D")),
    (TH2F*) (lepSF_ele_file_E->Get("EGamma_SF2D")),
    (TH2F*) (lepSF_ele_file_F->Get("EGamma_SF2D")),
};

TH2F* lepSFMapsMuon[4] = {

    (TH2F*)(recoToLoose_leptonSF_mu1->Get("SF")),
    (TH2F*)(recoToLoose_leptonSF_mu2->Get("SF")),
    (TH2F*)(recoToLoose_leptonSF_mu3->Get("SF")),
    (TH2F*) lepSF_mu_LeptonMVAfile->Get("sf"),

};

// btag SF
BTagCalibration calib_csvv2[2] {
  {"CSVv2BF", "data/btagSF/CSVv2_Moriond17_B_F.csv"},
  {"CSVv2GH", "data/btagSF/CSVv2_Moriond17_G_H.csv"},
};

const std::vector<std::string> otherSysTypes={"up_hf", "down_hf", "up_lf", "down_lf"};

BTagCalibrationReader readerBtag[2][3]{ { {BTagEntry::OP_RESHAPING, "central", otherSysTypes},
                                          {BTagEntry::OP_RESHAPING, "central", otherSysTypes},
                                          {BTagEntry::OP_RESHAPING, "central", otherSysTypes}},

                                        { {BTagEntry::OP_RESHAPING, "central", otherSysTypes},
                                          {BTagEntry::OP_RESHAPING, "central", otherSysTypes},
                                          {BTagEntry::OP_RESHAPING, "central", otherSysTypes}}

};

TFile *file_btagEff = TFile::Open("data/btagSF/btageff__ttbar_powheg_pythia8_25ns_Moriond17.root","READ"); // btagEff
TH2D* h_btagEff[3] = {
    (TH2D*)file_btagEff->Get("h2_BTaggingEff_csv_med_Eff_udsg"), 
    (TH2D*)file_btagEff->Get("h2_BTaggingEff_csv_med_Eff_c"),
    (TH2D*)file_btagEff->Get("h2_BTaggingEff_csv_med_Eff_b")
}; 

// trees for BDT
TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );   
TFile* fileDummy = new TFile("fileDummy.root", "RECREATE");
TTree* signalTree = new TTree("signalTree","signalTree");
TTree* bkgTree = new TTree("bkgTree","bkgTree");

double _weightEventInTree;
    
double minDeltaR;
double mtHighest;
double mtLowest;

double leadpt;
double trailpt;
double leadingJetPt;
double trailJetPt;

int nJLoc;
int nBLoc;
double HTLoc;
double MET;

Float_t userHTLoc, user_met, userele_mll, usermt, usermtlow, userleadpt, usertrailpt, userleadingjetpt, usertrailjetpt, userminDeltaR, usernJLoc, usernBLoc;


// For FR

const int nPt = 6;
const double ptBins[nPt] = {15., 20., 30., 45., 65., 100.};

const int nQt = 26;
const double qtBins[nQt] = {18., 24., 30., 38., 46., 52., 60, 68, 76, 84, 92, 100, 115, 130, 150, 175, 200, 225, 250, 275, 305, 335, 365, 400, 440, 500};

TH1D * histPhiZCorr[25][2][4];
TH1D * histPhiZUnCorr[25][2][4];

TH1D * histMetCorr[25][2][4];
TH1D * histMetUnCorr[25][2][4];

TH1D * sigmaParUnCorr[25][2][4];
TH1D * sigmaPerpUnCorr[25][2][4];

TH1D * sigmaParCorr[25][2][4];
TH1D * sigmaPerpCorr[25][2][4];

#endif
