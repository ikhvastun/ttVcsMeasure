#ifndef readTreeSync_H
#define readTreeSync_H

#include <iostream>
#include <fstream>

const int nFlavors = 2;

const int leptonSelectionAnalysis = 3;

const int nSamples = 100;
const int dataSample = 0;
//int uncertaintySample = distribsOrder.size() + 1; // +1 for data

TString flavorsString[2] = {"el", "mu"};
TString flavorComposString[4] = {"all", "b", "c", "light"};
//TString flavorComposString[4] = {"B_L", "B_C_L", "B_C_T_L", "B_T_L"};
TString additionalString[2] = {"_NC", ""};
TString histString[3] = {"all","tight","FR"};
TString rangeString[6] ={"10to15", "15to20", "20to30", "30to45", "45to65", "gr65"};
TString rangeEtaString[2] ={"barrel", "endcap"};

double leptonMVAcutAnalysis = leptonSelectionAnalysis == 2 ? 0.6 : 0.4; // 0.4 for ttZ
double magicFactorAnalysis = leptonSelectionAnalysis == 2 ? 0.9 : 0.85; // 0.85 for ttZ

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

const int nVars  = 16;

TString varN[nVars] = {
    "pt corr ele", "pt corr mu", "pt corr ele", "pt corr mu","pt corr ele", "pt corr mu","pt corr ele", "pt corr mu",
    "pt corr ele", "pt corr mu", "pt corr ele", "pt corr mu","pt corr ele", "pt corr mu","pt corr ele", "pt corr mu"
}; 

double varMin[nVars] = {
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0
};
    
double varMax[nVars] = {
    100, 100, 100, 100, 100, 100, 100, 100,
    100, 100, 100, 100, 100, 100, 100, 100
};
    
int nBins[nVars] = {
    20, 20, 20, 20, 20, 20, 20, 20,
    20, 20, 20, 20, 20, 20, 20, 20
};

const int nVars2D  = 4;

double varMin2D[2] = {
    0, 0
};
    
double varMax2D[2] = {
    200, 200
};
    
int nBins2D[2] = {
    100, 100
};



// Lepton SF
TFile *file_dataMC = TFile::Open("data/pileUpReweighing/puw_nTrueInt_Moriond2017_36p5fb_Summer16_central.root","READ"); // PU reweighing
TH1D *h_dataMC = (TH1D*)file_dataMC->Get("puw"); 

TFile *lepSF_ele_file = TFile::Open("data/leptonSF/scaleFactorsAll.root","READ");

TFile * recoToLoose_leptonSF_mu1 = TFile::Open("data/leptonSF/TnP_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root","read");
TFile * recoToLoose_leptonSF_mu2 = TFile::Open("data/leptonSF/TnP_NUM_MiniIsoLoose_DENOM_LooseID_VAR_map_pt_eta.root","read");
TFile * recoToLoose_leptonSF_mu3 = TFile::Open("data/leptonSF/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root","read");

TFile *lepSF_mu_trackBF = TFile::Open("data/leptonSF/Tracking_EfficienciesAndSF_BCDEF.root", "READ");
TFile *lepSF_mu_trackGH = TFile::Open("data/leptonSF/Tracking_EfficienciesAndSF_GH.root", "READ");

TFile *electronTrack = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D.root","READ");

TFile *lepSF_el_LeptonMVAfile = TFile::Open("data/leptonSF/scaleFactors.root","READ");
//TFile *lepSF_el_LeptonMVAfile = leptonSelectionAnalysis == 2 ? TFile::Open("leptonSF/lepMVAEffSF_el_2lss.root","READ") : TFile::Open("leptonSF/lepMVAEffSF_el_3l.root","READ");
TFile *lepSF_mu_LeptonMVAfile = leptonSelectionAnalysis == 2 ? TFile::Open("data/leptonSF/lepMVAEffSF_mu_2lss.root","READ") : TFile::Open("data/leptonSF/lepMVAEffSF_mu_3l.root","READ");


TGraphAsymmErrors* lepSFMaps1DMuon[2] = {
    (TGraphAsymmErrors*) lepSF_mu_trackBF->Get("ratio_eff_eta3_dr030e030_corr"),
    (TGraphAsymmErrors*) lepSF_mu_trackGH->Get("ratio_eff_eta3_dr030e030_corr"),
  };


TH2F* lepSFMapsElectron[5] = {
    (TH2F*) electronTrack->Get("EGamma_SF2D"),
    (TH2F*) (lepSF_ele_file->Get("GsfElectronToMVAVLooseFOIDEmuTightIP2D")),
    (TH2F*) (lepSF_ele_file->Get("MVAVLooseElectronToMini4")),
    (TH2F*) (lepSF_ele_file->Get("MVAVLooseElectronToConvVetoIHit1")),
    //(TH2F*) lepSF_el_LeptonMVAfile->Get("sf"),
    leptonSelectionAnalysis == 2 ? (TH2F*) lepSF_el_LeptonMVAfile->Get("GsfElectronToTTZ2017TightCharge") : (TH2F*) lepSF_el_LeptonMVAfile->Get("GsfElectronToTTZ2017")
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
TMVA::Reader *reader = new TMVA::Reader( "!Color:!SilesignalTreent" );   
TMVA::Reader *readerLeptonMVAele = new TMVA::Reader( "!Color:!Silent" );  
TMVA::Reader *readerLeptonMVAmu = new TMVA::Reader( "!Color:!Silent" );  

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

// For FR

const int nPt = 7;
const double ptBins[nPt] = {10., 15., 20., 30., 45., 65., 100.};

const int nEta = 4;
double etaBins[2][nEta] = {{0., 0.8, 1.442, 2.5}, {0., 1.2, 2.1, 2.4}};

const int rangePeriods = 6; // in lepton corrected pt 10-15, 15-20, 20-30, 30-45, 45-65, > 65 
const int rangeEtaPeriods = 2;
TH2D* fakeMapsCalc[nFlavors][3][rangePeriods][rangeEtaPeriods][nSamples+1]; 
TH1D* mtMaps[nFlavors][rangePeriods][rangeEtaPeriods][nSamples+1];
THStack* mtStack[nFlavors][rangePeriods][rangeEtaPeriods]; 

double passedPrompt_TTV[2][2];
double passedNonPrompt_TTV[2][2];
double passedPrompt[2][2][800];
double passedNonPrompt[2][2][800];
int nPoints = 800;

TH1D * hAveragePt[80];
double borderOfBarrelEndcap[2] = {1.479, 1.2};

Float_t user_pt, user_eta, user_trackMult,  user_miniIsoCharged, user_miniIsoNeutral, user_ptrel, user_ptratio, user_jetBtagCSV, user_sip3d, user_dxy, user_dz, user_segmComp, user_eleMVA, user_relIso;

#endif
