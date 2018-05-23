#ifndef readTreeSync_H
#define readTreeSync_H

#include <iostream>
#include <fstream>

const int nFlavors = 2;

const int leptonSelectionAnalysis = 4;

const int nSamples = 100;
const int dataSample = 0;

TString flavorsString[2] = {"el", "mu"};
TString additionalString[2] = {"_NC", ""};

double leptonMVAcutAnalysis = leptonSelectionAnalysis == 2 ? 0.6 : (leptonSelectionAnalysis == 3 ? 0.4 : -0.4);
double magicFactorAnalysis = leptonSelectionAnalysis == 2 ? 0.8 : 0.85;

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
      //{21, "2j"},
      //{22, "3j"},
      //{23, ">3j"},

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

std::vector<BinLabelOptions> flavourLabelOptionsFor4L = {
      
      {1, "#mu#mu#mu#mu"},
      {2, "#mu#mu ee"},
      {3, "eeee"},
      
    };


const int nVars  = 33;

 TString varN[nVars] = {
    "p_{T}^{leading}[GeV]", "p_{T}^{sub-leading} [GeV]", "p_{T}^{trailing} [GeV]",
    "M_{T}^{leading} [GeV]", "M_{T}^{trailing} [GeV]", 
    "R_{trail lep}^{jet}",
    "E_{T}^{miss}", "H_{T}", "N_{jets}", "N_{b jets}", "Jet p_{T}^{leading} [GeV]","Jet p_{T}^{trailing} [GeV]",
    "BDTG", "PU",
    "flavour",
    "SR",
    "m_{ll} [GeV]", "p_{T}^{Z} [GeV]", "Non-Z lepton p_{T} [GeV]",
    "m_{ll} [GeV]", "m_{ll} [GeV]", "m_{ll} [GeV]", "m_{ll} [GeV]",
    "met", "minDeltaR", "mt highest", "mt lowest", "leading jet pt", "trailing jet pt",
    "SR control region",
    "ele conv p_{T}", "ele conv eta",
    "number of PV interactions"
    }; 


double varMin[nVars] = {
    0, 0, 0,
    0, 0,
    0, 
    0, 0, -0.5, -0.5, 30, 30, 
    -1, 0, 
    0.5,
    -0.5,
    71., 0, 0,
    81., 81., 81., 81., 
    0, 0.4, 30, 30, 0, 0,
    19.5, 
    0, -2.5,
    -0.5
  };
    
double varMax[nVars] = {
    300,200, 200,
    300,200,
    3,
    400, 1000, 7.5, 5.5, 420, 300, 
    1, 50, 
    leptonSelectionAnalysis == 2 ? (static_cast<double>(flavourLabelOptionsFor2L.size()) + 0.5) : (leptonSelectionAnalysis == 3 ? (static_cast<double>(flavourLabelOptionsFor3L.size()) + 0.5) : (static_cast<double>(flavourLabelOptionsFor4L.size()) + 0.5)),
    leptonSelectionAnalysis == 2 ? (static_cast<double>(theSRLabelOptionsFor2L.size()) - 0.5) : (static_cast<double>(theSRLabelOptionsFor3L.size()) - 0.5),
    111., 400, 180,
    101, 101, 101, 101,
    300, 4, 300, 300, 200, 200, 
    22.5,
    200, 2.5,
    59.5
};
    
int nBins[nVars] = {
    20, 20, 20,
    30, 20,
    30,
    10, 10, 8, 6, 13, 18, 
    10, 25,
    leptonSelectionAnalysis == 2 ? (static_cast<int>(flavourLabelOptionsFor2L.size())) : (leptonSelectionAnalysis == 3 ? (static_cast<int>(flavourLabelOptionsFor3L.size())) : (static_cast<int>(flavourLabelOptionsFor4L.size()))) ,
    leptonSelectionAnalysis == 2 ? (static_cast<int>(theSRLabelOptionsFor2L.size())) : static_cast<int>(theSRLabelOptionsFor3L.size()),
    20, 16, 12,
    10, 10, 10, 10,
    15, 18, 18, 18, 20, 20,
    3,
    20, 20,
    20
};


// Lepton SF

TFile *file_dataMC = TFile::Open("data/pileUpReweighing/puWeights_36p8slashfb.root","READ"); // PU reweighing
TH1D *h_dataMC = (TH1D*)file_dataMC->Get("puw"); 

TFile *lepSF_ele_file_2016 = TFile::Open("data/leptonSF/scaleFactors2016.root","READ");
TFile *lepSF_ele_file_2017 = TFile::Open("data/leptonSF/scaleFactors2017.root","READ");

TFile * recoToLoose_leptonSF_mu1 = TFile::Open("data/leptonSF/TnP_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root","read");
TFile * recoToLoose_leptonSF_mu2 = TFile::Open("data/leptonSF/TnP_NUM_MiniIsoLoose_DENOM_LooseID_VAR_map_pt_eta.root","read");
TFile * recoToLoose_leptonSF_mu3 = TFile::Open("data/leptonSF/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root","read");

TFile *lepSF_mu_trackBF = TFile::Open("data/leptonSF/Tracking_EfficienciesAndSF_BCDEF.root", "READ");
TFile *lepSF_mu_trackGH = TFile::Open("data/leptonSF/Tracking_EfficienciesAndSF_GH.root", "READ");

TFile *electronTrack = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D.root","READ");

TFile *lepSF_el_LeptonMVAfile = leptonSelectionAnalysis == 2 ? TFile::Open("data/leptonSF/lepMVAEffSF_el_2lss.root","READ") : TFile::Open("data/leptonSF/lepMVAEffSF_el_3l.root","READ");
TFile *lepSF_mu_LeptonMVAfile = leptonSelectionAnalysis == 2 ? TFile::Open("data/leptonSF/lepMVAEffSF_mu_2lss.root","READ") : TFile::Open("data/leptonSF/lepMVAEffSF_mu_3l.root","READ");
//TFile *lepSF_el_LeptonMVAfile = TFile::Open("leptonSF/scaleFactors.root","READ");


TGraphAsymmErrors* lepSFMaps1DMuon[2] = {
    (TGraphAsymmErrors*) lepSF_mu_trackBF->Get("ratio_eff_eta3_dr030e030_corr"),
    (TGraphAsymmErrors*) lepSF_mu_trackGH->Get("ratio_eff_eta3_dr030e030_corr"),
  };

TH2F* lepSFMapsElectron[11] = {
    (TH2F*) electronTrack->Get("EGamma_SF2D"),
    (TH2F*) (lepSF_ele_file_2016->Get("EleToTTVLoose")),
    (TH2F*) (lepSF_ele_file_2016->Get("TTVLooseToTTVLeptonMvattW")),
    (TH2F*) (lepSF_ele_file_2016->Get("TTVLooseToTTVLeptonMvattZ3l")),
    (TH2F*) (lepSF_ele_file_2016->Get("TTVLooseToTTVLeptonMvattZ4l")),
    (TH2F*) (lepSF_ele_file_2016->Get("TTVLeptonMvattWToTightCharge")),
    (TH2F*) (lepSF_ele_file_2017->Get("EleToTTVLoose")),
    (TH2F*) (lepSF_ele_file_2017->Get("TTVLooseToTTVLeptonMvattW")),
    (TH2F*) (lepSF_ele_file_2017->Get("TTVLooseToTTVLeptonMvattZ3l")),
    (TH2F*) (lepSF_ele_file_2017->Get("TTVLooseToTTVLeptonMvattZ4l")),
    (TH2F*) (lepSF_ele_file_2017->Get("TTVLeptonMvattWToTightCharge")),
    //leptonSelectionAnalysis == 2 ? (TH2F*) lepSF_el_LeptonMVAfile->Get("GsfElectronToTTZ2017TightCharge") : (TH2F*) lepSF_el_LeptonMVAfile->Get("GsfElectronToTTZ2017")
};

TH2F* lepSFMapsMuon[4] = {

    (TH2F*)(recoToLoose_leptonSF_mu1->Get("SF")),
    (TH2F*)(recoToLoose_leptonSF_mu2->Get("SF")),
    (TH2F*)(recoToLoose_leptonSF_mu3->Get("SF")),
    (TH2F*) lepSF_mu_LeptonMVAfile->Get("sf"),

};


// used in ttV
/*
TFile *lepSF_ele_file = TFile::Open("/Users/illiakhvastunov/Desktop/CERN/ss2l_2016_fulldataset/LeptonSF/scaleFactors.root","READ");
TFile *lepSF_mu_ID = TFile::Open("/Users/illiakhvastunov/Desktop/CERN/ss2l_2016_fulldataset/LeptonSF/EfficienciesAndSF_BCDEF.root", "READ");
TFile *lepSF_mu_iso = TFile::Open("/Users/illiakhvastunov/Desktop/CERN/ss2l_2016_fulldataset/LeptonSF/EfficienciesAndSF_BCDEF_iso.root", "READ");
TFile *lepSF_mu_ID_GH = TFile::Open("/Users/illiakhvastunov/Desktop/CERN/ss2l_2016_fulldataset/LeptonSF/EfficienciesAndSF_GH.root", "READ");
TFile *lepSF_mu_iso_GH = TFile::Open("/Users/illiakhvastunov/Desktop/CERN/ss2l_2016_fulldataset/LeptonSF/EfficienciesAndSF_GH_iso.root", "READ");

TFile *lepSF_mu_trackBF = TFile::Open("/Users/illiakhvastunov/Desktop/CERN/ss2l_2016_fulldataset/LeptonSF/Tracking_EfficienciesAndSF_BCDEF.root", "READ");
TFile *lepSF_mu_trackGH = TFile::Open("/Users/illiakhvastunov/Desktop/CERN/ss2l_2016_fulldataset/LeptonSF/Tracking_EfficienciesAndSF_GH.root", "READ");

TFile *electronTrack = TFile::Open("/Users/illiakhvastunov/Desktop/CERN/ss2l_2016_fulldataset/LeptonSF/egammaEffi.txt_EGM2D.root","READ");



TH2D* lepSFMaps[8] = {
    (TH2D*) lepSF_ele_file->Get("GsfElectronToTTZ"),
    (TH2D*) lepSF_ele_file->Get("MVATightElectronToConvVetoIHit0"),
    (TH2D*) lepSF_ele_file->Get("MVATightConvIHit0ElectronToCharge"), 
    (TH2D*) lepSF_mu_ID->Get("MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio"),
    (TH2D*) lepSF_mu_iso->Get("LooseISO_MediumID_pt_eta/pt_abseta_ratio"),
    (TH2D*) lepSF_mu_ID_GH->Get("MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio"),
    (TH2D*) lepSF_mu_iso_GH->Get("LooseISO_MediumID_pt_eta/pt_abseta_ratio"),
    //(TH2D*) electronTrack->Get("EGamma_SF2D")

};

TGraphAsymmErrors* lepSFMaps1D[3] = {
    (TGraphAsymmErrors*) lepSF_mu_trackBF->Get("ratio_eff_eta3_dr030e030_corr"),
    (TGraphAsymmErrors*) lepSF_mu_trackGH->Get("ratio_eff_eta3_dr030e030_corr"),
    (TGraphAsymmErrors*) electronTrack->Get("grSF1D_0")
  };
  */

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
    
double minDeltaRlead;
double minDeltaR;
double mtHighest;
double mtLowest;

double leadpt;
double trailpt;
double leadeta;
double traileta;
double leadingJetPt;
double trailJetPt;

int nJLoc;
int nJLocNotB;
int nBLoc;
double HTLoc;
double MET;
int chargeOfLeptons;
double mll_ss;
double ll_deltaR;
double mt2ll_ss;

Float_t userHTLoc, user_met, userele_mll, usermt, usermtlow, userleadpt, usertrailpt, userleadeta, usertraileta, userleadingjetpt, usertrailjetpt, userminDeltaRlead, userminDeltaR, usernJLoc, usernBLoc, userchargeOfLeptons, usermll_ss, userll_deltaR, usermt2ll_ss;

TString eraRuns[2] = {"", "_GH"};

// For FR

const int nPt = 6;
const double ptBins[nPt] = {15., 20., 30., 45., 65., 100.};

#endif 
