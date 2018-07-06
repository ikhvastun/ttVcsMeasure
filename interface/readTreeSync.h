#ifndef readTreeSync_H
#define readTreeSync_H

#include <iostream>
#include <fstream>

const int nFlavors = 2;

const int leptonSelectionAnalysis = 3;

const int nSamples = 100;
const int dataSample = 0;

TString flavorsString[2] = {"el", "mu"};
TString additionalString[2] = {"_NC", ""};

double leptonMVAcutAnalysis = leptonSelectionAnalysis == 2 ? 0.6 : (leptonSelectionAnalysis == 3 ? 0.4 : 0.8); // should be 0.6 for 2L analysis
double magicFactorAnalysis = leptonSelectionAnalysis == 2 ? 0.9 : 0.85;

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
      /*
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
      */
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



/*
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
    "number of PV interactions",
    "p_{T}^{4th} [GeV]",
    "M_{lll} [GeV]",
    "eta^{leading}", "eta^{sub-leading}"
    "p_{T}^{leading}[GeV]", "p_{T}^{sub-leading}[GeV]",
    "eta^{leading}", "eta^{sub-leading}"
    "p_{T}^{leading}[GeV]", "p_{T}^{sub-leading}[GeV]",
    "eta^{leading}", "eta^{sub-leading}",
    "mt", "mt", "mt", "mt",
    "cosThetaStar",
    "R_{lead lep}^{jet}"
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
    -0.5,
    0,
    81,
    -2.5, -2.5,
    0, 0, -2.5, -2.5,
    0, 0, -2.5, -2.5,
    0, 0, 0, 0,
    -1,
    0.4
  };
    
double varMax[nVars] = {
    300,200, 100,
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
    59.5,
    200,
    101,
    2.5, 2.5,
    60, 60, 2.5, 2.5,
    60, 60, 2.5, 2.5,
    200, 200, 200, 200,
    1,
    3
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
    20,
    20,
    20,
    50, 50,
    6, 6, 50, 50,
    6, 6, 50, 50,
    20, 20, 20, 20,
    5,
    26 
};
*/


// Lepton SF

//TFile *file_dataMC_2016 = TFile::Open("data/pileUpReweighing/puWeights_ZZTo4L_13TeV_powheg_pythia8_Summer16.root","READ"); // PU reweighing
TFile *file_dataMC_2016 = TFile::Open("data/pileUpReweighing/puWeights_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8_Summer16.root","READ"); // PU reweighing
//TFile *file_dataMC_2016 = TFile::Open("data/pileUpReweighing/puWeights_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root","READ"); // PU reweighing
TH1D *h_dataMC_2016 = (TH1D*)file_dataMC_2016->Get("puw_Run2016Inclusive_central"); 

//TFile *file_dataMC_2016 = TFile::Open("data/pileUpReweighing/puw_nTrueInt_Moriond2017_36p5fb_Summer16_central.root","READ"); // PU reweighing
//TH1D *h_dataMC_2016 = (TH1D*)file_dataMC_2016->Get("puw"); 

TFile *file_dataMC_2017 = TFile::Open("data/pileUpReweighing/puWeights_ZZTo4L_13TeV_powheg_pythia8_Fall17.root","READ"); // PU reweighing
//TFile *file_dataMC_2017 = TFile::Open("data/pileUpReweighing/puWeights_WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_Fall17.root","READ"); // PU reweighing
//TFile *file_dataMC_2017 = TFile::Open("data/pileUpReweighing/puWeights_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_Fall17.root","READ"); // PU reweighing
TH1D *h_dataMC_2017 = (TH1D*)file_dataMC_2017->Get("puw_Run2017Inclusive_central"); 

//TFile *file_dataMC = TFile::Open("data/pileUpReweighing/puWeights_36p8slashfb.root","READ"); // PU reweighing
//TH1D *h_dataMC = (TH1D*)file_dataMC->Get("puw"); 

//TFile *lepSF_ele_file_2016 = TFile::Open("data/leptonSF/scaleFactors_electrons_2016_old.root","READ");
TFile *lepSF_ele_file_2016 = TFile::Open("data/leptonSF/scaleFactors_electrons_2016.root","READ");
TFile *lepSF_ele_file_2017 = TFile::Open("data/leptonSF/scaleFactors_electrons_2017.root","READ");

TFile *lepSF_mu_file_2016 = TFile::Open("data/leptonSF/scaleFactors_muons_2016.root","READ");
TFile *lepSF_mu_file_2017 = TFile::Open("data/leptonSF/scaleFactors_muons_2017.root","READ");

//TFile * recoToLoose_leptonSF_mu1 = TFile::Open("data/leptonSF/TnP_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root","read");
//TFile * recoToLoose_leptonSF_mu2 = TFile::Open("data/leptonSF/TnP_NUM_MiniIsoLoose_DENOM_LooseID_VAR_map_pt_eta.root","read");
//TFile * recoToLoose_leptonSF_mu3 = TFile::Open("data/leptonSF/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root","read");

TFile *lepSF_mu_trackBF = TFile::Open("data/leptonSF/Tracking_EfficienciesAndSF_BCDEF.root", "READ");
TFile *lepSF_mu_trackGH = TFile::Open("data/leptonSF/Tracking_EfficienciesAndSF_GH.root", "READ");
TFile *lepSF_mu_trackBH = TFile::Open("data/leptonSF/Tracking_EfficienciesAndSF_BCDEFGH.root", "READ");

TFile *electronTrack_2016 = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D.root","READ");
TFile *electronTrack_2017 = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_2017.root","READ");
TFile *electronTrack_lowEt_2017 = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt_2017.root","READ");

TFile *lepSF_el_LeptonMVAfile = leptonSelectionAnalysis == 2 ? TFile::Open("data/leptonSF/lepMVAEffSF_el_2lss.root","READ") : TFile::Open("data/leptonSF/lepMVAEffSF_el_3l.root","READ");
TFile *lepSF_mu_LeptonMVAfile = leptonSelectionAnalysis == 2 ? TFile::Open("data/leptonSF/lepMVAEffSF_mu_2lss.root","READ") : TFile::Open("data/leptonSF/lepMVAEffSF_mu_3l.root","READ");
//TFile *lepSF_el_LeptonMVAfile = TFile::Open("leptonSF/scaleFactors.root","READ");

TGraphAsymmErrors* lepSFMaps1DMuon[2] = {
    (TGraphAsymmErrors*) lepSF_mu_trackBH->Get("ratio_eff_eta3_dr030e030_corr"),
    // this is not used at the moment 
    (TGraphAsymmErrors*) lepSF_mu_trackGH->Get("ratio_eff_eta3_dr030e030_corr"),
  };

TH2F* lepSFMapsElectron[13] = {
    (TH2F*) electronTrack_2016->Get("EGamma_SF2D"),
    (TH2F*) (lepSF_ele_file_2016->Get("EleToTTVLoose")),
    (TH2F*) (lepSF_ele_file_2016->Get("TTVLooseToTTVLeptonMvattW")),
    (TH2F*) (lepSF_ele_file_2016->Get("TTVLooseToTTVLeptonMvattZ3l")),
    //(TH2F*) (lepSF_ele_file_2016->Get("TTVLooseToTTVLeptonMvattZ4l")),
    (TH2F*) (lepSF_ele_file_2016->Get("TTVLooseToTTVLeptonMvatZq")),
    (TH2F*) (lepSF_ele_file_2016->Get("TTVLeptonMvattWToTightCharge")),
    (TH2F*) electronTrack_2017->Get("EGamma_SF2D"),
    (TH2F*) electronTrack_lowEt_2017->Get("EGamma_SF2D"),
    (TH2F*) (lepSF_ele_file_2017->Get("EleToTTVLoose")),
    (TH2F*) (lepSF_ele_file_2017->Get("TTVLooseToTTVLeptonMvattW")),
    (TH2F*) (lepSF_ele_file_2017->Get("TTVLooseToTTVLeptonMvattZ3l")),
    (TH2F*) (lepSF_ele_file_2017->Get("TTVLooseToTTVLeptonMvattZ4l")),
    (TH2F*) (lepSF_ele_file_2017->Get("TTVLeptonMvattWToTightCharge")),
};

TH2F* lepSFMapsMuon[10] = {

    (TH2F*) (lepSF_mu_file_2016->Get("MuonToTTVLoose")),
    (TH2F*) (lepSF_mu_file_2016->Get("TTVLooseToTTVLeptonMvattW")),
    (TH2F*) (lepSF_mu_file_2016->Get("TTVLooseToTTVLeptonMvattZ3l")),
    //(TH2F*) (lepSF_mu_file_2016->Get("TTVLooseToTTVLeptonMvattZ4l")),
    (TH2F*) (lepSF_mu_file_2016->Get("TTVLooseToTTVLeptonMvatZq")),
    (TH2F*) (lepSF_mu_file_2016->Get("TTVLeptonMvattWTotkSigmaPtOverPtCut")),
    (TH2F*) (lepSF_mu_file_2017->Get("MuonToTTVLoose")),
    (TH2F*) (lepSF_mu_file_2017->Get("TTVLooseToTTVLeptonMvattW")),
    (TH2F*) (lepSF_mu_file_2017->Get("TTVLooseToTTVLeptonMvattZ3l")),
    (TH2F*) (lepSF_mu_file_2017->Get("TTVLooseToTTVLeptonMvattZ4l")),
    (TH2F*) (lepSF_mu_file_2017->Get("TTVLeptonMvattWTotkSigmaPtOverPtCut")),

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
  //{"CSVv2BF", "data/btagSF/CSVv2_Moriond17_B_F.csv"},
  //{"CSVv2GH", "data/btagSF/CSVv2_Moriond17_G_H.csv"},
  {"deepCSV_2016", "data/btagSF/DeepCSV_Moriond17_B_H.csv"},
  {"deepCSV_2017", "data/btagSF/DeepCSV_94XSF_V2_B_F.csv"},
};

// available for csv
//const std::vector<std::string> otherSysTypes={"up_hf", "down_hf", "up_lf", "down_lf"};
const std::vector<std::string> otherSysTypes={"up", "down"};

/*
BTagCalibrationReader readerBtag[2][3]{ { {BTagEntry::OP_RESHAPING, "central", otherSysTypes},
                                          {BTagEntry::OP_RESHAPING, "central", otherSysTypes},
                                          {BTagEntry::OP_RESHAPING, "central", otherSysTypes}},

                                        { {BTagEntry::OP_RESHAPING, "central", otherSysTypes},
                                          {BTagEntry::OP_RESHAPING, "central", otherSysTypes},
                                          {BTagEntry::OP_RESHAPING, "central", otherSysTypes}}

};
*/
BTagCalibrationReader readerBtag[2][3]{ { {BTagEntry::OP_MEDIUM, "central", otherSysTypes},
                                          {BTagEntry::OP_MEDIUM, "central", otherSysTypes},
                                          {BTagEntry::OP_MEDIUM, "central", otherSysTypes}},

                                        { {BTagEntry::OP_MEDIUM, "central", otherSysTypes},
                                          {BTagEntry::OP_MEDIUM, "central", otherSysTypes},
                                          {BTagEntry::OP_MEDIUM, "central", otherSysTypes}}
};

TFile *file_btagEff_ttZ4l_2016 = TFile::Open("data/btagSF/bTagEff_deepCSV_cleaned_ttZ4l_2016.root","READ"); // btagEff
TFile *file_btagEff_ttZ4l_2017 = TFile::Open("data/btagSF/bTagEff_deepCSV_cleaned_ttZ4l_2017.root","READ"); // btagEff
TFile *file_btagEff_ttZ3l_2016 = TFile::Open("data/btagSF/bTagEff_deepCSV_cleaned_ttZ3l_2016.root","READ"); // btagEff
TFile *file_btagEff_ttZ3l_2017 = TFile::Open("data/btagSF/bTagEff_deepCSV_cleaned_ttZ3l_2017.root","READ"); // btagEff
TFile *file_btagEff_ttW_2016 = TFile::Open("data/btagSF/bTagEff_deepCSV_cleaned_ttW_2016.root","READ"); // btagEff
TFile *file_btagEff_ttW_2017 = TFile::Open("data/btagSF/bTagEff_deepCSV_cleaned_ttW_2017.root","READ"); // btagEff
TH2D* h_btagEff[2][9] = {{ (TH2D*)file_btagEff_ttW_2016->Get("bTagEff_mediumudsg"), 
                           (TH2D*)file_btagEff_ttW_2016->Get("bTagEff_mediumcharm"),
                           (TH2D*)file_btagEff_ttW_2016->Get("bTagEff_mediumbeauty"), 
                           (TH2D*)file_btagEff_ttZ3l_2016->Get("bTagEff_mediumudsg"), 
                           (TH2D*)file_btagEff_ttZ3l_2016->Get("bTagEff_mediumcharm"),
                           (TH2D*)file_btagEff_ttZ3l_2016->Get("bTagEff_mediumbeauty"), 
                           (TH2D*)file_btagEff_ttZ4l_2016->Get("bTagEff_mediumudsg"), 
                           (TH2D*)file_btagEff_ttZ4l_2016->Get("bTagEff_mediumcharm"),
                           (TH2D*)file_btagEff_ttZ4l_2016->Get("bTagEff_mediumbeauty"),},
                         { (TH2D*)file_btagEff_ttW_2017->Get("bTagEff_mediumudsg"), 
                           (TH2D*)file_btagEff_ttW_2017->Get("bTagEff_mediumcharm"),
                           (TH2D*)file_btagEff_ttW_2017->Get("bTagEff_mediumbeauty"), 
                           (TH2D*)file_btagEff_ttZ3l_2017->Get("bTagEff_mediumudsg"), 
                           (TH2D*)file_btagEff_ttZ3l_2017->Get("bTagEff_mediumcharm"),
                           (TH2D*)file_btagEff_ttZ3l_2017->Get("bTagEff_mediumbeauty"), 
                           (TH2D*)file_btagEff_ttZ4l_2017->Get("bTagEff_mediumudsg"), 
                           (TH2D*)file_btagEff_ttZ4l_2017->Get("bTagEff_mediumcharm"),
                           (TH2D*)file_btagEff_ttZ4l_2017->Get("bTagEff_mediumbeauty")}
}; 

int jetFlavourNumber[6] = {0,0,0,0,1,2};
/*
  KEY: TH2D     bTagEff_mediumudsg;1    bTagEff_medium_udsg_numerator
  KEY: TH2D     bTagEff_mediumcharm;1   bTagEff_medium_charm_numerator
  KEY: TH2D     bTagEff_mediumbeauty;1  bTagEff_medium_beauty_numerator
*/

// trees for BDT
TMVA::Reader *readerTTWcsttbar = new TMVA::Reader( "!Color:!Silent" );   
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

//vector<int> orderForSaveFiles = {0, 1, 2, 8, 9, 15, 14, 12, 16, 17, 18, 4, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48};

struct histInfo{
  std::string fancyName;
  int index;
  //std::string usualName;
  double varMin;
  double varMax;
  int nBins;
  //bool isEnVar;
};

const int nVars  = 52;

std::map<TString, histInfo> figNames  =         {{"ptlead",  {"Leading lepton p_{T} [GeV]", 0, 0, 300, 30}},
                                                 {"sublead", {"Sub-leading lepton p_{T} [GeV]", 1, 0, 200, 20}},
                                                 {"trail",   {"Trailing lepton p_{T} [GeV]", 2, 0, 200, 20}},
                                                 {"pt4th",   {"4th lepton p_{T} [GeV]", 3, 0, 100, 20}},
                                                 {"mtW",     {"m_{T}^{W} [GeV]", 4, 0, 200, 20}},
                                                 {"njets",   {"N_{j}", 5, -0.5, 7.5, 8}},
                                                 {"nbjets",  {"N_{b}", 6, -0.5, 4.5, 5}},
                                                 {"BDT",     {"BDT", 7, -1, 1, 10}},
                                                 {"flavour", {"flavour", 8, 0.5, (leptonSelectionAnalysis == 2 ? (static_cast<double>(flavourLabelOptionsFor2L.size()) + 0.5) : (leptonSelectionAnalysis == 3 ? (static_cast<double>(flavourLabelOptionsFor3L.size()) + 0.5) : (static_cast<double>(flavourLabelOptionsFor4L.size()) + 0.5))), (leptonSelectionAnalysis == 2 ? (static_cast<int>(flavourLabelOptionsFor2L.size())) : (leptonSelectionAnalysis == 3 ? (static_cast<int>(flavourLabelOptionsFor3L.size())) : (static_cast<int>(flavourLabelOptionsFor4L.size()))))}},
                                                 {"SR",      {"", 9, -0.5, leptonSelectionAnalysis == 2 ? (static_cast<double>(theSRLabelOptionsFor2L.size()) - 0.5) : (static_cast<double>(theSRLabelOptionsFor3L.size()) - 0.5), leptonSelectionAnalysis == 2 ? (static_cast<int>(theSRLabelOptionsFor2L.size())) : static_cast<int>(theSRLabelOptionsFor3L.size())}},
                                                 {"mll",     {"M(ll) [GeV]", 10, 81., 101., 10}},
                                                 {"ptZ",     {"p_{T}^{Z} [GeV]", 11, 0, 400, 16}},
                                                 {"ptNonZ",  {"Non-Z lepton p_{T} [GeV]", 12, 0, 200, 20}},
                                                 {"mll3e",   {"M(ll) in 3e [GeV]", 13, 0, 200, 20}},
                                                 {"mll2e1mu",{"M(ll) in 2e1mu [GeV]", 14, 0, 200, 20}},
                                                 {"mll1e2mu",{"M(ll) in 1e2mu [GeV]", 15, 0, 200, 20}},
                                                 {"mll3mu",  {"M(ll) in 3mu [GeV]", 16, 0, 200, 20}},
                                                 {"met",     {"E_{T}^{miss} [GeV]", 17, 0, 300, 15}},
                                                 {"deltaR",  {"#Delta R(jet, trailing lepton)", 18, 0.4, 3., 13}},
                                                 {"deltaRlead",  {"#Delta R(jet, leading lepton)", 19, 0.4, 3., 13}},
                                                 {"mtLeading",{"Leading lepton M_{T} [GeV]", 20, 0, 300, 20}},
                                                 {"mtTrailing",{"Trailing lepton M_{T} [GeV]", 21, 0, 200, 20}},
                                                 {"leadJetPt", {"Leading non-b jet p_{T} [GeV]", 22, 30, 310, 14}},
                                                 {"trailJetPt", {"Trailing non-b jet p_{T} [GeV]", 23, 30, 310, 14}},
                                                 {"SRnpCR", {"", 24, -0.5, 2.5, 3}},
                                                 {"nPV", {"number of PV", 25, -0.5, 49.5, 25}},
                                                 {"mlll", {"M(lll) [GeV]", 26, 81., 101., 10}},
                                                 {"etaLead", {"Leading lepton #eta", 27, -2.5, 2.5, 20}},
                                                 {"etaSubl", {"Sub-leading lepton #eta", 28, -2.5, 2.5, 20}},
                                                 {"etaTrail", {"Trailing lepton #eta", 29, -2.5, 2.5, 20}},
                                                 {"eta4th", {"4th lepton #eta", 30, -2.5, 2.5, 20}},
                                                 {"mt_3m", {"m_{T}^{W} in 3#mu [GeV]", 31, 0, 200, 20}},
                                                 {"mt_2m1e", {"m_{T}^{W} in 2#mu e [GeV]", 32, 0, 200, 20}}, 
                                                 {"mt_1m2e", {"m_{T}^{W} in 1#mu 2e [GeV]", 33, 0, 200, 20}}, 
                                                 {"mt_3e", {"m_{T}^{W} in 3e [GeV]", 34, 0, 200, 20}},
                                                 {"cosThetaStar", {"cos(#Theta^{*})", 35, -1, 1, 5}},
                                                 {"mll_ss",  {"Invariant mass of ss 2l pair [GeV]", 36, 0, 300, 20}},
                                                 {"chargeOfLeptons",  {"Charge of the leptons in ss2l channel", 37, -1.5, 1.5, 3}},
                                                 {"ll_deltaR",  {"#Delta R(leading lepton, trailing lepton)", 38, 0, 7., 35}},
                                                 {"mt2ll_ss",  {"M_{T2}^{ll} [GeV]", 39, 0, 200., 20}}
                                           };
#endif 
