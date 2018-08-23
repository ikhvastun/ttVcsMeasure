#ifndef readTreeSync_H
#define readTreeSync_H

#include <iostream>
#include <fstream>

const int nFlavors = 2;

const int leptonSelectionAnalysis = 3;

const int nSamples = 11;
const int dataSample = 0;

const int numberOfSyst = 11;

TString flavorsString[2] = {"el", "mu"};
TString additionalString[2] = {"_NC", ""};

struct BinLabelOptions{
  int index;
  std::string labelSR;
};

std::vector<BinLabelOptions> theSRLabelOptionsForttZCleanPTZ = {
      {1, "0-75"},
      {2, "75-150"},
      {3, "150-250"},
      {4, ">250"}
};

std::vector<BinLabelOptions> theSRLabelOptionsForttZCleanCosTheta = {
      {1, "[-1,-0.5]"},
      {2, "[-0.5,0]"},
      {3, "[0, 0.5]"},
      {4, "[0.5, 1]]"}
};

std::vector<BinLabelOptions> theSRLabelOptionsForTTZ = {

    /*
      {1, "1"}, // WZ CR, first 4 bins, nbjets = 0, njets = 1, 2, 3, > 3
      {2, "2"},
      {3, "3"},
      {4, "4"},
      {5, "5"}, // ZZ CR, nbjets = 0, njets = 1, > 1; nbjets > 0, njets 1, > 1
      {6, "6"},
      {7, "7"},
      {8, "8"},
      {9, "9"}, // nonprompt CR, same as main selection, but off Z or noOSSF
      {10, "10"},
      {11, "11"},
      {12, "12"},
      {13, "13"},
      {14, "14"},
      {15, "15"},
      {16, "16"},
      {17, "17"},
      {18, "1"},// ttZ 3L, nbjets = 1, njets = 2,3,4,>4; nbjets >1, njets = 2,3,4,>4
      {19, "2"},
      {20, "3"},
      {21, "4"},
      {22, "5"},
      {23, "6"},
      {24, "7"},
      {25, "8"},
      {26, "9"}, // ttZ 4L, nbjets = 0, 1
      {27, "10"},
      */
      {1, "1"}, // WZ CR, first 4 bins, nbjets = 0, njets = 1, 2, 3, > 3
      {2, "2"},
      {3, "3"},
      {4, "> 3"},
      {5, "2"},// ttZ 3L, nbjets = 1, njets = 2,3,4,>4; nbjets >1, njets = 2,3,4,>4
      {6, "3"},
      {7, "4"},
      {8, "> 4"},
      {9, "2"},
      {10, "3"},
      {11, "4"},
      {12, "> 4"},
      {13, "0"}, // ttZ 4L, nbjets = 0, 1
      {14, "> 0"},
};

std::vector<BinLabelOptions> theSRLabelOptionsForWZCR = {

      {1, "1"}, // WZ CR, first 4 bins, nbjets = 0, njets = 1, 2, 3, > 3
      {2, "2"},
      {3, "3"},
      {4, "> 3"},
};

std::vector<BinLabelOptions> theSRLabelOptionsForZZCR = {
      {1, "nb=0, nj=1"}, // ZZ CR, nbjets = 0, njets = 1, > 1; nbjets > 0, njets 1, > 1
      {2, "nb=0, nj>1"},
      {3, "nb>0, nj=1"},
      {4, "nb>0, nj>1"},
};

std::vector<BinLabelOptions> theSRLabelOptionsForTTCR = {

      {1, "0-2"}, // nonprompt CR, same as main selection, but off Z or noOSSF
      {2, "3"},
      {3, ">3"},
      {4, "0-2"},
      {5, "3"},
      {6, ">3"},
      {7, "0-2"},
      {8, "3"},
      {9, ">3"},

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

    /*
      {1, "2"},
      {2, "3"},
      {3, ">3"},
      {4, "2"},
      {5, "3"},
      {6, ">3"},
      {7, "2"},
      {8, "3"},
      {9, ">3"},
      */
      {1, "2"},
      {2, "3"},
      {3, "4"},
      {4, ">4"},
      {5, "2"},
      {6, "3"},
      {7, "4"},
      {8, ">4"},
    };

std::vector<BinLabelOptions> theSRLabelOptionsFor4L = {

      {1, "0"},
      {2, ">0"},
      
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
      {2, "#mu#mu#mu e"},
      {3, "#mu#mu ee"},
      {4, "#mu eee"},
      {5, "eeee"},
      
    };

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

// here we define histos
struct histInfo{
  std::string fancyName;
  int index;
  //std::string usualName;
  double varMin;
  double varMax;
  int nBins;
  bool isEnVar;
};

std::map<TString, histInfo> figNames  =         {{"ptlead",  {"Leading lepton p_{T} [GeV]", 0, 0, 300, 30, true}},
                                                 {"sublead", {"Sub-leading lepton p_{T} [GeV]", 1, 0, 200, 20, true}},
                                                 {"trail",   {"Trailing lepton p_{T} [GeV]", 2, 0, 200, 20, true}},
                                                 {"pt4th",   {"4th lepton p_{T} [GeV]", 3, 0, 100, 20, true}},
                                                 {"mtW",     {"m_{T}^{W} [GeV]", 4, 0, 200, 20, true}},
                                                 {"njets",   {"N_{j}", 5, -0.5, 7.5, 8, false}},
                                                 {"nbjets",  {"N_{b}", 6, -0.5, 4.5, 5, false}},
                                                 {"BDTpp",   {"BDT in pp category", 7, -1, 1, 10, false}},
                                                 //{"flavour", {"", 8, 0.5, (leptonSelectionAnalysis == 2 ? (static_cast<double>(flavourLabelOptionsFor2L.size()) + 0.5) : (leptonSelectionAnalysis == 3 ? (static_cast<double>(flavourLabelOptionsFor3L.size()) + 0.5) : (static_cast<double>(flavourLabelOptionsFor4L.size()) + 0.5))), (leptonSelectionAnalysis == 2 ? (static_cast<int>(flavourLabelOptionsFor2L.size())) : (leptonSelectionAnalysis == 3 ? (static_cast<int>(flavourLabelOptionsFor3L.size())) : (static_cast<int>(flavourLabelOptionsFor4L.size())))), false}},
                                                 //{"SR",      {leptonSelectionAnalysis == 2 ? "" : (leptonSelectionAnalysis == 3 ? "N_{j}" : "N_{b}"), 9, -0.5, leptonSelectionAnalysis == 2 ? (static_cast<double>(theSRLabelOptionsFor2L.size()) - 0.5) : (leptonSelectionAnalysis == 3 ? static_cast<double>(theSRLabelOptionsFor3L.size()) - 0.5 : static_cast<double>(theSRLabelOptionsFor4L.size()) - 0.5), leptonSelectionAnalysis == 2 ? (static_cast<int>(theSRLabelOptionsFor2L.size())) : (leptonSelectionAnalysis == 3 ? static_cast<int>(theSRLabelOptionsFor3L.size()) : static_cast<int>(theSRLabelOptionsFor4L.size())), false}},
                                                 {"SR3L",    {"N_{j}", 8, -0.5, static_cast<double>(theSRLabelOptionsFor3L.size()) - 0.5, static_cast<int>(theSRLabelOptionsFor3L.size()), false}},
                                                 {"SR4L",    {"N_{b}", 9, -0.5, static_cast<double>(theSRLabelOptionsFor4L.size()) - 0.5, static_cast<int>(theSRLabelOptionsFor4L.size()), false}},
                                                 {"mll",     {"M(ll) [GeV]", 10, 81., 101., 10, true}},
                                                 {"ptZ",     {"p_{T}^{Z} [GeV]", 11, 0, 400, 16, true}},
                                                 {"ptNonZ",  {"Non-Z lepton p_{T} [GeV]", 12, 0, 200, 20, true}},
                                                 {"mll3e",   {"M(ll) in 3e [GeV]", 13, 0, 200, 20, true}},
                                                 {"mll2e1mu",{"M(ll) in 2e1mu [GeV]", 14, 0, 200, 20, true}},
                                                 {"mll1e2mu",{"M(ll) in 1e2mu [GeV]", 15, 0, 200, 20, true}},
                                                 {"mll3mu",  {"M(ll) in 3mu [GeV]", 16, 0, 200, 20, true}},
                                                 {"met",     {"E_{T}^{miss} [GeV]", 17, 0, 300, 15, true}},
                                                 {"deltaR",  {"#Delta R(jet, trailing lepton)", 18, 0.4, 3., 13, false}},
                                                 {"deltaRlead",  {"#Delta R(jet, leading lepton)", 19, 0.4, 3., 13, false}},
                                                 {"mtLeading",{"Leading lepton M_{T} [GeV]", 20, 0, 300, 20, true}},
                                                 {"mtTrailing",{"Trailing lepton M_{T} [GeV]", 21, 0, 200, 20, true}},
                                                 {"leadJetPt", {"Leading non-b jet p_{T} [GeV]", 22, 30, 310, 14, true}},
                                                 {"trailJetPt", {"Trailing non-b jet p_{T} [GeV]", 23, 30, 310, 14, true}},
                                                 {"SRnpCR", {"", 24, -0.5, 2.5, 3, false}},
                                                 {"nPV", {"number of PV", 25, -0.5, 49.5, 25, false}},
                                                 {"mlll", {"M(lll) [GeV]", 26, 81., 101., 10, true}},
                                                 {"etaLead", {"Leading lepton #eta", 27, -2.5, 2.5, 20, false}},
                                                 {"etaSubl", {"Sub-leading lepton #eta", 28, -2.5, 2.5, 20, false}},
                                                 {"etaTrail", {"Trailing lepton #eta", 29, -2.5, 2.5, 20, false}},
                                                 {"eta4th", {"4th lepton #eta", 30, -2.5, 2.5, 20, false}},
                                                 {"mt_3m", {"m_{T}^{W} in 3#mu [GeV]", 31, 0, 200, 20, true}},
                                                 {"mt_2m1e", {"m_{T}^{W} in 2#mu e [GeV]", 32, 0, 200, 20, true}}, 
                                                 {"mt_1m2e", {"m_{T}^{W} in 1#mu 2e [GeV]", 33, 0, 200, 20, true}}, 
                                                 {"mt_3e", {"m_{T}^{W} in 3e [GeV]", 34, 0, 200, 20, true}},
                                                 {"cosThetaStar", {"cos(#theta^{*})", 35, -1, 1, 5, false}},
                                                 {"mll_ss",  {"Invariant mass of ss 2l pair [GeV]", 36, 0, 300, 20, true}},
                                                 {"chargeOfLeptons",  {"Charge of the leptons in ss2l channel", 37, -1.5, 1.5, 3, false}},
                                                 {"ll_deltaR",  {"#Delta R(leading lepton, trailing lepton)", 38, 0, 7., 35, false}},
                                                 {"mt2ll_ss",  {"M_{T2}^{ll} [GeV]", 39, 0, 200., 20, true}},
                                                 {"BDTmm",   {"BDT in mm category", 40, -1, 1, 10, false}},
                                                 {"HT",     {"H_{T} [GeV]", 41, 0, 400, 20, true}},
                                                 {"SRallTTZ", {"", 42, -0.5, static_cast<double>(theSRLabelOptionsForTTZ.size()) - 0.5, static_cast<int>(theSRLabelOptionsForTTZ.size()), false}},
                                                 {"SRWZCR", {"N_{j}", 43, -0.5, static_cast<double>(theSRLabelOptionsForWZCR.size()) - 0.5, static_cast<int>(theSRLabelOptionsForWZCR.size()), false}},
                                                 {"SRZZCR", {"", 44, -0.5, static_cast<double>(theSRLabelOptionsForZZCR.size()) - 0.5, static_cast<int>(theSRLabelOptionsForZZCR.size()), false}},
                                                 {"SRTTCR", {"N_{j}", 45, -0.5, static_cast<double>(theSRLabelOptionsForTTCR.size()) - 0.5, static_cast<int>(theSRLabelOptionsForTTCR.size()), false}},
                                                 {"SRttZCleanPTZ", {"p_{T}^{Z} [GeV]", 46, -0.5, static_cast<double>(theSRLabelOptionsForttZCleanPTZ.size()) - 0.5, static_cast<int>(theSRLabelOptionsForttZCleanPTZ.size()), true}},
                                                 {"SRttZCleanCosTheta", {"cos(#theta^{*})", 47, -0.5, static_cast<double>(theSRLabelOptionsForttZCleanCosTheta.size()) - 0.5, static_cast<int>(theSRLabelOptionsForttZCleanCosTheta.size()), false}},
                                                 {"flavour3L", {"", 48, 0.5, static_cast<double>(flavourLabelOptionsFor3L.size()) + 0.5, static_cast<int>(flavourLabelOptionsFor3L.size()), false}},
                                                 {"flavour4L", {"", 49, 0.5, static_cast<double>(flavourLabelOptionsFor4L.size()) + 0.5, static_cast<int>(flavourLabelOptionsFor4L.size()), false}},
                                           };

const int nVars  = figNames.size() ;

const unsigned int indexSRttZcleanCosTheta = 47;
const unsigned int indexSRttZcleanPTZ = 46;
const unsigned int indexSRTTCR = 45;
const unsigned int indexSRZZCR = 44;
const unsigned int indexSRWZCR = 43;
const unsigned int indexSRTTZ = 42;
const unsigned int indexSR3L = 8;
const unsigned int indexSR4L = 9;

const unsigned int indexFlavour3L = 48;
const unsigned int indexFlavour4L = 49;

std::map<std::string, std::vector<TString>> listToPrint;

std::map<std::string, double> uncOnNorm = {{"ttZ", 0.0}, 
                                           {"ttH", 0.11}, 
                                           {"ttX", 0.11}, 
                                           {"WZ", 0.1}, 
                                           {"ZZ", 0.2}, 
                                           {"Xgamma", 0.2}, 
                                           {"rare", 0.5}, 
};

std::vector<double> ttZISRUpW = {1.03, 1.02, 1.005, 0.98, 1.02, 1.01, 0.99, 0.945, 1.035, 1.025, 1.00, 0.95, 0.98, 0.995}; // 14 SRs
std::vector<double> ttZISRDownW = {0.965, 0.975, 1.00, 1.025, 0.98, 0.99, 1.015, 1.07, 0.96, 0.97, 1.00, 1.065, 1.025, 1.005}; // 14 SRs

std::vector<double> ttZFSRUpW = {0.96, 0.965, 0.99, 0.985, 0.99, 1.00, 1.00, 1.00, 1.01, 1.02, 1.025, 1.03, 0.98, 1.015}; // 14 SRs
std::vector<double> ttZFSRDownW = {1.07, 1.06, 1.02, 1.005, 1.025, 1.01, 0.99, 0.975, 0.98, 0.97, 0.96, 0.97, 1.075, 0.98}; // 14 SRs

#endif 
