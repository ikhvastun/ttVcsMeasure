#ifndef readTreeSync_H
#define readTreeSync_H

#include <iostream>
#include <fstream>

#include "TString.h"

const int nFlavors = 2;

//const int leptonSelectionAnalysis = 3;

const int nProcesses = 11;
const int dataSample = 0;
const int ttWSample = 1;

const int numberOfSyst = 24;
//const int numberOfSyst = 3;
const int pdfUncIndex = 12;

TString flavorsString[2] = {"el", "mu"};
TString additionalString[2] = {"_NC", ""};
TString flavorComposString[4] = {"all", "b", "c", "light"};
//TString flavorComposString[4] = {"B_L", "B_C_L", "B_C_T_L", "B_T_L"};


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
      {4, "[0.5, 1]"}
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
      //{5, "1"},// ttZ 3L, nbjets = 1, njets = 2,3,4,>4; nbjets >1, njets = 2,3,4,>4
      {5, "2"},
      {6, "3"},
      {7, "4"},
      {8, "> 4"},
      {9, "2"},
      {10, "3"},
      {11, "4"},
      {12, "> 4"},
      //{14, "nj1nb1"}, // ttZ 4L, nbjets = 0, 1
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
      
      
std::vector<BinLabelOptions> theSRLabelOptionsForTTZ8SR3L = {

      {1, "2"},// ttZ 3L, nbjets = 1, njets = 2,3,4,>4; nbjets >1, njets = 2,3,4,>4
      {2, "3"},
      {3, "4"},
      {4, "> 4"},
      {5, "2"},
      {6, "3"},
      {7, "4"},
      {8, "> 4"},
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
    /*
      {1, "2"},
      {2, "3"},
      {3, "4"},
      {4, ">4"},
      {5, "2"},
      {6, "3"},
      {7, "4"},
      {8, ">4"},
      */

      {1, "1"}, // WZ CR, first 4 bins, nbjets = 0, njets = 1, 2, 3, > 3
      {2, "2"},
      {3, "3"},
      {4, "> 3"},
      //{5, "1"},
      {5, "2"},// ttZ 3L, nbjets = 1, njets = 2,3,4,>4; nbjets >1, njets = 2,3,4,>4
      {6, "3"},
      {7, "4"},
      {8, "> 4"},
      {9, "2"},
      {10, "3"},
      {11, "4"},
      {12, "> 4"},
    };

std::vector<BinLabelOptions> theSRLabelOptionsFor4L = {

      //{1, "nj1nb1"},
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

std::vector<BinLabelOptions> flavourLabelOptionsFor3L4L = {
      
      {1, "#mu#mu#mu(#mu)"},
      {2, "#mu#mue(#mu)"},
      {3, "#muee(#mu/e)"},
      {4, "eee(e)"},
      
    };

std::vector<BinLabelOptions> flavourLabelOptionsFor3L = {
      
      {1, "#mu#mu#mu"},
      {2, "#mu#mue"},
      {3, "#muee"},
      {4, "eee"},
      
    };

std::vector<BinLabelOptions> flavourLabelOptionsFor4L = {
      
      {1, "#mu#mu#mu#mu"},
      {2, "#mu#mu#mue"},
      {3, "#mu#muee"},
      {4, "#mueee"},
      {5, "eeee"},
      
    };

std::vector<BinLabelOptions> flavourLabelOptionsFor4LZZ = {
      
      {1, "#mu#mu#mu#mu"},
      {2, "#mu#muee"},
      {3, "eeee"},
      
    };

// trees for BDT
TMVA::Reader *readerTTWcsttbar = new TMVA::Reader( "!Color:!Silent" );   
TFile* fileDummy = new TFile("fileDummy.root", "RECREATE");
TTree* signalTree = new TTree("signalTree","signalTree");
TTree* bkgTree = new TTree("bkgTree","bkgTree");

double _jetPt1, _jetEta1, _jetPhi1, _jetE1, _jetCSV1;
double _jetPt2, _jetEta2, _jetPhi2, _jetE2, _jetCSV2;
double _jetPt3, _jetEta3, _jetPhi3, _jetE3, _jetCSV3;
double _jetPt4, _jetEta4, _jetPhi4, _jetE4, _jetCSV4;
double _jetPt5, _jetEta5, _jetPhi5, _jetE5, _jetCSV5;
double _jetPt6, _jetEta6, _jetPhi6, _jetE6, _jetCSV6;

double _lepPt1, _lepEta1, _lepPhi1, _lepE1, _lepCharge1;
double _lepPt2, _lepEta2, _lepPhi2, _lepE2, _lepCharge2;

double _metPt1, _metEta1, _metPhi1, _metE1;

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

// lepton + jet
double minMLeptonJet;
double maxMLeptonJet;
double minDeltaRLeptonJet;
double maxDeltaRLeptonJet;
double minDeltaPhiLeptonJet;
double maxDeltaPhiLeptonJet;
double minpTLeptonJet;
double maxpTLeptonJet;

//lepton bjet
double minMLeptonbJet;
double maxMLeptonbJet;
double minDeltaRLeptonbJet;
double maxDeltaRLeptonbJet;
double minDeltaPhiLeptonbJet;
double maxDeltaPhiLeptonbJet;
double minpTLeptonbJet;
double maxpTLeptonbJet;

//jet jet
double minMJetJet;
double maxMJetJet;
double minDeltaRJetJet;
double maxDeltaRJetJet;
double minDeltaPhiJetJet;
double maxDeltaPhiJetJet;
double minpTJetJet;
double maxpTJetJet;

//lepton + MET 
double minDeltaPhiLeptonMET;
double maxDeltaPhiLeptonMET;
double minmTLeptonMET;
double maxmTLeptonMET;
double minpTLeptonMET;
double maxpTLeptonMET;

//jet + MET
double minDeltaPhiJetMET;
double maxDeltaPhiJetMET;
double minmTJetMET;
double maxmTJetMET;
double minpTJetMET;
double maxpTJetMET;

//bjet + MET 
double minDeltaPhiBJetMET;
double maxDeltaPhiBJetMET;
double minmTBJetMET;
double maxmTBJetMET;
double minpTBJetMET;
double maxpTBJetMET; 

Float_t userHTLoc, user_met, userele_mll, usermt, usermtlow, userleadpt, usertrailpt, userleadeta, usertraileta, userleadingjetpt, usertrailjetpt, userminDeltaRlead, userminDeltaR, usernJLoc, usernBLoc, userchargeOfLeptons, usermll_ss, userll_deltaR, usermt2ll_ss, user_maxMJetJet, user_maxDeltaPhiJetJet, user_maxMLeptonbJet, user_maxMLeptonJet, user_maxpTLeptonJet, user_maxpTLeptonbJet, user_maxpTJetJet, user_minDeltaPhiJetJet, user_minDeltaRJetJet, user_minMJetJet;

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

std::map<TString, histInfo> figNames  =         {
                                                 {"ptlead",  {"Leading lepton p_{T} [GeV]", 0, 40, 300, 13, true}},
                                                 {"sublead", {"Sub-leading lepton p_{T} [GeV]", 1, 20, 200, 18, true}},
                                                 {"trail",   {"Trailing lepton p_{T} [GeV]", 2, 10, 120, 11, true}},
                                                 {"pt4th",   {"4th lepton p_{T} [GeV]", 3, 10, 100, 18, true}},
                                                 {"mtW",     {"m_{T}^{W} [GeV]", 4, 0, 200, 20, true}},
                                                 {"njets",   {"N_{j}", 5, -0.5, 7.5, 8, false}},
                                                 {"nbjets",  {"N_{b}", 6, -0.5, 3.5, 4, false}},
                                                 {"BDTpp",   {"BDT in pp category", 7, -1, 1, 10, false}},
                                                 {"SR3L",    {"N_{j}", 8, -0.5, static_cast<double>(theSRLabelOptionsFor3L.size()) - 0.5, static_cast<int>(theSRLabelOptionsFor3L.size()), false}},
                                                 {"SR4L",    {"N_{b}", 9, -0.5, static_cast<double>(theSRLabelOptionsFor4L.size()) - 0.5, static_cast<int>(theSRLabelOptionsFor4L.size()), false}},
                                                 {"mll",     {"M(ll) [GeV]", 10, 81., 101., 10, true}},
                                                 {"ptZ",     {"p_{T}^{Z} [GeV]", 11, 0, 400, 16, true}},
                                                 {"ptNonZ",  {"Non-Z lepton p_{T} [GeV]", 12, 10, 200, 19, true}},
                                                 {"mll3e",   {"M(ll) in 3e [GeV]", 13, 81, 101, 10, true}},
                                                 {"mll2e1mu",{"M(ll) in 2e1mu [GeV]", 14, 81, 101, 10, true}},
                                                 {"mll1e2mu",{"M(ll) in 1e2mu [GeV]", 15, 81, 101, 10, true}},
                                                 {"mll3mu",  {"M(ll) in 3mu [GeV]", 16, 81, 101, 10, true}},
                                                 {"met",     {"E_{T}^{miss} [GeV]", 17, 0, 300, 15, true}},
                                                 {"deltaR",  {"#Delta R(jet, trailing lepton)", 18, 0.4, 3., 13, false}},
                                                 {"deltaRlead",  {"#Delta R(jet, leading lepton)", 19, 0.4, 3., 13, false}},
                                                 {"mtLeading",{"Leading lepton M_{T} [GeV]", 20, 0, 300, 20, true}},
                                                 {"mtTrailing",{"Trailing lepton M_{T} [GeV]", 21, 0, 200, 20, true}},
                                                 {"leadJetPt", {"Leading non-b jet p_{T} [GeV]", 22, 30, 310, 14, true}},
                                                 {"trailJetPt", {"Trailing non-b jet p_{T} [GeV]", 23, 30, 310, 14, true}},
                                                 {"SRnpCR", {"", 24, -0.5, 2.5, 3, false}},
                                                 {"nPV", {"number of PV", 25, -0.5, 49.5, 25, false}},
                                                 {"mlll", {"M(lll) [GeV]", 26, 0, 200, 50, true}},
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
                                                 {"flavour4LZZ", {"", 50, 0.5, static_cast<double>(flavourLabelOptionsFor4LZZ.size()) + 0.5, static_cast<int>(flavourLabelOptionsFor4LZZ.size()), false}},
                                                 {"SR3L3m",    {"N_{j}", 51, -0.5, static_cast<double>(theSRLabelOptionsFor3L.size()) - 0.5, static_cast<int>(theSRLabelOptionsFor3L.size()), false}},
                                                 {"SR3L2m1e",    {"N_{j}", 52, -0.5, static_cast<double>(theSRLabelOptionsFor3L.size()) - 0.5, static_cast<int>(theSRLabelOptionsFor3L.size()), false}},
                                                 {"SR3L1m2e",    {"N_{j}", 53, -0.5, static_cast<double>(theSRLabelOptionsFor3L.size()) - 0.5, static_cast<int>(theSRLabelOptionsFor3L.size()), false}},
                                                 {"SR3L3e",    {"N_{j}", 54, -0.5, static_cast<double>(theSRLabelOptionsFor3L.size()) - 0.5, static_cast<int>(theSRLabelOptionsFor3L.size()), false}},
                                                 {"flavour3L4L", {"", 55, 0.5, static_cast<double>(flavourLabelOptionsFor3L4L.size()) + 0.5, static_cast<int>(flavourLabelOptionsFor3L4L.size()), false}},
                                                 {"SRTTZ8SR3L", {"", 56, -0.5, static_cast<double>(theSRLabelOptionsForTTZ8SR3L.size()) - 0.5, static_cast<int>(theSRLabelOptionsForTTZ8SR3L.size()), false}},
                                                 {"mllnoZcut",     {"M(ll) [GeV]", 57, 40., 140., 20, true}},

                                                 {"minMLeptonJet", {"M_{l + jet}^{min} (GeV)", 58, 0, 200, 30, true}},
                                                 {"maxMLeptonJet", {"M_{l + jet}^{max} (GeV)", 59, 0, 800, 30, true}},
                                                 {"minDeltaRLeptonJet", {"min(#DeltaR(l, jet))", 60, 0, 7, 30, false}},
                                                 {"maxDeltaRLeptonJet", {"max(#DeltaR(l, jet))", 61, 0, 10, 30, false}},
                                                 {"minDeltaPhiLeptonJet", {"min(#Delta#Phi(l, jet))", 62, 0, 3.15, 30, false}},
                                                 {"maxDeltaPhiLeptonJet", {"max(#Delta#Phi(l, jet))", 63, 0, 3.15, 30, false}},
                                                 {"minpTLeptonJet", {"P_{T}^{min}(l + jet) (GeV)", 64, 0, 200, 30, true}},
                                                 {"maxpTLeptonJet", {"P_{T}^{max}(l + jet) (GeV)", 65, 0, 300, 30, true}},

                                                 {"minMLeptonbJet", {"M_{l + bjet}^{min} (GeV)", 66, 0, 300, 30, true}},
                                                 {"maxMLeptonbJet", {"M_{l + bjet}^{max} (GeV)", 67, 0, 800, 30, true}},
                                                 {"minDeltaRLeptonbJet", {"min(#DeltaR(l, bjet))", 68, 0, 7, 30, false}},
                                                 {"maxDeltaRLeptonbJet", {"max(#DeltaR(l, bjet))", 69, 0, 10, 30, false}},
                                                 {"minDeltaPhiLeptonbJet", {"min(#Delta#Phi(l, bjet))", 70, 0, 3.15, 30, false}},
                                                 {"maxDeltaPhiLeptonbJet", {"max(#Delta#Phi(l, bjet))", 71, 0, 3.15, 30, false}},
                                                 {"minpTLeptonbJet", {"P_{T}^{min}(l + bjet) (GeV)", 72, 0, 200, 30, true}},
                                                 {"maxpTLeptonbJet", {"P_{T}^{max}(l + bjet) (GeV)", 73, 0, 300, 30, true}},

                                                 {"minMJetJet", {"M_{jet + jet}^{min} (GeV)", 74, 0, 600, 30, true}},
                                                 {"maxMJetJet", {"M_{jet + jet}^{max} (GeV)", 75, 0, 1200, 30, true}},
                                                 {"minDeltaRJetJet", {"min(#DeltaR(jet, jet))", 76, 0, 7, 30, false}},
                                                 {"maxDeltaRJetJet", {"max(#DeltaR(jet, jet))", 77, 0, 10, 30, false}},
                                                 {"minDeltaPhiJetJet", {"min(#Delta#Phi(jet, jet))", 78, 0, 3.15, 30, false}},
                                                 {"maxDeltaPhiJetJet", {"max(#Delta#Phi(jet, jet))", 79, 0, 3.15, 30, false}},
                                                 {"minpTJetJet", {"P_{T}^{min}(jet + jet) (GeV)", 80, 0, 200, 30, true}},
                                                 {"maxpTJetJet", {"P_{T}^{max}(jet + jet) (GeV)", 81, 0, 300, 30, true}},

                                                 {"minDeltaPhiLeptonMET", {"min(#Delta#Phi(l, MET))", 82, 0, 3.15, 30, false}},
                                                 {"maxDeltaPhiLeptonMET", {"max(#Delta#Phi(l, MET))", 83, 0, 3.15, 30, false}},
                                                 {"minmTLeptonMET", {"min(M_{T}(l + MET)) (GeV)", 84, 0, 300, 30, true}},
                                                 {"maxmTLeptonMET", {"max(M_{T}(l + MET)) (GeV)", 85, 0, 400, 30, true}},
                                                 {"minpTLeptonMET", {"min(P_{T}(l + MET))", 86, 0, 300, 30, true}},
                                                 {"maxpTLeptonMET", {"max(P_{T}(l + MET))", 87, 0, 400, 30, true}},

                                                 {"minDeltaPhiJetMET", {"min(#Delta#Phi(jet, MET))", 88, 0, 3.15, 30, false}},
                                                 {"maxDeltaPhiJetMET", {"max(#Delta#Phi(jet, MET))", 89, 0, 3.15, 30, false}},
                                                 {"minmTJetMET", {"min(M_{T}(jet + MET)) (GeV)", 90, 0, 300, 30, true}},
                                                 {"maxmTJetMET", {"max(M_{T}(jet + MET)) (GeV)", 91, 0, 400, 30, true}},
                                                 {"minpTJetMET", {"min(P_{T}(jet + MET))", 92, 0, 300, 30, true}},
                                                 {"maxpTJetMET", {"max(P_{T}(jet + MET))", 93, 0, 400, 30, true}},

                                                 {"minDeltaPhiBJetMET", {"min(#Delta#Phi(bjet, MET))", 94, 0, 3.15, 30, false}},
                                                 {"maxDeltaPhiBJetMET", {"max(#Delta#Phi(bjet, MET))", 95, 0, 3.15, 30, false}},
                                                 {"minmTBJetMET", {"min(M_{T}(bjet + MET)) (GeV)", 96, 0, 300, 30, true}},
                                                 {"maxmTBJetMET", {"max(M_{T}(bjet + MET)) (GeV)", 97, 0, 400, 30, true}},
                                                 {"minpTBJetMET", {"min(P_{T}(bjet + MET))", 98, 0, 300, 30, true}},
                                                 {"maxpTBJetMET", {"max(P_{T}(bjet + MET))", 99, 0, 400, 30, true}},
                                                 {"ptMuonPassedTight",  {"Leading lepton p_{T} [GeV]", 100, 10, 200, 38, true}},
                                                 {"etaMuonPassedTight",  {"Leading lepton #eta", 101, -2.5, 2.5, 20, false}},

                                                 {"ptLepPassedLooseForEff",  {"Leading lepton p_{T} [GeV]", 102, 10, 200, 38, true}},
                                                 {"etaLepPassedLooseForEff",  {"Leading lepton #eta", 103, -2.5, 2.5, 20, false}},
                                                 {"ptLepPassedTightForEff",  {"Leading lepton p_{T} [GeV]", 104, 10, 200, 38, true}},
                                                 {"etaLepPassedTightForEff",  {"Leading lepton #eta", 105, -2.5, 2.5, 20, false}},

                                                 {"etaLepPassedLooseForEffLowPt",  {"Muon #eta (p_{T} < 25 GeV)", 106, -2.5, 2.5, 20, false}},
                                                 {"etaLepPassedTightForEffLowPt",  {"Muon #eta (p_{T} < 25 GeV)", 107, -2.5, 2.5, 20, false}},
                                                 {"etaLepPassedLooseForEffHighPt",  {"Muon #eta (p_{T} > 25 GeV)", 108, -2.5, 2.5, 20, false}},
                                                 {"etaLepPassedTightForEffHighPt",  {"Muon #eta (p_{T} > 25 GeV)", 109, -2.5, 2.5, 20, false}},
                                           };


const int nVars  = figNames.size() ;

const unsigned int indexSR3L3m = 51;
const unsigned int indexSR3L2m1e = 52;
const unsigned int indexSR3L1m2e = 53;
const unsigned int indexSR3L3e = 54;
const unsigned int indexFlavour3L4L = 55;
const unsigned int indexSRTTZ8SR3L = 56;

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
const unsigned int indexFlavour4LZZ = 50;

const unsigned int indexLeadPt = 0;
const unsigned int indexTrailPt = 2;

std::map<std::string, std::vector<TString>> listToPrint;

std::map<std::string, double> uncOnNorm = {{"ttZ", 0.0}, 
                                           {"ttW", 0.11}, 
                                           {"ttH", 0.11}, 
                                           {"ttX", 0.11}, 
                                           {"WZ", 0.1}, 
                                           {"ZZ", 0.1}, 
                                           {"Xgamma", 0.2}, 
                                           {"rare", 0.5}, 
};

/*
// effect of ISR and FSR for ttZ to be varied by factor 2
std::vector<double> ttZISRUpW = {1.03, 1.02, 1.005, 0.98, 1.02, 1.01, 0.99, 0.945, 1.035, 1.025, 1.00, 0.95, 0.98, 0.995}; // 14 SRs
std::vector<double> ttZISRDownW = {0.965, 0.975, 1.00, 1.025, 0.98, 0.99, 1.015, 1.07, 0.96, 0.97, 1.00, 1.065, 1.025, 1.005}; // 14 SRs

std::vector<double> ttZFSRUpW = {0.96, 0.965, 0.99, 0.985, 0.99, 1.00, 1.00, 1.00, 1.01, 1.02, 1.025, 1.03, 0.98, 1.015}; // 14 SRs
std::vector<double> ttZFSRDownW = {1.07, 1.06, 1.02, 1.005, 1.025, 1.01, 0.99, 0.975, 0.98, 0.97, 0.96, 0.97, 1.075, 0.98}; // 14 SRs
*/

std::vector<double> ttZISRUpW = {1.03, 1.025, 1.01, 1.045, 1.03, 1.015, 1.015, 1.07, 1.035, 1.025, 1.015, 1.065, 1.025, 1.01}; // 14 SRs
std::vector<double> ttZISRDownW = {0.965, 0.975, 0.995, 0.965, 0.965, 0.985, 0.99, 0.945, 0.96, 0.97, 0.98, 0.95, 0.98, 0.995}; // 14 SRs

std::vector<double> ttZFSRUpW = {1.03, 1.025, 1.01, 1.045, 1.03, 1.015, 1.015, 1.07, 1.035, 1.025, 1.015, 1.065, 1.025, 1.01}; // 14 SRs
std::vector<double> ttZFSRDownW = {0.965, 0.975, 0.995, 0.965, 0.965, 0.985, 0.99, 0.945, 0.96, 0.97, 0.98, 0.95, 0.98, 0.995}; // 14 SRs

//std::vector<double> osToss = {1.28, 1.23, 1.24, 1.17, 1.1, 1.}; // for nbjets >= 1
std::vector<double> osToss = {1.26, 1.23, 1.25, 1.15, 1.08, 1.}; // for nbjets = 1
//std::vector<double> osToss = {1.37, 1.22, 1.2, 1.19, 1.12, 1.}; // for nbjets >= 2
//std::vector<double> osToss = {1.19, 1.13, 1.22, 1.17, 1.1, 1.}; // for nbjets = 0 for ttbar


// pt and eta ranges for FR measurement 
const int nPt = 7;
const double ptBins[nPt] = {10., 15., 20., 30., 45., 65., 100.};

const int nEta = 4;
double etaBins[2][nEta] = {{0., 0.8, 1.442, 2.5}, {0., 1.2, 2.1, 2.4}};

double borderOfBarrelEndcap[2] = {1.479, 1.2};

#endif 
