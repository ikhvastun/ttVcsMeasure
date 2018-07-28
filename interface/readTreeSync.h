#ifndef readTreeSync_H
#define readTreeSync_H

#include <iostream>
#include <fstream>

const int nFlavors = 2;

const int leptonSelectionAnalysis = 3;

const int nSamples = 50;
const int dataSample = 0;

const unsigned int indexSR = 9;
const unsigned int indexFlavour = 8;

const int numberOfSyst = 6;

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
      {1, "nB0,nJ2"},
      {2, "nB0,nJ3"},
      {3, "nB0,nJ#geq4"},
      {4, "nB1,nJ2"},
      {5, "nB1,nJ3"},
      {6, "nB1,nJ#geq4"},
      {7, "nB#geq2,nJ2"},
      {8, "nB#geq2,nJ3"},
      {9, "nB#geq2,nJ#geq4"},
      */
      {1, "2"},
      {2, "3"},
      {3, ">3"},
      {4, "2"},
      {5, "3"},
      {6, ">3"},
      {7, "2"},
      {8, "3"},
      {9, ">3"},
      
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
                                                 {"BDT",     {"BDT", 7, -1, 1, 10, false}},
                                                 {"flavour", {"", 8, 0.5, (leptonSelectionAnalysis == 2 ? (static_cast<double>(flavourLabelOptionsFor2L.size()) + 0.5) : (leptonSelectionAnalysis == 3 ? (static_cast<double>(flavourLabelOptionsFor3L.size()) + 0.5) : (static_cast<double>(flavourLabelOptionsFor4L.size()) + 0.5))), (leptonSelectionAnalysis == 2 ? (static_cast<int>(flavourLabelOptionsFor2L.size())) : (leptonSelectionAnalysis == 3 ? (static_cast<int>(flavourLabelOptionsFor3L.size())) : (static_cast<int>(flavourLabelOptionsFor4L.size())))), false}},
                                                 {"SR",      {leptonSelectionAnalysis == 2 ? "" : "N_{j}", 9, -0.5, leptonSelectionAnalysis == 2 ? (static_cast<double>(theSRLabelOptionsFor2L.size()) - 0.5) : (static_cast<double>(theSRLabelOptionsFor3L.size()) - 0.5), leptonSelectionAnalysis == 2 ? (static_cast<int>(theSRLabelOptionsFor2L.size())) : static_cast<int>(theSRLabelOptionsFor3L.size()), false}},
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
                                                 {"cosThetaStar", {"cos(#Theta^{*})", 35, -1, 1, 5, false}},
                                                 {"mll_ss",  {"Invariant mass of ss 2l pair [GeV]", 36, 0, 300, 20, true}},
                                                 {"chargeOfLeptons",  {"Charge of the leptons in ss2l channel", 37, -1.5, 1.5, 3, false}},
                                                 {"ll_deltaR",  {"#Delta R(leading lepton, trailing lepton)", 38, 0, 7., 35, false}},
                                                 {"mt2ll_ss",  {"M_{T2}^{ll} [GeV]", 39, 0, 200., 20, true}}
                                           };

//const int nVars  = 40;
const int nVars  = figNames.size() ;

std::map<std::string, std::vector<TString>> listToPrint;
/*
listToPrint["WZ"] = {"ptlead", "sublead", "trail", "njets", "nbjets", "flavour", "mll", "ptZ", "ptNonZ", "mtW", "mll3e", "mll2e1mu", "mll1e2mu", "mll3mu", "met", "nPV", "mt_3m", "mt_2m1e",  "mt_1m2e", "mt_3e", "cosThetaStar"};
listToPrint["Zgamma"] = {"ptlead", "sublead", "trail", "njets", "nbjets", "flavour", "met", "nPV", "mlll"};
listToPrint["ttbar"] = {"ptlead", "sublead", "trail", "njets", "nbjets", "flavour", "met", "nPV"};
listToPrint["DY"] = {"ptlead", "sublead", "trail", "njets", "nbjets", "flavour", "met", "nPV", "mll", "ptZ", "ptNonZ", "mtW", "mll3e", "mll2e1mu", "mll1e2mu", "mll3mu", "mt_3m", "mt_2m1e",  "mt_1m2e", "mt_3e"};
listToPrint["ttW"] = {"ptlead", "sublead", "njets", "nbjets", "flavour", "met", "nPV", "deltaR", "deltaRlead", "mtLeading", "mtTrailing", "leadJetPt", "trailJetPt", "etaLead", "etaSubl",     "mll_ss", "chargeOfLeptons", "ll_deltaR", "mt2ll_ss", "SR", "BDT"}; 
listToPrint["ZZ"] = {"ptlead", "sublead", "trail", "pt4th", "njets", "nbjets", "flavour", "met", "nPV", "mll", "ptZ", "etaLead", "etaSubl", "etaTrail", "eta4th"};
listToPrint["ttZ3L"] = {"ptlead", "sublead", "trail", "njets", "nbjets", "flavour", "mll", "ptZ", "ptNonZ", "SR", "met", "cosThetaStar"};
*/

#endif 
