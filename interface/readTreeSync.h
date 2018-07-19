#ifndef readTreeSync_H
#define readTreeSync_H

#include <iostream>
#include <fstream>

const int nFlavors = 2;

const int leptonSelectionAnalysis = 3;

const int nSamples = 50;
const int dataSample = 0;

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
