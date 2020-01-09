#ifndef OUTPUT_H
#define OUTPUT_H

#include "readTreeSync.h"

#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>

namespace Output{

  struct Unc{
    Unc():unc(numberOfSyst){}
    std::vector<TH1D> unc;
    void FillUnc(double value, int number, double lastBin, double weight){
      unc.at(number).Fill(TMath::Min(value, lastBin-0.1), weight);
    };
  };

  struct PDF{
    PDF():var(100){}
    std::vector<TH1D> var;
  };

  struct DistribsAll{
  	DistribsAll():vectorHisto(nProcesses+1), vectorHistoUncUp(nProcesses+1), vectorHistoUncDown(nProcesses+1), vectorHistoPDF(nProcesses+1){}
  	std::vector<TH1D> vectorHisto;
  	TH1D vectorHistoTotalUnc;
  	std::vector<Unc> vectorHistoUncUp;
  	std::vector<Unc> vectorHistoUncDown;
  	std::vector<PDF> vectorHistoPDF;
  	THStack stack;
  };

  struct DistribsAll2D{
  	// 2 eta categories \times pass or fail \times nonprompt stemming from all flavours, b, c and light = 16
    DistribsAll2D():vectorHisto(16){}
    std::vector<TH2D> vectorHisto;
  };

  struct DistribsAllForFR{
  	// 2 eta categories \times pass or fail \times nonprompt stemming from all flavours, b, c and light = 16
    DistribsAllForFR():vectorHisto(16){}
    std::vector<TH1D> vectorHisto;
  };

  struct DistribsAllForCT{
  	// only 2 categories: events that passed full selection and events used for nonprompt background estimation
    DistribsAllForCT():vectorHisto(2){}
    std::vector<TH1D> vectorHisto;
  };

  std::vector<DistribsAll> distribs(nVars);
  // separately FR for electrons and muons
  std::vector<DistribsAllForFR> distribs1DForFR(2);

  std::vector<DistribsAllForCT> distribs1DForCT(nVars);
  std::vector<DistribsAll2D> distribs2D(2);

  // Histograms needed to measure FR in data
  TH2D* fakeMapsCalc[nFlavors][3][nPt-1][nEta-2][nProcesses+1]; // nEta - 2 (4-2) -> consider only barrel and endcap, 3 here stands for 3 categories (status): passed, all, passed / all
  TH1D mtMaps[nFlavors][nPt-1][nEta-2][nProcesses+1];
  THStack* mtStack[nFlavors][nPt-1][nEta-2];
};

#endif
