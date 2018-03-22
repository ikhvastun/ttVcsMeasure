#ifndef OUTPUT_H
#define OUTPUT_H

#include "readTreeSync.h"

namespace Output{

  struct Unc
  {
    Unc():unc(3){}
    std::vector<TH1D> unc;
    void FillUnc(double uncJES, double uncJER, double uncUncl, double& lastBin, double & weight){
      unc.at(0).Fill(TMath::Min(uncJES, lastBin-0.1), weight);
      unc.at(1).Fill(TMath::Min(uncJER, lastBin-0.1), weight);
      unc.at(2).Fill(TMath::Min(uncUncl, lastBin-0.1), weight);
    };
  };

  struct Era
  {
    Era():runEras(3), runErasUncUp(3), runErasUncDown(3){}
    
    std::vector<TH1D> runEras;
    std::vector<Unc> runErasUncUp;
    std::vector<Unc> runErasUncDown;
  };


	struct DistribsAll{
  		DistribsAll():vectorHisto(nSamples+1), vectorHistoUp(nSamples+1), vectorHistoDown(nSamples+1), vectorHistoEras(nSamples+1){}
  		std::vector<TH1D> vectorHisto;
  		TH1D vectorHistoTotalUnc;

      std::vector<Era> vectorHistoEras;
      
      //std::vector<std::vector<Era>> vectorHistoErasUp;
      //std::vector<std::vector<Era>> vectorHistoErasDown;
  		
      std::vector<Unc> vectorHistoUp;
  		std::vector<Unc> vectorHistoDown;

  		TH1D histDataEras[3][2];
  		THStack stack;
      THStack stackJESUp;
      THStack stackJESDown;
      THStack stackJERUp;
      THStack stackJERDown;
      THStack stackUnclUp;
      THStack stackUnclDown;
	};


	std::vector<DistribsAll> distribs(nVars);
};

#endif
