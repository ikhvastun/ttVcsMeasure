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
    Era():runEras(4), runErasUncUp(4), runErasUncDown(4){}
    
    std::vector<TH1D> runEras;
    std::vector<Unc> runErasUncUp;
    std::vector<Unc> runErasUncDown;
    void Fill(double value, int era, double valueUncJESUp, double valueUncJESDown, double valueUncJERUp, double valueUncJERDown, double valueUncUnclUp, double valueUncUnclDown, double lastBin, double weight){
        runEras.at(era).Fill(TMath::Min(value, lastBin-0.001), weight);
        runErasUncUp.at(era).FillUnc(valueUncJESUp, valueUncJERUp, valueUncUnclUp, lastBin, weight);
        runErasUncDown.at(era).FillUnc(valueUncJESDown, valueUncJERDown, valueUncUnclDown, lastBin, weight);

        runEras.at(3).Fill(TMath::Min(value, lastBin-0.001), weight);
        runErasUncUp.at(3).FillUnc(valueUncJESUp, valueUncJERUp, valueUncUnclUp, lastBin, weight);
        runErasUncDown.at(3).FillUnc(valueUncJESDown, valueUncJERDown, valueUncUnclDown, lastBin, weight);
    }
  };


	struct DistribsAll{
  		DistribsAll():vectorHisto(nSamples+1), vectorHistoUp(nSamples+1), vectorHistoDown(nSamples+1), vectorHistoEras(nSamples+1){}
  		TH1D vectorHistoTotalUnc;

        std::vector<Era> vectorHistoEras;
      
        //std::vector<std::vector<Era>> vectorHistoErasUp;
        //std::vector<std::vector<Era>> vectorHistoErasDown;
  		
  		std::vector<TH1D> vectorHisto;
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
