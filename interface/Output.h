#ifndef OUTPUT_H
#define OUTPUT_H

#include "readTreeSync.h"

namespace Output{

    struct Unc{
      Unc():unc(numberOfSyst){}
      std::vector<TH1D> unc;
      void FillUnc(double value, int number, double lastBin, double weight){
        unc.at(number).Fill(TMath::Min(value, lastBin-0.1), weight);
      };
    };

	struct DistribsAll{
  		DistribsAll():vectorHisto(nSamples+1), vectorHistoUncUp(nSamples+1), vectorHistoUncDown(nSamples+1){}
  		std::vector<TH1D> vectorHisto;
  		TH1D vectorHistoTotalUnc;
  		std::vector<Unc> vectorHistoUncUp;
  		std::vector<Unc> vectorHistoUncDown;
  		THStack stack;
	};

	std::vector<DistribsAll> distribs(nVars);
};

#endif
