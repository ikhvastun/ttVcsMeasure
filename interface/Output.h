#ifndef OUTPUT_H
#define OUTPUT_H

#include "readTreeSync.h"

#include <TH1D.h>
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
        DistribsAll2D():vectorHisto(nSamples){}
        std::vector<TH2D> vectorHisto;
    };

    DistribsAll2D distribs2D;
	std::vector<DistribsAll> distribs(nVars);
};

#endif
