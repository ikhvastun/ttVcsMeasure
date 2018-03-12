#ifndef OUTPUT_H
#define OUTPUT_H

#include "readTreeSync.h"

namespace Output{
	struct DistribsAll{
  		DistribsAll():vectorHisto(nSamples+1), vectorHistoUp(nSamples+1), vectorHistoDown(nSamples+1), vectorHisto2D(nSamples+1){}
  		std::vector<TH1D> vectorHisto;
  		TH1D vectorHistoTotalUnc;
  		std::vector<TH1D> vectorHistoUp;
  		std::vector<TH1D> vectorHistoDown;
  		std::vector<TH2D> vectorHisto2D;
  		THStack stack;
	};
	struct DistribsAll2D{
  		DistribsAll2D():vectorHisto(nSamples){}
  		std::vector<TH2D> vectorHisto;
    };

	std::vector<DistribsAll> distribs(nVars);
    DistribsAll2D distribs2D;
	//std::vector<DistribsAll2D> distribsForFR(nVars2D);
};

#endif
