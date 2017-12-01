#ifndef OUTPUT_H
#define OUTPUT_H

#include "readTreeSync.h"

namespace Output{
	struct DistribsAll{
  		DistribsAll():vectorHisto(nSamples+1), vectorHistoUp(nSamples+1), vectorHistoDown(nSamples+1){}
  		std::vector<TH1D> vectorHisto;
  		TH1D vectorHistoTotalUnc;
  		std::vector<TH1D> vectorHistoUp;
  		std::vector<TH1D> vectorHistoDown;
  		THStack stack;
	};

	std::vector<DistribsAll> distribs(nVars);
};

#endif
