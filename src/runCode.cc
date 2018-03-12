#include <algorithm>
#include <vector>
#include <map>
#include <iomanip>

#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TApplication.h"
#include "TFrame.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "../interface/BTagCalibrationStandalone.h" 

#include "../interface/showHist.h"
#include "../interface/readTreeSync.h"
#include "../interface/Tools.h"
#include "../interface/analysisTools.h"
#include "../interface/Output.h"

#include "../interface/errors.h"
#include "../interface/treeReader.h"

#include "../interface/analysisTools.h"
#include "../interface/fillDatacards.h"

#include "tdrStyle.C"

using namespace std;
using namespace tools;

Errors LastError::lasterror = Errors::UNKNOWN;
using Output::distribs;
using Output::distribs2D;

void treeReader::Analyze(){

  leptonSelection = leptonSelectionAnalysis;
  //Set CMS plotting style
  setTDRStyle();
  gROOT->SetBatch(kTRUE);
  //read samples and cross sections from txt file
  readSamples("data/samples_FOtuning.txt"); // 
  
  std::vector<std::string> namesOfSamples = treeReader::getNamesOfTheSample();
  initdistribs(namesOfSamples);

  addVariablesToBDT();

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample();

      Color_t color = assignColor(std::get<0>(samples[sam]));
      setStackColors(color, sam);

      std::cout<<"Entries in "<< std::get<1>(samples[sam]) << " " << nEntries << std::endl;
      double progress = 0;  //for printing progress bar
      for(long unsigned it = 0; it < nEntries; ++it){
          //print progress bar  
          if(it%100 == 0 && it != 0){
            progress += (double) (100./nEntries);
            tools::printProgress(progress);
          } 
          else if(it == nEntries -1){
              progress = 1.;
              tools::printProgress(progress);
          }

          GetEntry(it);
          //if(it > 100000) break;
          
          std::vector<unsigned> ind, indLoose;
          const unsigned lCount = selectLep(ind);
          const unsigned lCountLoose = selectLooseLep(indLoose);

          if(lCountLoose < 1) continue;

          int featureCategory = -99;
          int nLocEle = getElectronNumber(indLoose);

          std::vector<unsigned> indJets;

          unsigned third = -9999;
          double mll = 99999;
          double pt_Z = 999999;
          double phi_Z = 999999;
          double ptNonZ = 999999;

          nJLoc = nJets(0, true, indJets, std::get<0>(samples[sam]) == "nonpromptData");
          nBLoc = nBJets(0, false, true, 1, std::get<0>(samples[sam]) == "nonpromptData");
          //double dMZ = deltaMZ(ind, third, mll, pt_Z, ptNonZ, phi_Z);

          for(auto & i : indLoose){
            double mvaVL = -999.;
            vector<Float_t> varForBDT = { (Float_t)_lPt[i], (Float_t)_lEta[i], (Float_t)_selectedTrackMult[i], (Float_t)_miniIsoCharged[i], (Float_t)(_miniIso[i] - _miniIsoCharged[i]), (Float_t)_ptRel[i], (Float_t)_ptRatio[i], (Float_t)_relIso[i], 
            (Float_t)(_closestJetDeepCsv_bb[i] + _closestJetDeepCsv_b[i]),
            //(Float_t)_closestJetCsvV2[i], 
            (Float_t)_3dIPSig[i], (Float_t)_dxy[i], (Float_t)_dz[i], _lFlavor[i] == 0 ? (Float_t)_lElectronMva[i] : (Float_t)_lMuonSegComp[i]};
            fillBDTvariables(varForBDT, _lFlavor[i]);
            mvaVL =  _leptonMvatZqTTV[i];
            //mvaVL =  _lFlavor[i] == 0 ? readerLeptonMVAele->EvaluateMVA( "BDTG method") : readerLeptonMVAmu->EvaluateMVA( "BDTG method");

            //if(_lFlavor[i] == 1 && !_lPOGMedium[i]) continue;
            if(_lProvenanceCompressed[i] == 0) continue;
            if (mvaVL > 0.85) featureCategory = 0;
            else if (mvaVL > -0.8) featureCategory = 1;
            else continue;

            //if(_gen_partonPt[i] == 0) continue;
            distribs[_lFlavor[i]].vectorHisto[featureCategory].Fill(TMath::Min(double(featureCategory == 0 ? _lPt[i] : _lPt[i] / _ptRatio[i]),varMax[0]-0.001), weight);
            distribs[_lFlavor[i]].vectorHisto2D[featureCategory].Fill(TMath::Min(double(featureCategory == 0 ? _lPt[i] : _lPt[i] / _ptRatio[i]),varMax[0]-0.001), TMath::Abs(_lEta[i]), weight);
            if(featureCategory == 1){
                //distribs2D.vectorHisto[0].Fill(double(_lPt[i] * (1 + std::max(_relIso[i] - 0.1, 0.))), double(_lPt[i] / _ptRatio[i]));
                //distribs2D.vectorHisto[0].Fill(_gen_partonPt[i], double(_lPt[i] / _ptRatio[i]));
                distribs2D.vectorHisto[0].Fill(_gen_partonPt[i], double(_lPt[i] * (1 + std::max(_relIso[i] - 0.1, 0.))));
                //std::cout << "ptWithRelIso/ptWithClosJetPt: " << double(_lPt[i] * (1 + std::max(_relIso[i] - 0.1, 0.))) << " " << double(_lPt[i] / _ptRatio[i]) << std::endl;
            }
          }
          
      }

      std::cout << std::endl;
      cout << "Total number of events: " << distribs[0].vectorHisto[sam].Integral() << endl;
      std::cout << std::endl;
  }

  TLegend* mtleg = new TLegend(0.77,0.89,0.95,0.62);
  //mtleg->SetNColumns(1);
  mtleg->SetFillColor(0);
  mtleg->SetFillStyle(0);
  mtleg->SetBorderSize(0);
  mtleg->SetTextFont(42);


  //mtleg->AddEntry(&distribs[0].vectorHisto[dataSample],"Data","lep"); //data
  int count = 0;

  for (std::vector<std::string>::iterator it = samplesOrderNames.begin(); it != samplesOrderNames.end(); it++) {

        cout << "count and sample name: " << count << " " << *it << " " << samplesOrder.at(count) << endl;

        mtleg->AddEntry(&distribs[0].vectorHisto[samplesOrder.at(count)],(*it).c_str(),"f");
        count++;
  }

  for (int i=0; i!=nVars; ++i)  {

    //cout << "The sample size is " << samples.size() << endl;
    for (int sam=0; sam != samples.size(); ++sam){
      if(std::get<0>(samples[sam]) == "data") continue;

      // is used in MC CT
      //if(sam == 0) continue;

      //cout << "the sample is added: " << std::get<0>(samples[sam]) << endl;
      distribs[i].vectorHistoTotalUnc.Add(&distribs[i].vectorHisto[sam]);
    }

    for (int ibin = 1; ibin!=nBins[i]+1; ++ibin) {
      //get syst. uncertainty band:
      double err = 0.;
      for (int sam=0; sam != samples.size(); ++sam) {
        if(std::get<0>(samples[sam]) == "data") continue;
        // is used in MC CT
        //if(sam == 0) continue;
        if(distribs[i].vectorHisto[sam].GetBinContent(ibin) != 0){
          err += TMath::Power(distribs[i].vectorHisto[sam].GetBinError(ibin), 2); //  + TMath::Power(distribs[i].vectorHisto[sam].GetBinContent(ibin) * systematicsLeptonID[sam], 2)
        }
        else{
          err += 0;
        }
      }

      err = sqrt(err);
      distribs[i].vectorHistoTotalUnc.SetBinError(ibin, err);
    }

  }

  double scale_num = 1.6;

  TCanvas* plot[nVars];
  TCanvas* plot2D[nVars];

  for(int i = 0; i < nVars; i++) plot[i] = new TCanvas(Form("plot_%d", i),"",500,450);
  for(int i = 0; i < 2; i++) plot2D[i] = new TCanvas(Form("plot_2D_%d", i),"",500,450);

  vector<std::string> figNames = {"ptCorrEle", "ptCorrMuon"};
  vector<TString> namesForSaveFiles = {"ptCorrEle", "ptCorrMuon"};

  for(int varPlot = 0; varPlot < nVars; varPlot++){
    plot[varPlot]->cd();
    showHist(plot[varPlot],distribs[varPlot],"",figNames.at(varPlot),"Events", scale_num, mtleg, false, false); // + std::to_string(int((varMax[varPlot] - varMin[varPlot])/nBins[varPlot]))
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + ".pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + ".png");
    plot[varPlot]->cd();
    showHist(plot[varPlot],distribs[varPlot],"",figNames.at(varPlot),"Events", scale_num, mtleg, true, false); //  + std::to_string(int((varMax[varPlot] - varMin[varPlot])/nBins[varPlot]))
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + "Log.pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + "Log.png");
  }

  /*
  showHist2D(plot2D[0],distribs2D); 
  plot2D[0]->SaveAs("plotsForSave/correlation.pdf");
  plot2D[0]->SaveAs("plotsForSave/correlation.png");
  plot2D[0]->SaveAs("plotsForSave/correlation.root");
  */
  showHist2D(plot2D[0],distribs[0], 0); 
  //plot2D[0]->SaveAs("plotsForSave/eleFR.root");
  showHist2D(plot2D[1],distribs[1], 1); 
  //plot2D[1]->SaveAs("plotsForSave/muFR.root");

  return;

}

int main(int argc, char *argv[]){

    TApplication *rootapp = new TApplication("example", &argc, argv);
    treeReader reader;
    reader.Analyze();
    rootapp->Run();

    return 0;
}

