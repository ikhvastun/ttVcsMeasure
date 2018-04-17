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
  //readSamples("data/samples_FOtuning.txt"); // 
  //readSamples("data/samples_FOtuning_ttbar_2017.txt"); // 
  readSamples("data/samples_FOtuning_ttbar_2017_3L.txt"); // 
  
  std::vector<std::string> namesOfSamples = treeReader::getNamesOfTheSample();
  initdistribs(namesOfSamples);

  LastError::lasterror = Errors::UNKNOWN;

  vector<TH2D> fakeMaps;
  getFRmaps(fakeMaps);

  if(LastError::lasterror != Errors::OK){
     cout << "FR maps not found" << endl;
     return;
  }

  addVariablesToBDT();

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample();

      Color_t color = assignColor(std::get<0>(samples[sam]));
      setStackColors(color, sam);

      //if(std::get<0>(samples[sam]) == "loose") continue; // valid only for 2017 in ss2L

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
          const unsigned lCountLoose = selectLooseLep(indLoose);

          //if(lCountLoose != 2) continue;
          //if(!passPtCuts2L(indLoose)) continue;
          //if(_lCharge[indLoose.at(0)] * _lCharge[indLoose.at(1)] < 0) continue;

          if(lCountLoose != 3) continue;
          if(!passPtCuts3L(indLoose)) continue;
          
          int nLocEle = getElectronNumber(indLoose);
          //if(nLocEle != 3) continue;

          std::vector<unsigned> indJets;

          unsigned third = -9999;
          double mll = 99999;
          double pt_Z = 999999;
          double phi_Z = 999999;
          double ptNonZ = 999999;

          nJLoc = nJets(0, true, indJets, std::get<0>(samples[sam]) == "nonpromptData");
          nBLoc = nBJets(0, true, true, 1, std::get<0>(samples[sam]) == "nonpromptData");
          double dMZ = deltaMZ(indLoose, third, mll, pt_Z, ptNonZ, phi_Z);
          HTLoc = HTCalc(indJets);

          if (nJLoc < 2) continue;
          //if(_met < 30) continue;
          //if(nBLoc != 2) continue;
          //if(dMZ > 20) continue;
          if(mll < 12) continue;

          int featureCategory = 0;
          int promptCategory = 0;
          for(auto & i : indLoose){
            if (_leptonMvatZqTTV[i] > 0.4) featureCategory += 1;
            if(_lIsPrompt[i] && _lMatchPdgId[i] != 22) promptCategory += 1;
            //if(leptonIsPrompt(i)) promptCategory += 1;
            /*
            cout << "lepton info: " << _lPt[i] << " " << _lEta[i] << " " << _lPhi[i] << " " << _lFlavor[i] << " " << _lIsPrompt[i] << " " << _lProvenance[i] << " " << _lProvenanceCompressed[i] << endl;
            for(int j = 0; j < _gen_nL; j++){
                cout << "all gen leptons in the event: " << _gen_lPt[j] << " " << _gen_lEta[j] << " " << _gen_lPhi[j] << " " << _gen_lFlavor[j] << " " << _gen_lIsPrompt[j] << endl;
            }
            */
          }
        
          //cout << "prompt cat: " <<promptCategory << endl;
          //if(promptCategory == 2)
          //    cout << _runNb << " " << _lumiBlock << " " << _eventNb << endl; 
          //if(promptCategory > 1) continue;
          double FRloc = 1.;
          //if(featureCategory < 2){ 
          if(featureCategory < 3){ 

            int nFakeLepCounter = 0;

            for(auto & i : indLoose){
              if (_leptonMvatZqTTV[i] > 0.4) continue;

              // used in ttV
              const double magicNumber = 0.9;
              double leptFakePtCorr = magicNumber * _lPt[i] / _ptRatio[i];
              double FRloc_loc = fakeMaps.at(_lFlavor[i]).GetBinContent(fakeMaps.at(_lFlavor[i]).FindBin(TMath::Min(leptFakePtCorr,100-1.), fabs(_lEta[i])));

              //FRloc *= FRloc_loc/(1.-FRloc_loc);
              FRloc *= FRloc_loc;
              nFakeLepCounter++;

            }

            FRloc *= TMath::Power(-1, nFakeLepCounter + 1);
          }

          //cout << "FRloc is: " << FRloc << endl;

          weight = weight * FRloc;  

          //distribs[0].vectorHisto[featureCategory >= 2 ? 0 : 1].Fill(TMath::Min(mtHighest,varMax[0]-0.1), weight);
          distribs[1].vectorHisto[featureCategory >= leptonSelectionAnalysis ? 0 : 1].Fill(TMath::Min(_met,varMax[1]-0.1), weight);
          distribs[2].vectorHisto[featureCategory >= leptonSelectionAnalysis ? 0 : 1].Fill(TMath::Min(HTLoc,varMax[2]-0.1), weight);
          distribs[3].vectorHisto[featureCategory >= leptonSelectionAnalysis ? 0 : 1].Fill(TMath::Min(double(nJLoc),varMax[3]-0.1), weight);
          distribs[4].vectorHisto[featureCategory >= leptonSelectionAnalysis ? 0 : 1].Fill(TMath::Min(double(nBLoc),varMax[4]-0.1), weight);
          distribs[5].vectorHisto[featureCategory >= leptonSelectionAnalysis ? 0 : 1].Fill(TMath::Min(ptCorrV[0].first,varMax[5]-0.1), weight);
          distribs[6].vectorHisto[featureCategory >= leptonSelectionAnalysis ? 0 : 1].Fill(TMath::Min(ptCorrV[1].first,varMax[6]-0.1), weight);
          if(leptonSelectionAnalysis == 3)
            distribs[7].vectorHisto[featureCategory >= leptonSelectionAnalysis ? 0 : 1].Fill(TMath::Min(ptCorrV[2].first,varMax[7]-0.1), weight);
          distribs[8].vectorHisto[featureCategory >= leptonSelectionAnalysis ? 0 : 1].Fill(TMath::Min(_lEta[ptCorrV[0].second],varMax[8]-0.01), weight);
          distribs[9].vectorHisto[featureCategory >= leptonSelectionAnalysis ? 0 : 1].Fill(TMath::Min(_lEta[ptCorrV[1].second],varMax[9]-0.01), weight);
          if(leptonSelectionAnalysis == 3)
            distribs[10].vectorHisto[featureCategory >= leptonSelectionAnalysis ? 0 : 1].Fill(TMath::Min(_lEta[ptCorrV[2].second],varMax[10]-0.01), weight);
          distribs[11].vectorHisto[featureCategory >= leptonSelectionAnalysis ? 0 : 1].Fill(TMath::Min(mll,varMax[11]-0.1), weight);
          //distribs[12].vectorHisto[featureCategory >= 2 ? 0 : 1].Fill(SRID(nJLoc, nBLoc, mvaValueRegion, _charges[maxPtInd]), FRloc*scale*_weight);
          if(leptonSelectionAnalysis == 2)
            distribs[13].vectorHisto[featureCategory >= leptonSelectionAnalysis ? 0 : 1].Fill(flavourCategory2L(nLocEle, _lCharge[indLoose.at(0)]), weight);
          if(leptonSelectionAnalysis == 3)
            distribs[13].vectorHisto[featureCategory >= leptonSelectionAnalysis ? 0 : 1].Fill(flavourCategory3L(nLocEle), weight);

          
      }

      std::cout << std::endl;
      //cout << "Total number of events: " << distribs[1].vectorHisto[sam].Integral() << endl;
      //std::cout << std::endl;
  }

  cout << "Total number expected: " << distribs[1].vectorHisto[0].Integral() << endl;
  cout << "Total number predicted: " << distribs[1].vectorHisto[1].Integral() << endl;

  TLegend* mtleg = new TLegend(0.67,0.89,0.95,0.62);
  //mtleg->SetNColumns(1);
  mtleg->SetFillColor(0);
  mtleg->SetFillStyle(0);
  mtleg->SetBorderSize(0);
  mtleg->SetTextFont(42);


  //mtleg->AddEntry(&distribs[0].vectorHisto[dataSample],"Data","lep"); //data
  int count = 0;

  /*
  for (std::vector<std::string>::iterator it = samplesOrderNames.begin(); it != samplesOrderNames.end(); it++) {

        cout << "count and sample name: " << count << " " << *it << " " << samplesOrder.at(count) << endl;

        mtleg->AddEntry(&distribs[0].vectorHisto[samplesOrder.at(count)],(*it).c_str(),"f");
        count++;
  }
  */
  mtleg->AddEntry(&distribs[0].vectorHisto[0],"MC observed","lep");
  mtleg->AddEntry(&distribs[0].vectorHisto[1],"Tight-to-loose prediction","f");

  for (int i=0; i!=nVars; ++i)  {

    //cout << "The sample size is " << samples.size() << endl;
    for (int sam=0; sam != samples.size(); ++sam){
      if(std::get<0>(samples[sam]) == "data") continue;

      // is used in MC CT
      if(sam == 0) continue;

      //cout << "the sample is added: " << std::get<0>(samples[sam]) << endl;
      distribs[i].vectorHistoTotalUnc.Add(&distribs[i].vectorHisto[sam]);
    }

    for (int ibin = 1; ibin!=nBins[i]+1; ++ibin) {
      //get syst. uncertainty band:
      double err = 0.;
      for (int sam=0; sam != samples.size(); ++sam) {
        if(std::get<0>(samples[sam]) == "data") continue;
        // is used in MC CT
        if(sam == 0) continue;
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

  vector<std::string> figNames = {"Leading lepton M_{T} [GeV]", "Missing E_{T} [GeV]", "H_{T} [GeV]", "N_{jets}", "N_{b jets}", "Leading lepton p_{T}^{corr} [GeV]", "Sub-leading p_{T}^{corr} [GeV]", "Trailing p_{T}^{corr} [GeV]", "Leading lepton #eta", "Sub-leading lepton #eta", "Trailing lepton #eta", "M_{ll} [GeV]", "SR", "", "Sub-leading lepton M_{T} [GeV]", "Leading jet p_{T} [GeV]", "Sub-leading jet p_{T} [GeV]", "Min #Delta R(jet, trailing lepton)", "BDT score"};
  vector<TString> namesForSaveFiles = {"mtlead", "met", "HT", "njets", "nbjets", "ptCorrLead", "ptCorrSubLead", "ptCorrTrail", "etaCorrLead", "etaCorrSubLead", "etaCorrTrail", "mll", "SR", "flavour", "mtsublead", "jetptlead", "jetpttrail", "mindeltaR", "mvaVL"};

  for(int varPlot = 0; varPlot < nVars; varPlot++){
    //if(varPlot == 0 || varPlot == 7 || varPlot == 10 || varPlot == 12 || varPlot > 13) continue;
    if(varPlot == 0 || varPlot == 12 || varPlot > 13) continue;
    plot[varPlot]->cd();
    showHist(plot[varPlot],distribs[varPlot],"",figNames.at(varPlot),"Events", scale_num, mtleg, false, false); // + std::to_string(int((varMax[varPlot] - varMin[varPlot])/nBins[varPlot]))
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + ".pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + ".png");
    //plot[varPlot]->cd();
    //showHist(plot[varPlot],distribs[varPlot],"",figNames.at(varPlot),"Events", scale_num, mtleg, true, false); //  + std::to_string(int((varMax[varPlot] - varMin[varPlot])/nBins[varPlot]))
    //plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + "Log.pdf");
    //plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + "Log.png");
  }

  return;

}

int main(int argc, char *argv[]){

    TApplication *rootapp = new TApplication("example", &argc, argv);
    treeReader reader;
    reader.Analyze();
    rootapp->Run();

    return 0;
}

