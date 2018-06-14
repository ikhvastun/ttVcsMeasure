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
  magicFactor = magicFactorAnalysis;
  leptonMVAcut = leptonMVAcutAnalysis;
  //Set CMS plotting style
  setTDRStyle();
  gROOT->SetBatch(kTRUE);
  //read samples and cross sections from txt file
  readSamples("data/samples_chargeMisId.txt"); 
  //readSamples("data/samples_chargeMisId_2017.txt"); 
  
  std::vector<std::string> namesOfSamples = treeReader::getNamesOfTheSample();
  initdistribs(namesOfSamples);

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample();

      //Color_t color = assignColor(std::get<0>(samples[sam]));
      //setStackColors(color, sam);

      if(std::get<0>(samples[sam]) == "loose") continue;

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
          //if(it > 10000) break;
          //if(it > nEntries / 20) break;
          
          std::vector<unsigned> ind;
          const unsigned lCount = selectLep(ind);

          int nLocEle =  getElectronNumber(ind);

          //cout << "number of leptons: " << lCount << " " << nLocEle << endl;

          if(lCount != 2) continue;
          //if(nLocEle != 2) continue;

          TLorentzVector l1p4, l2p4;
          l1p4.SetPtEtaPhiE(_lPt[ind.at(0)], _lEta[ind.at(0)], _lPhi[ind.at(0)], _lE[ind.at(0)]);
          l2p4.SetPtEtaPhiE(_lPt[ind.at(1)], _lEta[ind.at(1)], _lPhi[ind.at(1)], _lE[ind.at(1)]);

          //if((l1p4+l2p4).M() < 71 || (l1p4+l2p4).M() > 111) continue;

          //cout << "found two leptons with mass close to Z" << endl;

          /*
          std::vector<unsigned> indJets;

          unsigned third = -9999;
          double mll = 99999;
          double pt_Z = 999999;
          double phi_Z = 999999;
          double ptNonZ = 999999;

          nJLoc = nJets(0, true, indJets, std::get<0>(samples[sam]) == "nonpromptData");
          nBLoc = nBJets(0, false, true, 1, std::get<0>(samples[sam]) == "nonpromptData");
          HTLoc = HTCalc(indJets);
          */

          for(auto & lep : ind){
            if(_lFlavor[lep] != 0) continue;
            if(!_lIsPrompt[lep]) continue;
            if(leptonGenCharge(lep) == 99) continue;
            //cout << "lepton with pt: " << _lPt[lep] << "; charge reco: " << _lCharge[lep] << "; charge gen: " << leptonGenCharge(lep) << endl;
            chargeMisIDMapsCalc[_lCharge[lep] != leptonGenCharge(lep)]->Fill(TMath::Min(_lPt[lep], ptBins[nPt-1]-1.), fabs(_lEta[lep]),weight);
          }

      }

      std::cout << std::endl;
      //cout << "Total number of events: " << distribs[0].vectorHisto[sam].Integral() << endl;
      std::cout << std::endl;
  }

  TCanvas* c1 = new TCanvas("chargeMisIDMaps","chargeMisIDMaps",400,400);

  chargeMisIDMapsCalc[1]->Divide(chargeMisIDMapsCalc[0]);
  chargeMisIDMapsCalc[1]->SetName("passed");
  chargeMisIDMapsCalc[1]->GetZaxis()->SetRangeUser(1e-5,1e-1);
  chargeMisIDMapsCalc[1]->SetMarkerSize(1.5);
  chargeMisIDMapsCalc[1]->Draw("col");
  chargeMisIDMapsCalc[1]->Draw("text e same");
  chargeMisIDMapsCalc[1]->SaveAs("maps/chargeMisID_MC.root");

  return;

}

int main(int argc, char *argv[]){

    TApplication *rootapp = new TApplication("example", &argc, argv);
    treeReader reader;
    reader.Analyze();
    rootapp->Run();

    return 0;
}

