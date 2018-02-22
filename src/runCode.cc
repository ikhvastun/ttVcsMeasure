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

void treeReader::Analyze(){

  leptonSelection = leptonSelectionAnalysis;
  //Set CMS plotting style
  setTDRStyle();
  gROOT->SetBatch(kTRUE);
  //read samples and cross sections from txt file
  //readSamples("samples_ttbar_emu_2017data.txt");
  //readSamples("samples_Zll_2017data.txt");
  readSamples("data/samples_LeptonMVAtraining.txt"); // 
  
  std::vector<std::string> namesOfSamples = treeReader::getNamesOfTheSample();
  initdistribs(namesOfSamples);
  /*
  cout << "string before quiting the program" << endl;
  return ;
  */

  /*
  readerBtag[0][0].load(calib_csvv2[0], BTagEntry::FLAV_B, "iterativefit");
  readerBtag[0][1].load(calib_csvv2[0], BTagEntry::FLAV_C, "iterativefit");
  readerBtag[0][2].load(calib_csvv2[0], BTagEntry::FLAV_UDSG, "iterativefit");

  readerBtag[1][0].load(calib_csvv2[1], BTagEntry::FLAV_B, "iterativefit");
  readerBtag[1][1].load(calib_csvv2[1], BTagEntry::FLAV_C, "iterativefit");
  readerBtag[1][2].load(calib_csvv2[1], BTagEntry::FLAV_UDSG, "iterativefit");
  */

  addVariablesToBDT();

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample();

      Color_t color = assignColor(std::get<0>(samples[sam]));
      setStackColors(color, sam);

      //if(leptonSelectionAnalysis == 3)
      //  if(std::get<0>(samples[sam]) == "chargeMisID") continue;

      if(std::get<0>(samples[sam]) == "data") continue;
    
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
          //if(it > 30000) break;
          
          std::vector<unsigned> ind, indLoose;
          const unsigned lCount = selectLep(ind);
          const unsigned lCountLoose = selectLooseLep(indLoose);

          if(lCountLoose < 1) continue;

          int samCategory = sam;
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
            mvaVL =  _lFlavor[i] == 0 ? readerLeptonMVAele->EvaluateMVA( "BDTG method") : readerLeptonMVAmu->EvaluateMVA( "BDTG method");

            //bool isTruePrompt = leptonIsPrompt(i);
            bool isTruePrompt = _lIsPrompt[i] && _lProvenance[i] != 1;

            if(std::get<0>(samples[sam]) == "signal" && isTruePrompt){
                
                for(int cut = 0; cut < nPoints; cut++){
                    //if(_leptonMvaTTH[i] > -1 + 2. / nPoints * cut){
                    if(mvaVL > -1 + 2. / nPoints * cut){
                        if(_lPt[i] < 25)
                            passedPrompt[_lFlavor[i]][0][cut] += 1;
                        if(_lPt[i] > 25)
                            passedPrompt[_lFlavor[i]][1][cut] += 1;
                    }
                }
                if(lepIsGood_TTV(i)){
                    if(_lPt[i] < 25)
                        passedPrompt_TTV[_lFlavor[i]][0] += 1;
                    if(_lPt[i] > 25)
                         passedPrompt_TTV[_lFlavor[i]][1] += 1;
                }
            }

            //isTruePrompt = leptonIsPrompt(i);
            isTruePrompt = _lProvenanceCompressed[i] == 0;
            if(std::get<0>(samples[sam]) == "background" && !isTruePrompt){ // && !leptonIsPrompt(i) 
            
                for(int cut = 0; cut < nPoints; cut++){
                    //if(_leptonMvaTTH[i] > -1 + 2. / nPoints * cut){
                    if(mvaVL > -1 + 2. / nPoints * cut){
                        if(_lPt[i] < 25)
                            passedNonPrompt[_lFlavor[i]][0][cut] += 1;
                        if(_lPt[i] > 25)
                            passedNonPrompt[_lFlavor[i]][1][cut] += 1;
                    }
                }
                if(lepIsGood_TTV(i)){
                    if(_lPt[i] < 25)
                        passedNonPrompt_TTV[_lFlavor[i]][0] += 1;
                    if(_lPt[i] > 25)
                        passedNonPrompt_TTV[_lFlavor[i]][1] += 1;
                }
            }
          }
          

          

      }

      std::cout << std::endl;
      cout << "Total number of events: " << distribs[0].vectorHisto[sam].Integral() << endl;
      std::cout << std::endl;
  }

  double pointsToDrawPrompt[2][2][nPoints-1];
  double pointsToDrawNonPrompt[2][2][nPoints-1];
  std::vector<int> flavorsNumber = {0, 1};
  std::vector<int> ptRanges = {0, 1};
  for(auto & flavour : flavorsNumber){
     for(auto & ptRange : ptRanges){
        for(int point = 0; point < nPoints; point++){
            if(point == 0) continue;
            pointsToDrawNonPrompt[flavour][ptRange][point-1] = passedNonPrompt[flavour][ptRange][point] / passedNonPrompt[flavour][ptRange][0];
            pointsToDrawPrompt[flavour][ptRange][point-1] = passedPrompt[flavour][ptRange][point] / passedPrompt[flavour][ptRange][0];
        }
        passedNonPrompt_TTV[flavour][ptRange] = passedNonPrompt_TTV[flavour][ptRange] / passedNonPrompt[flavour][ptRange][0];
        passedPrompt_TTV[flavour][ptRange] = passedPrompt_TTV[flavour][ptRange] /  passedPrompt[flavour][ptRange][0];
    }
  }

  TGraph* gr[2][2]; // = new TGraph(nPoints-1,pointsToDrawNonPrompt[fl][ptr],pointsToDrawPrompt[fl][ptr]);
  TCanvas *cROC = new TCanvas("cROC","cROC");
  TFile f("plotsForSave/eleCurve.root","RECREATE");
  cROC->Divide(2,2);
  for(int fl = 0; fl < 2; fl++){
    for(int ptr = 0; ptr < 2; ptr++){
        cROC->cd(2 * fl + ptr + 1);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TString name = Form("graph_%d_%d",fl,ptr);
        gr[fl][ptr] = new TGraph(nPoints-1,pointsToDrawNonPrompt[fl][ptr],pointsToDrawPrompt[fl][ptr]);
        gr[fl][ptr]->SetTitle(name);
        gr[fl][ptr]->SetName(name);
        gr[fl][ptr]->SetLineWidth(3);
        gr[fl][ptr]->Draw("AC");
        gr[fl][ptr]->Write();
    }
  }
  f.Close();
  //cROC->SaveAs("plotsForSave/eleCurve.root");


  cout << "Efficiency for POG ID: " << endl;
  cout << "     electrons: pt < 25 -> (" << passedNonPrompt_TTV[0][0] << "," << passedPrompt_TTV[0][0] << ")" << endl;
  cout << "     electrons: pt > 25 -> (" << passedNonPrompt_TTV[0][1] << "," << passedPrompt_TTV[0][1] << ")" << endl;
  cout << "     muons: pt < 25 -> (" << passedNonPrompt_TTV[1][0] << "," << passedPrompt_TTV[1][0] << ")" << endl;
  cout << "     mouns: pt > 25 -> (" << passedNonPrompt_TTV[1][1] << "," << passedPrompt_TTV[1][1] << ")" << endl;

  return;

}

int main(int argc, char *argv[]){

    TApplication *rootapp = new TApplication("example", &argc, argv);
    treeReader reader;
    reader.Analyze();
    rootapp->Run();

    return 0;
}

