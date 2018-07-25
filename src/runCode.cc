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
  //readSamples("data/samples_QCD_data.txt"); // 
  readSamples("data/samples_QCD_data_2017.txt"); // 
  
  std::vector<std::string> namesOfSamples = treeReader::getNamesOfTheSample();
  initdistribs(namesOfSamples);

  addVariablesToBDT();

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample();

      Color_t color = assignColor(std::get<0>(samples[sam]));
      setStackColors(color, sam);

      //if(std::get<1>(samples[sam]).find("MuEnriched") != std::string::npos ) continue;

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
          //if(it > nEntries / 50) break;
          
          std::vector<unsigned> indFO;
          const unsigned lCountFO = selectFakeLep(indFO);

          if(lCountFO != 1) continue;

          bool allLeptonsArePrompt = true;
          if(std::get<0>(samples[sam]) != "data")
            allLeptonsArePrompt = promptLeptons(indFO);

          if((std::get<0>(samples[sam]) == "DYjets" || std::get<0>(samples[sam]) == "Wjets" || std::get<0>(samples[sam]) == "ttbar") && !allLeptonsArePrompt) continue;

          int nLocEle = getElectronNumber(indFO);

          std::vector<unsigned> indJets;

          unsigned third = -9999;
          double mll = 99999;
          double pt_Z = 999999;
          double phi_Z = 999999;
          double ptNonZ = 999999;

          nJLoc = nJets(0, true, indJets, std::get<0>(samples[sam]) == "nonpromptData");
          nBLoc = nBJets(0, false, true, 1, std::get<0>(samples[sam]) == "nonpromptData");
          HTLoc = HTCalc(indJets);

          double leptFakePtCorr = lepIsGood(indFO.at(0)) ? _lPt[indFO.at(0)] : magicFactorAnalysis * _lPt[indFO.at(0)] / _ptRatio[indFO.at(0)];
          TLorentzVector l0p4;
          //l0p4.SetPtEtaPhiE(_lPt[indFO.at(0)], _lEta[indFO.at(0)], _lPhi[indFO.at(0)], _lE[indFO.at(0)]);
          // let's implement the formula from ttH AN
          l0p4.SetPtEtaPhiE(35., _lEta[indFO.at(0)], _lPhi[indFO.at(0)], _lE[indFO.at(0)]);
          double mtL = mtCalc(l0p4, _met, _metPhi);

          int fakeCS = 0;
          if (_met > 20) {
            fakeCS = 1; //region for EWK subtraction
          } else if (!(_met < 20 && mtL < 20)) continue;

          //double ptCorrCutEle = magicFactorAnalysis * 17 / 0.3; // 17 is pt threshold on high pt prescalled trigger, 0.3 stands for pt ratio cut in FO object
          //double ptCorrCutMu  = magicFactorAnalysis * 17 / 0.3;

          // for electrons there are 3 triggers: ele8, ele17 and ele23 
          // here let's consider ele17 will start from 30 GeV, mu27 from 50
          // for muons 4 triggers: mu3, mu8, mu17 and mu27
          // here let's consider mu3 from 10, mu8 from 15, mu17 from 30 and mu27 from 50
          /*
          int rangePeriod = 0;
          if(_lFlavor[indFO.at(0)] == 0)
            rangePeriod = (leptFakePtCorr > magicFactorAnalysis * 17 / (is2017 ? 0.3 : 0.4)) + (leptFakePtCorr > magicFactorAnalysis * 23 / (is2017 ? 0.3 : 0.4));
          if(_lFlavor[indFO.at(0)] == 1)
            rangePeriod = (leptFakePtCorr > magicFactorAnalysis * 17 / (is2017 ? 0.3 : 0.4)) + (leptFakePtCorr > magicFactorAnalysis * 27 / (is2017 ? 0.3 : 0.4));

          bool eleTrigDecision = (_HLT_Ele8_CaloIdM_TrackIdM_PFJet30 && leptFakePtCorr < magicFactorAnalysis * 17 / (is2017 ? 0.3 : 0.4)) || (_HLT_Ele17_CaloIdM_TrackIdM_PFJet30 && leptFakePtCorr > magicFactorAnalysis * 17 / (is2017 ? 0.3 : 0.4) && leptFakePtCorr < magicFactorAnalysis * 23 / (is2017 ? 0.3 : 0.4)) || (_HLT_Ele23_CaloIdM_TrackIdM_PFJet30 && leptFakePtCorr > magicFactorAnalysis * 23 / (is2017 ? 0.3 : 0.4)); 
          bool muTrigDecision = ((_HLT_Mu3_PFJet40 || _HLT_Mu8) && leptFakePtCorr < magicFactorAnalysis * 17 / (is2017 ? 0.3 : 0.4)) || (_HLT_Mu17 && leptFakePtCorr > magicFactorAnalysis * 17 / (is2017 ? 0.3 : 0.4) && leptFakePtCorr < magicFactorAnalysis * 27 / (is2017 ? 0.3 : 0.4)) || (_HLT_Mu27 && leptFakePtCorr > magicFactorAnalysis * 27 / (is2017 ? 0.3 : 0.4));

          */
          /*
          int rangePeriod = (leptFakePtCorr > 30) + (leptFakePtCorr > 50);
          bool eleTrigDecision = _HLT_Ele8_CaloIdM_TrackIdM_PFJet30 || 
                                 (_HLT_Ele17_CaloIdM_TrackIdM_PFJet30 && leptFakePtCorr > 30) || 
                                 (_HLT_Ele23_CaloIdM_TrackIdM_PFJet30 && leptFakePtCorr > 50); 
           
          bool muTrigDecision = (_HLT_Mu3_PFJet40 || _HLT_Mu8) || 
                                (_HLT_Mu17 && leptFakePtCorr > 30) || 
                                (_HLT_Mu27 && leptFakePtCorr > 50);
          */

          int rangePeriod = (leptFakePtCorr > 15) + (leptFakePtCorr > 20) + (leptFakePtCorr > 30) + (leptFakePtCorr > 45) + (leptFakePtCorr > 65);
          int rangeEtaPeriod = _lFlavor[indFO.at(0)] ? fabs(_lEta[indFO.at(0)]) > 1.2 : fabs(_lEta[indFO.at(0)]) > 1.479;
          bool eleTrigDecision = _HLT_Ele8_CaloIdM_TrackIdM_PFJet30 || 
                                 (_HLT_Ele17_CaloIdM_TrackIdM_PFJet30 && leptFakePtCorr > 25) || 
                                 (_HLT_Ele23_CaloIdM_TrackIdM_PFJet30 && leptFakePtCorr > 32); 
           
          bool muTrigDecision = (_HLT_Mu3_PFJet40 || _HLT_Mu8) || 
                                (_HLT_Mu17 && leptFakePtCorr > 32) || 
                                (_HLT_Mu27 && leptFakePtCorr > 45);

          bool triggerDecision[2] = {eleTrigDecision, muTrigDecision};

          //if(sam == 0)
          if(!triggerDecision[_lFlavor[indFO.at(0)]]) continue;

          if (fakeCS == 0) { //FR measurement region
            if(lepIsGood(indFO.at(0))){
               fakeMapsCalc[_lFlavor[indFO.at(0)]][1][rangePeriod][rangeEtaPeriod][sam]->Fill(TMath::Min(leptFakePtCorr, ptBins[nPt-1]-1.), fabs(_lEta[indFO.at(0)]),weight);
            }
            fakeMapsCalc[_lFlavor[indFO.at(0)]][0][rangePeriod][rangeEtaPeriod][sam]->Fill(TMath::Min(leptFakePtCorr, ptBins[nPt-1]-1.), fabs(_lEta[indFO.at(0)]),weight);
          }
          else if (lepIsGood(indFO.at(0))){
            mtMaps[_lFlavor[indFO.at(0)]][rangePeriod][rangeEtaPeriod][sam]->Fill(mtL,weight);
          }

      }

      std::cout << std::endl;
      //cout << "Total number of events: " << distribs[0].vectorHisto[sam].Integral() << endl;
      std::cout << std::endl;
  }

  TCanvas* c1 = new TCanvas("fakeMaps","fakeMaps",800,400);
  c1->Divide(2,1);

  TCanvas* c2 = new TCanvas("promptCont","promptCont",800,400);
  c2->Divide(rangePeriods,2*2);
  TPad* c2pads[2][rangePeriods][rangeEtaPeriods][2];

  //gStyle->SetPaintTextFormat("1.2f");
  gStyle->SetOptTitle(1);
  gStyle->SetPadTopMargin(0.10);

  //TGaxis::SetMaxDigits(3);

  TLegend* leg[2];
  for(int i=0; i!=nFlavors; ++i) {
     for(int rangeEta = 0; rangeEta < rangeEtaPeriods; rangeEta++){
        for(int range = 0; range < rangePeriods; range++){

            //derive normalization for prompt contamination:
            double datayields = mtMaps[i][range][rangeEta][0]->Integral(80/10+1,150/10+1); // 70 to 120

            double MCyields = 0;
            for (int sam=1; sam!=nSamples; ++sam) {
                if(sam > (is2017 ? 4 : 5)) continue;
                MCyields += mtMaps[i][range][rangeEta][sam]->Integral(80/10+1,150/10+1);
                mtMaps[i][range][rangeEta][nSamples]->Add(mtMaps[i][range][rangeEta][sam]);
            }
            double EWKsf = datayields/MCyields;
            std::cout<<"SF for "<<flavorsString[i]<<" is "<<EWKsf<<std::endl;

            mtMaps[i][range][rangeEta][0]->Scale(1./EWKsf);
            //plot histograms showing MT distributions
            c2->cd(1+rangePeriods*rangeEtaPeriods*i+rangeEta*rangePeriods+range);
            c2pads[i][range][rangeEta][0] = new TPad(Form("pad_%d_%d_%d_0",i,range,rangeEta),"",0,0.3,1,1);
            c2pads[i][range][rangeEta][0]->Draw();
            c2pads[i][range][rangeEta][0]->cd();
            mtMaps[i][range][rangeEta][0]->Draw("pe");
            mtStack[i][range][rangeEta]->Draw("hist same");
            mtMaps[i][range][rangeEta][0]->Draw("pe same");
            mtMaps[i][range][rangeEta][0]->Draw("axis same");

            //leg[i] = new TLegend(0.6,0.85,0.9,0.5,Form("SF=%3.2f",EWKsf),"NDC");
            leg[i] = new TLegend(0.6,0.85,0.9,0.5,"pt:" + rangeString[range] + ", eta:" + rangeEtaString[rangeEta],"NDC");
            leg[i]->AddEntry(mtMaps[i][0][0][0],"data","pel");
            for (int sam=1; sam!=nSamples; ++sam) {
                if(sam > 3) continue;
                leg[i]->AddEntry(mtMaps[i][0][0][sam],samplesOrderNames.at(sam).c_str(),"f");
            }
            leg[i]->Draw("same");
            c2->cd(1+rangePeriods*rangeEtaPeriods*i+rangeEta*rangePeriods+range);
            c2pads[i][range][rangeEta][1] = new TPad(Form("pad_%d_%d_%d_1",i,range,rangeEta),"",0,0.0,1,0.3);
            c2pads[i][range][rangeEta][1]->Draw();
            c2pads[i][range][rangeEta][1]->cd();
            mtMaps[i][range][rangeEta][nSamples]->Divide(mtMaps[i][range][rangeEta][0],mtMaps[i][range][rangeEta][nSamples]);
            mtMaps[i][range][rangeEta][nSamples]->GetYaxis()->SetTitle("data/MC");
            mtMaps[i][range][rangeEta][nSamples]->GetYaxis()->SetRangeUser(0,6);
            mtMaps[i][range][rangeEta][nSamples]->GetYaxis()->SetNdivisions(505);
            mtMaps[i][range][rangeEta][nSamples]->GetXaxis()->SetTitle("");
            mtMaps[i][range][rangeEta][nSamples]->GetYaxis()->SetTitleOffset(0.3/0.7*1.25);
            mtMaps[i][range][rangeEta][nSamples]->GetYaxis()->SetTitleSize(0.7/0.3*0.06);
            mtMaps[i][range][rangeEta][nSamples]->GetXaxis()->SetTitleSize(0.7/0.3*0.06);
            mtMaps[i][range][rangeEta][nSamples]->GetYaxis()->SetLabelSize(0.7/0.3*0.05);
            mtMaps[i][range][rangeEta][nSamples]->GetXaxis()->SetLabelSize(0.7/0.3*0.05);
            mtMaps[i][range][rangeEta][nSamples]->Draw("pe");
            mtMaps[i][range][rangeEta][nSamples]->Fit("pol0","","",80,150);
  
            //obtain and plot FR
            for (int sam=1; sam!=nSamples; ++sam) {
                if(sam > (is2017 ? 4 : 5)) continue;
                fakeMapsCalc[i][1][range][rangeEta][0]->Add(fakeMapsCalc[i][1][range][rangeEta][sam],-EWKsf);
                fakeMapsCalc[i][0][range][rangeEta][0]->Add(fakeMapsCalc[i][0][range][rangeEta][sam],-EWKsf);
            }
        }
    }

    for(int rangeEta = 0; rangeEta < rangeEtaPeriods; rangeEta++){
        for(int range = 1; range < rangePeriods; range++){
            fakeMapsCalc[i][1][0][rangeEta][0]->Add(fakeMapsCalc[i][1][range][rangeEta][0]);
            fakeMapsCalc[i][0][0][rangeEta][0]->Add(fakeMapsCalc[i][0][range][rangeEta][0]);
        }
    }
    fakeMapsCalc[i][1][0][0][0]->Add(fakeMapsCalc[i][1][0][1][0]);
    fakeMapsCalc[i][0][0][0][0]->Add(fakeMapsCalc[i][0][0][1][0]);

    //fakeMapsCalc[i][2][0][0]->Divide(fakeMapsCalc[i][1][0][0],fakeMapsCalc[i][0][0][0]);
    TH2D *cloneFirst = (TH2D*) fakeMapsCalc[i][1][0][0][0]->Clone("passed");
    cloneFirst->Divide(fakeMapsCalc[i][0][0][0][0]);

    c1->cd(1+i);

    //fakeMapsCalc[i][2][0][0]->GetZaxis()->SetRangeUser(0.,1.0);
    //fakeMapsCalc[i][2][0][0]->SetMarkerSize(1.5);

    //fakeMapsCalc[i][2][0][0]->Draw("col");
    //fakeMapsCalc[i][2][0][0]->Draw("text e same");

    //fakeMapsCalc[i][2][0][0]->SaveAs("maps/fakerate_"+flavorsString[i] + "_data.root");

    cloneFirst->GetZaxis()->SetRangeUser(0.,1.0);
    cloneFirst->SetMarkerSize(1.5);
    cloneFirst->Draw("col");
    cloneFirst->Draw("text e same");
    cloneFirst->SaveAs("maps/fakerate_"+flavorsString[i] + "_data.root");

 }

 //c1->SaveAs("maps/data_fake_maps.pdf");
 c1->SaveAs("maps/data_fake_maps.root");
 //c2->SaveAs("maps/data_fake_EWK.pdf");
 c2->SaveAs("maps/data_fake_EWK.root");
 return;

}

int main(int argc, char *argv[]){

    TApplication *rootapp = new TApplication("example", &argc, argv);
    treeReader reader;
    reader.Analyze();
    rootapp->Run();

    return 0;
}

