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
using Output::fakeMapsCalc;
using Output::mtMaps;
using Output::mtStack;

void treeReader::Analyze(){

  double leptonMVAcut = 0.4; // 0.4 for ttZ, 0.8 for 3L tZq
  double magicFactor = 0.85; // 0.85 for ttZ, 0.95 for tZq

  // ptRatio cut that is used to select FO object in the ttZ analysis
  // be careful if you wanna use this for another measurement (for example TOP-18-008)
  // there the ptRatioCut is different
  double ptRatioCut = 0.4; 

  leptonSelection = 3;
  //Set CMS plotting style
  setTDRStyle();
  gROOT->SetBatch(kTRUE);
  //read samples and cross sections from txt file
  //readSamples("data/samples/FRmeasurement/FRInData/samples_QCD_data.txt"); // 
  readSamples("data/samples/FRmeasurement/FRInData/samples_QCD_data_2017.txt"); // 
  
  initdistribsForFRInData();

  //std::ofstream myfile;
  //myfile.open("myevents.txt");

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample("FRInData");

      Color_t color = assignColor(samples[sam].getProcessName());
      setStackColors(color, sam);

      //if(!((samples[sam].getProcessName()).find("Wjets") != std::string::npos)) continue;

      std::cout<<"Entries in "<< samples[sam].getFileName() << " " << nEntries << std::endl;
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
          //if(it > 100) break;
          //if(it > nEntries / 10) break;
          
          std::vector<unsigned> indFO;
          const unsigned lCountFO = selectFakeLep(indFO, leptonSelection);

          if(lCountFO != 1) continue;

          bool allLeptonsArePrompt = true;
          if(samples[sam].getProcessName() != "data")
            allLeptonsArePrompt = promptLeptons(indFO);

          if((samples[sam].getProcessName() == "DYjets" || samples[sam].getProcessName() == "Wjets" || samples[sam].getProcessName() == "ttbar") && !allLeptonsArePrompt) continue;

          int nLocEle = getElectronNumber(indFO);

          std::vector<unsigned> indJets;
          std::vector<unsigned> indBJets;

          unsigned third = -9999;
          double mll = 99999;
          double pt_Z = 999999;
          double phi_Z = 999999;
          double ptNonZ = 999999;

          nJLoc = nJets(0, true, indJets, samples[sam].is2017());
          nBLoc = nBJets(0, true, true, indBJets, 1, samples[sam].is2017());
          HTLoc = HTCalc(indJets);

          double leptFakePtCorr = lepIsGood(indFO.at(0), leptonSelection) ? _lPt[indFO.at(0)] : magicFactor * _lPt[indFO.at(0)] / _ptRatio[indFO.at(0)];
          TLorentzVector l0p4;
          //l0p4.SetPtEtaPhiE(_lPt[indFO.at(0)], _lEta[indFO.at(0)], _lPhi[indFO.at(0)], _lE[indFO.at(0)]);
          // let's implement the formula from ttH AN
          l0p4.SetPtEtaPhiE(35., _lEta[indFO.at(0)], _lPhi[indFO.at(0)], _lE[indFO.at(0)]);
          double mtL = mtCalc(l0p4, _met, _metPhi);

          int fakeCS = 0;
          if (_met > 20) {
            fakeCS = 1; //region for EWK subtraction
          } else if (!(_met < 20 && mtL < 20)) continue;

          int rangePeriod = (leptFakePtCorr > 15) + (leptFakePtCorr > 20) + (leptFakePtCorr > 30) + (leptFakePtCorr > 45) + (leptFakePtCorr > 65);
          int rangeEtaPeriod = _lFlavor[indFO.at(0)] ? fabs(_lEta[indFO.at(0)]) > 1.2 : fabs(_lEta[indFO.at(0)]) > 1.479;

          // for electrons there are 3 triggers: ele8, ele17 and ele23 
          // here let's consider ele17 will start from 25 GeV, ele23 from 32 
          // for muons 4 triggers: mu3, mu8, mu17 and mu27, but mu8 is not used in TOP-18-009 (should be checked why in the next iteration of the analysis)
          // here let's consider mu3 from 10, mu17 from 32 and mu27 from 45 

          // pt cut should be in principle magic factor / ptratio cut for FO not tight object
          // here it's done as it's used in TOP-18-009
          // in next iteration of the analysis this should be fixed
          bool eleTrigDecision = _HLT_Ele8_CaloIdM_TrackIdM_PFJet30 || 
                                 //(_HLT_Ele17_CaloIdM_TrackIdM_PFJet30 && leptFakePtCorr > (17 * magicFactor / ptRatioCut)) || 
                                 //(_HLT_Ele23_CaloIdM_TrackIdM_PFJet30 && leptFakePtCorr > (23 * magicFactor / ptRatioCut)); 
                                 (_HLT_Ele17_CaloIdM_TrackIdM_PFJet30 && leptFakePtCorr > 25) ||
                                 (_HLT_Ele23_CaloIdM_TrackIdM_PFJet30 && leptFakePtCorr > 32);

           
          bool muTrigDecision = (_HLT_Mu3_PFJet40 || _HLT_Mu8) || 
                                //(_HLT_Mu17 && leptFakePtCorr > (17 * magicFactor / ptRatioCut)) || 
                                //(_HLT_Mu27 && leptFakePtCorr > (27 * magicFactor / ptRatioCut));
                                (_HLT_Mu17 && leptFakePtCorr > 32) ||
                                (_HLT_Mu27 && leptFakePtCorr > 45);


          bool triggerDecision[2] = {eleTrigDecision, muTrigDecision};

          if(!triggerDecision[_lFlavor[indFO.at(0)]]) continue;

          if (fakeCS == 0) { //FR measurement region
            if(lepIsGood(indFO.at(0), leptonSelection)){
               fakeMapsCalc[_lFlavor[indFO.at(0)]][1][rangePeriod][rangeEtaPeriod][sam]->Fill(TMath::Min(leptFakePtCorr, ptBins[nPt-1]-1.), fabs(_lEta[indFO.at(0)]),weight);
            }
            fakeMapsCalc[_lFlavor[indFO.at(0)]][0][rangePeriod][rangeEtaPeriod][sam]->Fill(TMath::Min(leptFakePtCorr, ptBins[nPt-1]-1.), fabs(_lEta[indFO.at(0)]),weight);
          }
          else if (lepIsGood(indFO.at(0), leptonSelection)){
            mtMaps[_lFlavor[indFO.at(0)]][rangePeriod][rangeEtaPeriod][sam].Fill(mtL,weight);
          }

          //myfile << _runNb << " " << _lumiBlock << " " << _eventNb << endl;

      }

      std::cout << std::endl;
      std::cout << std::endl;
  }


  TCanvas* c1 = new TCanvas("fakeMaps","fakeMaps",800,400);
  c1->Divide(2,1);

  TCanvas* c2; // = new TCanvas("promptCont","promptCont",400,400);
  TPad* c2pads[2][nPt-1][nEta-2]; // for eta keep 2 regions: barrel and endcap

  std::map<std::string, std::string> listForLegend; 
  listForLegend["data"] = "data";
  listForLegend["Wjets"] = "W + jets";
  listForLegend["DYjets"] = "DY + jets";
  listForLegend["ttbar"] = "t#bar{t}";

  std::map<std::string, std::string> listFlavourForLegend; 
  listFlavourForLegend["mu"] = "#mu";
  listFlavourForLegend["el"] = "e";

  gStyle->SetOptTitle(1);
  gStyle->SetPadTopMargin(0.10);

  TLegend* leg[2];
  for(int fl=0; fl!=nFlavors; ++fl) {
     for(int rangeEta = 0; rangeEta < nEta-2; rangeEta++){ // consider barrel and endcap only 
        for(int rangePt = 0; rangePt < nPt-1; rangePt++){

            c2 = new TCanvas("promptCont","promptCont",400,400);
            c2->cd();
            //derive normalization for prompt contamination:
            double datayields = mtMaps[fl][rangePt][rangeEta][0].Integral(80/10+1,150/10+1); // 70 to 120

            double MCyields = 0.;
            for (int sam=1; sam!=nProcesses; ++sam) {
                if(sam > (is2017() ? 4 : 5)) continue;
                MCyields += mtMaps[fl][rangePt][rangeEta][sam].Integral(80/10+1,150/10+1);
            }
            double EWKsf = datayields/MCyields;
            std::cout<<"SF for "<<flavorsString[fl]<<" is "<<EWKsf<<std::endl;

            // scale data to get same number of events as in MC
            mtMaps[fl][rangePt][rangeEta][0].Scale(1./EWKsf);
            //plot histograms showing MT distributions
            c2pads[fl][rangePt][rangeEta] = new TPad(Form("pad_%d_%d_%d",fl,rangePt,rangeEta),"",0,0.,1,1);
            c2pads[fl][rangePt][rangeEta]->SetTopMargin(0.07);
            c2pads[fl][rangePt][rangeEta]->Draw();
            c2pads[fl][rangePt][rangeEta]->cd();
            mtMaps[fl][rangePt][rangeEta][0].Draw("pe");
            mtStack[fl][rangePt][rangeEta]->Draw("hist same");
            mtMaps[fl][rangePt][rangeEta][0].Draw("pe same");
            mtMaps[fl][rangePt][rangeEta][0].Draw("axis same");

            // 0 and 1 in second index correspond to denominator and numerator
            // here we subtract the prompt contamination both in denominator and numerator
            for (int sam=1; sam!=nProcesses; ++sam) {
                if(sam > (is2017() ? 4 : 5)) continue;
                fakeMapsCalc[fl][1][rangePt][rangeEta][0]->Add(fakeMapsCalc[fl][1][rangePt][rangeEta][sam],-EWKsf);
                fakeMapsCalc[fl][0][rangePt][rangeEta][0]->Add(fakeMapsCalc[fl][0][rangePt][rangeEta][sam],-EWKsf);
            }

            double lumi = 35.9;
            if(is2017()) lumi = 41.5;
            CMS_lumi( c2pads[fl][rangePt][rangeEta], 4, 0, lumi);

            leg[fl] = new TLegend(0.6,0.85,0.9,0.5,Form("SF=%3.2f",EWKsf),"NDC");
            leg[fl] = new TLegend(0.6,0.85,0.9,0.5,flavorsString[fl] + ", p_{T}:" + rangeString[rangePt] + ", #eta:" + rangeEtaString[rangeEta],"NDC");
            leg[fl]->AddEntry(&mtMaps[fl][0][0][0],"data","pel");
            for (int sam=1; sam!=nProcesses; ++sam) {
                if(sam > 3) continue;
                leg[fl]->AddEntry(&mtMaps[fl][0][0][sam],listForLegend[samplesOrderNames.at(sam)].c_str(),"f");
            }
            leg[fl]->Draw("same");
            gSystem->Exec("mkdir -p plotsForSave/maps/split");
            c2->SaveAs(Form("plotsForSave/maps/split/data_fake_EWK%d%d%d.root", fl, rangePt, rangeEta));
            delete c2;
        }
    }

    // here we sum over all eta and pt regions
    for(int rangeEta = 0; rangeEta < nEta-2; rangeEta++){
        for(int rangePt = 1; rangePt < nPt-1; rangePt++){
            fakeMapsCalc[fl][1][0][rangeEta][0]->Add(fakeMapsCalc[fl][1][rangePt][rangeEta][0]);
            fakeMapsCalc[fl][0][0][rangeEta][0]->Add(fakeMapsCalc[fl][0][rangePt][rangeEta][0]);
        }
    }
    fakeMapsCalc[fl][1][0][0][0]->Add(fakeMapsCalc[fl][1][0][1][0]);
    fakeMapsCalc[fl][0][0][0][0]->Add(fakeMapsCalc[fl][0][0][1][0]);

    // finally obtain FR maps
    TH2D *cloneFirst = (TH2D*) fakeMapsCalc[fl][1][0][0][0]->Clone("passed");
    cloneFirst->Divide(fakeMapsCalc[fl][0][0][0][0]);

    c1->cd(1+fl);

    cloneFirst->GetZaxis()->SetRangeUser(0.,1.0);
    cloneFirst->SetMarkerSize(1.5);
    cloneFirst->Draw("col");
    cloneFirst->Draw("text e same");
    cloneFirst->SaveAs("plotsForSave/maps/fakerate_"+flavorsString[fl] + "_data.root");

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

