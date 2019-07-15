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
#include "../interface/kinematicTools.h"
#include "../interface/fillDatacards.h"

#include "tdrStyle.C"

using namespace std;
using namespace tools;

Errors LastError::lasterror = Errors::UNKNOWN;
using Output::distribs1DForCT;

void treeReader::Analyze(){

  // here define parameters that are needed to be optimized
  // leptonSelection -> either dilepton or trilepton 
  // magicFactor - factor A in the definition of the corrected lepton pt, see AN-18-025 for the definition
  // leptonMVAcut is defined by the analysis itself
  leptonSelection = 3;
  magicFactor = 0.85;
  leptonMVAcut = 0.4;
  //Set CMS plotting style
  setTDRStyle();
  gROOT->SetBatch(kTRUE);
  //read samples and cross sections from txt file
  readSamples("data/samples/ClosureTestInMC/samples_FOtuning_ttbar_2017.txt"); // 
  
  std::string selection = "CTInMC";
  initListToPrint(selection);

  std::vector<std::string> namesOfProcesses = treeReader::getNamesOfTheProcesses();
  initdistribsForCT(namesOfProcesses, selection);

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample("ttZ");

      // only for 3L
      //if(std::get<1>(samples[sam]).find("DY") != std::string::npos ) continue;
      
      cout << "the sample is 2017: " << samples[sam].is2017() << endl;
      cout << "sample name is " << samples[sam].getProcessName() << endl;
      if(samples[sam].is2017() && samples[sam].getProcessName().find("TTToSemiLeptonic") != std::string::npos ) continue;
      if(samples[sam].is2017() && samples[sam].getProcessName().find("TTTo2L2Nu") != std::string::npos ) continue;

      if(samples[sam].getProcessName().find("TTJets_SingleLeptFromT") != std::string::npos ) continue;
      if(samples[sam].getProcessName().find("TTJets_DiLept") != std::string::npos ) continue;

      std::cout<<"Entries in "<< samples[sam].getProcessName() << " " << nEntries << std::endl;
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
          
          // select number of events passing FO selection
          std::vector<unsigned> indFO;
          const unsigned lCountFO = selectFakeLep(indFO, leptonSelection);

          if(leptonSelection == 3){
            if(lCountFO < 3) continue;
            if(!passPtCuts3L(indFO)) continue;
          }
          if(leptonSelection == 2){
            if(lCountFO != 2) continue;
            if(!passPtCuts2L(indFO)) continue;
            if(_lCharge[indFO.at(0)] * _lCharge[indFO.at(1)] < 0) continue;
          }
          

          // here we define number of prompt leptons in the event using promptCategory variable
          // featureCategory is responsible for event clasification into TTT, TTF, etc. categories
          int featureCategory = 0;
          int promptCategory = 0;
          for(auto & i : indFO){
            if (lepIsGood(i, leptonSelection)) featureCategory += 1;
            if(_lIsPrompt[i]) promptCategory += 1;
          }
        
          // don't consider leptons steming from the photon conversions, they are also considered as prompt leptons 
          if(leptonSelection == 3 && promptCategory > 2) continue;
          if(leptonSelection == 2 && promptCategory > 1) continue;

          // here we estimate the FR weight for events in which at least 1 lepton doesn't pass tight selection
          // also check if there is a prompt lepton in the sideband, if yes, discard this event
          //double FRloc = 1.;
          bool promptInSideband = false;
          if(featureCategory < leptonSelection){ 
            for(auto & i : indFO){
              if(_lIsPrompt[i]) promptInSideband = true;
            }
            weight *= -1 * fakeRateWeight();
          }
          if(promptInSideband) continue;


          std::vector<unsigned> indJets;
          std::vector<unsigned> indBJets;
          nJLoc = nJets(0, true, indJets, samples[sam].is2017());
          nBLoc = nBJets(0, true, true, indBJets, 1, samples[sam].is2017());
          unsigned third = -9999;
          double mll = 99999;
          double mlll = 99999;
          double ptZ = 999999;
          double ptNonZ = -999999;
          TLorentzVector Zboson, lnegative;
          std::vector<unsigned> indOf2LonZ;
          double dMZ = deltaMZ(indFO, third, mll, ptZ, ptNonZ, mlll, indOf2LonZ, Zboson, lnegative);

          if(leptonSelection == 2){
            //if(nJLoc < 2) continue;
            //if(nBLoc < 1) continue;

            TLorentzVector l0p4, l1p4;
            l0p4.SetPtEtaPhiE(ptCorrV[0].first, _lEta[indFO.at(0)], _lPhi[indFO.at(0)], _lE[indFO.at(0)] * ptCorrV[0].first / _lPt[indFO.at(0)]);
            l1p4.SetPtEtaPhiE(ptCorrV[1].first, _lEta[indFO.at(1)], _lPhi[indFO.at(1)], _lE[indFO.at(1)] * ptCorrV[1].first / _lPt[indFO.at(1)]);

            double ele_mll = (l0p4+l1p4).M();

            if(ele_mll < 12) continue;
            int nLocEle = getElectronNumber(indFO);
            if(ele_mll > 76 && ele_mll < 106 && nLocEle == 2) continue;
            if(_met < 30) continue;
          }

          double maxMJetJet = -999.;
          TLorentzVector recoilingJet(0,0,0,0);
          TLorentzVector taggedBJet(0,0,0,0);
          unsigned jetCount = 0;
          unsigned bJetCount = 0;
          double minDeltaRLeptonbJet = 0.;
          if(leptonSelection == 3){
             if(dMZ > 10) continue;
             //make lorentzvectors for leptons
             TLorentzVector lepV[lCountFO];
             for(unsigned l = 0; l < lCountFO; ++l) lepV[l].SetPtEtaPhiE(ptCorrV[l].first, _lEta[indFO.at(l)], _lPhi[indFO.at(l)], _lE[indFO.at(l)] * ptCorrV[l].first / _lPt[indFO.at(l)]);

             // met value
             TLorentzVector met;
             met.SetPtEtaPhiE(_met, 0, _metPhi, _met);

             // fill njets and nbjets variables
             std::vector<unsigned> jetInd, bJetInd;
             jetCount = nJets(0, true, indJets, samples[sam].is2017());
             bJetCount = nBJets(0, true, true, indBJets, 1, samples[sam].is2017());

             //if(jetCount < 2) continue;
             // fill vector of all jets
             TLorentzVector jetV[(const unsigned) _nJets];
             for(unsigned j = 0; j < _nJets; ++j) jetV[j].SetPtEtaPhiE(_jetPt[j], _jetEta[j], _jetPhi[j], _jetE[j]);

             std::vector<unsigned> taggedJetI; //0 -> b jet from tZq, 1 -> forward recoiling jet
             TLorentzVector neutrino = findBestNeutrinoAndTop(lepV[ptCorrV[0].second], met, taggedJetI, jetInd, bJetInd, jetV);

             if(taggedJetI[0] != 99) taggedBJet = jetV[taggedJetI[0]];
             if(taggedJetI[1] != 99) recoilingJet = jetV[taggedJetI[1]];
             maxMJetJet = kinematics::maxMass(jetV, jetInd);

             std::vector<unsigned> bjetVecInd;
             //for(unsigned l = 0; l < bJetCount; ++l) bjetVecInd.push_back(bJetInd[l]);

             std::vector<unsigned> lepVecInd;
             for(unsigned l = 0; l < lCountFO; ++l) lepVecInd.push_back(indFO[l]);

             minDeltaRLeptonbJet = kinematics::minDeltaR(lepV, lepVecInd, jetV, bjetVecInd);
             HTLoc = HTCalc(jetInd);

             //cout << "new values are " << fabs(recoilingJet.Eta()) << " " << std::max(maxMJetJet, 0.) << endl;
          }

          double mt1 = 9999;
          double mvaVL = 0;
          double mll1stpair = -999., mll2ndpair = -999.;
          double cosTSt = -999;

          map<TString, double> fillVar; 
          fillVar["ptlead"] = ptCorrV[0].first;
          fillVar["sublead"] = ptCorrV[1].first;
          fillVar["trail"] = leptonSelection > 2 ? ptCorrV[2].first : 0.;
          fillVar["njets"] = double(nJLoc);
          fillVar["nbjets"] = double(nBLoc);
          fillVar["mll"] = mll;
          fillVar["met"] = _met;
          fillVar["HT"] = HTLoc;

          //cout << "Number of tight : " << featureCategory << "; number of FO: " << lCountFO << endl; 
          //cout << "The weight is: " << weight << endl;

          for ( const auto & varPair : fillVar) {
              TString name = varPair.first;
              double var = varPair.second;
              int index = figNames[name].index;
              double varMaxOfVar = figNames[name].varMax;
              distribs1DForCT[index].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(TMath::Min(var,varMaxOfVar-0.1),weight);
          }

      }

      std::cout << std::endl;
      int index = figNames["ptlead"].index;
      cout << "Total number of events: " << distribs1DForCT[index].vectorHisto[sam].Integral() << endl;
      std::cout << std::endl;
  }

  int indexOfFirstKinVar = figNames[listToPrint[selection].at(0)].index;

  cout << "Total number expected: " << distribs1DForCT[indexOfFirstKinVar].vectorHisto[0].Integral() << endl;
  cout << "Total number predicted: " << distribs1DForCT[indexOfFirstKinVar].vectorHisto[1].Integral() << endl;

  TLegend* mtleg = new TLegend(0.17,0.89,0.95,0.72);
  mtleg->SetNColumns(2);
  mtleg->SetFillColor(0);
  mtleg->SetFillStyle(0);
  mtleg->SetBorderSize(0);
  mtleg->SetTextFont(42);

  mtleg->AddEntry(&distribs[indexOfFirstKinVar].vectorHisto[0],"Monte Carlo","lep");
  mtleg->AddEntry(&distribs[indexOfFirstKinVar].vectorHisto[1],"Tight-to-loose prediction","f");

  double scale_num = 1.4;

  TCanvas* plot[nVars];
  TCanvas* plot2D[nVars];

  for(int i = 0; i < nVars; i++) plot[i] = new TCanvas(Form("plot_%d", i),"",500,450);
  for(int i = 0; i < 2; i++) plot2D[i] = new TCanvas(Form("plot_2D_%d", i),"",500,450);

  for(int varPlot = 0; varPlot < listToPrint[selection].size(); varPlot++){

    plot[varPlot]->cd();
    showHistCT(plot[varPlot],distribs1DForCT[figNames[listToPrint[selection].at(varPlot)].index],figNames[listToPrint[selection].at(varPlot)], scale_num, mtleg, false);
    
    plot[varPlot]->SaveAs("plotsForSave/" + listToPrint[selection].at(varPlot) + ".pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + listToPrint[selection].at(varPlot) + ".png");
    plot[varPlot]->SaveAs("plotsForSave/" + listToPrint[selection].at(varPlot) + ".root");
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
    //rootapp->Run();

    return 0;
}

