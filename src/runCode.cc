#include <algorithm>
#include <vector>
#include <map>
#include <iomanip>

#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
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
  //gROOT->SetBatch(kTRUE);
  //read samples and cross sections from txt file
  //readSamples("samples_ttWtraining.txt");
  //readSamples("samples_ttW.txt");
  //readSamples("nonpromptMCForCT.txt");
  readSamples("samples_ttbar_emu.txt");
  setTDRStyle(); 
  

  std::vector<std::string> namesOfSamples = treeReader::getNamesOfTheSample();
  initdistribs(namesOfSamples);

  readerBtag[0][0].load(calib_csvv2[0], BTagEntry::FLAV_B, "iterativefit");
  readerBtag[0][1].load(calib_csvv2[0], BTagEntry::FLAV_C, "iterativefit");
  readerBtag[0][2].load(calib_csvv2[0], BTagEntry::FLAV_UDSG, "iterativefit");

  readerBtag[1][0].load(calib_csvv2[1], BTagEntry::FLAV_B, "iterativefit");
  readerBtag[1][1].load(calib_csvv2[1], BTagEntry::FLAV_C, "iterativefit");
  readerBtag[1][2].load(calib_csvv2[1], BTagEntry::FLAV_UDSG, "iterativefit");

  
  LastError::lasterror = Errors::UNKNOWN;

  vector<TH2D> fakeMaps;
  getFRmaps(fakeMaps);

  if(LastError::lasterror != Errors::OK){
     cout << "FR maps not found" << endl;
     return;
  }
  

  if(leptonSelectionAnalysis == 2){

    addVariablesToBDT();
    //addBranchToBDTTreeVariables();
  }



  std::ofstream myfile;
  myfile.open("myevents.txt");

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample();

      Color_t color = assignColor(std::get<0>(samples[sam]));
      setStackColors(color, sam);

      //if(leptonSelectionAnalysis == 3)
      //  if(std::get<0>(samples[sam]) == "chargeMisID") continue;

      //if(std::get<0>(samples[sam]) != "nonpromptData") continue;
    
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
          //if(it > 5000) break;

          bool printAddInfo = false;
          //if(_eventNb != 2901696869) continue;

          if(!(_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || _HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ)) continue;

          std::vector<unsigned> indTight, indFake;
          //select leptons
          const unsigned lCount = selectLep(indTight);
          const unsigned lCountFake = selectFakeLep(indFake);

          if(printAddInfo)
            cout << "number of leptons: " << lCount << " " << lCountFake << endl;

          std::vector<unsigned> & ind = indTight;
          if(lCount != 2) continue;

          int samCategory = sam;
          //if(_lFlavor[ind.at(0)] + _lFlavor[ind.at(1)] != 1) continue;

          //cout << "new event#################" << endl;

          bool allLeptonsArePrompt = true;
          
          if(std::get<0>(samples[sam]) != "data" && std::get<0>(samples[sam]) != "nonpromptData")
            allLeptonsArePrompt = promptLeptons(ind);
          
          if(std::get<0>(samples[sam]) == "chargeMisID" && !allLeptonsArePrompt) continue;
          //if(std::get<0>(samples[sam]) == "nonprompt" && allLeptonsArePrompt) continue; // works just for MC

          if((std::get<0>(samples[sam]) == "ttZ" || std::get<0>(samples[sam]) == "ttX" || std::get<0>(samples[sam]) == "WZ" || std::get<0>(samples[sam]) == "rare") && !allLeptonsArePrompt) continue;

          int nLocEle = getElectronNumber(ind);
          if(nLocEle != 1) continue;

          if(!passPtCuts2L(ind)) continue;

          std::vector<unsigned> indJets;

          unsigned third = -9999;
          double mll = 99999;
          double ptZ = 999999;
          double ptNonZ = 999999;

          nJLoc = nJets(0, true, indJets);
          nBLoc = nBJets(0, true, true, 1);
          double dMZ = deltaMZ(ind, third, mll, ptZ, ptNonZ);

          HTLoc = HTCalc(indJets);
          
          if(printAddInfo){
            cout << "number of jets/bjets/deltaMZ: " << nJLoc << " " << nBLoc << " " << dMZ << endl;
            cout << endl;
          }

          if(nJLoc < 2) continue;
          //if(nBLoc < 1) continue;
          
          double dataMCSF = 1.;
          double lepSF = 1.;
          double lepSFUp = 1.;
          double lepSFDown = 1.;

          if(printAddInfo){
            cout << "total weight is: " << weight << endl;
          }
          

          double mvaVL = 0;

          distribs[0].vectorHisto[samCategory].Fill(TMath::Min(ptCorrV[0].first,varMax[0]-0.1),weight);
          distribs[1].vectorHisto[samCategory].Fill(TMath::Min(ptCorrV[1].first,varMax[1]-0.1),weight);
          if(leptonSelectionAnalysis > 2)
            distribs[2].vectorHisto[samCategory].Fill(TMath::Min(ptCorrV[2].first,varMax[2]-0.1),weight);

          distribs[8].vectorHisto[samCategory].Fill(TMath::Min(double(nJLoc),varMax[8]-0.1),weight);
          distribs[9].vectorHisto[samCategory].Fill(TMath::Min(double(nBLoc),varMax[9]-0.1),weight);

          distribs[16].vectorHisto[samCategory].Fill(TMath::Min(mll,varMax[16]-0.1),weight);
          distribs[17].vectorHisto[samCategory].Fill(TMath::Min(ptZ,varMax[17]-0.1),weight);
          distribs[18].vectorHisto[samCategory].Fill(TMath::Min(ptNonZ,varMax[18]-0.1),weight);

          //distribs[19].vectorHisto[samCategory].Fill(TMath::Min(_jetCsvV2[indJets.at(0)],varMax[19]-0.001),weight);
          //distribs[20].vectorHisto[samCategory].Fill(TMath::Min(_jetCsvV2[indJets.at(1)],varMax[20]-0.001),weight);

          distribs[19].vectorHisto[samCategory].Fill(TMath::Min(_jetDeepCsv_b[indJets.at(0)] + _jetDeepCsv_bb[indJets.at(0)],varMax[19]-0.001),weight);
          distribs[20].vectorHisto[samCategory].Fill(TMath::Min(_jetDeepCsv_b[indJets.at(1)] + _jetDeepCsv_bb[indJets.at(1)],varMax[20]-0.001),weight);
          if(nJLoc > 2)
            distribs[21].vectorHisto[samCategory].Fill(TMath::Min(_jetDeepCsv_b[indJets.at(2)] + _jetDeepCsv_bb[indJets.at(2)],varMax[21]-0.001),weight);


          if(SRID3L(nJLoc, nBLoc) == 8)
            myfile << _runNb << " " << _lumiBlock << " " << _eventNb << endl;
          
          //cout << "Lepton events: " << _lPt[ind.at(0)] << " " << _lPt[ind.at(1)] << " " << _lPt[ind.at(2)] << endl;
      }

      std::cout << std::endl;
      cout << "Total number of events: " << distribs[0].vectorHisto[sam].Integral() << endl;
      std::cout << std::endl;
  }

  fileDummy->cd();
  signalTree->Write();
  bkgTree->Write();

  fileDummy->Close();
  //return;


  TLegend* mtleg = new TLegend(0.25,0.89,0.95,0.77); 
  mtleg->SetNColumns(4);
  mtleg->SetFillColor(0);
  mtleg->SetFillStyle(0);
  mtleg->SetBorderSize(0);
  mtleg->SetTextFont(42);
  
  
  mtleg->AddEntry(&distribs[0].vectorHisto[dataSample],"Data","lep"); //data
  int count = 0;
  for (std::vector<std::string>::iterator it = samplesOrderNames.begin()+1; it != samplesOrderNames.end(); it++) {
        count++;
        if(samplesOrderNames.at(count) == "ttH") continue;
        if(samplesOrderNames.at(count) == "ZZ") continue;

        cout << "count and sample name: " << count << " " << *it << " " << samplesOrder.at(count) << endl;

        mtleg->AddEntry(&distribs[0].vectorHisto[samplesOrder.at(count)],(*it).c_str(),"f");
  }
  
  /*
  Color_t color = kBlack;
  setStackColors(color, 0);
  color = kBlue-9;
  setStackColors(color, 1);
  mtleg->AddEntry(&distribs[0].vectorHisto[0],"MC observed","lep"); 
  mtleg->AddEntry(&distribs[0].vectorHisto[1],"Tight-to-loose prediction","f"); //data
  */

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
  
  TCanvas* plot[14];
      
  for(int i = 0; i < 14; i++){
      plot[i] = new TCanvas(Form("plot_%d", i),"",500,450);
  }

  
  plot[0]->cd();
  showHist(plot[0],distribs[0],"","Leading lepton p_{T} [GeV]","Events / " + std::to_string(int((varMax[0] - varMin[0])/nBins[0])) + " GeV",scale_num, mtleg);

  plot[1]->cd();
  showHist(plot[1],distribs[1],"","Sub-leading lepton p_{T} [GeV]","Events / " + std::to_string(int((varMax[1] - varMin[1])/nBins[1])) + " GeV",scale_num, mtleg);

  plot[2]->cd();
  showHist(plot[2],distribs[2],"","Trailing lepton p_{T} [GeV]","Events / " + std::to_string(int((varMax[2] - varMin[2])/nBins[2])) + " GeV",scale_num, mtleg);

  plot[3]->cd();
  showHist(plot[3],distribs[8],"","N_{j}","Events",scale_num, mtleg);

  plot[4]->cd();
  showHist(plot[4],distribs[9],"","N_{b}","Events",scale_num, mtleg);

  plot[5]->cd();
  showHist(plot[5],distribs[15],"","SR","Events",scale_num, mtleg);

  plot[6]->cd();
  showHist(plot[6],distribs[14],"","flavour","Events",scale_num, mtleg);

  plot[7]->cd();
  showHist(plot[7],distribs[12],"","BDT","Events",scale_num, mtleg);

  plot[8]->cd();
  showHist(plot[8],distribs[16],"","M(ll) [GeV]","Events",scale_num, mtleg);

  plot[9]->cd();
  showHist(plot[9],distribs[17],"","p_{T}^{Z} [GeV]","Events",scale_num, mtleg);

  plot[10]->cd();
  showHist(plot[10],distribs[18],"","Non-Z lepton p_{T} [GeV]","Events",scale_num, mtleg);

  plot[11]->cd();
  showHist(plot[11],distribs[19],"","deepcsv of leading jet","Events",scale_num, mtleg);

  plot[12]->cd();
  showHist(plot[12],distribs[20],"","deepcsv of subleading jet","Events",scale_num, mtleg);

  plot[13]->cd();
  showHist(plot[13],distribs[21],"","deepcsv of 3rd jet","Events",scale_num, mtleg);

  vector<TString> namesForSaveFiles = {"ptlead", "sublead", "trail", "njets", "nbjets", "SR", "flavour", "BDT", "mll", "ptZ", "ptNonZ", "csv1st", "csv2nd", "csv3rd"};
  int countPlot = 0;
  for(auto & i : namesForSaveFiles){
    plot[countPlot]->SaveAs("plotsForSave/" + i + ".pdf");
    countPlot++;
  }
  //plot[8]->cd();
  //drawSystUnc(plot[8], distribs[15], 8);
  

  //fillDatacards(distribs[15], samplesOrderNames, samplesOrder);
}

int main(int argc, char *argv[]){

    TApplication *rootapp = new TApplication("example", &argc, argv);
    treeReader reader;
    reader.Analyze();
    rootapp->Run();

    return 0;
}

