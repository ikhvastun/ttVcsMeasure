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
  //readSamples("data/samples_ttbar_emu_2016data.txt");
  readSamples("data/samples_Zll_2016data.txt");
  //readSamples("data/samples_LeptonMVAtraining_2017MC.txt");
  
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

      //if(std::get<0>(samples[sam]) == "data") continue;
      //if(std::get<0>(samples[sam]) != "background") continue;
    
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
          //cout << "NEW EVENT" << endl;
          
          /*
          const bool is_in = eventsList.find(_eventNb) != eventsList.end();
          if(!is_in) continue;
          */

          //bool printAddInfo = false;
          
          //if(_eventNb != 1492195608) continue;
          /*
          if(printAddInfo){
            cout << "################################################################" << endl;
            cout << "Event number: " << _runNb << " " << _lumiBlock << " " << _eventNb << endl;
          }
          */

          //
          
          //if(isData){
              // should be used for Muon Prompt Reco
              //_passMETFilters = _Flag_HBHENoiseFilter && _Flag_HBHENoiseIsoFilter && _Flag_EcalDeadCellTriggerPrimitiveFilter && _Flag_goodVertices && _Flag_eeBadScFilter && _Flag_global;
              
              //if(!_passMETFilters) continue;
              //if(!_2017_mm) continue;
          //}

          /*
          if(printAddInfo){
            cout << "Filter info: " << _passMETFilters << " " << _Flag_HBHENoiseFilter << " " << _Flag_HBHENoiseIsoFilter << " " << _Flag_EcalDeadCellTriggerPrimitiveFilter << " " << _Flag_goodVertices << " " << _Flag_eeBadScFilter << " " << _Flag_globalTightHalo2016Filter << " " << _Flag_BadPFMuonFilter << " " << _Flag_BadChargedCandidateFilter << endl;
          }
          */

          //myfile << _runNb << " " << _lumiBlock << " " << _eventNb << " " << _met << " " << _passMETFilters << endl;

          std::vector<unsigned> ind, indLoose;
          const unsigned lCount = selectLep(ind);
          const unsigned lCountLoose = selectLooseLep(indLoose);

          /*
          if(printAddInfo)
            cout << "number of leptons: " << lCount << endl;
            */

          //if(lCountLoose != 2) continue;
          //if(_lCharge[indLoose.at(0)] * _lCharge[indLoose.at(1)] < 0) continue;
          if(lCountLoose != 2) continue;

          int samCategory = sam;
          int nLocEle = getElectronNumber(indLoose);
          int nLocMu = lCountLoose - nLocEle;
          if(nLocEle != 2) continue;
          //if(nLocMu != 2) continue;

          /*
          if(nLocEle < 1) continue;
          int nLocMu = lCountLoose - nLocEle;
          if(nLocMu < 1) continue;
          */

          if(!passPtCuts2L(indLoose)) continue;

          std::vector<unsigned> indJets;

          unsigned third = -9999;
          double mll = 99999;
          double pt_Z = 999999;
          double phi_Z = 999999;
          double ptNonZ = 999999;

          nJLoc = nJets(0, true, indJets, std::get<0>(samples[sam]) == "nonpromptData");
          nBLoc = nBJets(0, false, true, 1, std::get<0>(samples[sam]) == "nonpromptData");
          double dMZ = deltaMZ(indLoose, third, mll, pt_Z, ptNonZ, phi_Z);

          //HTLoc = HTCalc(indJets);
          
          /*
          if(printAddInfo){
            cout << "number of jets/bjets/deltaMZ: " << nJLoc << " " << nBLoc << " " << dMZ << endl;
          }
          */

          if(dMZ > 10) continue;
          
          /*
          if(nJLoc < 2) continue;
          if(nBLoc < 1) continue;
          if(_met < 50) continue;
          */

          /*
          if(nJLoc < 3) continue;
          if(nJLoc > 4) continue;
          if(nBLoc != 1) continue;
          if(_met < 50) continue;
          */

          double dataMCSF = 1.;
          double lepSF = 1.;
          double lepSFUp = 1.;
          double lepSFDown = 1.;

          /*
          double tempValue = (double) rand() / (RAND_MAX);
          int leptonFileDicision = -99;

          if(tempValue < 20./35.9) leptonFileDicision = 0;
          else leptonFileDicision = 1;  
          */

          if(std::get<0>(samples[sam]) != "nonpromptData" && std::get<0>(samples[sam]) != "data"){

            dataMCSF = h_dataMC->GetBinContent(h_dataMC->GetXaxis()->FindBin(_nTrueInt));
              
          /*
            for(unsigned int leptonInd = 0; leptonInd < leptonSelectionAnalysis; leptonInd++){

              lepSF *= getLeptonSF(_lFlavor[ind.at(leptonInd)], _lPt[ind.at(leptonInd)], _lEta[ind.at(leptonInd)], 0, leptonFileDicision);

              lepSFUp *= getLeptonSF(_lFlavor[ind.at(leptonInd)], _lPt[ind.at(leptonInd)], _lEta[ind.at(leptonInd)], 1, leptonFileDicision);
              lepSFDown *= getLeptonSF(_lFlavor[ind.at(leptonInd)], _lPt[ind.at(leptonInd)], _lEta[ind.at(leptonInd)], -1, leptonFileDicision);
            }
          */
          }
          
          dataMCSF *= lepSF;
          weight = weight * dataMCSF;
          
          /*
          double btagSF_event = 1.;
          double btagSF_event_Up = 1.;
          double btagSF_event_Down = 1.;

          double btagSF_light_event = 1.;
          double btagSF_light_event_Up = 1.;
          double btagSF_light_event_Down = 1.;

          int btagFileDicision = 0;
             
          if(tempValue < 20./35.9)
            btagFileDicision = 0;
          else
            btagFileDicision = 1;  

          for(unsigned int i = 0; i < nJLoc; ++i){

            int j = indJets.at(i);

            double csv = _jetCsvV2[j];
            if( csv < 0.0 ) csv = -0.05;
            if( csv > 1.0 ) csv = 1.0;

            double pt = _jetPt[j];
            if( pt > 1000 ) pt = 999.;

            if(std::get<0>(samples[sam]) != "nonpromptData" && std::get<0>(samples[sam]) != "data"){

              double centralValue = getBTagSF(btagFileDicision, 0, _jetHadronFlavor[j], TMath::Abs(_jetEta[j]), pt, csv);
              btagSF_event *= centralValue;
                    
              if(fabs(_jetHadronFlavor[j]) != 0){
                btagSF_event_Up *= centralValue;
                btagSF_event_Down *= centralValue;
              } 
              else{
                btagSF_event_Up *= getBTagSF(btagFileDicision, 1, _jetHadronFlavor[j], TMath::Abs(_jetEta[j]), pt, csv);
                btagSF_event_Down *= getBTagSF(btagFileDicision, -1, _jetHadronFlavor[j], TMath::Abs(_jetEta[j]), pt, csv);
              }

              if(fabs(_jetHadronFlavor[j]) != 5){
                btagSF_light_event_Up *= centralValue;
                btagSF_light_event_Down *= centralValue;
              }
              else{
                btagSF_light_event_Up *= getBTagSF(btagFileDicision, 2, _jetHadronFlavor[j], TMath::Abs(_jetEta[j]), pt, csv);
                btagSF_light_event_Down *= getBTagSF(btagFileDicision, -2, _jetHadronFlavor[j], TMath::Abs(_jetEta[j]), pt, csv);
              }

            }

          }

          weight = weight * btagSF_event;
          */

          for(auto & i : indLoose){
            double mvaVL = -999.;
            vector<Float_t> varForBDT = { (Float_t)_lPt[i], (Float_t)_lEta[i], (Float_t)_selectedTrackMult[i], (Float_t)_miniIsoCharged[i], (Float_t)(_miniIso[i] - _miniIsoCharged[i]), (Float_t)_ptRel[i], (Float_t)_relIso[i], (Float_t)_ptRatio[i], (Float_t)_closestJetCsvV2[i], (Float_t)_3dIPSig[i], (Float_t)_dxy[i], (Float_t)_dz[i], _lFlavor[i] == 0 ? (Float_t)_lElectronMva[i] : (Float_t)_lMuonSegComp[i]};
            fillBDTvariables(varForBDT, _lFlavor[i]);
            mvaVL =  _lFlavor[i] == 0 ? readerLeptonMVAele->EvaluateMVA( "BDTG method") : readerLeptonMVAmu->EvaluateMVA( "BDTG method");
            _leptonMvaSUSY[i] = mvaVL;
          }


          distribs[0].vectorHisto[samCategory].Fill(TMath::Min(_met,varMax[0]-0.1),weight);
          
          distribs[1].vectorHisto[samCategory].Fill(TMath::Min(_lPt[indLoose.at(0)],varMax[1]-0.001), weight);
          distribs[2].vectorHisto[samCategory].Fill(TMath::Min(_lPt[indLoose.at(1)],varMax[2]-0.001), weight);

          distribs[3].vectorHisto[samCategory].Fill(TMath::Min(_lEta[indLoose.at(0)],varMax[3]-0.001), weight);
          distribs[4].vectorHisto[samCategory].Fill(TMath::Min(_lEta[indLoose.at(1)],varMax[4]-0.001), weight);

          distribs[5].vectorHisto[samCategory].Fill(TMath::Min(double(_nVertex),varMax[5]-0.001), weight);

          distribs[6].vectorHisto[samCategory].Fill(TMath::Min(mll,varMax[6]-0.001), weight);
        
          for(auto & i : indLoose){
            //if(lepIsGood(i)) continue;
            //if(i == indLoose.at(0)) continue;
            //if(_lFlavor[i] == 1) continue;
            distribs[7].vectorHisto[samCategory].Fill(TMath::Min(double(_closestJetCsvV2[i]),varMax[7]-0.001), weight);
            distribs[8].vectorHisto[samCategory].Fill(TMath::Min(double(_dxy[i]),varMax[8]-0.001), weight);
            distribs[9].vectorHisto[samCategory].Fill(TMath::Min(double(_dz[i]),varMax[9]-0.001), weight);
            distribs[10].vectorHisto[samCategory].Fill(TMath::Min(double(_3dIPSig[i]),varMax[10]-0.001), weight);
            distribs[11].vectorHisto[samCategory].Fill(TMath::Min(double(_ptRatio[i]),varMax[11]-0.001), weight);
            distribs[12].vectorHisto[samCategory].Fill(TMath::Min(double(_ptRel[i]),varMax[12]-0.001), weight);
            distribs[23].vectorHisto[samCategory].Fill(TMath::Min(double(_relIso[i]),varMax[23]-0.001), weight);

            if(_lFlavor[i] == 0){
                distribs[13].vectorHisto[samCategory].Fill(TMath::Min(double(_lElectronMvaHZZ[i]),varMax[13]-0.001), weight);
                distribs[14].vectorHisto[samCategory].Fill(TMath::Min(double(_lElectronMva[i]),varMax[14]-0.001), weight);
            }
            distribs[15].vectorHisto[samCategory].Fill(TMath::Min(double(_leptonMvaSUSY[i]),varMax[15]-0.001), weight);
            distribs[16].vectorHisto[samCategory].Fill(TMath::Min(double(_leptonMvaTTH[i]),varMax[16]-0.001), weight);
            distribs[17].vectorHisto[samCategory].Fill(TMath::Min(double(_miniIso[i]),varMax[17]-0.001), weight);
            distribs[18].vectorHisto[samCategory].Fill(TMath::Min(double(_miniIsoCharged[i]),varMax[18]-0.001), weight);
            if(_lFlavor[i] == 1){
                distribs[21].vectorHisto[samCategory].Fill(TMath::Min(double(_lMuonSegComp[i]),varMax[21]-0.001), weight);
                distribs[22].vectorHisto[samCategory].Fill(TMath::Min(double(_relIso0p4Mu[i]),varMax[22]-0.001), weight);
            }
         }
         distribs[19].vectorHisto[samCategory].Fill(TMath::Min(double(nJLoc),varMax[19]-0.001), weight);
         distribs[20].vectorHisto[samCategory].Fill(TMath::Min(double(nBLoc),varMax[20]-0.001), weight);

      }

      std::cout << std::endl;
      cout << "Total number of events: " << distribs[0].vectorHisto[sam].Integral() << endl;
      std::cout << std::endl;
  }

  /*
  fileDummy->cd();
  signalTree.Write();
  bkgTree.Write();

  fileDummy->Close();
  */

  TLegend* mtleg = new TLegend(0.77,0.89,0.95,0.62); 
  //mtleg->SetNColumns(1);
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
      
  for(int i = 0; i < nVars; i++){
      plot[i] = new TCanvas(Form("plot_%d", i),"",500,450);
  }

  vector<std::string> figNames = {"Type I E_{T}^{miss} [GeV]", "p_{T}^{leading} [GeV]", "p_{T}^{trailing} [GeV]", "#eta_{T}^{leading} [GeV]", "#eta_{T}^{trailing} [GeV]", "NPV", "M_{ll} [GeV]", "closestJetCSVv2", "dxy", "dz", "SIP3D", "ptratio", "ptrel","electron HZZ MVA","electron GP MVA", "SUSY lepton MVA", "TTH lepton MVA", "miniIso", "miniIsoCharged", "nJets", "nBJets", "muon segment comp", "relIso0p4", "relIso0p3"};  
  vector<TString> namesForSaveFiles = {"met", "ptlead", "pttrail", "etalead", "etatrail", "npv", "mll", "closejetcsv", "dxy", "dz", "sip3d", "ptratio", "ptrel", "electronMVAHZZ", "electronMVA", "leptonMVASUSY", "leptonMVATTH", "miniIso", "miniIsoChar", "njets", "nbjets", "muonSegmComp", "relIso0p4", "relIso0p3"};

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

  return;

}

int main(int argc, char *argv[]){

    TApplication *rootapp = new TApplication("example", &argc, argv);
    treeReader reader;
    reader.Analyze();
    rootapp->Run();

    return 0;
}

