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
  readSamples("samples_ttZ_2017data.txt");
  //readSamples("nonpromptMCForCT.txt");
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

          if(!isData)
            if(!(_HLT_Ele27_WPTight_Gsf || _HLT_IsoMu24 || _HLT_IsoTkMu24)) continue;
          else
            if(!(_HLT_Ele35_WPTight_Gsf || _HLT_IsoMu27)) continue;

          std::vector<unsigned> indTight, indFake;
          //select leptons
          const unsigned lCount = selectLep(indTight);
          const unsigned lCountFake = selectFakeLep(indFake);

          if(printAddInfo)
            cout << "number of leptons: " << lCount << " " << lCountFake << endl;

          std::vector<unsigned> ind;
          //used for main analysis
          if(std::get<0>(samples[sam]) != "nonpromptData"){
            if(lCount != leptonSelectionAnalysis) continue;
            for(auto & i : indTight)
            ind.push_back(i);
          }
          else{
            if(lCountFake != leptonSelectionAnalysis) continue;
            if(lCount >= leptonSelectionAnalysis) continue;
            for(auto & i : indFake)
            ind.push_back(i);
          }

          int samCategory = sam;

          // for CT
          /*
          int samCategory = -999;
          if(lCount == leptonSelectionAnalysis && lCountFake == leptonSelection)
            samCategory = 0;
          else if(lCount < leptonSelectionAnalysis && lCountFake == leptonSelection)
            samCategory = 1;
          else
            continue;
          */
          
          //std::vector<unsigned> & ind = indTight; // indFake is used for MC


          if(leptonSelectionAnalysis == 2)
            if(_lCharge[ind.at(0)] * _lCharge[ind.at(1)] < 0) continue;

          //cout << "new event#################" << endl;

          bool allLeptonsArePrompt = true;

          
          if(std::get<0>(samples[sam]) != "data" && std::get<0>(samples[sam]) != "nonpromptData")
            allLeptonsArePrompt = promptLeptons(ind);
          
          if(std::get<0>(samples[sam]) == "chargeMisID" && !allLeptonsArePrompt) continue;
          //if(std::get<0>(samples[sam]) == "nonprompt" && allLeptonsArePrompt) continue; // works just for MC

          if((std::get<0>(samples[sam]) == "ttZ" || std::get<0>(samples[sam]) == "ttX" || std::get<0>(samples[sam]) == "WZ" || std::get<0>(samples[sam]) == "rare") && !allLeptonsArePrompt) continue;

          int nLocEle = getElectronNumber(ind);


          if(leptonSelectionAnalysis == 3)
            if(!passPtCuts3L(ind)) continue;
          if(leptonSelectionAnalysis == 2)
            if(!passPtCuts2L(ind)) continue;

          std::vector<unsigned> indJets;

          unsigned third = -9999;
          double mll = 99999;
          double ptZ = 999999;
          double ptNonZ = 999999;

          nJLoc = nJets(0, true, indJets, std::get<0>(samples[sam]) == "nonpromptData");
          nBLoc = nBJets(0, false, true, 1, std::get<0>(samples[sam]) == "nonpromptData");
          double dMZ = deltaMZ(ind, third, mll, ptZ, ptNonZ);

          HTLoc = HTCalc(indJets);
          
          if(printAddInfo){
            cout << "number of jets/bjets/deltaMZ: " << nJLoc << " " << nBLoc << " " << dMZ << endl;
            cout << endl;
          }

          
          double mt1 = 9999;
          if(leptonSelectionAnalysis == 3){
            
              /*
            if(nJLoc < 3) continue;
            if(nBLoc < 1) continue;
            if(dMZ > 10) continue;
            */
            

            // WZ CR
            if(dMZ > 10) continue;

            if(nBLoc != 0) continue;
            if(_met < 30) continue;

            int leptThirdIndex = -999;

            for(int leptThird = 0; leptThird < ptCorrV.size(); leptThird++){
              if(ptCorrV[leptThird].second == third)
                leptThirdIndex = leptThird;
            }

            TLorentzVector l0p4;
            l0p4.SetPtEtaPhiE(ptCorrV[leptThirdIndex].first, _lEta[third], _lPhi[third], _lE[third] * ptCorrV[leptThirdIndex].first / _lPt[third]);

            mt1 = mtCalc(l0p4, _met, _metPhi);

            //if(mt1 < 50) continue;

            // ttbar CR
            //if(!(dMZ == 999999 || !((dMZ < 10) || (nBLoc < 1)))) continue;
          }
          
          
          if(leptonSelectionAnalysis == 2){
            if(nJLoc < 2) continue;
            if(nBLoc < 1) continue;
            if(nLocEle == 2 && dMZ < 10) continue;
            if(_met < 30) continue;
          }

          //if(nLocEle != 3) continue;
          
          double dataMCSF = 1.;
          double lepSF = 1.;
          double lepSFUp = 1.;
          double lepSFDown = 1.;

          
          double tempValue = (double) rand() / (RAND_MAX);
          int leptonFileDicision = -99;

          if(tempValue < 20./35.9) leptonFileDicision = 0;
          else leptonFileDicision = 1;  

          if(std::get<0>(samples[sam]) != "nonpromptData" && std::get<0>(samples[sam]) != "data"){

            dataMCSF = h_dataMC->GetBinContent(h_dataMC->GetXaxis()->FindBin(_nTrueInt));
              
            for(unsigned int leptonInd = 0; leptonInd < leptonSelectionAnalysis; leptonInd++){

              lepSF *= getLeptonSF(_lFlavor[ind.at(leptonInd)], _lPt[ind.at(leptonInd)], _lEta[ind.at(leptonInd)], 0, leptonFileDicision);

              lepSFUp *= getLeptonSF(_lFlavor[ind.at(leptonInd)], _lPt[ind.at(leptonInd)], _lEta[ind.at(leptonInd)], 1, leptonFileDicision);
              lepSFDown *= getLeptonSF(_lFlavor[ind.at(leptonInd)], _lPt[ind.at(leptonInd)], _lEta[ind.at(leptonInd)], -1, leptonFileDicision);
            }
          }
          
          dataMCSF *= lepSF;
          weight = weight * dataMCSF;
          

          
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
          

          
          double FRloc = 1.;

          if(std::get<0>(samples[sam]) == "nonpromptData"){ 

            int nFakeLepCounter = 0;

            for(int i = 0; i < lCountFake; i++){
              if(lepIsGood(ind.at(i))) continue;    
              double leptFakePtCorr = ptFake(_lPt[ind.at(i)], _ptRatio[ind.at(i)], _lFlavor[ind.at(i)], _leptonMvaTTH[ind.at(i)], _lPOGMedium[ind.at(i)]);
              double FRloc_loc = fakeMaps.at(_lFlavor[ind.at(i)]).GetBinContent(fakeMaps.at(_lFlavor[ind.at(i)]).FindBin(TMath::Min(leptFakePtCorr,ptBins[nPt-1]-1.), fabs(_lEta[ind.at(i)])));

              FRloc *= FRloc_loc/(1.-FRloc_loc);
              //FRloc *= FRloc_loc;
              nFakeLepCounter++;

            }

            FRloc *= TMath::Power(-1, nFakeLepCounter + 1);
          }

          if(printAddInfo){
            cout << "the FR is: " << FRloc << endl;
          }
          
          weight = weight * FRloc;  

          if(printAddInfo){
            cout << "total weight is: " << weight << endl;
          }
          

          double mvaVL = 0;

          if(leptonSelectionAnalysis == 2){

            TLorentzVector l0p4, l1p4;

            l0p4.SetPtEtaPhiE(ptCorrV[0].first, _lEta[ind.at(0)], _lPhi[ind.at(0)], _lE[ind.at(0)] * ptCorrV[0].first / _lPt[ind.at(0)]);
            l1p4.SetPtEtaPhiE(ptCorrV[1].first, _lEta[ind.at(1)], _lPhi[ind.at(1)], _lE[ind.at(1)] * ptCorrV[1].first / _lPt[ind.at(1)]);

            minDeltaR = deltaRCalc(indJets, ind.at(1));

            double mt1, mt2;
            mt1 = mtCalc(l0p4, _met, _metPhi);
            mt2 = mtCalc(l1p4, _met, _metPhi);
            mtHighest = mt1 > mt2 ? mt1 : mt2;
            mtLowest  = mt1 > mt2 ? mt2 : mt1;

            leadpt = ptCorrV[0].first;
            trailpt = ptCorrV[1].first;
            leadingJetPt = _jetPt[indJets.at(0)];
            trailJetPt = _jetPt[indJets.at(1)];

            //nJLoc = nJets;
            //nBLoc = nBJets;
            //HTLoc = HTCalc(indJets);
            MET = _met;

            /*
            _weightEventInTree = TMath::Abs(weight);
            if(std::get<0>(samples[sam]) == "TTW")
              signalTree->Fill();
            else if(std::get<0>(samples[sam]) == "nonprompt")
              bkgTree->Fill();
            */

            vector<Float_t> varForBDT = {(Float_t)nJLoc, (Float_t)nBLoc, (Float_t)HTLoc, (Float_t)MET, (Float_t)minDeltaR, (Float_t)leadpt, (Float_t)trailpt, (Float_t)mtHighest, (Float_t)mtLowest, (Float_t) leadingJetPt, (Float_t) trailJetPt};

            fillBDTvariables(varForBDT);

            mvaVL = reader->EvaluateMVA( "BDTG method"           );
          }

          /*
          if(leptonSelectionAnalysis == 2){
            if(mvaVL > 0) continue;
          }
          */

          distribs[4].vectorHisto[samCategory].Fill(TMath::Min(mt1,varMax[4]-0.1),weight);

          if(leptonSelectionAnalysis == 3)
            if(mt1 < 50) continue;

          distribs[12].vectorHisto[samCategory].Fill(TMath::Min(mvaVL,varMax[12]-0.001),weight);

          int mvaValueRegion = -999;
          if(mvaVL > 0 && mvaVL < 0.6)
            mvaValueRegion = 0;
          if(mvaVL > 0.6)
            mvaValueRegion = 1;

          if(mvaVL < 0)
            mvaValueRegion = 2;

          distribs[0].vectorHisto[samCategory].Fill(TMath::Min(ptCorrV[0].first,varMax[0]-0.1),weight);
          distribs[1].vectorHisto[samCategory].Fill(TMath::Min(ptCorrV[1].first,varMax[1]-0.1),weight);
          if(leptonSelectionAnalysis > 2)
            distribs[2].vectorHisto[samCategory].Fill(TMath::Min(ptCorrV[2].first,varMax[2]-0.1),weight);

          distribs[8].vectorHisto[samCategory].Fill(TMath::Min(double(nJLoc),varMax[8]-0.1),weight);
          distribs[9].vectorHisto[samCategory].Fill(TMath::Min(double(nBLoc),varMax[9]-0.1),weight);

          if(leptonSelectionAnalysis == 3){
            distribs[14].vectorHisto[samCategory].Fill(flavourCategory3L(nLocEle),weight);
            distribs[15].vectorHisto[samCategory].Fill(SRID3L(nJLoc, nBLoc),weight);

            distribs[15].vectorHistoUp[samCategory].Fill(SRID3L(nJLoc, nBLoc),weight / lepSF * lepSFUp);
            distribs[15].vectorHistoDown[samCategory].Fill(SRID3L(nJLoc, nBLoc),weight / lepSF * lepSFDown);
          }

          if(leptonSelectionAnalysis == 2){
            distribs[14].vectorHisto[samCategory].Fill(flavourCategory2L(nLocEle, _lCharge[ind.at(0)]),weight);
            distribs[15].vectorHisto[samCategory].Fill(SRID2L(nJLoc, nBLoc, mvaValueRegion, _lCharge[ind.at(0)]),weight);

            distribs[15].vectorHistoUp[samCategory].Fill(SRID2L(nJLoc, nBLoc, mvaValueRegion, _lCharge[ind.at(0)]),weight / lepSF * lepSFUp);
            distribs[15].vectorHistoDown[samCategory].Fill(SRID2L(nJLoc, nBLoc, mvaValueRegion, _lCharge[ind.at(0)]),weight / lepSF * lepSFDown);
          }

          distribs[16].vectorHisto[samCategory].Fill(TMath::Min(mll,varMax[16]-0.1),weight);
          distribs[17].vectorHisto[samCategory].Fill(TMath::Min(ptZ,varMax[17]-0.1),weight);
          distribs[18].vectorHisto[samCategory].Fill(TMath::Min(ptNonZ,varMax[18]-0.1),weight);

          if(nLocEle == 3)
            distribs[19].vectorHisto[samCategory].Fill(TMath::Min(mll,varMax[19]-0.1),weight);
          if(nLocEle == 2)
            distribs[20].vectorHisto[samCategory].Fill(TMath::Min(mll,varMax[20]-0.1),weight);
          if(nLocEle == 1)
            distribs[21].vectorHisto[samCategory].Fill(TMath::Min(mll,varMax[21]-0.1),weight);
          if(nLocEle == 0)
            distribs[22].vectorHisto[samCategory].Fill(TMath::Min(mll,varMax[22]-0.1),weight);

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
  
  TCanvas* plot[16];
      
  for(int i = 0; i < 16; i++){
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
  showHist(plot[11],distribs[4],"","m_{T}^{W} [GeV]","Events",scale_num, mtleg);

  plot[12]->cd();
  showHist(plot[12],distribs[19],"","M(ll) in 3e [GeV]","Events",scale_num, mtleg);

  plot[13]->cd();
  showHist(plot[13],distribs[20],"","M(ll) in 2e1mu [GeV]","Events",scale_num, mtleg);

  plot[14]->cd();
  showHist(plot[14],distribs[21],"","M(ll) in 1e2mu [GeV]","Events",scale_num, mtleg);

  plot[15]->cd();
  showHist(plot[15],distribs[22],"","M(ll) in 3mu [GeV]","Events",scale_num, mtleg);

  vector<TString> namesForSaveFiles = {"ptlead", "sublead", "trail", "njets", "nbjets", "SR", "flavour", "BDT", "mll", "ptZ", "ptNonZ", "mtW", "mll3e", "mll2e1mu", "mll1e2mu", "mll3mu"};
  int countPlot = 0;
  for(auto & i : namesForSaveFiles){
    plot[countPlot]->SaveAs("plotsForSave/" + i + ".pdf");
    countPlot++;
  }
  //plot[8]->cd();
  //drawSystUnc(plot[8], distribs[15], 8);
  

  fillDatacards(distribs[15], samplesOrderNames, samplesOrder);
}

int main(int argc, char *argv[]){

    TApplication *rootapp = new TApplication("example", &argc, argv);
    treeReader reader;
    reader.Analyze();
    rootapp->Run();

    return 0;
}

