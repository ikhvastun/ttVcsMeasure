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
  //gROOT->SetBatch(kTRUE);
  //read samples and cross sections from txt file
  //readSamples("samples_Zll_2017data.txt");
  readSamples("samples_Zll_2017data_2018MoriondMC.txt");
  //readSamples("test.txt");
  
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

  
  /*
  LastError::lasterror = Errors::UNKNOWN;

  vector<TH2D> fakeMaps;
  getFRmaps(fakeMaps);

  if(LastError::lasterror != Errors::OK){
     cout << "FR maps not found" << endl;
     return;
  }
  */

  //std::set<ULong64_t> eventsList = {1492195608, 146025307, 508933991, 321197298, 277740241, 892590768, 26652538, 315966641, 320474674, 59566750, 66424217, 144766203, 252985782, 65013805, 133493915, 102651335, 455363823, 111273860, 373187240, 239729165, 239078987, 195340703, 1856389, 461581059, 506460369, 656450946, 184079568, 963204027, 963441575, 336829296, 115879231, 41506348, 640566084, 78165014, 1255495340, 64807712, 2666356, 107674345, 529788719, 293132512, 39548847, 267438002, 55582650, 111863351, 1061848933, 217235889, 562481829, 704561987, 480406039, 744805798, 202616779, 327142299, 167976892, 391402062, 327599224, 503087687, 98671267, 198248356, 90666260, 406200233, 272203871, 163786055, 253703343, 142520294, 382842138, 405922138, 111982906, 210867535, 135955294, 252681643, 1150234279, 1103452740, 1262721932, 36893531, 312773489, 253682178, 450349644, 37685662, 38382410, 400846856, 552906843, 82761276, 664083895, 325380020, 585216792, 102898869, 18194142, 637834750, 573681949, 318118288, 529995837, 144910725, 14517387, 226434855, 192484572, 1617293625, 196222566, 75686784, 1147847462, 108502231, 1534207907, 1306527199, 131111757, 48742638, 131076926, 385195332, 694375811, 384522728, 894287433, 370558190, 21779924, 1099098547, 898202951, 1216158608, 253240568, 963202886, 367361713, 587481963, 115793233, 322860364, 211654778, 1206592648, 1655550800, 1350511364, 1695619142, 1754680495, 267327098, 1110306499, 5545114, 27071751, 211503876, 3377288405, 55721552, 618991682, 1117123544, 142703839, 517934154, 834180370, 1371104318, 1269000161, 2763422568, 1083197265, 1444157466, 1421307868, 1027796765, 1861225520, 1372234311, 172620241, 1040650747, 342353721, 234086389, 236104091, 76206227, 204990413, 39938667, 191786789, 325353201, 183963204, 158559514, 96976446, 157532648, 73280873, 97907529, 197833912, 697733424, 418136939, 305060877, 1273299716, 778786389, 131142593, 842492785, 332241563, 760618473, 244050843, 675289609, 694369509, 1053432589, 1266575101, 598168272, 638247408, 14814758, 1151944380, 1301114924, 101525287, 274827271, 1041988330, 321245229, 390563015, 508574825, 344972376, 522431310, 221056243, 40355959, 409853777, 200561756, 74081286, 89251178, 119661967, 207168464, 58509301, 251183962, 28591491, 30036062, 1124290562, 289965738, 1090947268, 699472740, 75939013, 411070835, 536170736, 144462225, 528847863, 260371307, 202531851, 1106760047, 82655785, 773526112, 248826088, 83999606, 297716316, 18276366, 1085186745, 372661440, 729677885, 1318455870};             

  //std::ofstream myfile;
  //myfile.open("myevents.txt");

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample();

      Color_t color = assignColor(std::get<0>(samples[sam]));
      setStackColors(color, sam);

      //if(std::get<0>(samples[sam]) != "data" && std::get<0>(samples[sam]) != "DY") continue; // PU reco reweigning 

      if(std::get<0>(samples[sam]) != "data") continue;
    
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
          //if(it > 1000000) break;
          
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
          
          if(isData){
              // should be used for Muon Prompt Reco
              //_passMETFilters = _Flag_HBHENoiseFilter && _Flag_HBHENoiseIsoFilter && _Flag_EcalDeadCellTriggerPrimitiveFilter && _Flag_goodVertices && _Flag_eeBadScFilter && _Flag_globalTightHalo2016Filter && _Flag_BadPFMuonFilter && _Flag_BadChargedCandidateFilter;
              
              if(!_passMETFilters) continue;
              if(!_2017_mm) continue;
          }

          /*
          if(printAddInfo){
            cout << "Filter info: " << _passMETFilters << " " << _Flag_HBHENoiseFilter << " " << _Flag_HBHENoiseIsoFilter << " " << _Flag_EcalDeadCellTriggerPrimitiveFilter << " " << _Flag_goodVertices << " " << _Flag_eeBadScFilter << " " << _Flag_globalTightHalo2016Filter << " " << _Flag_BadPFMuonFilter << " " << _Flag_BadChargedCandidateFilter << endl;
          }
          */

          //myfile << _runNb << " " << _lumiBlock << " " << _eventNb << " " << _met << " " << _passMETFilters << endl;

          std::vector<unsigned> ind;
          const unsigned lCount = selectLep(ind);

          /*
          if(printAddInfo)
            cout << "number of leptons: " << lCount << endl;
            */

          if(lCount != 2) continue;

          int samCategory = sam;
          int nLocEle = getElectronNumber(ind);

          if(nLocEle != 0) continue;

          if(!passPtCuts2L(ind)) continue;

          std::vector<unsigned> indJets;

          unsigned third = -9999;
          double mll = 99999;
          double pt_Z = 999999;
          double phi_Z = 999999;
          double ptNonZ = 999999;

          nJLoc = nJets(0, true, indJets, std::get<0>(samples[sam]) == "nonpromptData");
          //nBLoc = nBJets(0, false, true, 1, std::get<0>(samples[sam]) == "nonpromptData");
          double dMZ = deltaMZ(ind, third, mll, pt_Z, ptNonZ, phi_Z);

          //HTLoc = HTCalc(indJets);
          
          /*
          if(printAddInfo){
            cout << "number of jets/bjets/deltaMZ: " << nJLoc << " " << nBLoc << " " << dMZ << endl;
          }
          */

          //if(dMZ > 10) continue;
          if(nJLoc == 0) continue;

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

            dataMCSF = h_dataMC->GetBinContent(h_dataMC->GetXaxis()->FindBin(_nVertex));
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
          

          /*
          if(printAddInfo){
            cout << "This event passes all selection criteria" << endl;
          }
          */
          
          /*
          if(_met > 200){
            myfile << _runNb << " " << _lumiBlock << " " << _eventNb << " " << _met << " " << _passMETFilters << endl;
          }
          */

          double uPara =     ((( -_met*TMath::Cos(_metPhi )  - pt_Z * TMath::Cos( phi_Z))* pt_Z*TMath::Cos( phi_Z )+(- _met* TMath::Sin(_metPhi )- pt_Z*TMath::Sin(phi_Z ))* pt_Z*TMath::Sin( phi_Z ))/pt_Z + pt_Z);
          double uPara_raw = ((( -_rawmet*TMath::Cos(_rawmetPhi )  - pt_Z * TMath::Cos( phi_Z))* pt_Z*TMath::Cos( phi_Z )+(- _rawmet* TMath::Sin(_rawmetPhi )- pt_Z*TMath::Sin(phi_Z ))* pt_Z*TMath::Sin( phi_Z ))/pt_Z + pt_Z);
            
          double uPerp =     ((( -_met*TMath::Cos(_metPhi )  - pt_Z * TMath::Cos( phi_Z))* pt_Z*TMath::Sin( phi_Z )-(- _met* TMath::Sin(_metPhi )- pt_Z*TMath::Sin(phi_Z ))* pt_Z*TMath::Cos( phi_Z ))/pt_Z);
          double uPerp_raw = ((( -_rawmet*TMath::Cos(_rawmetPhi )  - pt_Z * TMath::Cos( phi_Z))* pt_Z*TMath::Sin( phi_Z )-(- _rawmet* TMath::Sin(_rawmetPhi )- pt_Z*TMath::Sin(phi_Z ))* pt_Z*TMath::Cos( phi_Z ))/pt_Z);
            

          distribs[0].vectorHisto[samCategory].Fill(TMath::Min(_rawmet,varMax[0]-0.1),weight);
          distribs[1].vectorHisto[samCategory].Fill(TMath::Min(_met,varMax[1]-0.1),weight);
          
          distribs[2].vectorHisto[samCategory].Fill(TMath::Min(uPara,varMax[2]-0.001), weight);
          distribs[3].vectorHisto[samCategory].Fill(TMath::Min(uPerp,varMax[3]-0.001), weight);

          distribs[4].vectorHisto[samCategory].Fill(TMath::Min(_lPt[ind.at(0)],varMax[4]-0.001), weight);
          distribs[5].vectorHisto[samCategory].Fill(TMath::Min(_lPt[ind.at(1)],varMax[5]-0.001), weight);

          distribs[6].vectorHisto[samCategory].Fill(TMath::Min(_lEta[ind.at(0)],varMax[6]-0.001), weight);
          distribs[7].vectorHisto[samCategory].Fill(TMath::Min(_lEta[ind.at(1)],varMax[7]-0.001), weight);

          distribs[8].vectorHisto[samCategory].Fill(TMath::Min(double(_nVertex),varMax[8]-0.001), weight);

          distribs[9].vectorHisto[samCategory].Fill(TMath::Min(mll,varMax[9]-0.001), weight);

          int runBin = -999;
          if(_runNb < 299329){
            runBin = 0;
          }
          else if(_runNb < 304826){
            runBin = 1;
          }
          else{
            runBin = 2;
          }


          distribs[0].histDataEras[runBin].Fill(TMath::Min(_rawmet,varMax[0]-0.1),weight);
          distribs[1].histDataEras[runBin].Fill(TMath::Min(_met,varMax[1]-0.1),weight);
          
          distribs[2].histDataEras[runBin].Fill(TMath::Min(uPara,varMax[2]-0.001), weight);
          distribs[3].histDataEras[runBin].Fill(TMath::Min(uPerp,varMax[3]-0.001), weight);

          distribs[4].histDataEras[runBin].Fill(TMath::Min(_lPt[ind.at(0)],varMax[4]-0.001), weight);
          distribs[5].histDataEras[runBin].Fill(TMath::Min(_lPt[ind.at(1)],varMax[5]-0.001), weight);

          distribs[6].histDataEras[runBin].Fill(TMath::Min(_lEta[ind.at(0)],varMax[6]-0.001), weight);
          distribs[7].histDataEras[runBin].Fill(TMath::Min(_lEta[ind.at(1)],varMax[7]-0.001), weight);

          distribs[8].histDataEras[runBin].Fill(TMath::Min(double(_nVertex),varMax[8]-0.001), weight);

          distribs[9].histDataEras[runBin].Fill(TMath::Min(mll,varMax[9]-0.001), weight);

          
          if(pt_Z < 18) continue;
          //continue;
          int binForPtZ = 0;
          for(int i = 0; i < nQt; i++){
            if(pt_Z < qtBins[i]){
              binForPtZ = i - 1;
              break;
            }
          }

          histMetCorr[binForPtZ]->Fill((uPara - pt_Z) / pt_Z);
          histMetUnCorr[binForPtZ]->Fill((uPara_raw - pt_Z) / pt_Z);

          sigmaParUnCorr[binForPtZ]->Fill(uPara_raw);
          sigmaPerpUnCorr[binForPtZ]->Fill(uPerp_raw);

          sigmaParCorr[binForPtZ]->Fill(uPara);
          sigmaPerpCorr[binForPtZ]->Fill(uPerp);

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


  vector<std::string> figNames = {"Raw E_{T}^{miss} [GeV]", "Type I E_{T}^{miss} [GeV]", "u_{||} + q_{T} [GeV]", "u_{#perp}   [GeV]", "p_{T}^{leading} [GeV]", "p_{T}^{trailing} [GeV]", "#eta_{T}^{leading} [GeV]", "#eta_{T}^{trailing} [GeV]", "NPV", "M_{ll} [GeV]"};
  vector<TString> namesForSaveFiles = {"rawmet", "met", "upara", "uperp", "ptlead", "pttrail", "etalead", "etatrail", "npv", "mll"};

  /*
  for(int varPlot = 0; varPlot < nVars; varPlot++){
    plot[varPlot]->cd();
    showHist(plot[varPlot],distribs[varPlot],"",figNames.at(varPlot),"Events", scale_num, mtleg); //  + std::to_string(int((varMax[varPlot] - varMin[varPlot])/nBins[varPlot]))
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + ".pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + ".png");
    plot[varPlot]->cd();
    showHist(plot[varPlot],distribs[varPlot],"",figNames.at(varPlot),"Events", scale_num, mtleg, true); //  + std::to_string(int((varMax[varPlot] - varMin[varPlot])/nBins[varPlot]))
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + "Log.pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + "Log.png");
  }
  */

  TLegend* mtlegData = new TLegend(0.77,0.89,0.95,0.62); 
  mtlegData->SetFillColor(0);
  mtlegData->SetFillStyle(0);
  mtlegData->SetBorderSize(0);
  mtlegData->SetTextFont(42);
  
  
  mtlegData->AddEntry(&distribs[0].histDataEras[runB],"run B","lep"); 
  mtlegData->AddEntry(&distribs[0].histDataEras[runCDE],"run CDE","lep"); //data
  mtlegData->AddEntry(&distribs[0].histDataEras[runF],"run F","lep"); //data


  for(int varPlot = 0; varPlot < nVars; varPlot++){
    plot[varPlot]->cd();
    showDataComp(plot[varPlot],distribs[varPlot],"",figNames.at(varPlot),"Events", scale_num, mtlegData); //  + std::to_string(int((varMax[varPlot] - varMin[varPlot])/nBins[varPlot]))
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + ".pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + ".png");
    plot[varPlot]->cd();
    showDataComp(plot[varPlot],distribs[varPlot],"",figNames.at(varPlot),"Events", scale_num, mtlegData, true); //  + std::to_string(int((varMax[varPlot] - varMin[varPlot])/nBins[varPlot]))
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + "Log.pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + "Log.png");
  }


  /*  
  //distribs[8].vectorHisto[1].Scale(1./distribs[8].vectorHisto[1].Integral());

  //distribs[8].vectorHisto[dataSample].Scale(1./distribs[8].vectorHisto[dataSample].Integral());

  TH1F *puReweign = (TH1F*)distribs[8].vectorHisto[dataSample].Clone("puw");
  TH1F *puReweignMC = (TH1F*)distribs[8].vectorHisto[1].Clone("MC"); // 1 stands for DY sample

  puReweign->Divide(puReweignMC);

  TCanvas *cRewCanvs = new TCanvas("cRewCanvs", "cRewCanvs");
  TFile *file = TFile::Open("puWeights.root","RECREATE");
  puReweign->Draw();
  puReweign->Write();
  file->Close();
*/
  return;


  TCanvas * c1 = new TCanvas("c1", "c1");
  c1->Divide(5,5);

  TH1D * scaleChoice[2];
  TH1D * resChoice[2][2];
  for(int i = 0; i < 2; i++){
    scaleChoice[i] = new TH1D(Form("scaleChoice_%d", i), Form("scaleChoice_%d", i), nQt - 1, qtBins);
    for(int j = 0; j < 2; j++)
      resChoice[i][j] = new TH1D(Form("resChoice_%d_%d", i, j), Form("resChoice_%d_%d", i, j), nQt - 1, qtBins);
  }


    TF1 * f1 = new TF1("f1", "[0] * TMath::Voigt(x - [1], [2], [3], 4) + [4] * x + [5]"   );

    f1->SetParameters(histMetCorr[0]->Integral(), histMetCorr[0]->GetMean(), histMetCorr[0]->GetRMS(), 0.25);

    for(int i = 0; i < 25; i++){

      c1->cd(i+1);
      histMetCorr[i]->Draw();
      histMetCorr[i]->Fit(f1, "", "", histMetCorr[i]->GetMean() - 4 * histMetCorr[i]->GetRMS(), histMetCorr[i]->GetMean() + 4 * histMetCorr[i]->GetRMS());
      histMetCorr[i]->GetXaxis()->SetRangeUser(histMetCorr[i]->GetMean() - 7 * histMetCorr[i]->GetRMS(), histMetCorr[i]->GetMean() + 7 * histMetCorr[i]->GetRMS());
     
      f1->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3), f1->GetParameter(4), f1->GetParameter(5));

      TLatex latex;
      latex.DrawLatex(histMetCorr[i]->GetMean() + 3 * histMetCorr[i]->GetRMS(), histMetCorr[i]->GetMaximum() / 2, Form("\\chi^{2} / ndf = %.1f/%i", f1->GetChisquare(), f1->GetNDF()));

      TLatex latex2;
      latex2.DrawLatex(histMetCorr[i]->GetMean() + 3 * histMetCorr[i]->GetRMS(), histMetCorr[i]->GetMaximum() / 3, Form("mean = %.2f", f1->GetParameter(1)));
      
      scaleChoice[0]->SetBinContent(i+1, -1 * f1->GetParameter(1));
      scaleChoice[0]->SetBinError(i+1, -1 * f1->GetParError(1));
    }

    TCanvas * c2 = new TCanvas("c2", "c2");
    c2->Divide(5,5);


    f1->SetParameters(histMetUnCorr[0]->Integral(), histMetUnCorr[0]->GetMean(), histMetUnCorr[0]->GetRMS(), 0.25);

    for(int i = 0; i < 25; i++){
      
      //f1->SetParameters(histMetUnCorr[i]->Integral(), histMetUnCorr[i]->GetMean(), histMetUnCorr[i]->GetRMS(), 0.1);
      c2->cd(i+1);
      histMetUnCorr[i]->Draw();
      histMetUnCorr[i]->Fit(f1, "", "", histMetUnCorr[i]->GetMean() - 4 * histMetUnCorr[i]->GetRMS(), histMetUnCorr[i]->GetMean() + 4 * histMetUnCorr[i]->GetRMS());
      histMetUnCorr[i]->GetXaxis()->SetRangeUser(histMetUnCorr[i]->GetMean() - 7 * histMetCorr[i]->GetRMS(), histMetUnCorr[i]->GetMean() + 7 * histMetCorr[i]->GetRMS());
      
      TLatex latex;
      latex.DrawLatex(histMetUnCorr[i]->GetMean() + 3 * histMetUnCorr[i]->GetRMS(), histMetUnCorr[i]->GetMaximum() / 2, Form("\\chi^{2} / ndf = %.1f/%i", f1->GetChisquare(), f1->GetNDF()));

      TLatex latex2;
      latex2.DrawLatex(histMetUnCorr[i]->GetMean() + 3 * histMetUnCorr[i]->GetRMS(), histMetUnCorr[i]->GetMaximum() / 3, Form("mean = %.2f", f1->GetParameter(1)));
      
      scaleChoice[1]->SetBinContent(i+1, -1 * f1->GetParameter(1));
      scaleChoice[1]->SetBinError(i+1, -1 * f1->GetParError(1));

      f1->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3), f1->GetParameter(4), f1->GetParameter(5));

    }

    TCanvas * c3 = new TCanvas("c3", "c3");
    c3->Divide(5,5);


    f1->SetParameters(sigmaPerpCorr[0]->Integral(), sigmaPerpCorr[0]->GetMean(), sigmaPerpCorr[0]->GetRMS());
    for(int i = 0; i < 25; i++){
      c3->cd(i+1);
      sigmaPerpCorr[i]->Draw();
      sigmaPerpCorr[i]->Fit(f1, "", "", sigmaPerpCorr[i]->GetMean() - 4 * sigmaPerpCorr[i]->GetRMS(), sigmaPerpCorr[i]->GetMean() + 4 * sigmaPerpCorr[i]->GetRMS());

      TLatex latex;
      latex.DrawLatex(sigmaPerpCorr[i]->GetMean() + 3 * sigmaPerpCorr[i]->GetRMS(), sigmaPerpCorr[i]->GetMaximum() / 2, Form("\\chi^{2} / ndf = %.1f/%i", f1->GetChisquare(), f1->GetNDF()));

      TLatex latex2;
      latex2.DrawLatex(sigmaPerpCorr[i]->GetMean() + 3 * sigmaPerpCorr[i]->GetRMS(), sigmaPerpCorr[i]->GetMaximum() / 3, Form("sigma = %.1f", f1->GetParameter(2)));
      

      resChoice[0][0]->SetBinContent(i+1, f1->GetParameter(2));
      resChoice[0][0]->SetBinError(i+1, f1->GetParError(2));

      f1->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3), f1->GetParameter(4), f1->GetParameter(5));

    }

    TCanvas * c4 = new TCanvas("c4", "c4");
    c4->Divide(5,5);


    f1->SetParameters(sigmaParCorr[0]->Integral(), sigmaParCorr[0]->GetMean(), sigmaParCorr[0]->GetRMS());
    for(int i = 0; i < 25; i++){
      c4->cd(i+1);
      sigmaParCorr[i]->Draw();
      sigmaParCorr[i]->Fit(f1, "", "", sigmaParCorr[i]->GetMean() - 4 * sigmaParCorr[i]->GetRMS(), sigmaParCorr[i]->GetMean() + 4 * sigmaParCorr[i]->GetRMS());
      
      TLatex latex;
      latex.DrawLatex(sigmaParCorr[i]->GetMean() + 3 * sigmaParCorr[i]->GetRMS(), sigmaParCorr[i]->GetMaximum() / 2, Form("\\chi^{2} / ndf = %.1f/%i", f1->GetChisquare(), f1->GetNDF()));

      TLatex latex2;
      latex2.DrawLatex(sigmaParCorr[i]->GetMean() + 3 * sigmaParCorr[i]->GetRMS(), sigmaParCorr[i]->GetMaximum() / 3, Form("sigma = %.1f", f1->GetParameter(2)));
      

      resChoice[0][1]->SetBinContent(i+1, f1->GetParameter(2));
      resChoice[0][1]->SetBinError(i+1, f1->GetParError(2));

      f1->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3), f1->GetParameter(4), f1->GetParameter(5));

    }

    TCanvas * c5 = new TCanvas("c5", "c5");
    c5->Divide(5,5);


    f1->SetParameters(sigmaPerpUnCorr[0]->Integral(), sigmaPerpUnCorr[0]->GetMean(), sigmaPerpUnCorr[0]->GetRMS());
    for(int i = 0; i < 25; i++){

      c5->cd(i+1);
      sigmaPerpUnCorr[i]->Draw();
      sigmaPerpUnCorr[i]->Fit(f1, "", "", sigmaPerpUnCorr[i]->GetMean() - 4 * sigmaPerpUnCorr[i]->GetRMS(), sigmaPerpUnCorr[i]->GetMean() + 4 * sigmaPerpUnCorr[i]->GetRMS());

      TLatex latex;
      latex.DrawLatex(sigmaPerpUnCorr[i]->GetMean() + 3 * sigmaPerpUnCorr[i]->GetRMS(), sigmaPerpUnCorr[i]->GetMaximum() / 2, Form("\\chi^{2} / ndf = %.1f/%i", f1->GetChisquare(), f1->GetNDF()));

      TLatex latex2;
      latex2.DrawLatex(sigmaPerpUnCorr[i]->GetMean() + 3 * sigmaPerpUnCorr[i]->GetRMS(), sigmaPerpUnCorr[i]->GetMaximum() / 3, Form("sigma = %.1f", f1->GetParameter(2)));
      
      resChoice[1][0]->SetBinContent(i+1, f1->GetParameter(2));
      resChoice[1][0]->SetBinError(i+1, f1->GetParError(2));

      f1->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3), f1->GetParameter(4), f1->GetParameter(5));

    }

    TCanvas * c6 = new TCanvas("c6", "c6");
    c6->Divide(5,5);


    f1->SetParameters(sigmaParUnCorr[0]->Integral(), sigmaParUnCorr[0]->GetMean(), sigmaParUnCorr[0]->GetRMS());
    for(int i = 0; i < 25; i++){

      c6->cd(i+1);
      sigmaParUnCorr[i]->Draw();
      sigmaParUnCorr[i]->Fit(f1, "", "", sigmaParUnCorr[i]->GetMean() - 4 * sigmaParUnCorr[i]->GetRMS(), sigmaParUnCorr[i]->GetMean() + 4 * sigmaParUnCorr[i]->GetRMS());
      
      TLatex latex;
      latex.DrawLatex(sigmaParUnCorr[i]->GetMean() + 3 * sigmaParUnCorr[i]->GetRMS(), sigmaParUnCorr[i]->GetMaximum() / 2, Form("\\chi^{2} / ndf = %.1f/%i", f1->GetChisquare(), f1->GetNDF()));

      TLatex latex2;
      latex2.DrawLatex(sigmaParUnCorr[i]->GetMean() + 3 * sigmaParUnCorr[i]->GetRMS(), sigmaParUnCorr[i]->GetMaximum() / 3, Form("sigma = %.1f", f1->GetParameter(2)));
      

      resChoice[1][1]->SetBinContent(i+1, f1->GetParameter(2));
      resChoice[1][1]->SetBinError(i+1, f1->GetParError(2));

      f1->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3), f1->GetParameter(4), f1->GetParameter(5));

    }

    TCanvas * c7 = new TCanvas("c7", "c7");
    scaleChoice[0]->Draw();
    scaleChoice[0]->SaveAs("scale/cor/scaleEl.root");

    TCanvas * c8 = new TCanvas("c8", "c8");
    scaleChoice[1]->Draw();
    
    scaleChoice[1]->SaveAs("scale/uncor/scaleEl.root");

    resChoice[0][0]->SaveAs("scale/cor/sigmaPerpEl.root");
    resChoice[0][1]->SaveAs("scale/cor/sigmaParEl.root");

    resChoice[1][0]->SaveAs("scale/uncor/sigmaPerpEl.root");
    resChoice[1][1]->SaveAs("scale/uncor/sigmaParEl.root");

    return;

}

int main(int argc, char *argv[]){

    TApplication *rootapp = new TApplication("example", &argc, argv);
    treeReader reader;
    reader.Analyze();
    rootapp->Run();

    return 0;
}

