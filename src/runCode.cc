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
  readSamples("samples_Zll_2017data.txt");
  //readSamples("test.txt");
  setTDRStyle(); 
  
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
  std::set<ULong64_t> eventsList = {1492195608, 146025307, 508933991, 321197298, 277740241, 892590768, 26652538, 315966641, 320474674, 59566750, 66424217, 144766203, 252985782, 65013805, 133493915, 102651335, 455363823, 111273860, 373187240, 239729165, 239078987, 195340703, 1856389, 461581059, 506460369, 656450946, 184079568, 963204027, 963441575, 336829296, 115879231, 41506348, 640566084, 78165014, 1255495340, 64807712, 2666356, 107674345, 529788719, 293132512, 39548847, 267438002, 55582650, 111863351, 1061848933, 217235889, 562481829, 704561987, 480406039, 744805798, 202616779, 327142299, 167976892, 391402062, 327599224, 503087687, 98671267, 198248356, 90666260, 406200233, 272203871, 163786055, 253703343, 142520294, 382842138, 405922138, 111982906, 210867535, 135955294, 252681643, 1150234279, 1103452740, 1262721932, 36893531, 312773489, 253682178, 450349644, 37685662, 38382410, 400846856, 552906843, 82761276, 664083895, 325380020, 585216792, 102898869, 18194142, 637834750, 573681949, 318118288, 529995837, 144910725, 14517387, 226434855, 192484572, 1617293625, 196222566, 75686784, 1147847462, 108502231, 1534207907, 1306527199, 131111757, 48742638, 131076926, 385195332, 694375811, 384522728, 894287433, 370558190, 21779924, 1099098547, 898202951, 1216158608, 253240568, 963202886, 367361713, 587481963, 115793233, 322860364, 211654778, 1206592648, 1655550800, 1350511364, 1695619142, 1754680495, 267327098, 1110306499, 5545114, 27071751, 211503876, 3377288405, 55721552, 618991682, 1117123544, 142703839, 517934154, 834180370, 1371104318, 1269000161, 2763422568, 1083197265, 1444157466, 1421307868, 1027796765, 1861225520, 1372234311, 172620241, 1040650747, 342353721, 234086389, 236104091, 76206227, 204990413, 39938667, 191786789, 325353201, 183963204, 158559514, 96976446, 157532648, 73280873, 97907529, 197833912, 697733424, 418136939, 305060877, 1273299716, 778786389, 131142593, 842492785, 332241563, 760618473, 244050843, 675289609, 694369509, 1053432589, 1266575101, 598168272, 638247408, 14814758, 1151944380, 1301114924, 101525287, 274827271, 1041988330, 321245229, 390563015, 508574825, 344972376, 522431310, 221056243, 40355959, 409853777, 200561756, 74081286, 89251178, 119661967, 207168464, 58509301, 251183962, 28591491, 30036062, 1124290562, 289965738, 1090947268, 699472740, 75939013, 411070835, 536170736, 144462225, 528847863, 260371307, 202531851, 1106760047, 82655785, 773526112, 248826088, 83999606, 297716316, 18276366, 1085186745, 372661440, 729677885, 1318455870};             

  std::ofstream myfile;
  myfile.open("myevents.txt");

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample();

      Color_t color = assignColor(std::get<0>(samples[sam]));
      setStackColors(color, sam);

      //if(leptonSelectionAnalysis == 3)
      //  if(std::get<0>(samples[sam]) == "chargeMisID") continue;

      //if(std::get<0>(samples[sam]) != "data") continue;
    
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
          
          /*
          const bool is_in = eventsList.find(_eventNb) != eventsList.end();
          if(!is_in) continue;
          */

          bool printAddInfo = false;
          //if(_eventNb != 1492195608) continue;
          if(printAddInfo){
            cout << "################################################################" << endl;
            cout << "Event number: " << _runNb << " " << _lumiBlock << " " << _eventNb << endl;
          }

          /*
          if(!isData)
            if(!(_HLT_Ele27_WPTight_Gsf || _HLT_IsoMu24 || _HLT_IsoTkMu24)) continue;
          else
            if(!(_HLT_Ele35_WPTight_Gsf || _HLT_IsoMu27)) continue;
            */
          if(!_2017_mm) continue;
          _passMETFilters = _Flag_HBHENoiseFilter && _Flag_HBHENoiseIsoFilter && _Flag_EcalDeadCellTriggerPrimitiveFilter && _Flag_goodVertices && _Flag_eeBadScFilter && _Flag_globalTightHalo2016Filter && _Flag_BadPFMuonFilter && _Flag_BadChargedCandidateFilter;
          if(printAddInfo){
            cout << "Filter info: " << _passMETFilters << " " << _Flag_HBHENoiseFilter << " " << _Flag_HBHENoiseIsoFilter << " " << _Flag_EcalDeadCellTriggerPrimitiveFilter << " " << _Flag_goodVertices << " " << _Flag_eeBadScFilter << " " << _Flag_globalTightHalo2016Filter << " " << _Flag_BadPFMuonFilter << " " << _Flag_BadChargedCandidateFilter << endl;
          }

          myfile << _runNb << " " << _lumiBlock << " " << _eventNb << " " << _met << " " << _passMETFilters << endl;

          std::vector<unsigned> ind;
          //select leptons
          const unsigned lCount = selectLep(ind);

          if(printAddInfo)
            cout << "number of leptons: " << lCount << endl;

          //used for main analysis
          if(lCount != 2) continue;

          int samCategory = sam;

          int nLocEle = getElectronNumber(ind);

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
          }

          if(dMZ > 10) continue;
          if(nJLoc != 0) continue;

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
          

          if(printAddInfo){
            cout << "This event passes all selection criteria" << endl;
            //cout << "total weight is: " << weight << endl;
          }
          
          /*
          if(_met > 200){
            myfile << _runNb << " " << _lumiBlock << " " << _eventNb << " " << _met << " " << _passMETFilters << endl;
          }
          */

          distribs[0].vectorHisto[samCategory].Fill(TMath::Min(_rawmet,varMax[0]-0.1),weight);
          distribs[1].vectorHisto[samCategory].Fill(TMath::Min(_met,varMax[1]-0.1),weight);

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

  TLegend* mtleg = new TLegend(0.45,0.89,0.95,0.77); 
  mtleg->SetNColumns(3);
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
  
  TCanvas* plot[2];
      
  for(int i = 0; i < 2; i++){
      plot[i] = new TCanvas(Form("plot_%d", i),"",500,450);
  }

  plot[0]->cd();
  showHist(plot[0],distribs[0],"","Raw E_{T}^{miss} [GeV]","Events / " + std::to_string(int((varMax[0] - varMin[0])/nBins[0])) + " GeV",scale_num, mtleg);

  plot[1]->cd();
  showHist(plot[1],distribs[1],"","Type I E_{T}^{miss} [GeV]","Events / " + std::to_string(int((varMax[1] - varMin[1])/nBins[1])) + " GeV",scale_num, mtleg);

  vector<TString> namesForSaveFiles = {"rawmet", "met"};
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

