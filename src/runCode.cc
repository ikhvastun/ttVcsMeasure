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
  readSamples("data/samples_Zll_2017data_2018MoriondMC.txt");
  //readSamples("test.txt");
  
  std::vector<std::string> namesOfSamples = treeReader::getNamesOfTheSample();
  initdistribs(namesOfSamples);

  // b tag SF are not applied here
  /*
  readerBtag[0][0].load(calib_csvv2[0], BTagEntry::FLAV_B, "iterativefit");
  readerBtag[0][1].load(calib_csvv2[0], BTagEntry::FLAV_C, "iterativefit");
  readerBtag[0][2].load(calib_csvv2[0], BTagEntry::FLAV_UDSG, "iterativefit");

  readerBtag[1][0].load(calib_csvv2[1], BTagEntry::FLAV_B, "iterativefit");
  readerBtag[1][1].load(calib_csvv2[1], BTagEntry::FLAV_C, "iterativefit");
  readerBtag[1][2].load(calib_csvv2[1], BTagEntry::FLAV_UDSG, "iterativefit");
  */

  
  //std::set<ULong64_t> eventsList = {1492195608, 146025307, 508933991, 321197298, 277740241, 892590768, 26652538, 315966641, 320474674, 59566750, 66424217, 144766203, 252985782, 65013805, 133493915, 102651335, 455363823, 111273860, 373187240, 239729165, 239078987, 195340703, 1856389, 461581059, 506460369, 656450946, 184079568, 963204027, 963441575, 336829296, 115879231, 41506348, 640566084, 78165014, 1255495340, 64807712, 2666356, 107674345, 529788719, 293132512, 39548847, 267438002, 55582650, 111863351, 1061848933, 217235889, 562481829, 704561987, 480406039, 744805798, 202616779, 327142299, 167976892, 391402062, 327599224, 503087687, 98671267, 198248356, 90666260, 406200233, 272203871, 163786055, 253703343, 142520294, 382842138, 405922138, 111982906, 210867535, 135955294, 252681643, 1150234279, 1103452740, 1262721932, 36893531, 312773489, 253682178, 450349644, 37685662, 38382410, 400846856, 552906843, 82761276, 664083895, 325380020, 585216792, 102898869, 18194142, 637834750, 573681949, 318118288, 529995837, 144910725, 14517387, 226434855, 192484572, 1617293625, 196222566, 75686784, 1147847462, 108502231, 1534207907, 1306527199, 131111757, 48742638, 131076926, 385195332, 694375811, 384522728, 894287433, 370558190, 21779924, 1099098547, 898202951, 1216158608, 253240568, 963202886, 367361713, 587481963, 115793233, 322860364, 211654778, 1206592648, 1655550800, 1350511364, 1695619142, 1754680495, 267327098, 1110306499, 5545114, 27071751, 211503876, 3377288405, 55721552, 618991682, 1117123544, 142703839, 517934154, 834180370, 1371104318, 1269000161, 2763422568, 1083197265, 1444157466, 1421307868, 1027796765, 1861225520, 1372234311, 172620241, 1040650747, 342353721, 234086389, 236104091, 76206227, 204990413, 39938667, 191786789, 325353201, 183963204, 158559514, 96976446, 157532648, 73280873, 97907529, 197833912, 697733424, 418136939, 305060877, 1273299716, 778786389, 131142593, 842492785, 332241563, 760618473, 244050843, 675289609, 694369509, 1053432589, 1266575101, 598168272, 638247408, 14814758, 1151944380, 1301114924, 101525287, 274827271, 1041988330, 321245229, 390563015, 508574825, 344972376, 522431310, 221056243, 40355959, 409853777, 200561756, 74081286, 89251178, 119661967, 207168464, 58509301, 251183962, 28591491, 30036062, 1124290562, 289965738, 1090947268, 699472740, 75939013, 411070835, 536170736, 144462225, 528847863, 260371307, 202531851, 1106760047, 82655785, 773526112, 248826088, 83999606, 297716316, 18276366, 1085186745, 372661440, 729677885, 1318455870};             

  //std::ofstream myfile;
  //myfile.open("myevents.txt");

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample();

      Color_t color = assignColor(std::get<0>(samples[sam]));
      setStackColors(color, sam);

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
          
          // check if event is in the list and print info for it
          /*
          const bool is_in = eventsList.find(_eventNb) != eventsList.end();
          if(!is_in) continue;
          */

          //bool printAddInfo = false;
          //if(_eventNb != 1492195608) continue;

          if(isData){
              if(!_passMETFilters) continue;
              // check all possible dilepton unprescaled triggers
              if(!(_passTrigger_mm || _passTrigger_ee)) continue;
          }

          std::vector<unsigned> ind;
          const unsigned lCount = selectLep(ind);

          if(lCount != 2) continue;

          int samCategory = sam;
          int nLocEle = getElectronNumber(ind);

          if(nLocEle != 0) continue;
          if(!passPtCuts2L(ind)) continue;

          int runBin = -999;
          int runBinForPU = -99;
          if(isData){
            if(_runNb < 299329) runBin = 0;
            else if(_runNb < 304826) runBin = 1;
            else runBin = 2;
          }
          else{
            double tempValue = ((double) rand() / (RAND_MAX));
            if(tempValue < 0.115424){ // lumi Run B / full 2017 dataset
              runBin = 0;
              runBinForPU = 0;
            }
            else if (tempValue < 0.673971){ // lumi Run CDE / full 2017 dataset
              runBin = 1;
              if(tempValue < 0.347458) //  lumi Run C / full 2017 dataset
                  runBinForPU = 1;
              else if (tempValue < 0.449734) // lumi Run D / full 2017 dataset
                  runBinForPU = 2;
              else // lumi Run E / full 2017 dataset
                  runBinForPU = 3;
            }
            else { // lumi Run F / full 2017 dataset
              runBin = 2;
              runBinForPU = 4;
            }
          }

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
          
          // onZ is already asked in skimming, no need here
          //if(dMZ > 10) continue;
          //if(nJLoc != 0) continue;

          double dataMCSF = 1.;
          double lepSF = 1.;
          double lepSFUp = 1.;
          double lepSFDown = 1.;

          if(std::get<0>(samples[sam]) != "nonpromptData" && std::get<0>(samples[sam]) != "data" && std::get<0>(samples[sam]) != "dataPrompt"){

            //dataMCSF = h_dataMC->GetBinContent(h_dataMC->GetXaxis()->FindBin(_nVertex));
            
            //if(_nTrueInt < 10) dataMCSF = 0.;
            dataMCSF = h_dataMC_2017[runBinForPU]->GetBinContent(h_dataMC_2017[runBinForPU]->GetXaxis()->FindBin(TMath::Min(double(_nTrueInt), 99.)));
            for(unsigned int leptonInd = 0; leptonInd < 2; leptonInd++){

              lepSF *= getLeptonSF(_lFlavor[ind.at(leptonInd)], _lPt[ind.at(leptonInd)], (_lFlavor[ind.at(leptonInd)] ? _lEta[ind.at(leptonInd)] : _lEtaSC[ind.at(leptonInd)]), 0, runBinForPU);

              //lepSFUp *= getLeptonSF(_lFlavor[ind.at(leptonInd)], _lPt[ind.at(leptonInd)], _lEta[ind.at(leptonInd)], 1, );
              //lepSFDown *= getLeptonSF(_lFlavor[ind.at(leptonInd)], _lPt[ind.at(leptonInd)], _lEta[ind.at(leptonInd)], -1, );
            }
            
          }
          
          
          weight = weight * dataMCSF;
          weight = weight * lepSF;
          
          double uPara =     paraCalc(_met, _metPhi, phi_Z, pt_Z);
          double uPara_raw = paraCalc(_rawmet, _rawmetPhi, phi_Z, pt_Z);
            
          double uPerp =     perpCalc(_met, _metPhi, phi_Z, pt_Z);
          double uPerp_raw = perpCalc(_rawmet, _rawmetPhi, phi_Z, pt_Z);
            
          
          // used for total number of events
          distribs[0].vectorHisto[samCategory].Fill(TMath::Min(_rawmet,varMax[0]-0.1),weight);
          
          distribs[0].vectorHistoEras[samCategory].Fill(_rawmet, runBin, _rawmetJECUp, _rawmetJECDown, _rawmetJetResUp, _rawmetJetResDown, _rawmetUnclUp, _rawmetUnclDown, varMax[0], weight);
          distribs[1].vectorHistoEras[samCategory].Fill(_met, runBin, _metJECUp, _metJECDown, _metJetResUp, _metJetResDown, _metUnclUp, _metUnclDown, varMax[1], weight);
          distribs[2].vectorHistoEras[samCategory].Fill(uPara, runBin, paraCalc(_metJECUp, _metPhiJECUp, phi_Z, pt_Z), paraCalc(_metJECDown, _metPhiJECDown, phi_Z, pt_Z), paraCalc(_metJetResUp, _metPhiJetResUp, phi_Z, pt_Z), paraCalc(_metJetResDown, _metPhiJetResDown, phi_Z, pt_Z), paraCalc(_metUnclUp, _metPhiUnclUp, phi_Z, pt_Z), paraCalc(_metUnclDown, _metPhiUnclDown, phi_Z, pt_Z), varMax[2], weight);
          distribs[3].vectorHistoEras[samCategory].Fill(uPerp, runBin, perpCalc(_metJECUp, _metPhiJECUp, phi_Z, pt_Z), perpCalc(_metJECDown, _metPhiJECDown, phi_Z, pt_Z), perpCalc(_metJetResUp, _metPhiJetResUp, phi_Z, pt_Z), perpCalc(_metJetResDown, _metPhiJetResDown, phi_Z, pt_Z), perpCalc(_metUnclUp, _metPhiUnclUp, phi_Z, pt_Z), perpCalc(_metUnclDown, _metPhiUnclDown, phi_Z, pt_Z), varMax[3], weight);
          distribs[4].vectorHistoEras[samCategory].Fill(_lPt[ind.at(0)], runBin, 0, 0, 0, 0, 0, 0, varMax[4], weight);
          distribs[5].vectorHistoEras[samCategory].Fill(_lPt[ind.at(1)], runBin, 0, 0, 0, 0, 0, 0, varMax[5], weight);
          distribs[6].vectorHistoEras[samCategory].Fill(_lEta[ind.at(0)], runBin, 0, 0, 0, 0, 0, 0, varMax[6], weight);
          distribs[7].vectorHistoEras[samCategory].Fill(_lEta[ind.at(1)], runBin, 0, 0, 0, 0, 0, 0, varMax[7], weight);
          distribs[8].vectorHistoEras[samCategory].Fill(double(_nVertex), runBin, 0, 0, 0, 0, 0, 0, varMax[8], weight);
          distribs[9].vectorHistoEras[samCategory].Fill(mll, runBin, 0, 0, 0, 0, 0, 0, varMax[9], weight);
          distribs[10].vectorHistoEras[samCategory].Fill(_met_sm, runBin, _metJECUp_sm, _metJECDown_sm, _metJetResUp_sm, _metJetResDown_sm, _metUnclUp_sm, _metUnclDown_sm, varMax[10], weight);

          if(pt_Z < 18) continue;
          //continue;
          int binForPtZ = 0;
          for(int i = 0; i < nQt; i++){
            if(pt_Z < qtBins[i]){
              binForPtZ = i - 1;
              break;
            }
          }

          histMetCorr[binForPtZ][samCategory == 0 ? 0 : 1][runBin]->Fill((uPara - pt_Z) / pt_Z, weight);
          histMetUnCorr[binForPtZ][samCategory == 0 ? 0 : 1][runBin]->Fill((uPara_raw - pt_Z) / pt_Z, weight);
          histMetCorr[binForPtZ][samCategory == 0 ? 0 : 1][3]->Fill((uPara - pt_Z) / pt_Z, weight);
          histMetUnCorr[binForPtZ][samCategory == 0 ? 0 : 1][3]->Fill((uPara_raw - pt_Z) / pt_Z, weight);

          sigmaParUnCorr[binForPtZ][samCategory == 0 ? 0 : 1][runBin]->Fill(uPara_raw, weight);
          sigmaPerpUnCorr[binForPtZ][samCategory == 0 ? 0 : 1][runBin]->Fill(uPerp_raw, weight);
          sigmaParUnCorr[binForPtZ][samCategory == 0 ? 0 : 1][3]->Fill(uPara_raw, weight);
          sigmaPerpUnCorr[binForPtZ][samCategory == 0 ? 0 : 1][3]->Fill(uPerp_raw, weight);

          sigmaParCorr[binForPtZ][samCategory == 0 ? 0 : 1][runBin]->Fill(uPara, weight);
          sigmaPerpCorr[binForPtZ][samCategory == 0 ? 0 : 1][runBin]->Fill(uPerp, weight);
          sigmaParCorr[binForPtZ][samCategory == 0 ? 0 : 1][3]->Fill(uPara, weight);
          sigmaPerpCorr[binForPtZ][samCategory == 0 ? 0 : 1][3]->Fill(uPerp, weight);
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

  vector<std::string> figNames = {"Raw E_{T}^{miss} [GeV]", "Type I E_{T}^{miss} [GeV]", "u_{||} + q_{T} [GeV]", "u_{#perp}   [GeV]", "p_{T}^{leading} [GeV]", "p_{T}^{trailing} [GeV]", "#eta_{T}^{leading} [GeV]", "#eta_{T}^{trailing} [GeV]", "NPV", "M_{ll} [GeV]", "E_{T}^{miss} smeared [GeV]"};
  vector<TString> namesForSaveFiles = {"rawmet", "met", "upara", "uperp", "ptlead", "pttrail", "etalead", "etatrail", "npv", "mll", "metSM"};
  vector<TString> runEraNames = {"RunB", "RunCDE", "RunF", "fullDataset"};

  /*
  for(int varPlot = 0; varPlot < nVars; varPlot++){

    for(int runEra = 0; runEra < 4; runEra++){
        plot[varPlot]->cd();
        showHist(plot[varPlot],distribs[varPlot],"",figNames.at(varPlot),"Events", scale_num, mtleg, false, runEra); //  + std::to_string(int((varMax[varPlot] - varMin[varPlot])/nBins[varPlot]))
        plot[varPlot]->SaveAs("plotsForSave/" + runEraNames[runEra] + "/nonlog/" + namesForSaveFiles.at(varPlot) + ".pdf");
        plot[varPlot]->SaveAs("plotsForSave/" + runEraNames[runEra] + "/nonlog/" + namesForSaveFiles.at(varPlot) + ".png");
        plot[varPlot]->SaveAs("plotsForSave/" + runEraNames[runEra] + "/nonlog/" + namesForSaveFiles.at(varPlot) + ".root");
        plot[varPlot]->cd();
        showHist(plot[varPlot],distribs[varPlot],"",figNames.at(varPlot),"Events", scale_num, mtleg, true, runEra); //  + std::to_string(int((varMax[varPlot] - varMin[varPlot])/nBins[varPlot]))
        plot[varPlot]->SaveAs("plotsForSave/" + runEraNames[runEra] + "/log/" + namesForSaveFiles.at(varPlot) + "Log.pdf");
        plot[varPlot]->SaveAs("plotsForSave/" + runEraNames[runEra] + "/log/" + namesForSaveFiles.at(varPlot) + "Log.png");
        plot[varPlot]->SaveAs("plotsForSave/" + runEraNames[runEra] + "/log/" + namesForSaveFiles.at(varPlot) + "Log.root");
    }
    
  }
  */

  TH1D * scaleChoice[2][2][4];
  TH1D * resChoice[2][2][2][4];
  for(int i = 0; i < 2; i++){
    for(int k = 0; k < 2; k++){
      for(int l = 0; l < 4; l++){
        scaleChoice[i][k][l] = new TH1D(Form("scaleChoice_%d_%d_%d", i, k, l), Form("scaleChoice_%d_%d_%d", i, k, l), nQt - 1, qtBins);
        for(int j = 0; j < 2; j++)
          resChoice[i][j][k][l] = new TH1D(Form("resChoice_%d_%d_%d_%d", i, j, k, l), Form("resChoice_%d_%d_%d_%d", i, j, k, l), nQt - 1, qtBins);
      }
    }
  }

  TF1 * f1 = new TF1("f1", "[0] * TMath::Voigt(x - [1], [2], [3], 4) + [4] * x + [5]"   );
  for(int runEra = 0; runEra < 4; runEra++){

    TCanvas * c1[2];
    for(int j = 0; j < 2; j++){
      c1[j] = new TCanvas(Form("c1_%d", j), Form("c1_%d", j));
      f1->SetParameters(histMetCorr[0][j][runEra]->Integral(), histMetCorr[0][j][runEra]->GetMean(), histMetCorr[0][j][runEra]->GetRMS(), 0.1);
      c1[j]->Divide(5,5);
      for(int i = 0; i < 25; i++){
        c1[j]->cd(i+1);
      
        histMetCorr[i][j][runEra]->Draw("hist");
        f1->SetParLimits(3, 0, 1);
        histMetCorr[i][j][runEra]->Fit(f1, "", "", histMetCorr[i][j][runEra]->GetMean() - 4 * histMetCorr[i][j][runEra]->GetRMS(), histMetCorr[i][j][runEra]->GetMean() + 4 * histMetCorr[i][j][runEra]->GetRMS());
        histMetCorr[i][j][runEra]->GetXaxis()->SetRangeUser(histMetCorr[i][j][runEra]->GetMean() - 7 * histMetCorr[i][j][runEra]->GetRMS(), histMetCorr[i][j][runEra]->GetMean() + 7 * histMetCorr[i][j][runEra]->GetRMS());
        //f1->Draw("same");
        f1->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3), f1->GetParameter(4), f1->GetParameter(5));

        TLatex latex;
        latex.DrawLatex(histMetCorr[i][j][runEra]->GetMean() + 3 * histMetCorr[i][j][runEra]->GetRMS(), histMetCorr[i][j][runEra]->GetMaximum() / 2, Form("\\chi^{2} / ndf = %.1f/%i", f1->GetChisquare(), f1->GetNDF()));

        TLatex latex2;
        latex2.DrawLatex(histMetCorr[i][j][runEra]->GetMean() + 3 * histMetCorr[i][j][runEra]->GetRMS(), histMetCorr[i][j][runEra]->GetMaximum() / 3, Form("mean = %.2f", f1->GetParameter(1)));
        
        //scaleChoice[0][j][runEra]->SetBinContent(i+1, -1 * f1->GetParameter(1));
        //scaleChoice[0][j][runEra]->SetBinError(i+1, -1 * f1->GetParError(1));

        scaleChoice[0][j][runEra]->SetBinContent(i+1, -1 * histMetCorr[i][j][runEra]->GetMean());
        scaleChoice[0][j][runEra]->SetBinError(i+1, 0.0001);
      }
      c1[j]->SaveAs("scale/fits/" + (TString) runEraNames[runEra] + "/" + (TString) (j == 0 ? "data" : "MC") + "scaleCorr.pdf");
      c1[j]->SaveAs("scale/fits/" + (TString) runEraNames[runEra] + "/" + (TString) (j == 0 ? "data" : "MC") + "scaleCorr.png");
    }

    TCanvas * c2[2];

    for(int j = 0; j < 2; j++){

      c2[j] = new TCanvas(Form("c2_%d", j), Form("c2_%d", j));
      f1->SetParameters(histMetUnCorr[0][j][runEra]->Integral(), histMetUnCorr[0][j][runEra]->GetMean(), histMetUnCorr[0][j][runEra]->GetRMS(), 0.1);
      c2[j]->Divide(5,5);

      for(int i = 0; i < 25; i++){
        c2[j]->cd(i+1);

        histMetUnCorr[i][j][runEra]->Draw("hist");
        f1->SetParLimits(3, 0, 1);
        histMetUnCorr[i][j][runEra]->Fit(f1, "", "", histMetUnCorr[i][j][runEra]->GetMean() - 4 * histMetUnCorr[i][j][runEra]->GetRMS(), histMetUnCorr[i][j][runEra]->GetMean() + 4 * histMetUnCorr[i][j][runEra]->GetRMS());
        histMetUnCorr[i][j][runEra]->GetXaxis()->SetRangeUser(histMetUnCorr[i][j][runEra]->GetMean() - 7 * histMetCorr[i][j][runEra]->GetRMS(), histMetUnCorr[i][j][runEra]->GetMean() + 7 * histMetCorr[i][j][runEra]->GetRMS());
        //f1->Draw("same");
        f1->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3), f1->GetParameter(4), f1->GetParameter(5));
      
        TLatex latex;
        latex.DrawLatex(histMetUnCorr[i][j][runEra]->GetMean() + 3 * histMetUnCorr[i][j][runEra]->GetRMS(), histMetUnCorr[i][j][runEra]->GetMaximum() / 2, Form("\\chi^{2} / ndf = %.1f/%i", f1->GetChisquare(), f1->GetNDF()));

        TLatex latex2;
        latex2.DrawLatex(histMetUnCorr[i][j][runEra]->GetMean() + 3 * histMetUnCorr[i][j][runEra]->GetRMS(), histMetUnCorr[i][j][runEra]->GetMaximum() / 3, Form("mean = %.2f", f1->GetParameter(1)));
        
        //scaleChoice[1][j][runEra]->SetBinContent(i+1, -1 * f1->GetParameter(1));
        //scaleChoice[1][j][runEra]->SetBinError(i+1, -1 * f1->GetParError(1));

        scaleChoice[1][j][runEra]->SetBinContent(i+1, -1 * histMetUnCorr[i][j][runEra]->GetMean());
        scaleChoice[1][j][runEra]->SetBinError(i+1, 0.0001);

      }
      c2[j]->SaveAs("scale/fits/" + (TString) runEraNames[runEra] + "/" + (TString) (j == 0 ? "data" : "MC") + "scaleUnCorr.pdf");
      c2[j]->SaveAs("scale/fits/" + (TString) runEraNames[runEra] + "/" + (TString) (j == 0 ? "data" : "MC") + "scaleUnCorr.png");
    }

    
    TCanvas * c3[2];
    
    for(int j = 0; j < 2; j++){

      c3[j] = new TCanvas(Form("c3_%d", j), Form("c3_%d", j));
      f1->SetParameters(sigmaPerpCorr[0][j][runEra]->Integral(), sigmaPerpCorr[0][j][runEra]->GetMean(), sigmaPerpCorr[0][j][runEra]->GetRMS());
      c3[j]->Divide(5,5);

      for(int i = 0; i < 25; i++){
        c3[j]->cd(i+1);
        sigmaPerpCorr[i][j][runEra]->Draw();
        f1->SetParLimits(3, 0, 1);
        sigmaPerpCorr[i][j][runEra]->Fit(f1, "", "", sigmaPerpCorr[i][j][runEra]->GetMean() - 4 * sigmaPerpCorr[i][j][runEra]->GetRMS(), sigmaPerpCorr[i][j][runEra]->GetMean() + 4 * sigmaPerpCorr[i][j][runEra]->GetRMS());

        TLatex latex;
        latex.DrawLatex(sigmaPerpCorr[i][j][runEra]->GetMean() + 3 * sigmaPerpCorr[i][j][runEra]->GetRMS(), sigmaPerpCorr[i][j][runEra]->GetMaximum() / 2, Form("\\chi^{2} / ndf = %.1f/%i", f1->GetChisquare(), f1->GetNDF()));

        TLatex latex2;
        latex2.DrawLatex(sigmaPerpCorr[i][j][runEra]->GetMean() + 3 * sigmaPerpCorr[i][j][runEra]->GetRMS(), sigmaPerpCorr[i][j][runEra]->GetMaximum() / 3, Form("sigma = %.1f", f1->GetParameter(2)));
        

        //resChoice[0][0][j][runEra]->SetBinContent(i+1, f1->GetParameter(2));
        //resChoice[0][0][j][runEra]->SetBinError(i+1, f1->GetParError(2));

        resChoice[0][0][j][runEra]->SetBinContent(i+1, sigmaPerpCorr[i][j][runEra]->GetRMS());
        resChoice[0][0][j][runEra]->SetBinError(i+1, 0.0001);

        f1->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3), f1->GetParameter(4), f1->GetParameter(5));
      }
      c3[j]->SaveAs("scale/fits/" + (TString) runEraNames[runEra] + "/" + (TString) (j == 0 ? "data" : "MC") + "sigmaPerpCorr.pdf");
      c3[j]->SaveAs("scale/fits/" + (TString) runEraNames[runEra] + "/" + (TString) (j == 0 ? "data" : "MC") + "sigmaPerpCorr.png");
    }

    TCanvas * c4[2];

    for(int j = 0; j < 2; j++){

      c4[j] = new TCanvas(Form("c4_%d", j), Form("c4_%d", j));
      f1->SetParameters(sigmaParCorr[0][j][runEra]->Integral(), sigmaParCorr[0][j][runEra]->GetMean(), sigmaParCorr[0][j][runEra]->GetRMS());
      c4[j]->Divide(5,5);

    
      for(int i = 0; i < 25; i++){
        c4[j]->cd(i+1);
        sigmaParCorr[i][j][runEra]->Draw();
        f1->SetParLimits(3, 0, 1);
        sigmaParCorr[i][j][runEra]->Fit(f1, "", "", sigmaParCorr[i][j][runEra]->GetMean() - 4 * sigmaParCorr[i][j][runEra]->GetRMS(), sigmaParCorr[i][j][runEra]->GetMean() + 4 * sigmaParCorr[i][j][runEra]->GetRMS());
        
        TLatex latex;
        latex.DrawLatex(sigmaParCorr[i][j][runEra]->GetMean() + 3 * sigmaParCorr[i][j][runEra]->GetRMS(), sigmaParCorr[i][j][runEra]->GetMaximum() / 2, Form("\\chi^{2} / ndf = %.1f/%i", f1->GetChisquare(), f1->GetNDF()));

        TLatex latex2;
        latex2.DrawLatex(sigmaParCorr[i][j][runEra]->GetMean() + 3 * sigmaParCorr[i][j][runEra]->GetRMS(), sigmaParCorr[i][j][runEra]->GetMaximum() / 3, Form("sigma = %.1f", f1->GetParameter(2)));
        

        //resChoice[0][1][j][runEra]->SetBinContent(i+1, f1->GetParameter(2));
        //resChoice[0][1][j][runEra]->SetBinError(i+1, f1->GetParError(2));

        resChoice[0][1][j][runEra]->SetBinContent(i+1, sigmaParCorr[i][j][runEra]->GetRMS());
        resChoice[0][1][j][runEra]->SetBinError(i+1, 0.0001);

        f1->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3), f1->GetParameter(4), f1->GetParameter(5));
      }
      c4[j]->SaveAs("scale/fits/" + (TString) runEraNames[runEra] + "/" + (TString) (j == 0 ? "data" : "MC") + "sigmaParCorr.pdf");
      c4[j]->SaveAs("scale/fits/" + (TString) runEraNames[runEra] + "/" + (TString) (j == 0 ? "data" : "MC") + "sigmaParCorr.png");
    }

    TCanvas * c5[2];

    for(int j = 0; j < 2; j++){

      c5[j] = new TCanvas(Form("c5_%d", j), Form("c5_%d", j));
      f1->SetParameters(sigmaPerpUnCorr[0][j][runEra]->Integral(), sigmaPerpUnCorr[0][j][runEra]->GetMean(), sigmaPerpUnCorr[0][j][runEra]->GetRMS());
      c5[j]->Divide(5,5);
    
      for(int i = 0; i < 25; i++){

        c5[j]->cd(i+1);
        sigmaPerpUnCorr[i][j][runEra]->Draw();
        f1->SetParLimits(3, 0, 1);
        sigmaPerpUnCorr[i][j][runEra]->Fit(f1, "", "", sigmaPerpUnCorr[i][j][runEra]->GetMean() - 4 * sigmaPerpUnCorr[i][j][runEra]->GetRMS(), sigmaPerpUnCorr[i][j][runEra]->GetMean() + 4 * sigmaPerpUnCorr[i][j][runEra]->GetRMS());

        TLatex latex;
        latex.DrawLatex(sigmaPerpUnCorr[i][j][runEra]->GetMean() + 3 * sigmaPerpUnCorr[i][j][runEra]->GetRMS(), sigmaPerpUnCorr[i][j][runEra]->GetMaximum() / 2, Form("\\chi^{2} / ndf = %.1f/%i", f1->GetChisquare(), f1->GetNDF()));

        TLatex latex2;
        latex2.DrawLatex(sigmaPerpUnCorr[i][j][runEra]->GetMean() + 3 * sigmaPerpUnCorr[i][j][runEra]->GetRMS(), sigmaPerpUnCorr[i][j][runEra]->GetMaximum() / 3, Form("sigma = %.1f", f1->GetParameter(2)));
        
        //resChoice[1][0][j][runEra]->SetBinContent(i+1, f1->GetParameter(2));
        //resChoice[1][0][j][runEra]->SetBinError(i+1, f1->GetParError(2));

        resChoice[1][0][j][runEra]->SetBinContent(i+1, sigmaPerpUnCorr[i][j][runEra]->GetRMS());
        resChoice[1][0][j][runEra]->SetBinError(i+1, 0.0001);

        f1->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3), f1->GetParameter(4), f1->GetParameter(5));
      }
      c5[j]->SaveAs("scale/fits/" + (TString) runEraNames[runEra] + "/" + (TString) (j == 0 ? "data" : "MC") + "sigmaPerpUnCorr.pdf");
      c5[j]->SaveAs("scale/fits/" + (TString) runEraNames[runEra] + "/" + (TString) (j == 0 ? "data" : "MC") + "sigmaPerpUnCorr.png");
    }

    TCanvas * c6[2];
    for(int j = 0; j < 2; j++){

      c6[j] = new TCanvas(Form("c6_%d", j), Form("c6_%d", j));
      f1->SetParameters(sigmaParUnCorr[0][j][runEra]->Integral(), sigmaParUnCorr[0][j][runEra]->GetMean(), sigmaParUnCorr[0][j][runEra]->GetRMS());
      c6[j]->Divide(5,5);
    
      for(int i = 0; i < 25; i++){

        c6[j]->cd(i+1);
        sigmaParUnCorr[i][j][runEra]->Draw();
        f1->SetParLimits(3, 0, 1);
        sigmaParUnCorr[i][j][runEra]->Fit(f1, "", "", sigmaParUnCorr[i][j][runEra]->GetMean() - 4 * sigmaParUnCorr[i][j][runEra]->GetRMS(), sigmaParUnCorr[i][j][runEra]->GetMean() + 4 * sigmaParUnCorr[i][j][runEra]->GetRMS());
        
        TLatex latex;
        latex.DrawLatex(sigmaParUnCorr[i][j][runEra]->GetMean() + 3 * sigmaParUnCorr[i][j][runEra]->GetRMS(), sigmaParUnCorr[i][j][runEra]->GetMaximum() / 2, Form("\\chi^{2} / ndf = %.1f/%i", f1->GetChisquare(), f1->GetNDF()));

        TLatex latex2;
        latex2.DrawLatex(sigmaParUnCorr[i][j][runEra]->GetMean() + 3 * sigmaParUnCorr[i][j][runEra]->GetRMS(), sigmaParUnCorr[i][j][runEra]->GetMaximum() / 3, Form("sigma = %.1f", f1->GetParameter(2)));
        

        //resChoice[1][1][j][runEra]->SetBinContent(i+1, f1->GetParameter(2));
        //resChoice[1][1][j][runEra]->SetBinError(i+1, f1->GetParError(2));

        resChoice[1][1][j][runEra]->SetBinContent(i+1, sigmaParUnCorr[i][j][runEra]->GetRMS());
        resChoice[1][1][j][runEra]->SetBinError(i+1, 0.0001);

        f1->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3), f1->GetParameter(4), f1->GetParameter(5));
      }

      c6[j]->SaveAs("scale/fits/" + (TString) runEraNames[runEra] + "/" + (TString) (j == 0 ? "data" : "MC") + "sigmaParUnCorr.pdf");
      c6[j]->SaveAs("scale/fits/" + (TString) runEraNames[runEra] + "/" + (TString) (j == 0 ? "data" : "MC") + "sigmaParUnCorr.png");
    }
  }

  
  TCanvas * cScale[2][4];
  for(int runEra = 0; runEra < 4; runEra++){
    
    TLegend * leg = new TLegend(0.2, 0.7, 0.5, 0.9);
    leg->AddEntry(scaleChoice[0][0][runEra], "data", "lep");
    leg->AddEntry(scaleChoice[0][1][runEra], "MC", "lep");
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);

    for(int i = 0; i < 2; i++){
      cScale[i][runEra] = new TCanvas(Form("cScale_%d_%d", i, runEra), Form("cScale_%d_%d", i, runEra));
      scaleChoice[i][0][runEra]->GetXaxis()->SetTitle("q_{T} [GeV]");
      scaleChoice[i][0][runEra]->GetYaxis()->SetTitle("-<u_{||}>/<q_{T}>");
      scaleChoice[i][0][runEra]->GetYaxis()->SetRangeUser((i == 0 ? 0.95 : 0.8), (i == 0 ? 1.05 : 1.0));
      scaleChoice[i][0][runEra]->Draw();
      scaleChoice[i][1][runEra]->SetLineColor(kRed);
      scaleChoice[i][1][runEra]->SetMarkerColor(kRed);
      scaleChoice[i][1][runEra]->Draw("same");
      leg->Draw("same");
      scaleChoice[i][0][runEra]->SaveAs("scale/" + (TString) runEraNames[runEra] + "/" + "scale" + (TString)(i == 0 ? "Cor" : "UnCor") + "Data.root");
      scaleChoice[i][1][runEra]->SaveAs("scale/" + (TString) runEraNames[runEra] + "/" + "scale" + (TString)(i == 0 ? "Cor" : "UnCor") + "MC.root");
      //cScale[i][runEra]->SaveAs("scale/" + (runEra == 3 ? (TString) "fullDataset/" : (TString) "Run" + runEraNames[runEra] + "/") + "scale" + (TString)(i == 0 ? "Cor" : "UnCor") +".pdf");
      //cScale[i][runEra]->SaveAs("scale/" + (runEra == 3 ? (TString) "fullDataset/" : (TString) "Run" + runEraNames[runEra] + "/") + "scale" + (TString)(i == 0 ? "Cor" : "UnCor") +".png");
    
    }
    
    TCanvas * cRes[2][2][4];
    for(int i = 0; i < 2; i++){
      for(int j = 0; j < 2; j++){
        cRes[i][j][runEra] = new TCanvas(Form("cRes_%d_%d_%d", i, j,runEra), Form("cRes_%d_%d_%d", i, j,runEra));
        resChoice[i][j][0][runEra]->GetXaxis()->SetTitle("q_{T} [GeV]");
        resChoice[i][j][0][runEra]->GetYaxis()->SetTitle("#sigma(u_{" + (TString)(j == 0 ? "#perp}  " : "||}") + ") [GeV]");
        resChoice[i][j][0][runEra]->Draw();
        resChoice[i][j][1][runEra]->SetLineColor(kRed);
        resChoice[i][j][1][runEra]->SetMarkerColor(kRed);
        resChoice[i][j][1][runEra]->Draw("same");
        leg->Draw("same");
        resChoice[i][j][0][runEra]->SaveAs("scale/" + (TString) runEraNames[runEra] + "/" + "sigma" + (TString)(j == 0 ? "Perp" : "Par") + (i == 0 ? "Cor" : "UnCor") + "Data.root");
        resChoice[i][j][1][runEra]->SaveAs("scale/" + (TString) runEraNames[runEra] + "/" + "sigma" + (TString)(j == 0 ? "Perp" : "Par") + (i == 0 ? "Cor" : "UnCor") + "MC.root");
        //cRes[i][j][runEra]->SaveAs("scale/" + (runEra == 3 ? (TString) "fullDataset/" : (TString) "Run" + runEraNames[runEra] + "/") + "sigma" + (TString)(j == 0 ? "Perp" : "Par")  + (TString)(i == 0 ? "Cor" : "UnCor") +".pdf");
        //cRes[i][j][runEra]->SaveAs("scale/" + (runEra == 3 ? (TString) "fullDataset/" : (TString) "Run" + runEraNames[runEra] + "/") + "sigma" + (TString)(j == 0 ? "Perp" : "Par")  + (TString)(i == 0 ? "Cor" : "UnCor") +".png");
    
      }
    }
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

