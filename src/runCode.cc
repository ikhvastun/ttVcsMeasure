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

void treeReader::Analyze(const string& fileToAnalyse, const std::string option, const string& sampleToDebug, long evNb){
//void treeReader::Analyze(){

  debug = (option == "debug" ? true : false);
  leptonSelection = leptonSelectionAnalysis;
  //Set CMS plotting style
  setTDRStyle();
  //gROOT->SetBatch(kTRUE);
  //read samples and cross sections from txt file
  cout << "reading sample file...." << endl;
  readSamples(fileToAnalyse);
  cout << "finished with reading"<< endl;
  setTDRStyle(); 

  std::vector<std::string> namesOfSamples = treeReader::getNamesOfTheSample();
  cout << "initiating histos...." << endl;
  initdistribs(namesOfSamples);
  cout << "finished with initiating of histos"<< endl;

  
  if(leptonSelectionAnalysis == 2){

    // for reading values from xml files
    addVariablesToBDT(is2017);

    // for adding variables to the trees
    //addBranchToBDTTreeVariables();
  }

  std::ofstream myfile;
  myfile.open("myevents.txt");

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample();

      Color_t color = assignColor(samples[sam].getProcessName());
      setStackColors(color, sam);

      if((option == "runOnOneProcess" || debug) && (samples[sam].getProcessName()) != sampleToDebug) continue;
      if(samples[sam].getProcessName() == "nonpromptData"){
          cout << "Total number of events: " << distribs[0].vectorHisto[sam].Integral() << endl;
          continue;
      }

      if(leptonSelectionAnalysis == 3 || leptonSelectionAnalysis == 4)
        if((samples[sam].getProcessName()) == "chargeMisID") continue;

      std::cout<<"Entries in "<< (samples[sam].getFileName()) << " " << nEntries << std::endl;
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
          if(debug && _eventNb != evNb) continue;
          if(debug) cout << "event weight is " << weight << " after get entry" << endl; 
          
          if(!(_passTrigger_e || _passTrigger_m || _passTrigger_ee || _passTrigger_em || _passTrigger_mm || _passTrigger_eee || _passTrigger_eem || _passTrigger_emm || _passTrigger_mmm)) continue;
          if(!_passMETFilters) continue;
          
          //if(it > 10000) break;
          //if(it > nEntries / 20) break;

          std::vector<unsigned> indTight, indFake, indOf2LonZ;
          //select leptons
          const unsigned lCount = selectLep(indTight, leptonSelection);
          const unsigned lCountFake = selectFakeLep(indFake, leptonSelection);

          if(invMassOfAny2Lbelow12GeV(indFake)) continue;

          int samCategory = sam;

          std::vector<unsigned> ind;
          //used for main analysis
          // remove additional FO lepton from ss2l selection
          if(leptonSelection == 2 && lCountFake != leptonSelection) continue;
          if(leptonSelection == 3){
            std::vector<unsigned> indTight4L;
            const unsigned lCount4L = selectLep(indTight4L, 4);
            if(lCount4L >= 4) continue;
          }
          if(leptonSelection == 2){
            std::vector<unsigned> indTight3L;
            const unsigned lCount3L = selectLep(indTight3L, 3);
            if(lCount3L >= 3) continue;
          }
          if(lCount > leptonSelection) continue;
          if(lCountFake < leptonSelection) continue;
          if(lCount == leptonSelection) {
              samCategory = sam;
              ind = indTight;
          }
          if(lCount < leptonSelection){
              if(leptonSelectionAnalysis == 4) continue;
              samCategory = nonPromptSample;
              ind = indFake;
          }

          if(leptonSelectionAnalysis == 2)
            if(_lCharge[ind.at(0)] * _lCharge[ind.at(1)] < 0) continue;

          bool allLeptonsArePrompt = true;
          
          if((samples[sam].getProcessName()) != "data" && (samples[sam].getProcessName()) != "nonpromptData")
            allLeptonsArePrompt = promptLeptons(ind);
          
          if((samples[sam].getProcessName()) == "chargeMisID" && !allLeptonsArePrompt) continue;
          if((samples[sam].getProcessName()) == "Nonprompt" && allLeptonsArePrompt) continue; // works just for MC

          if(((samples[sam].getProcessName()) == "ttW" || (samples[sam].getProcessName()) == "ttH" || (samples[sam].getProcessName()) == "ttZ" || (samples[sam].getProcessName()) == "ttX" || (samples[sam].getProcessName()) == "WZ" || (samples[sam].getProcessName()) == "Z#gamma"  || (samples[sam].getProcessName()) == "ZZ" || (samples[sam].getProcessName()) == "rare") && !allLeptonsArePrompt) continue;

          int nLocEle = getElectronNumber(ind);
          if(nLocEle != 3) continue;

          if(leptonSelectionAnalysis == 4)
            if(!passPtCuts4L(ind)) continue;
          if(leptonSelectionAnalysis == 3)
            if(!passPtCuts3L(ind)) continue;
          if(leptonSelectionAnalysis == 2)
            if(!passPtCuts2L(ind)) continue;

          std::vector<unsigned> indJets;
          std::vector<unsigned> indJetsNotB;

          unsigned third = -9999;
          double mll = 99999;
          double mlll = 99999;
          double ptZ = 999999;
          double ptNonZ = 999999;

          nJLoc = nJets(0, true, indJets, is2017);
          nJLocNotB = nJetsNotB(0, true, indJetsNotB, 2, is2017);
          nBLoc = nBJets(0, true, true, 1, is2017);

          TLorentzVector Zboson, lnegative;
          double dMZ = deltaMZ(ind, third, mll, ptZ, ptNonZ, mlll, indOf2LonZ, Zboson, lnegative);
          double mll1stpair, mll2ndpair;
          double cosTSt = -999;

          if(debug) cout << "number of jets/bjets/dMZ: " << nJLoc << " " << nBLoc << " " << dMZ << endl;

          HTLoc = HTCalc(indJets);
          
          double mt1 = 9999;
          if(leptonSelectionAnalysis == 4){
            // ZZ control region
            if(dMZ > 20) continue; // 20 GeV originally in TOP-17-005
            mll1stpair = mll;
            std::vector<unsigned> vectorRemoved;
            vectorRemoved = ind;
            vectorRemoved.erase(std::remove(vectorRemoved.begin(), vectorRemoved.end(), indOf2LonZ.at(0)), vectorRemoved.end());
            vectorRemoved.erase(std::remove(vectorRemoved.begin(), vectorRemoved.end(), indOf2LonZ.at(1)), vectorRemoved.end());
            dMZ = deltaMZ(vectorRemoved, third, mll, ptZ, ptNonZ, mlll, indOf2LonZ, Zboson, lnegative);
            if(debug) cout << "second deltaMZ is : " << dMZ << endl;
            if(dMZ > 20) continue;
            mll2ndpair = mll;
          }

          if(leptonSelectionAnalysis == 3){
            
            if(!passTTZSelection(nJLoc, dMZ)) continue;
            //if(!passTTZCleanSelection(nJLoc, nBLoc, dMZ)) continue;
            //if(!passWZCRSelection(nBLoc, dMZ)) continue;
            //if(!passDYCRSelection(dMZ, ptNonZ, third, _met, _metPhi, nJLoc, nBLoc)) continue;
            //if(!passttbarCRSelection(nBLoc, dMZ)) continue;
            //passZGCRSelection(mlll, dMZ);
            cosTSt = cosThetaStar(Zboson, lnegative);

          }
          
          if(leptonSelectionAnalysis == 2){
            pass2Lpreselection(nJLoc, nBLoc, ind, _met, nLocEle);
          }
          
          double triggerSF = 1.;
          if(!isData) triggerSF = getTriggerSF(leptonSelection, ptCorrV[0].first);
          if(debug) cout << "trigger SF: " << triggerSF << endl;
          weight = weight * triggerSF;

          if(debug) cout << "weight is " << weight <<  " after trigger SF" << endl;
          
          double mvaVL = 0;

          if(leptonSelectionAnalysis == 2){

            TLorentzVector l0p4, l1p4, metTL;

            l0p4.SetPtEtaPhiE(ptCorrV[0].first, _lEta[ind.at(0)], _lPhi[ind.at(0)], _lE[ind.at(0)] * ptCorrV[0].first / _lPt[ind.at(0)]);
            l1p4.SetPtEtaPhiE(ptCorrV[1].first, _lEta[ind.at(1)], _lPhi[ind.at(1)], _lE[ind.at(1)] * ptCorrV[1].first / _lPt[ind.at(1)]);
            metTL.SetPtEtaPhiE(_met, 0, _metPhi, 0);

            minDeltaRlead = deltaRCalc(indJets, ind.at(0), true);
            minDeltaR = deltaRCalc(indJets, ind.at(1), true);

            double mt1, mt2;
            mt1 = mtCalc(l0p4, _met, _metPhi);
            mt2 = mtCalc(l1p4, _met, _metPhi);
            mtHighest = mt1 > mt2 ? mt1 : mt2;
            mtLowest  = mt1 > mt2 ? mt2 : mt1;

            leadpt = ptCorrV[0].first;
            trailpt = ptCorrV[1].first;
            leadeta = _lEta[ptCorrV[0].second];
            traileta = _lEta[ptCorrV[1].second];
            if(nJLocNotB > 0)
                leadingJetPt = _jetPt[indJetsNotB.at(0)];
            else
                leadingJetPt = -999.;
            if(nJLocNotB > 1)
                trailJetPt = _jetPt[indJetsNotB.at(1)];
            else
                trailJetPt = -999.;

            //nJLoc = nJets;
            //nBLoc = nBJets;
            //HTLoc = HTCalc(indJets);
            MET = _met;
            chargeOfLeptons = _lCharge[ind.at(0)];
            ll_deltaR = l0p4.DeltaR(l1p4);
            mll_ss = (l0p4+l1p4).M();
            mt2ll_ss = mt2ll(l0p4, l1p4, metTL);

            _weightEventInTree = weight;
            //_weightEventInTree = 1;
            /*
            if(std::get<0>(samples[sam]) == "ttW")
              signalTree->Fill();
            else if(std::get<0>(samples[sam]) == "Nonprompt" || std::get<0>(samples[sam]) == "chargeMisID" || std::get<0>(samples[sam]) == "ttZ" || std::get<0>(samples[sam]) == "ttH" || std::get<0>(samples[sam]) == "WZ" || std::get<0>(samples[sam]) == "nonpromptData")
              bkgTree->Fill();
            */
            
            vector<Float_t> varForBDT = {(Float_t)nJLoc, (Float_t)nBLoc, (Float_t)HTLoc, (Float_t)MET, (Float_t)minDeltaRlead, (Float_t)minDeltaR, (Float_t)l0p4.DeltaR(l1p4), (Float_t)mtHighest, (Float_t)mtLowest, (Float_t)leadpt, (Float_t)trailpt, (Float_t) leadingJetPt, (Float_t) trailJetPt, (Float_t)leadeta, (Float_t)traileta, (Float_t)chargeOfLeptons, (Float_t)(l0p4+l1p4).M(), (Float_t)mt2ll_ss};
            fillBDTvariables(varForBDT);
            mvaVL = readerTTWcsttbar->EvaluateMVA( "BDTG method");
            /*
            cout << "all variables: ";
            for(auto & value : varForBDT)
                cout << value << " ";
            cout << endl;
            cout << "mva value is " << mvaVL << endl;
            */

            //if(mvaVL > 0) continue;
            
          }

          //if(chargeOfLeptons != -1) continue;
          if((samples[sam].getProcessName()) != "data" && (samples[sam].getProcessName()) != "nonpromptData")
            weight *= sfWeight();
          if(samples[sam].getProcessName() == "data" && samCategory == nonPromptSample)
            weight *= fakeRateWeight();

          if(leptonSelectionAnalysis == 3){

            // WZ CR
            TLorentzVector l0p4;
            l0p4.SetPtEtaPhiE(ptNonZ, _lEta[third], _lPhi[third], _lE[third] * ptNonZ / _lPt[third]);

            mt1 = mtCalc(l0p4, _met, _metPhi);
            //if(mt1 < 50) continue;
          }

          if(debug) cout << "weight of event is " << weight << endl;

          int mvaValueRegion = 0;
          /*
          int mvaValueRegion = -999;
          if(mvaVL > -0.2 && mvaVL < 0.4)
            mvaValueRegion = 0;
          if(mvaVL > 0.4)
            mvaValueRegion = 1;

          if(mvaVL < -0.2)
            mvaValueRegion = 2;
          */

          if(samCategory == dataSample && nLocEle == 0)
            myfile << _runNb << " " << _lumiBlock << " " << _eventNb << endl; //  << FRloc << endl; 

          vector<double> fillVar = {ptCorrV[0].first, ptCorrV[1].first, leptonSelectionAnalysis > 2 ? ptCorrV[2].first : 0., leptonSelectionAnalysis > 3 ? ptCorrV[3].first : 0.,
                                   mt1, double(nJLoc), double(nBLoc), mvaVL,
                                   (leptonSelectionAnalysis == 4 ? flavourCategory4L(nLocEle) : (leptonSelectionAnalysis == 3 ? flavourCategory3L(nLocEle) : flavourCategory2L(nLocEle,_lCharge[ind.at(0)]))), 
                                   leptonSelectionAnalysis == 4 ? 0. : (leptonSelectionAnalysis == 3 ? SRID3L(nJLoc, nBLoc) : SRID2L(nJLoc, nBLoc, mvaValueRegion, _lCharge[ind.at(0)])),
                                   (leptonSelectionAnalysis != 4 ? mll:mll1stpair),ptZ,ptNonZ, (nLocEle == 3?mll:-999.), (nLocEle==2?mll:-999.), (nLocEle==1?mll:-999.), (nLocEle==0? mll: -999.),
                                   _met, minDeltaR, minDeltaRlead, mtHighest, mtLowest, leadingJetPt, trailJetPt, 0., double(_nVertex), mlll,
                                   _lEta[ptCorrV[0].second], _lEta[ptCorrV[1].second], (leptonSelectionAnalysis > 2 ? ptCorrV[2].second : -999.), (leptonSelectionAnalysis > 3 ? ptCorrV[3].second : -999.),
                                   (nLocEle == 3?mt1:-999.), (nLocEle==2?mt1:-999.), (nLocEle==1?mt1:-999.), (nLocEle==0? mt1: -999.),
                                   cosTSt, mll_ss, double(chargeOfLeptons), ll_deltaR, mt2ll_ss};
          vector<TString> fncName = {"ptlead", "sublead", "trail", "pt4th", 
                                     "mtW", "njets", "nbjets", "BDT", 
                                     "flavour", 
                                     "SR",
                                     "mll", "ptZ", "ptNonZ", "mll3e", "mll2e1mu", "mll1e2mu", "mll3mu",
                                     "met", "deltaR", "deltaRlead", "mtLeading", "mtTrailing", "leadJetPt", "trailJetPt", "SRnpCR", "nPV", "mlll",
                                     "etaLead", "etaSubl", "etaTrail", "eta4th", 
                                     "mt_3m", "mt_2m1e", "mt_1m2e", "mt_3e", 
                                     "cosThetaStar", "mll_ss", "chargeOfLeptons", "ll_deltaR", "mt2ll_ss"};
                                   
          for(int dist = 0; dist < fillVar.size(); dist++){
            distribs[dist].vectorHisto[samCategory].Fill(TMath::Min(fillVar.at(dist),figNames[fncName.at(dist)].varMax-0.1),weight);

            /*
            if((samples[sam].getProcessName()) != "data"){
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * leptonWeight(1) / leptonWeight(0));
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * leptonWeight(2) / leptonWeight(0));

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 1, figNames[fncName.at(dist)].varMax-0.1, weight*puWeight(1)/puWeight(0));
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 1, figNames[fncName.at(dist)].varMax-0.1, weight*puWeight(2)/puWeight(0));

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 2, figNames[fncName.at(dist)].varMax-0.1, weight*_lheWeight[8]);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 2,figNames[fncName.at(dist)].varMax-0.1, weight*_lheWeight[4]);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 3, figNames[fncName.at(dist)].varMax-0.1, weight*bTagWeight_udsg(1)/bTagWeight_udsg(0));
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 3, figNames[fncName.at(dist)].varMax-0.1, weight*bTagWeight_udsg(2)/bTagWeight_udsg(0));

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 4, figNames[fncName.at(dist)].varMax-0.1, weight*bTagWeight_c(1)*bTagWeight_b(1)/bTagWeight_c(0)/bTagWeight_b(0));
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 4, figNames[fncName.at(dist)].varMax-0.1, weight*bTagWeight_c(2)*bTagWeight_b(2)/bTagWeight_c(0)/bTagWeight_b(0));

            }
            else if(samCategory == nonPromptSample){
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * 1.3);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * 0.7);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 1, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 1, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 2, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 2,figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 3, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 3, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 4, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 4, figNames[fncName.at(dist)].varMax-0.1, weight);

            }
            */
          }
      }

      cout << endl;
      cout << "Total number of events: " << distribs[0].vectorHisto[sam].Integral() << endl;
      if(leptonSelectionAnalysis != 4)
        cout << "Total number of events in non prompt category: " << distribs[0].vectorHisto[nonPromptSample].Integral() << endl;
      cout << endl;
  }

  fileDummy->cd();
  signalTree->Write();
  bkgTree->Write();

  fileDummy->Close();

  TLegend* mtleg = new TLegend(0.15,0.89,0.95,0.72); 
  mtleg->SetNColumns(5);
  mtleg->SetFillColor(0);
  mtleg->SetFillStyle(0);
  mtleg->SetBorderSize(0);
  mtleg->SetTextFont(42);
  
  /*
  mtleg->AddEntry(&distribs[0].vectorHisto[dataSample],"Nonprompt DD","lep"); //data
  mtleg->AddEntry(&distribs[0].vectorHisto[1],"Nonprompt MC","f"); //data
  */

  mtleg->AddEntry(&distribs[0].vectorHisto[dataSample],"Data","lep"); //data
  int count = 0;
  for (std::vector<std::string>::iterator it = samplesOrderNames.begin()+1; it != samplesOrderNames.end(); it++) {
        count++;
        if(samplesOrderNames.at(count) == "ttH") continue;
        //if(samplesOrderNames.at(count) == "ZZ") continue;

        cout << "count and sample name: " << count << " " << *it << " " << samplesOrder.at(count) << endl;

        if(leptonSelectionAnalysis == 3 && samplesOrderNames.at(count) == "chargeMisID") continue;
        if(samplesOrderNames.at(count) != "nonpromptData") // && samplesOrderNames.at(count) != "Zgamma")
            mtleg->AddEntry(&distribs[0].vectorHisto[samplesOrder.at(count)],(*it).c_str(),"f");
        else if(samplesOrderNames.at(count) == "nonpromptData")
            mtleg->AddEntry(&distribs[0].vectorHisto[samplesOrder.at(count)],"Nonprompt","f");
        //else if(samplesOrderNames.at(count) == "Zgamma")
        //    mtleg->AddEntry(&distribs[0].vectorHisto[samplesOrder.at(count)],"Z#gamma","f");
  }

  count = 0;
  for (std::vector<std::string>::iterator it = samplesOrderNames.begin(); it != samplesOrderNames.end(); it++) {
      // here let's draw the yields output for each flavour category for each component
      std::cout << samplesOrderNames.at(count).c_str() << " ";
      for(int binN = 1; binN < 5; binN++){
        double outputValue = 0;
        for(int begOfSam = samplesOrder.at(count); begOfSam < samplesOrder.at(count+1); begOfSam++)
            outputValue += distribs[8].vectorHisto[begOfSam].GetBinContent(binN) ;
        cout << outputValue << " "; 
      }
      cout << endl;
      count++;
  }

  gSystem->Exec("rm plotsForSave/*");
  double scale_num = 1.6;
  
  TCanvas* plot[nVars];
  for(int i = 0; i < nVars; i++){
      plot[i] = new TCanvas(Form("plot_%d", i),"",500,450);
  }

//  general list, don't erase
//  vector<TString> namesForSaveFiles = {"ptlead", "sublead", "trail", "njets", "nbjets", "SR", "flavour", "BDT", "mll", "ptZ", "ptNonZ", "mtW", "mll3e", "mll2e1mu", "mll1e2mu", "mll3mu", "met", "deltaR", "mtLeading", "mtTrailing", "leadJetPt", "trailJetPt", "SRnpCR", "pteleconv", "etaeleconv", "nPV", "pt4th", "mlll", "etaLead", "etaSubl", "ptleadForw", "subleadForw", "etaLeadForw", "etaSublForw", "ptleadMiddle", "subleadMiddle", "etaLeadMiddle", "etaSublMiddle", "mt_3m", "mt_2m1e", "mt_1m2e", "mt_3e"};

  // list for WZ
  //vector<TString> listToPrint = {"ptlead", "sublead", "trail", "njets", "nbjets", "flavour", "mll", "ptZ", "ptNonZ", "mtW", "mll3e", "mll2e1mu", "mll1e2mu", "mll3mu", "met", "nPV", "mt_3m", "mt_2m1e", "mt_1m2e", "mt_3e", "cosThetaStar"};
  // list for Zgamma
  //vector<TString> listToPrint = {"ptlead", "sublead", "trail", "njets", "nbjets", "flavour", "met", "nPV", "mlll"};
  // list for ttbar CR
  //vector<TString> listToPrint = {"ptlead", "sublead", "trail", "njets", "nbjets", "flavour", "met", "nPV"};
  // list for DY nonprompt CR
  //vector<TString> listToPrint = {"ptlead", "sublead", "trail", "njets", "nbjets", "flavour", "met", "nPV", "mll", "ptZ", "ptNonZ", "mtW", "mll3e", "mll2e1mu", "mll1e2mu", "mll3mu", "mt_3m", "mt_2m1e", "mt_1m2e", "mt_3e"};
  // list for ss2l ttW
  //vector<TString> listToPrint = {"ptlead", "sublead", "njets", "nbjets", "flavour", "met", "nPV", "deltaR", "deltaRlead", "mtLeading", "mtTrailing", "leadJetPt", "trailJetPt", "etaLead", "etaSubl", "mll_ss", "chargeOfLeptons", "ll_deltaR", "mt2ll_ss", "SR", "BDT"}; // 
  // list for 4L ZZ CR
  //vector<TString> listToPrint = {"ptlead", "sublead", "trail", "pt4th", "njets", "nbjets", "flavour", "met", "nPV", "mll", "ptZ", "etaLead", "etaSubl", "etaTrail", "eta4th"};
  // list for ttZ 3l signal seleciton
  vector<TString> listToPrint = {"ptlead", "sublead", "trail", "njets", "nbjets", "flavour", "mll", "ptZ", "ptNonZ", "SR", "met", "cosThetaStar"};
  //vector<TString> listToPrint = {"SR"};
  for(int varPlot = 0; varPlot < listToPrint.size(); varPlot++){
    plot[varPlot]->cd();
    showHist(plot[varPlot],distribs[figNames[listToPrint.at(varPlot)].index],"",figNames[listToPrint.at(varPlot)].fancyName,"Events", scale_num, mtleg, false, false, dataLumi);
    plot[varPlot]->SaveAs("plotsForSave/" + listToPrint.at(varPlot) + ".pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + listToPrint.at(varPlot) + ".png");
    plot[varPlot]->SaveAs("plotsForSave/" + listToPrint.at(varPlot) + ".root");
    plot[varPlot]->cd();
    showHist(plot[varPlot],distribs[figNames[listToPrint.at(varPlot)].index],"",figNames[listToPrint.at(varPlot)].fancyName,"Events", scale_num, mtleg, true, false, dataLumi);
    plot[varPlot]->SaveAs("plotsForSave/" + listToPrint.at(varPlot) + "Log.pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + listToPrint.at(varPlot) + "Log.png");
    plot[varPlot]->SaveAs("plotsForSave/" + listToPrint.at(varPlot) + "Log.root");
  }
  
  //fillDatacards(distribs[indexSR], samplesOrderNames, samplesOrder);
}

int main(int argc, const char **argv)
{
    int rargc = 1; char *rargv[1] = {""};
    cout << "Number of input arguments " << argc << endl;
    for(int i = 0; i < argc; ++i){
        cout << "Argument " << i << " " << argv[i] << endl;
    }
    TApplication *rootapp = new TApplication("example", &rargc, rargv);
    treeReader reader;
    if(argc == 1){
        std::cerr << "please specify input file with samples from data/samples directory" << std::endl;
        return 1;
    }
    if(argc == 2) reader.Analyze(std::string(argv[1]));
    if(argc > 2){
        if(argc == 3) {
            if(string(argv[2]) == "debug"){
                std::cerr << "please specify sample to debug" << std::endl;
                return 1;
            }
            if(string(argv[2]) == "runOnOneProcess"){
                std::cerr << "please specify sample to run on" << std::endl;
                return 1;
            }
        }
        if(argc == 4) reader.Analyze(std::string(argv[1]), std::string(argv[2]), std::string(argv[3]));
        if(argc == 5) reader.Analyze(std::string(argv[1]), std::string(argv[2]), std::string(argv[3]), atol(argv[4]));
    }
    rootapp->Run();
    return 0;
}
