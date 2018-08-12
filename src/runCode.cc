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

void treeReader::Analyze(const string& fileToAnalyse, const std::string option, const std::string selection, const string& sampleToDebug, long evNb){

  debug = (option == "debug" ? true : false);
  leptonSelection = leptonSelectionAnalysis;
  initListsToPrint(selection);
  //Set CMS plotting style
  setTDRStyle();
  gROOT->SetBatch(kTRUE);
  //read samples and cross sections from txt file
  cout << "reading sample file...." << endl;
  readSamples(fileToAnalyse);
  cout << "finished with reading"<< endl;
  setTDRStyle(); 

  std::vector<std::string> namesOfSamples = treeReader::getNamesOfTheSample();
  cout << "initiating histos...." << endl;
  initdistribs(namesOfSamples);
  cout << "finished with initiating of histos"<< endl;
  setLabelsForHistos();

  /*
  if(leptonSelection == 2){

    // for reading values from xml files
    addVariablesToBDT(is2017);

    // for adding variables to the trees
    //addBranchToBDTTreeVariables();
  }
  */

  //std::ofstream myfile;
  //myfile.open("myevents.txt");

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample();

      Color_t color = assignColor(samples[sam].getProcessName());
      setStackColors(color, sam);

      //if(!((samples[sam].getFileName()).find("ST_tW_") != std::string::npos)) continue;

      if((option == "runOnOneProcess" || debug) && (samples[sam].getProcessName()) != sampleToDebug) continue;
      if(samples[sam].getProcessName() == "nonpromptData"){
          cout << "Total number of events: " << distribs[0].vectorHisto[sam].Integral() << endl;
          continue;
      }

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
          if(debug && (_eventNb != evNb && evNb != -999)) continue;
          if(debug) cout << "######################### New Event ###########################################" << endl;
          if(debug) cout << "event " << _eventNb << " was found" << endl;
          
          // trigger and met filters
          if(debug) cout << "trigger decision: " << _passTrigger_e << " " << _passTrigger_m << " " << _passTrigger_ee << " " << _passTrigger_em << " " << _passTrigger_mm << " " << _passTrigger_eee << " " << _passTrigger_eem << " " << _passTrigger_emm << " " << _passTrigger_mmm << endl;
          if(!(_passTrigger_e || _passTrigger_m || _passTrigger_ee || _passTrigger_em || _passTrigger_mm || _passTrigger_eee || _passTrigger_eem || _passTrigger_emm || _passTrigger_mmm)) continue;
          if(debug) cout << "met filers flag: " << _passMETFilters << endl;
          if(!_passMETFilters) continue;
          
          //if(it > 10000) break;
          if(it > nEntries / 100) break;

          std::vector<unsigned> indTight, indFake, indOf2LonZ;
          //select leptons relative to the analysis
          leptonSelection = 3;
          const unsigned lCount = selectLep(indTight, leptonSelection);
          const unsigned lCountFake = selectFakeLep(indFake, leptonSelection);

          std::vector<unsigned> indTight4L;
          const unsigned lCount4L = selectLep(indTight4L, 4);

          // discard heavy flavour resonances
          if(debug) cout << "inva mass of any fake pair is below 12 GeV: " << invMassOfAny2Lbelow12GeV(indFake) << endl;
          if(invMassOfAny2Lbelow12GeV(indFake)) continue; 

          if(debug) cout << "number of ttZ3L tight and fo leptons: " << lCount << " " << lCountFake << endl;
          if(debug) cout << "number of ttZ4L tight leptons: " << lCount4L << endl;

          // selection of category for the event
          // 2L: possible contribution from TT, TF and FF; TTF is vetoed
          // 3L: TTT, TTF, TFF, FFF; for TTTF should consider if event pass 4L TTTT criteria
          // 4L: TTTT only combination is possible
          int samCategory = sam;
          std::vector<unsigned> ind;
          
          if(lCount4L == 4){
              leptonSelection = 4;
              samCategory = sam;
              ind = indTight4L;
          }
          else if (lCount == 3){
              if(lCountFake != 3) continue;
              samCategory = sam;
              ind = indTight;
          }
          else if (lCount < 3) {
              if(lCountFake != 3) continue;
              samCategory = nonPromptSample;
              ind = indFake;
          }
          else continue;
          /*
          //used for main analysis
          // remove additional FO lepton from ss2l selection
          //if(leptonSelection == 2 && lCountFake != leptonSelection) continue;
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
          if(lCountFake != leptonSelection) continue; // here used to be less than, but for the moment I switch to veto FO object
          if(lCount == leptonSelection) {
              samCategory = sam;
              ind = indTight;
          }
          if(lCount < leptonSelection){
              if(leptonSelection == 4) continue;
              if(leptonSelection == 2 && samCategory == CMIDSample) continue;
              samCategory = nonPromptSample;
              ind = indFake;
          }

          if(leptonSelection == 2 && samCategory != CMIDSample)
            if(_lCharge[ind.at(0)] * _lCharge[ind.at(1)] < 0) continue;
          */

          // consider only prompt leptons from the MC, all nonprompt should be taken into account by DD estimation
          bool allLeptonsArePrompt = true;
          
          if((samples[sam].getProcessName()) != "data" && (samples[sam].getProcessName()) != "nonpromptData" && (samples[sam].getProcessName()) != "chargeMisIDData")
            allLeptonsArePrompt = promptLeptons(ind);

          if(debug) cout << "all leptons are prompt ? " << allLeptonsArePrompt << endl;
          
          if((samples[sam].getProcessName()) == "chargeMisID" && !allLeptonsArePrompt) continue;
          if((samples[sam].getProcessName()) == "Nonprompt" && allLeptonsArePrompt) continue; // works just for MC

          if(((samples[sam].getProcessName()) == "ttW" || (samples[sam].getProcessName()) == "ttH" || (samples[sam].getProcessName()) == "ttZ" || (samples[sam].getProcessName()) == "ttX" || (samples[sam].getProcessName()) == "WZ" || (samples[sam].getProcessName()) == "Xgamma"  || (samples[sam].getProcessName()) == "ZZ" || (samples[sam].getProcessName()) == "rare") && !allLeptonsArePrompt) continue;

          int nLocEle = getElectronNumber(ind);
          //if(nLocEle != 3) continue;

          // lepton pt criteria
          if(leptonSelection == 4)
            if(!passPtCuts4L(ind)) continue;
          if(leptonSelection == 3)
            if(!passPtCuts3L(ind)) continue;
          if(leptonSelection == 2)
            if(!passPtCuts2L(ind)) continue;

          // select here jets, bjets, delta from M of Z boson, HT
          std::vector<unsigned> indJets;
          std::vector<unsigned> indJetsJECUp;
          std::vector<unsigned> indJetsJECDown;
          std::vector<unsigned> indJetsNotB;

          unsigned third = -9999;
          double mll = 99999;
          double mlll = 99999;
          double ptZ = 999999;
          double ptNonZ = 999999;

          nJLoc = nJets(0, true, indJets, is2017);
          int nJLocDown = nJets(1, true, indJetsJECDown, is2017);
          int nJLocUp = nJets(2, true, indJetsJECUp, is2017);
          nJLocNotB = nJetsNotB(0, true, indJetsNotB, 2, is2017);
          nBLoc = nBJets(0, true, true, 1, is2017);
          int nBLocDown = nBJets(1, true, true, 1, is2017);
          int nBLocUp = nBJets(2, true, true, 1, is2017);

          TLorentzVector Zboson, lnegative;
          double dMZ = deltaMZ(ind, third, mll, ptZ, ptNonZ, mlll, indOf2LonZ, Zboson, lnegative);
          double mll1stpair, mll2ndpair;
          double cosTSt = -999;

          if(debug) cout << "number of jets/bjets/dMZ: " << nJLoc << " " << nBLoc << " " << dMZ << endl;

          HTLoc = HTCalc(indJets);
          double HTLocJECUp = HTCalc(indJetsJECUp);
          double HTLocJECDown  = HTCalc(indJetsJECDown);
          
          double mt1 = 9999;
          if(leptonSelection == 4){
            // ZZ control region
            if(dMZ > 20) continue; 
            mll1stpair = mll;
            cosTSt = cosThetaStar(Zboson, lnegative);
            if(selection == "ttZ4L" && !passTTZ4LSelection(ind, indOf2LonZ, nJLoc)) continue;
            if(selection == "ZZ" && !passZZCRSelection(ind, indOf2LonZ)) continue;
            if(selection == "ttZ" && !(passTTZ4LSelection(ind, indOf2LonZ, nJLoc) || passZZCRSelection(ind, indOf2LonZ))) continue;
          }

          if(leptonSelection == 3){
            
            if(selection == "ttZ" && !(passTTZSelection(nJLoc, dMZ) || passWZCRSelection(nBLoc, dMZ) || passttbarCRSelection(nBLoc, dMZ))) continue;
            if(selection == "tZq" &&!passttbarCRintZqSelection(nJLoc, nBLoc, dMZ)) continue;
            if(selection == "ttZ3L" && !passTTZSelection(nJLoc, dMZ)) continue;
            if(selection == "ttZ3Lclean" && !passTTZCleanSelection(nJLoc, nBLoc, dMZ)) continue;
            if(selection == "WZ" && !passWZCRSelection(nBLoc, dMZ)) continue;
            if(selection == "DY" && !passDYCRSelection(dMZ, ptNonZ, third, _met, _metPhi, nJLoc, nBLoc)) continue;

            // if Z boson is reconstructed then we can calculate mt for 3 rd lepton and cos theta star
            if(dMZ < 10){
                TLorentzVector l0p4;
                l0p4.SetPtEtaPhiE(ptNonZ, _lEta[third], _lPhi[third], _lE[third] * ptNonZ / _lPt[third]);
                mt1 = mtCalc(l0p4, _met, _metPhi);
                cosTSt = cosThetaStar(Zboson, lnegative);
            }
            if(selection == "ttbar" &&!passttbarCRSelection(nBLoc, dMZ)) continue;
            if(selection == "Xgamma" && !passZGCRSelection(mlll, dMZ)) continue;

          }
          
          if(leptonSelection == 2){
            if(selection == "ttW" && !pass2Lpreselection(nJLoc, nBLoc, ind, _met, nLocEle)) continue;
            if(selection == "ttWclean" && !pass2Lcleanpreselection(nJLoc, nBLoc, ind, _met, nLocEle)) continue;
          }
          

          double mvaVL = 0;
          double mvaVLJECUp = 0;
          double mvaVLJECDown = 0;

          if(leptonSelection == 2){

            TLorentzVector l0p4, l1p4, metTL;

            l0p4.SetPtEtaPhiE(ptCorrV[0].first, _lEta[ind.at(0)], _lPhi[ind.at(0)], _lE[ind.at(0)] * ptCorrV[0].first / _lPt[ind.at(0)]);
            l1p4.SetPtEtaPhiE(ptCorrV[1].first, _lEta[ind.at(1)], _lPhi[ind.at(1)], _lE[ind.at(1)] * ptCorrV[1].first / _lPt[ind.at(1)]);
            metTL.SetPtEtaPhiM(_met, 0, _metPhi, 0);

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

            /*
            _weightEventInTree = weight;
            if((samples[sam].getProcessName()) == "ttW" && samCategory != nonPromptSample)
              signalTree->Fill();
            //else if((samples[sam].getProcessName()) == "Nonprompt" || (samples[sam].getProcessName()) == "chargeMisID" || (samples[sam].getProcessName()) == "ttZ" || (samples[sam].getProcessName()) == "ttH" || (samples[sam].getProcessName()) == "WZ" || (samples[sam].getProcessName()) == "nonpromptData")
            else if(((samples[sam].getProcessName() == "Nonprompt" || samples[sam].getProcessName() == "chargeMisID") && samCategory != nonPromptSample) || ((samples[sam].getProcessName() == "data") && samCategory == nonPromptSample))
              bkgTree->Fill();
              */
            
            vector<Float_t> varForBDT = {(Float_t)nJLoc, (Float_t)nBLoc, (Float_t)HTLoc, (Float_t)MET, (Float_t)minDeltaRlead, (Float_t)minDeltaR, (Float_t)l0p4.DeltaR(l1p4), (Float_t)mtHighest, (Float_t)mtLowest, (Float_t)leadpt, (Float_t)trailpt, (Float_t) leadingJetPt, (Float_t) trailJetPt, (Float_t)leadeta, (Float_t)traileta, (Float_t)chargeOfLeptons, (Float_t)(l0p4+l1p4).M(), (Float_t)mt2ll_ss};
            fillBDTvariables(varForBDT);
            mvaVL = readerTTWcsttbar->EvaluateMVA( "BDTG method");

            vector<Float_t> varForBDTJECUp = {(Float_t)nJLocUp, (Float_t)nBLocUp, (Float_t)HTLocJECUp, (Float_t)MET, (Float_t)minDeltaRlead, (Float_t)minDeltaR, (Float_t)l0p4.DeltaR(l1p4), (Float_t)mtHighest, (Float_t)mtLowest, (Float_t)leadpt, (Float_t)trailpt, (Float_t) leadingJetPt, (Float_t) trailJetPt, (Float_t)leadeta, (Float_t)traileta, (Float_t)chargeOfLeptons, (Float_t)(l0p4+l1p4).M(), (Float_t)mt2ll_ss};
            fillBDTvariables(varForBDTJECUp);
            mvaVLJECUp = readerTTWcsttbar->EvaluateMVA( "BDTG method");

            vector<Float_t> varForBDTJECDown = {(Float_t)nJLocDown, (Float_t)nBLocDown, (Float_t)HTLocJECDown, (Float_t)MET, (Float_t)minDeltaRlead, (Float_t)minDeltaR, (Float_t)l0p4.DeltaR(l1p4), (Float_t)mtHighest, (Float_t)mtLowest, (Float_t)leadpt, (Float_t)trailpt, (Float_t) leadingJetPt, (Float_t) trailJetPt, (Float_t)leadeta, (Float_t)traileta, (Float_t)chargeOfLeptons, (Float_t)(l0p4+l1p4).M(), (Float_t)mt2ll_ss};
            fillBDTvariables(varForBDTJECDown);
            mvaVLJECDown = readerTTWcsttbar->EvaluateMVA( "BDTG method");

          }

          // weight estimation for event
          //auto start = std::chrono::high_resolution_clock::now();
          if((samples[sam].getProcessName()) != "data" && (samples[sam].getProcessName()) != "nonpromptData" && (samples[sam].getProcessName()) != "chargeMisIDData")
            weight *= sfWeight();
          //auto finish = std::chrono::high_resolution_clock::now();
          //std::chrono::duration<double> elapsed = finish - start;
          //std::cout << "Elapsed time: " << elapsed.count() << std::endl;

          if(samples[sam].getProcessName() == "data" && samCategory == nonPromptSample && leptonSelection != 4)
            weight *= fakeRateWeight();
          if(samCategory == CMIDSample)
            weight *= CMIDRateWeight();

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

          if(debug) cout << "lepton selection is " << leptonSelection << " total SR: " << SRIDTTZ(ind, indOf2LonZ, nJLoc, nBLoc, dMZ) << endl; 

          vector<double> fillVar = {ptCorrV[0].first, ptCorrV[1].first, leptonSelection > 2 ? ptCorrV[2].first : 0., leptonSelection > 3 ? ptCorrV[3].first : 0.,
                                   mt1, double(nJLoc), double(nBLoc), (_lCharge[ind.at(0)] == 1 ?  mvaVL : -999),
                                   // currently here we will have ttZ3L and ttZ4L categories
                                   //(leptonSelection == 4 ? flavourCategory4L(nLocEle) : (leptonSelection == 3 ? flavourCategory3L(nLocEle) : flavourCategory2L(nLocEle,_lCharge[ind.at(0)]))), 
                                   //(leptonSelection == 4 ? SRID4L(nBLoc) : (leptonSelection == 3 ? SRID3L(nJLoc, nBLoc) : SRID2L(mvaVL, _lCharge[ind.at(0)]))),
                                   (leptonSelection == 3 && passTTZSelection(nJLoc, dMZ) ? SRID3L(nJLoc, nBLoc) : -999),
                                   (leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLoc)? SRID4L(nBLoc) : -999),
                                   (leptonSelection != 4 ? mll:mll1stpair),ptZ,ptNonZ, (nLocEle == 3?mll:-999.), (nLocEle==2?mll:-999.), (nLocEle==1?mll:-999.), (nLocEle==0? mll: -999.),
                                   _met, minDeltaR, minDeltaRlead, mtHighest, mtLowest, leadingJetPt, trailJetPt, 0., double(_nVertex), mlll,
                                   _lEta[ptCorrV[0].second], _lEta[ptCorrV[1].second], (leptonSelection > 2 ? _lEta[ptCorrV[2].second] : -999.), (leptonSelection > 3 ? _lEta[ptCorrV[3].second] : -999.),
                                   (nLocEle == 3?mt1:-999.), (nLocEle==2?mt1:-999.), (nLocEle==1?mt1:-999.), (nLocEle==0? mt1: -999.),
                                   cosTSt, mll_ss, double(chargeOfLeptons), ll_deltaR, mt2ll_ss,
                                   (_lCharge[ind.at(0)] == -1 ?  mvaVL : -999), HTLoc,
                                   SRIDTTZ(ind, indOf2LonZ, nJLoc, nBLoc, dMZ), SRIDWZCR(nJLoc, nBLoc, dMZ), SRIDZZCR(ind, indOf2LonZ, nJLoc, nBLoc), SRIDTTCR(nJLoc, nBLoc, dMZ)
                                   };

          vector<double> fillVarJecUp = {ptCorrV[0].first, ptCorrV[1].first, leptonSelection > 2 ? ptCorrV[2].first : 0., leptonSelection > 3 ? ptCorrV[3].first : 0.,
                                   mt1, double(nJLocUp), double(nBLocUp), (_lCharge[ind.at(0)] == 1 ?  mvaVLJECUp : -999),
                                   //(leptonSelection == 4 ? flavourCategory4L(nLocEle) : (leptonSelection == 3 ? flavourCategory3L(nLocEle) : flavourCategory2L(nLocEle,_lCharge[ind.at(0)]))), 
                                   //(leptonSelection == 4 ? SRID4L(nBLocUp) : (leptonSelection == 3 ? SRID3L(nJLocUp, nBLocUp) : SRID2L(mvaVLJECUp, _lCharge[ind.at(0)]))),
                                   (leptonSelection == 3 && passTTZSelection(nJLocUp, dMZ) ? SRID3L(nJLocUp, nBLocUp) : -999),
                                   (leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLocUp)? SRID4L(nBLocUp) : -999),
                                   (leptonSelection != 4 ? mll:mll1stpair),ptZ,ptNonZ, (nLocEle == 3?mll:-999.), (nLocEle==2?mll:-999.), (nLocEle==1?mll:-999.), (nLocEle==0? mll: -999.),
                                   _met, minDeltaR, minDeltaRlead, mtHighest, mtLowest, leadingJetPt, trailJetPt, 0., double(_nVertex), mlll,
                                   _lEta[ptCorrV[0].second], _lEta[ptCorrV[1].second], (leptonSelection > 2 ? _lEta[ptCorrV[2].second] : -999.), (leptonSelection > 3 ? _lEta[ptCorrV[3].second] : -999.),
                                   (nLocEle == 3?mt1:-999.), (nLocEle==2?mt1:-999.), (nLocEle==1?mt1:-999.), (nLocEle==0? mt1: -999.),
                                   cosTSt, mll_ss, double(chargeOfLeptons), ll_deltaR, mt2ll_ss,
                                   (_lCharge[ind.at(0)] == -1 ?  mvaVLJECUp : -999), HTLocJECUp,
                                   SRIDTTZ(ind, indOf2LonZ, nJLocUp, nBLocUp, dMZ), SRIDWZCR(nJLocUp, nBLocUp, dMZ), SRIDZZCR(ind, indOf2LonZ, nJLocUp, nBLocUp), SRIDTTCR(nJLocUp, nBLocUp, dMZ)
                                   };

          vector<double> fillVarJecDw = {ptCorrV[0].first, ptCorrV[1].first, leptonSelection > 2 ? ptCorrV[2].first : 0., leptonSelection > 3 ? ptCorrV[3].first : 0.,
                                   mt1, double(nJLocDown), double(nBLocDown), (_lCharge[ind.at(0)] == 1 ?  mvaVLJECDown : -999),
                                   //(leptonSelection == 4 ? flavourCategory4L(nLocEle) : (leptonSelection == 3 ? flavourCategory3L(nLocEle) : flavourCategory2L(nLocEle,_lCharge[ind.at(0)]))), 
                                   //(leptonSelection == 4 ? SRID4L(nBLocDown) : (leptonSelection == 3 ? SRID3L(nJLocDown, nBLocDown) : SRID2L(mvaVLJECDown, _lCharge[ind.at(0)]))),
                                   (leptonSelection == 3 && passTTZSelection(nJLocDown, dMZ) ? SRID3L(nJLocDown, nBLocDown) : -999),
                                   (leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLocDown)? SRID4L(nBLocDown) : -999),
                                   (leptonSelection != 4 ? mll:mll1stpair),ptZ,ptNonZ, (nLocEle == 3?mll:-999.), (nLocEle==2?mll:-999.), (nLocEle==1?mll:-999.), (nLocEle==0? mll: -999.),
                                   _met, minDeltaR, minDeltaRlead, mtHighest, mtLowest, leadingJetPt, trailJetPt, 0., double(_nVertex), mlll,
                                   _lEta[ptCorrV[0].second], _lEta[ptCorrV[1].second], (leptonSelection > 2 ? _lEta[ptCorrV[2].second] : -999.), (leptonSelection > 3 ? _lEta[ptCorrV[3].second] : -999.),
                                   (nLocEle == 3?mt1:-999.), (nLocEle==2?mt1:-999.), (nLocEle==1?mt1:-999.), (nLocEle==0? mt1: -999.),
                                   cosTSt, mll_ss, double(chargeOfLeptons), ll_deltaR, mt2ll_ss,
                                   (_lCharge[ind.at(0)] == -1 ?  mvaVLJECDown : -999), HTLocJECDown,
                                   SRIDTTZ(ind, indOf2LonZ, nJLocDown, nBLocDown, dMZ), SRIDWZCR(nJLocDown, nBLocDown, dMZ), SRIDZZCR(ind, indOf2LonZ, nJLocDown, nBLocDown), SRIDTTCR(nJLocDown, nBLocDown, dMZ)
                                   };
          vector<TString> fncName = {"ptlead", "sublead", "trail", "pt4th", 
                                     "mtW", "njets", "nbjets", "BDTpp", 
                                     //"flavour", 
                                     //"SR",
                                     "SR3L",
                                     "SR4L",
                                     "mll", "ptZ", "ptNonZ", "mll3e", "mll2e1mu", "mll1e2mu", "mll3mu",
                                     "met", "deltaR", "deltaRlead", "mtLeading", "mtTrailing", "leadJetPt", "trailJetPt", "SRnpCR", "nPV", "mlll",
                                     "etaLead", "etaSubl", "etaTrail", "eta4th", 
                                     "mt_3m", "mt_2m1e", "mt_1m2e", "mt_3e", 
                                     "cosThetaStar", "mll_ss", "chargeOfLeptons", "ll_deltaR", "mt2ll_ss", "BDTmm", "HT",
                                     "SRallTTZ", "SRWZCR", "SRZZCR", "SRTTCR"};
                                   
          if(debug) cout << "lep sf: " << leptonWeight(0) << " " << leptonWeight(1) << " " << leptonWeight(2) << endl;
          for(int dist = 0; dist < fillVar.size(); dist++){
            distribs[dist].vectorHisto[samCategory].Fill(TMath::Min(fillVar.at(dist),figNames[fncName.at(dist)].varMax-0.1),weight);

            if((samples[sam].getProcessName()) != "data" && (samples[sam].getProcessName()) != "chargeMisIDData"){

                //start = std::chrono::high_resolution_clock::now();
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

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 5, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 5, figNames[fncName.at(dist)].varMax-0.1, weight);

                //finish = std::chrono::high_resolution_clock::now();
                //elapsed = finish - start;
                //std::cout << "Elapsed time for bullshit fill: " << elapsed.count() << std::endl;

                //start = std::chrono::high_resolution_clock::now();

                double lepSF = leptonWeight(0);double lepSFUp = leptonWeight(1);double lepSFDown = leptonWeight(2);

                double puW = puWeight(0);double puWUp = puWeight(1);double puWDown = puWeight(2);

                double btagL = bTagWeight_udsg(0); double btagLUp = bTagWeight_udsg(1); double btagLDown = bTagWeight_udsg(2);
                double btagC = bTagWeight_c(0); double btagCUp = bTagWeight_c(1); double btagCDown = bTagWeight_c(2);
                double btagB = bTagWeight_b(0); double btagBUp = bTagWeight_b(1); double btagBDown = bTagWeight_b(2);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * lepSFUp / lepSF);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * lepSFDown / lepSF);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 1, figNames[fncName.at(dist)].varMax-0.1, weight*puWUp/puW);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 1, figNames[fncName.at(dist)].varMax-0.1, weight*puWDown/puW);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 2, figNames[fncName.at(dist)].varMax-0.1, weight*btagLUp/btagL);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 2, figNames[fncName.at(dist)].varMax-0.1, weight*btagLDown/btagL);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 3, figNames[fncName.at(dist)].varMax-0.1, weight*btagCUp*btagBUp/btagC/btagB);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 3, figNames[fncName.at(dist)].varMax-0.1, weight*btagCDown*btagBDown/btagC/btagB);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVarJecUp.at(dist), 4, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVarJecDw.at(dist), 4, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 5, figNames[fncName.at(dist)].varMax-0.1, weight*_lheWeight[8]*sumSimulatedEventWeights/sumSimulatedEventWeightsScaleUp);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 5,figNames[fncName.at(dist)].varMax-0.1, weight*_lheWeight[4]*sumSimulatedEventWeights/sumSimulatedEventWeightsScaleDown);

                /*
                // no need to run
                if(dist == indexSR3L || dist == indexSR4L || dist == indexSRTTZ){
                    for(int varPDF = 0; varPDF < 100; varPDF++){
                        distribs[dist].vectorHistoPDF[samCategory].var[varPDF].Fill(fillVar.at(dist), weight*_lheWeight[9+varPDF]);
                    }
                }
                */

                /*
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight);
                */

                //finish = std::chrono::high_resolution_clock::now();
                //elapsed = finish - start;
                //std::cout << "Elapsed time for needed fill: " << elapsed.count() << std::endl;

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

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 5, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 5, figNames[fncName.at(dist)].varMax-0.1, weight);

                /*
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight);
                */

            }
            else if(samCategory == CMIDSample){
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * 1.2);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * 0.8);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 1, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 1, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 2, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 2,figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 3, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 3, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 4, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 4, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 5, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 5, figNames[fncName.at(dist)].varMax-0.1, weight);

                /*
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight);
                */

            }
          }
      }

      cout << endl;
      cout << "Total number of events: " << distribs[0].vectorHisto[sam].Integral() << endl;
      if(leptonSelection != 4)
        cout << "Total number of events in non prompt category: " << distribs[0].vectorHisto[nonPromptSample].Integral() << endl;
      cout << endl;
  }

  // this should be done to be fully correct in PDF, takes a lot of time, an effect estimated with ttZ sample, unceratinty on acceptance is under 1%, simply assign flat uncertainty of 1 % to all signal and bkg
  /*
  for(unsigned sam = 0; sam < samples.size() + 1; ++sam){
    if(sam == dataSample) continue;
    if(sam == nonPromptSample) continue;
    for(unsigned bin = 1; bin < (unsigned) distribs[indexSRTTZ].vectorHistoUncUp[sam].unc.at(6).GetNbinsX() + 1; ++bin){
        double pdfVarRms = 0.;
        for(unsigned pdf = 0; pdf < 100; ++pdf){
            double variedBin = distribs[indexSRTTZ].vectorHistoPDF[sam].var[pdf].GetBinContent(bin);
            variedBin *= crossSectionRatio[sam][pdf];
            double diff = (  variedBin - distribs[indexSRTTZ].vectorHisto[sam].GetBinContent(bin) );
            pdfVarRms += diff * diff;
        }
        pdfVarRms = sqrt( 0.01 * pdfVarRms );
        //cout << "pdf rms for bin " << bin << " is equal to " << pdfVarRms << endl;
        distribs[indexSRTTZ].vectorHistoUncUp[sam].unc.at(6).SetBinContent(bin, distribs[indexSRTTZ].vectorHisto[sam].GetBinContent(bin) + pdfVarRms);
        distribs[indexSRTTZ].vectorHistoUncDown[sam].unc.at(6).SetBinContent(bin, distribs[indexSRTTZ].vectorHisto[sam].GetBinContent(bin) - pdfVarRms);
    }
  }
  */


  /*
  fileDummy->cd();
  signalTree->Write();
  bkgTree->Write();

  fileDummy->Close();
  */

  TLegend* mtleg = new TLegend(0.15,0.89,0.95,0.72); 
  mtleg->SetNColumns(4);
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

        if(leptonSelection == 3 && samplesOrderNames.at(count) == "chargeMisID") continue;
        if(samplesOrderNames.at(count) != "chargeMisIDData" && samplesOrderNames.at(count) != "nonpromptData" && samplesOrderNames.at(count) != "Xgamma")
            mtleg->AddEntry(&distribs[0].vectorHisto[samplesOrder.at(count)],(*it).c_str(),"f");
        else if(samplesOrderNames.at(count) == "nonpromptData")
            mtleg->AddEntry(&distribs[0].vectorHisto[samplesOrder.at(count)],"Nonprompt","f");
        else if(samplesOrderNames.at(count) == "Xgamma")
            mtleg->AddEntry(&distribs[0].vectorHisto[samplesOrder.at(count)],"X#gamma","f");
        else if(samplesOrderNames.at(count) == "chargeMisIDData")
            mtleg->AddEntry(&distribs[0].vectorHisto[samplesOrder.at(count)],"chargeMisID","f");
  }

  count = 0;
  for (std::vector<std::string>::iterator it = samplesOrderNames.begin(); it != samplesOrderNames.end(); it++) {
      // here let's draw the yields output for each flavour category for each component
      std::cout << samplesOrderNames.at(count).c_str() << " ";
      for(int binN = 1; binN < distribs[8].vectorHisto[0].GetNbinsX()+1; binN++){
        double outputValue = 0;
        for(int begOfSam = samplesOrder.at(count); begOfSam < samplesOrder.at(count+1); begOfSam++)
            outputValue += distribs[8].vectorHisto[begOfSam].GetBinContent(binN) ;
        cout << outputValue << " "; 
      }
      cout << endl;
      count++;
  }

  std::string processToStore = selection;
  gSystem->Exec("rm plotsForSave/" + (TString)(is2017 ? "2017/" : "2016/") + processToStore + "/*.{pdf,png,root}");
  gSystem->Exec("rmdir plotsForSave/" + (TString)(is2017 ? "2017/" : "2016/") + processToStore);
  gSystem->Exec("mkdir plotsForSave/" + (TString)(is2017 ? "2017/" : "2016/") + processToStore);
  double scale_num = 1.6;
  
  TCanvas* plot[nVars];
  for(int i = 0; i < nVars; i++){
      plot[i] = new TCanvas(Form("plot_%d", i),"",500,450);
  }

  std::string crToPrint = selection;

  for(int varPlot = 0; varPlot < listToPrint[crToPrint].size(); varPlot++){
    plot[varPlot]->cd();
    showHist(plot[varPlot],distribs[figNames[listToPrint[crToPrint].at(varPlot)].index],figNames[listToPrint[crToPrint].at(varPlot)], scale_num, mtleg, false, false, is2017);
    plot[varPlot]->SaveAs("plotsForSave/" + (TString)(is2017 ? "2017/" : "2016/") + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + ".pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + (TString)(is2017 ? "2017/" : "2016/") + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + ".png");
    plot[varPlot]->SaveAs("plotsForSave/" + (TString)(is2017 ? "2017/" : "2016/") + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + ".root");
    plot[varPlot]->cd();
    showHist(plot[varPlot],distribs[figNames[listToPrint[crToPrint].at(varPlot)].index],figNames[listToPrint[crToPrint].at(varPlot)], scale_num, mtleg, true, false, is2017);
    plot[varPlot]->SaveAs("plotsForSave/" + (TString)(is2017 ? "2017/" : "2016/") + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + "Log.pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + (TString)(is2017 ? "2017/" : "2016/") + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + "Log.png");
    plot[varPlot]->SaveAs("plotsForSave/" + (TString)(is2017 ? "2017/" : "2016/") + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + "Log.root");
  }
  
  //if(crToPrint == "ttW" || crToPrint == "ttZ3L" || crToPrint == "ttZ4L")
  //    fillDatacards(distribs[indexSR], samplesOrderNames, samplesOrder, is2017);
  if(crToPrint == "ttZ"){
    fillDatacards(distribs[indexSR3L], samplesOrderNames, samplesOrder, "SR3L", is2017); 
    fillDatacards(distribs[indexSR4L], samplesOrderNames, samplesOrder, "SR4L", is2017); 
    fillDatacards(distribs[indexSRTTZ], samplesOrderNames, samplesOrder, "SRallTTZ", is2017); 
    fillDatacards(distribs[indexSRWZCR], samplesOrderNames, samplesOrder, "SRWZCR", is2017); 
    fillDatacards(distribs[indexSRZZCR], samplesOrderNames, samplesOrder, "SRZZCR", is2017); 
    fillDatacards(distribs[indexSRTTCR], samplesOrderNames, samplesOrder, "SRTTCR", is2017); 
  }
}

int main(int argc, const char **argv)
{
    int rargc = 1; char *rargv[1] = {""};
    cout << "Number of input arguments " << argc << endl;
    for(int i = 0; i < argc; ++i){
        cout << "Argument " << i << " " << argv[i] << endl;
    }
    //TApplication *rootapp = new TApplication("example", &rargc, rargv);
    treeReader reader;
    if(argc == 1){
        std::cerr << "please specify input file with samples from data/samples directory" << std::endl;
        return 1;
    }
    if(argc == 2){
        std::cerr << "please option (runFullSelection, runOnOneProcess, debug)" << std::endl;
        return 1;
        //reader.Analyze(std::string(argv[1]));
    }    
    if(argc > 2){
        if(argc == 3) {
            if(string(argv[2]) == "runFullSelection"){
                std::cerr << "please specify control region (ttZ3L, WZ, ttbar, DY, Xgamma), use \'CR:\' before control region" << std::endl;
                return 1;
            }
            else if(string(argv[2]) == "debug"){
                std::cerr << "please specify process to debug" << std::endl;
                return 1;
            }
            else if(string(argv[2]) == "runOnOneProcess"){
                std::cerr << "please specify process to run on" << std::endl;
                return 1;
            }
        }
        if(argc == 4){
            if(string(argv[2]) == "runFullSelection"){
                if(string(argv[3]).find("selection:") != std::string::npos){
                    std::string selection = string(argv[3]);
                    selection.erase (selection.begin(), selection.begin()+10);
                    std::cout << "output folder is set to: " << selection<< std::endl;
                    reader.Analyze(std::string(argv[1]), std::string(argv[2]), selection); 
                }
            }
            else if(string(argv[2]) == "runBDTtraining"){
                if(string(argv[3]).find("selection:") != std::string::npos){
                    std::string selection = string(argv[3]);
                    selection.erase (selection.begin(), selection.begin()+10);
                    std::cout << "output folder is set to: " << selection<< std::endl;
                    reader.Analyze(std::string(argv[1]), std::string(argv[2]), selection); 
                }
            }
            else if(string(argv[2]) == "runOnOneProcess"){
                std::cerr << "please specify process to run on" << std::endl;
                return 1;
            }
            else{
                std::cerr << "option is unknown, please specify option (runFullSelection, runBDTtraining, runOnOneProcess, debug)" << std::endl;
                return 1;
            }
        }
        if(argc == 5){ 
            if(string(argv[2]) == "runOnOneProcess"){
                if(string(argv[3]).find("selection:") != std::string::npos){
                    std::string selection = string(argv[3]);
                    selection.erase (selection.begin(), selection.begin()+10);
                    std::cout << "output folder is set to: " << selection<< std::endl;
                    reader.Analyze(std::string(argv[1]), std::string(argv[2]), selection, std::string(argv[4])); 
                }
            }
            else if(string(argv[2]) == "debug"){
                if(string(argv[3]).find("selection:") != std::string::npos){
                    std::string selection = string(argv[3]);
                    selection.erase (selection.begin(), selection.begin()+10);
                    std::cout << "output folder is set to: " << selection<< std::endl;
                    reader.Analyze(std::string(argv[1]), std::string(argv[2]), selection, std::string(argv[4])); 
                }
            }
            else{
                std::cerr << "option is unknown, please specify option (runFullSelection, runBDTtraining, runOnOneProcess, debug)" << std::endl;
                return 1;
            }
        }
        if(argc == 6){ 
            if(string(argv[2]) == "debug"){
                if(string(argv[3]).find("selection:") != std::string::npos){
                    std::string selection = string(argv[3]);
                    selection.erase (selection.begin(), selection.begin()+10);
                    std::cout << "output folder is set to: " << selection<< std::endl;
                    reader.Analyze(std::string(argv[1]), std::string(argv[2]), selection, std::string(argv[4]), atol(argv[5]));
                }
            }
        }
    }
    //rootapp->Run();
    return 0;
}
