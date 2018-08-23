#include <algorithm>
#include <vector>
#include <map>
#include <iomanip>
#include <cstring>

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

void treeReader::Analyze(const vector<std::string> & filesToAnalyse, const std::string option, const std::string selection, const string& sampleToDebug, long evNb){

  debug = (option == "debug" ? true : false);
  leptonSelection = leptonSelectionAnalysis;
  initListsToPrint(selection);
  //Set CMS plotting style
  setTDRStyle();
  gROOT->SetBatch(kTRUE);
  //read samples and cross sections from txt file
  cout << "reading sample file...." << endl;
  samples.clear();
  for(auto & fileToAnalyse : filesToAnalyse)
    readSamples(fileToAnalyse);
  for(auto& sample : samples){
    std::cout << sample << std::endl;
  }
  cout << "finished with reading"<< endl;
  setTDRStyle(); 

  std::vector<std::string> namesOfFiles = treeReader::getNamesOfTheFiles();
  std::vector<std::string> namesOfProcesses = treeReader::getNamesOfTheProcesses();

  cout << "initiating histos...." << endl;
  initdistribs(namesOfProcesses);
  cout << "finished with initiating of histos"<< endl;
  setLabelsForHistos();

  std::ofstream myfile;
  myfile.open("myevents.txt");

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample();
      int samCategory = processIndex.at(samples[sam].getProcessName());

      Color_t color = assignColor(samples[sam].getProcessName());
      setStackColors(color, samCategory);

      //if(!(samples[sam].getFileName().find("ST_tWll_") != std::string::npos || samples[sam].getFileName().find("TTWJetsToLNu") != std::string::npos || samples[sam].getFileName().find("tZq_ll") != std::string::npos)) continue;
      //if(samples[sam].getProcessName() != "data" && samples[sam].getProcessName() != "nonpromptData" && samples[sam].getProcessName() != "Nonprompt") continue;

      if((option == "runOnOneProcess" || debug) && (samples[sam].getProcessName()) != sampleToDebug) continue;
      if(samples[sam].getProcessName() == "nonpromptData"){
          cout << "Total number of events: " << distribs[0].vectorHisto[samCategory].Integral() << endl;
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

          // in case during previous event run sam category was changed to nonprompt category 
          samCategory = processIndex.at(samples[sam].getProcessName());

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
          //if(it > nEntries / 20) break;

          std::vector<unsigned> indTight, indFake, indOf2LonZ;
          //select leptons relative to the analysis
          leptonSelection = 3;
          const unsigned lCount = selectLep(indTight, leptonSelection);
          const unsigned lCountFake = selectFakeLep(indFake, leptonSelection);

          std::vector<unsigned> indTight4L, indLoose4L;
          const unsigned lCount4L = selectLep(indTight4L, 4);
          const unsigned lCount4LLoose = selectFakeLep(indLoose4L, 4);

          // discard heavy flavour resonances
          if(debug) cout << "number of ttZ3L tight and fo leptons: " << lCount << " " << lCountFake << endl;
          if(debug) cout << "number of ttZ4L tight leptons: " << lCount4L << endl;

          // selection of category for the event
          // 2L: possible contribution from TT, TF and FF; TTF is vetoed
          // 3L: TTT, TTF, TFF, FFF; for TTTF should consider if event pass 4L TTTT criteria
          // 4L: TTTT only combination is possible
          
          std::vector<unsigned> ind;
          
          if(lCount4L == 4){
              if(lCount4LLoose != 4) continue;
              if(selection == "ttZ3L" || selection == "ttZ3Lclean" || selection == "DY" || selection == "Xgamma" || selection == "WZ" || selection == "ttbar") continue;
              leptonSelection = 4;
              //samCategory = sam;
              ind = indTight4L;
          }
          else if (lCount == 3){
              if(lCountFake != 3) continue;
              if(selection == "ZZ" || selection == "ttZ4L") continue;
              //samCategory = sam;
              ind = indTight;
          }
          else if (lCount < 3) {
              if(lCountFake != 3) continue;
              if(selection == "ZZ" || selection == "ttZ4L") continue;
              samCategory = nonPromptSample;
              ind = indFake;
          }
          else continue;

          if(debug) cout << "invariant mass of any fake pair is below 12 GeV: " << (leptonSelection == 3 ? invMassOfAny2Lbelow12GeV(indFake) : invMassOfAny2Lbelow12GeV(indLoose4L)) << endl;
          if(leptonSelection == 3 ? invMassOfAny2Lbelow12GeV(indFake) : invMassOfAny2Lbelow12GeV(indLoose4L)) continue; 

          if(debug) cout << "sum of all lepton charges: " << sumAllLeptonsCharge(ind) << endl;
          if(leptonSelection == 4 && sumAllLeptonsCharge(ind) != 0) continue;

          // consider only prompt leptons from the MC, all nonprompt should be taken into account by DD estimation
          bool allLeptonsArePrompt = true;
          
          if((samples[sam].getProcessName()) != "data" && (samples[sam].getProcessName()) != "nonpromptData" && (samples[sam].getProcessName()) != "chargeMisIDData")
            allLeptonsArePrompt = promptLeptons(ind);

          if(debug) cout << "all leptons are prompt ? " << allLeptonsArePrompt << endl;
          
          if((samples[sam].getProcessName()) == "chargeMisID" && !allLeptonsArePrompt) continue;
          if((samples[sam].getProcessName()) == "Nonprompt" && allLeptonsArePrompt) continue; // works just for MC

          if(leptonSelection == 3){
            if(((samples[sam].getProcessName()) == "ttW" || (samples[sam].getProcessName()) == "ttH" || (samples[sam].getProcessName()) == "ttZ" || (samples[sam].getProcessName()) == "ttX" 
                                                       || (samples[sam].getProcessName()) == "WZ" || (samples[sam].getProcessName()) == "Xgamma"  || (samples[sam].getProcessName()) == "ZZ" 
                                                       || (samples[sam].getProcessName()) == "rare") && !allLeptonsArePrompt) continue;
          }
          if(leptonSelection == 4){
            // WZ goes to nonprompt here
            if((samples[sam].getProcessName()) == "WZ") samCategory = nonPromptSample;
            if(!allLeptonsArePrompt) samCategory = nonPromptSample;
          }
          int nLocEle = getElectronNumber(ind);
          //if(nLocEle != 3) continue;

          // lepton pt criteria
          if(leptonSelection == 4)
            if(!passPtCuts4L(ind)) continue;
          if(leptonSelection == 3)
            if(!passPtCuts3L(ind)) continue;

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

          nJLoc = nJets(0, true, indJets, samples[sam].is2017());
          int nJLocDown = nJets(1, true, indJetsJECDown, samples[sam].is2017());
          int nJLocUp = nJets(2, true, indJetsJECUp, samples[sam].is2017());
          nJLocNotB = nJetsNotB(0, true, indJetsNotB, 2, samples[sam].is2017());
          nBLoc = nBJets(0, true, true, 1, samples[sam].is2017());
          int nBLocDown = nBJets(1, true, true, 1, samples[sam].is2017());
          int nBLocUp = nBJets(2, true, true, 1, samples[sam].is2017());

          TLorentzVector Zboson, lnegative;
          double dMZ = deltaMZ(ind, third, mll, ptZ, ptNonZ, mlll, indOf2LonZ, Zboson, lnegative);
          double mll1stpair, mll2ndpair;
          double cosTSt = -999;

          if(debug) cout << "number of jets/bjets/dMZ: " << nJLoc << " " << nBLoc << " " << dMZ << endl;
          if(debug && dMZ != 999999.) cout << "index of 2 leptons that makes 1st OSSF pair: " << indOf2LonZ.at(0) << " " << indOf2LonZ.at(1) << endl;

          HTLoc = HTCalc(indJets);
          double HTLocJECUp = HTCalc(indJetsJECUp);
          double HTLocJECDown  = HTCalc(indJetsJECDown);
          
          double mt1 = 9999;
          if(leptonSelection == 4){
            // used both in ttZ 4L and ZZ control region
            if(dMZ > 20) continue; 
            mll1stpair = mll;
            cosTSt = cosThetaStar(Zboson, lnegative);
            if(selection == "ttZ4L" && !passTTZ4LSelection(ind, indOf2LonZ, nJLoc)) continue;
            if(selection == "ZZ" && !passZZCRSelection(ind, indOf2LonZ, nJLoc)) continue;
            if(selection == "ttZ" && !(passTTZ4LSelection(ind, indOf2LonZ, nJLoc) || passZZCRSelection(ind, indOf2LonZ, nJLoc))) continue;
          }

          if(leptonSelection == 3){
            
            if(selection == "ttZ" && !(passTTZSelection(nJLoc, dMZ) || passWZCRSelection(nBLoc, dMZ) || passttbarCRSelection(nBLoc, dMZ, mlll))) continue;
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
            if(selection == "ttbar" &&!passttbarCRSelection(nBLoc, dMZ, mlll)) continue;
            if(selection == "Xgamma" && !passZGCRSelection(mlll, dMZ)) continue;

          }

          double mvaVL = 0;
          double mvaVLJECUp = 0;
          double mvaVLJECDown = 0;

          // weight estimation for event
          //auto start = std::chrono::high_resolution_clock::now();
          
          if((samples[sam].getProcessName()) != "data" && (samples[sam].getProcessName()) != "nonpromptData")
            weight *= sfWeight();
          if(samples[sam].getProcessName() == "data" && samCategory == nonPromptSample && leptonSelection == 3)
            weight *= fakeRateWeight();

          //auto finish = std::chrono::high_resolution_clock::now();
          //std::chrono::duration<double> elapsed = finish - start;
          //std::cout << "time needed to estimate event weight: " << elapsed.count() << std::endl;

          if(debug) cout << "weight of event is " << weight << endl;

          int mvaValueRegion = 0;

          if(debug) cout << "lepton selection is " << leptonSelection << " total SR: " << SRIDTTZ(ind, indOf2LonZ, nJLoc, nBLoc, dMZ, mlll) << endl; 
          if(leptonSelection == 4 && passZZCRSelection(ind, indOf2LonZ, nJLoc)) myfile << _runNb << " " << _lumiBlock << " " << _eventNb << endl;

          vector<double> fillVar = {ptCorrV[0].first, ptCorrV[1].first, leptonSelection > 2 ? ptCorrV[2].first : 0., leptonSelection > 3 ? ptCorrV[3].first : 0.,
                                   mt1, double(nJLoc), double(nBLoc), (_lCharge[ind.at(0)] == 1 ?  mvaVL : -999),
                                   // currently here we will have ttZ3L and ttZ4L categories
                                   //(leptonSelection == 4 ? flavourCategory4L(nLocEle) : (leptonSelection == 3 ? flavourCategory3L(nLocEle) : flavourCategory2L(nLocEle,_lCharge[ind.at(0)]))), 
                                   //(leptonSelection == 4 ? SRID4L(nBLoc) : (leptonSelection == 3 ? SRID3L(nJLoc, nBLoc) : SRID2L(mvaVL, _lCharge[ind.at(0)]))),
                                   (leptonSelection == 3 && passTTZSelection(nJLoc, dMZ) ? SRID3L(nJLoc, nBLoc) : -999),
                                   (leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLoc) ? SRID4L(nJLoc, nBLoc) : -999),
                                   (leptonSelection != 4 ? mll:mll1stpair),ptZ,ptNonZ, (nLocEle == 3?mll:-999.), (nLocEle==2?mll:-999.), (nLocEle==1?mll:-999.), (nLocEle==0? mll: -999.),
                                   _met, minDeltaR, minDeltaRlead, mtHighest, mtLowest, leadingJetPt, trailJetPt, 0., double(_nVertex), mlll,
                                   _lEta[ptCorrV[0].second], _lEta[ptCorrV[1].second], (leptonSelection > 2 ? _lEta[ptCorrV[2].second] : -999.), (leptonSelection > 3 ? _lEta[ptCorrV[3].second] : -999.),
                                   (nLocEle == 3?mt1:-999.), (nLocEle==2?mt1:-999.), (nLocEle==1?mt1:-999.), (nLocEle==0? mt1: -999.),
                                   cosTSt, mll_ss, double(chargeOfLeptons), ll_deltaR, mt2ll_ss,
                                   (_lCharge[ind.at(0)] == -1 ?  mvaVL : -999), HTLoc,
                                   SRIDTTZ(ind, indOf2LonZ, nJLoc, nBLoc, dMZ, mlll), SRIDWZCR(nJLoc, nBLoc, dMZ), SRIDZZCR(ind, indOf2LonZ, nJLoc, nBLoc), SRIDTTCR(nJLoc, nBLoc, dMZ, mlll),
                                   SRIDPTZ(ptZ), SRIDCosTheta(cosTSt),
                                   (leptonSelection == 3 ? flavourCategory3L(nLocEle) : -999),
                                   (leptonSelection == 4 ? flavourCategory4L(nLocEle) : -999),
                                   };

          vector<double> fillVarJecUp = {ptCorrV[0].first, ptCorrV[1].first, leptonSelection > 2 ? ptCorrV[2].first : 0., leptonSelection > 3 ? ptCorrV[3].first : 0.,
                                   mt1, double(nJLocUp), double(nBLocUp), (_lCharge[ind.at(0)] == 1 ?  mvaVLJECUp : -999),
                                   //(leptonSelection == 4 ? flavourCategory4L(nLocEle) : (leptonSelection == 3 ? flavourCategory3L(nLocEle) : flavourCategory2L(nLocEle,_lCharge[ind.at(0)]))), 
                                   //(leptonSelection == 4 ? SRID4L(nBLocUp) : (leptonSelection == 3 ? SRID3L(nJLocUp, nBLocUp) : SRID2L(mvaVLJECUp, _lCharge[ind.at(0)]))),
                                   (leptonSelection == 3 && passTTZSelection(nJLocUp, dMZ) ? SRID3L(nJLocUp, nBLocUp) : -999),
                                   (leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLocUp)? SRID4L(nJLocUp, nBLocUp) : -999),
                                   (leptonSelection != 4 ? mll:mll1stpair),ptZ,ptNonZ, (nLocEle == 3?mll:-999.), (nLocEle==2?mll:-999.), (nLocEle==1?mll:-999.), (nLocEle==0? mll: -999.),
                                   _met, minDeltaR, minDeltaRlead, mtHighest, mtLowest, leadingJetPt, trailJetPt, 0., double(_nVertex), mlll,
                                   _lEta[ptCorrV[0].second], _lEta[ptCorrV[1].second], (leptonSelection > 2 ? _lEta[ptCorrV[2].second] : -999.), (leptonSelection > 3 ? _lEta[ptCorrV[3].second] : -999.),
                                   (nLocEle == 3?mt1:-999.), (nLocEle==2?mt1:-999.), (nLocEle==1?mt1:-999.), (nLocEle==0? mt1: -999.),
                                   cosTSt, mll_ss, double(chargeOfLeptons), ll_deltaR, mt2ll_ss,
                                   (_lCharge[ind.at(0)] == -1 ?  mvaVLJECUp : -999), HTLocJECUp,
                                   SRIDTTZ(ind, indOf2LonZ, nJLocUp, nBLocUp, dMZ, mlll), SRIDWZCR(nJLocUp, nBLocUp, dMZ), SRIDZZCR(ind, indOf2LonZ, nJLocUp, nBLocUp), SRIDTTCR(nJLocUp, nBLocUp, dMZ, mlll),
                                   SRIDPTZ(ptZ), SRIDCosTheta(cosTSt),
                                   (leptonSelection == 3 ? flavourCategory3L(nLocEle) : -999),
                                   (leptonSelection == 4 ? flavourCategory4L(nLocEle) : -999),
                                   };

          vector<double> fillVarJecDw = {ptCorrV[0].first, ptCorrV[1].first, leptonSelection > 2 ? ptCorrV[2].first : 0., leptonSelection > 3 ? ptCorrV[3].first : 0.,
                                   mt1, double(nJLocDown), double(nBLocDown), (_lCharge[ind.at(0)] == 1 ?  mvaVLJECDown : -999),
                                   //(leptonSelection == 4 ? flavourCategory4L(nLocEle) : (leptonSelection == 3 ? flavourCategory3L(nLocEle) : flavourCategory2L(nLocEle,_lCharge[ind.at(0)]))), 
                                   //(leptonSelection == 4 ? SRID4L(nBLocDown) : (leptonSelection == 3 ? SRID3L(nJLocDown, nBLocDown) : SRID2L(mvaVLJECDown, _lCharge[ind.at(0)]))),
                                   (leptonSelection == 3 && passTTZSelection(nJLocDown, dMZ) ? SRID3L(nJLocDown, nBLocDown) : -999),
                                   (leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLocDown)? SRID4L(nJLocDown, nBLocDown) : -999),
                                   (leptonSelection != 4 ? mll:mll1stpair),ptZ,ptNonZ, (nLocEle == 3?mll:-999.), (nLocEle==2?mll:-999.), (nLocEle==1?mll:-999.), (nLocEle==0? mll: -999.),
                                   _met, minDeltaR, minDeltaRlead, mtHighest, mtLowest, leadingJetPt, trailJetPt, 0., double(_nVertex), mlll,
                                   _lEta[ptCorrV[0].second], _lEta[ptCorrV[1].second], (leptonSelection > 2 ? _lEta[ptCorrV[2].second] : -999.), (leptonSelection > 3 ? _lEta[ptCorrV[3].second] : -999.),
                                   (nLocEle == 3?mt1:-999.), (nLocEle==2?mt1:-999.), (nLocEle==1?mt1:-999.), (nLocEle==0? mt1: -999.),
                                   cosTSt, mll_ss, double(chargeOfLeptons), ll_deltaR, mt2ll_ss,
                                   (_lCharge[ind.at(0)] == -1 ?  mvaVLJECDown : -999), HTLocJECDown,
                                   SRIDTTZ(ind, indOf2LonZ, nJLocDown, nBLocDown, dMZ, mlll), SRIDWZCR(nJLocDown, nBLocDown, dMZ), SRIDZZCR(ind, indOf2LonZ, nJLocDown, nBLocDown), SRIDTTCR(nJLocDown, nBLocDown, dMZ, mlll),
                                   SRIDPTZ(ptZ), SRIDCosTheta(cosTSt),
                                   (leptonSelection == 3 ? flavourCategory3L(nLocEle) : -999),
                                   (leptonSelection == 4 ? flavourCategory4L(nLocEle) : -999),
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
                                     "SRallTTZ", "SRWZCR", "SRZZCR", "SRTTCR",
                                     "SRttZCleanPTZ", "SRttZCleanCosTheta",
                                     "flavour3L", "flavour4L"
                                   };
                                   
          //if(debug) cout << "lep sf: " << leptonWeight(0) << " " << leptonWeight(1) << " " << leptonWeight(2) << endl;

          //start = std::chrono::high_resolution_clock::now();
          double lepSF = 1.; double lepSFUp = 1.; double lepSFDown = 1.;

          double puW = 1.; double puWUp = 1.; double puWDown = 1.;

          double btagL = 1.; double btagLUp = 1.; double btagLDown = 1.;
          double btagC = 1.; double btagCUp = 1.; double btagCDown = 1.;
          double btagB = 1.; double btagBUp = 1.; double btagBDown = 1.;

          if((samples[sam].getProcessName()) != "data"){
            lepSF = leptonWeight(0); lepSFUp = leptonWeight(1); lepSFDown = leptonWeight(2);

            puW = puWeight(0); puWUp = puWeight(1); puWDown = puWeight(2);

            btagL = bTagWeight_udsg(0); btagLUp = bTagWeight_udsg(1); btagLDown = bTagWeight_udsg(2);
            btagC = bTagWeight_c(0); btagCUp = bTagWeight_c(1); btagCDown = bTagWeight_c(2);
            btagB = bTagWeight_b(0); btagBUp = bTagWeight_b(1); btagBDown = bTagWeight_b(2);

          }
          if(debug){
              if(leptonSelection == 3 && nJLoc == 2 && nBLoc  == 1)
                cout << "weights are (lepSF/pu/btagL/btagBC) : " << lepSF << " " << puW << " " << btagL << " " << btagC*btagB << endl;
          }

          //finish = std::chrono::high_resolution_clock::now();
          //elapsed = finish - start;
          //std::cout << "time needed to estimate all sf and deviations: " << elapsed.count() << std::endl;

          for(int dist = 0; dist < fillVar.size(); dist++){
            if(std::find(listToPrint[selection].begin(), listToPrint[selection].end(), fncName[dist]) == listToPrint[selection].end()) continue;
            //if(listToPrint[selection].find(fncName[dist]) == listToPrint[selection].end()) continue;
            distribs[dist].vectorHisto[samCategory].Fill(TMath::Min(fillVar.at(dist),figNames[fncName.at(dist)].varMax-0.1),weight);

            if((samples[sam].getProcessName()) != "data"){

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

                if(samples[sam].getProcessName() == "WZ" && nBLoc > 0){ // 8 % uncertainty for WZ + bb background in high nbjets categories
                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 5, figNames[fncName.at(dist)].varMax-0.1, weight * 1.08);
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 5, figNames[fncName.at(dist)].varMax-0.1, weight * 0.92);
                }
                else{
                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 5, figNames[fncName.at(dist)].varMax-0.1, weight);
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 5, figNames[fncName.at(dist)].varMax-0.1, weight);
                }

                if((samples[sam].getFileName().find("ST_tWll_") != std::string::npos || samples[sam].getFileName().find("TTWJetsToLNu") != std::string::npos || samples[sam].getFileName().find("tZq_ll") != std::string::npos || samples[sam].getFileName().find("TTZToLLNuNu_M-10_") != std::string::npos) && samples[sam].is2017()){ 
                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight*_psWeight[8]); // 10, 12 - factor 4; 6 and 8 - factor 2; 8, 9 - up factor 2, 6, 7 - down factor 2
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight*_psWeight[6]);

                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 7, figNames[fncName.at(dist)].varMax-0.1, weight*_psWeight[9]);
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 7, figNames[fncName.at(dist)].varMax-0.1, weight*_psWeight[7]);

                    //if(debug) cout << "ISR Up/Down: " << _psWeight[8] << " " << _psWeight[6] << "; FSR Up/Down: " << _psWeight[9] << " " << _psWeight[7] << endl;
                }
                else if(samples[sam].getFileName().find("TTZToLLNuNu_M-10_") != std::string::npos && samples[sam].is2016() && dist == indexSRTTZ){ // apply same weights as in 2017
                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight*ttZISRUpW[fillVar.at(dist)]); 
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight*ttZISRDownW[fillVar.at(dist)]);

                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 7, figNames[fncName.at(dist)].varMax-0.1, weight*ttZFSRUpW[fillVar.at(dist)]);
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 7, figNames[fncName.at(dist)].varMax-0.1, weight*ttZFSRDownW[fillVar.at(dist)]);
                }
                else{
                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight);
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight);

                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 7, figNames[fncName.at(dist)].varMax-0.1, weight);
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 7, figNames[fncName.at(dist)].varMax-0.1, weight);
                }

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 8, figNames[fncName.at(dist)].varMax-0.1, weight*_lheWeight[8]*sumSimulatedEventWeights/sumSimulatedEventWeightsScaleUp);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 8,figNames[fncName.at(dist)].varMax-0.1, weight*_lheWeight[4]*sumSimulatedEventWeights/sumSimulatedEventWeightsScaleDown);

                if(dist == indexSR3L || dist == indexSR4L || dist == indexSRTTZ || dist == indexSRWZCR || dist == indexSRZZCR || dist == indexSRTTCR){
                    for(int varPDF = 0; varPDF < 100; varPDF++){
                        distribs[dist].vectorHistoPDF[samCategory].var[varPDF].Fill(fillVar.at(dist), weight*_lheWeight[9+varPDF]);
                    }
                }
                else{
                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight);
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight);
                }


                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight*(1+uncOnNorm[samples[sam].getProcessName()]));
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight*(1-uncOnNorm[samples[sam].getProcessName()]));

                
            }
            else if(samCategory == nonPromptSample && leptonSelection != 4){
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

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 7, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 7, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 8, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 8, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight);

            }
          }
      }

      cout << endl;
      samCategory = processIndex.at(samples[sam].getProcessName());
      cout << "Total number of events: " << distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[samCategory].Integral() << endl;
      if(leptonSelection != 4)
        cout << "Total number of events in non prompt category: " << distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[nonPromptSample].Integral() << endl;
      cout << endl;
  }


  // this should be done to be fully correct in PDF, takes a lot of time, an effect estimated with ttZ sample, uncertainty on acceptance is under 1%, simply assign flat uncertainty of 1 % to all signal and bkg
  // for the moment draw this uncertainty only for SR and will propagate them to datacards 
  for(int dist = 0; dist < figNames.size(); dist++){
    if(!(dist == indexSR3L || dist == indexSR4L || dist == indexSRTTZ || dist == indexSRWZCR || dist == indexSRZZCR || dist == indexSRTTCR)) continue;

    for(unsigned sam = 0; sam < processIndex.size(); ++sam){
      if(sam == dataSample) continue;
      if(sam == nonPromptSample) continue;
      for(unsigned bin = 1; bin < (unsigned) distribs[dist].vectorHistoUncUp[sam].unc.at(9).GetNbinsX() + 1; ++bin){
          double pdfVarRms = 0.;
          for(unsigned pdf = 0; pdf < 100; ++pdf){
              double variedBin = distribs[dist].vectorHistoPDF[sam].var[pdf].GetBinContent(bin);
              variedBin *= crossSectionRatio[sam][pdf];
              double diff = (  variedBin - distribs[dist].vectorHisto[sam].GetBinContent(bin) );
              pdfVarRms += diff * diff;
          }
          pdfVarRms = sqrt( 0.01 * pdfVarRms );
          //cout << "pdf rms for bin " << bin << " is equal to " << pdfVarRms << endl;
          distribs[dist].vectorHistoUncUp[sam].unc.at(9).SetBinContent(bin, distribs[dist].vectorHisto[sam].GetBinContent(bin) + pdfVarRms);
          distribs[dist].vectorHistoUncDown[sam].unc.at(9).SetBinContent(bin, distribs[dist].vectorHisto[sam].GetBinContent(bin) - pdfVarRms);
      }
    }
  }

  // legend to print
  TLegend* mtleg = new TLegend(0.15,0.89,0.95,0.72); 
  mtleg->SetNColumns(4);
  mtleg->SetFillColor(0);
  mtleg->SetFillStyle(0);
  mtleg->SetBorderSize(0);
  mtleg->SetTextFont(42);

  mtleg->AddEntry(&distribs[0].vectorHisto[dataSample],"Data","lep"); //data

  std::map<int, std::string> processIndexReversed;
  std::vector<std::string> processOrder;
  for(map<std::string,int>::const_iterator it = processIndex.begin();it != processIndex.end(); ++it){
    processIndexReversed.insert(std::pair<int,std::string>(it->second, it->first));
  }

  // correct order according to increase
  for(map<int,std::string>::const_iterator it = processIndexReversed.begin();it != processIndexReversed.end(); ++it){
    processOrder.push_back(it->second);
  }
  
  for(map<int, std::string>::const_iterator it = processIndexReversed.begin();it != processIndexReversed.end(); ++it){
    //std::cout << it->first << " " << it->second << std::endl;
    if(it->second == "data") continue;
    if(it->second == "ttH") continue;

    if(selection == "ZZ" || selection == "ttZ4L"){
      if(it->second == "WZ") continue;
    }
    if(it->second != "nonpromptData" && it->second != "Xgamma")
      mtleg->AddEntry(&distribs[0].vectorHisto[it->first],it->second.c_str(),"f");
    else if(it->second == "nonpromptData")
      mtleg->AddEntry(&distribs[0].vectorHisto[it->first],"Nonprompt","f");
    else if(it->second == "Xgamma")
      mtleg->AddEntry(&distribs[0].vectorHisto[it->first],"X#gamma","f");
  }

  // plots to make with systematics and stat uncertainty on them
  std::string processToStore = selection;
  TString folderToStorePlots;
  int showLegendOption = 0; // 0 - 2016, 1 - 2017, 2 - 2016+2017
  // not implemented at the moment for more than 2 files
  if(filesToAnalyse.size() > 2)
    return;
  else if(filesToAnalyse.size() == 2){
    folderToStorePlots = "comb/";
    showLegendOption = 2;
  }
  else{
    if(filesToAnalyse.at(0).find("2017") != std::string::npos){
      cout << "file to analysis is " << filesToAnalyse.at(0) << endl;
      folderToStorePlots = "2017/";
      showLegendOption = 1;
    }
    else{
      folderToStorePlots = "2016/";
       showLegendOption = 0;
    }
  }
  gSystem->Exec("rm plotsForSave/" + folderToStorePlots + processToStore + "/*.{pdf,png,root}");
  gSystem->Exec("rmdir plotsForSave/" + folderToStorePlots + processToStore);
  gSystem->Exec("mkdir plotsForSave/" + folderToStorePlots + processToStore);
  double scale_num = 1.6;
  
  TCanvas* plot[nVars];
  for(int i = 0; i < nVars; i++){
      plot[i] = new TCanvas(Form("plot_%d", i),"",500,450);
  }

  std::string crToPrint = selection;

  for(int varPlot = 0; varPlot < listToPrint[crToPrint].size(); varPlot++){
    plot[varPlot]->cd();
    showHist(plot[varPlot],distribs[figNames[listToPrint[crToPrint].at(varPlot)].index],figNames[listToPrint[crToPrint].at(varPlot)], scale_num, mtleg, false, false, showLegendOption);
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + ".pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + ".png");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + ".root");
    plot[varPlot]->cd();
    showHist(plot[varPlot],distribs[figNames[listToPrint[crToPrint].at(varPlot)].index],figNames[listToPrint[crToPrint].at(varPlot)], scale_num, mtleg, true, false, showLegendOption);
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + "Log.pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + "Log.png");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + "Log.root");
  }
  
  /*
  if(crToPrint == "ttZ3Lclean"){
    fillDatacards(distribs[indexSRttZcleanPTZ], samplesOrderNames, samplesOrder, "SRttZCleanPTZ", is2017);
    fillDatacards(distribs[indexSRttZcleanCosTheta], samplesOrderNames, samplesOrder, "SRttZCleanCosTheta", is2017);
  }
  */
  // no need to fill datacards when running over 2016 and 2017 together
  if(showLegendOption == 2) return;
  if(crToPrint == "ttZ"){
    fillDatacards(distribs[indexSR3L], processOrder, "SR3L", (bool)showLegendOption); 
    fillDatacards(distribs[indexSR4L], processOrder, "SR4L", (bool)showLegendOption); 
    fillDatacards(distribs[indexSRTTZ], processOrder, "SRallTTZ", (bool)showLegendOption); 
    fillDatacards(distribs[indexSRWZCR], processOrder, "SRWZCR", (bool)showLegendOption); 
    fillDatacards(distribs[indexSRZZCR], processOrder, "SRZZCR", (bool)showLegendOption); 
    fillDatacards(distribs[indexSRTTCR], processOrder, "SRTTCR", (bool)showLegendOption); 
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
        std::cerr << "please specify one of the options (runFullSelection, runOnOneProcess, debug)" << std::endl;
        return 1;
        //reader.Analyze(std::string(argv[1]));
    }    
    if(argc > 2){
        if(argc == 3) {
            if(string(argv[2]) == "runFullSelection"){
                std::cerr << "please specify control region (ttZ3L, ttZ4L, ttZ, WZ, ZZ, ttbar, DY, Xgamma), use \'selection:\' before control region" << std::endl;
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
                    std::vector<std::string> inputFiles;
                    char * inF = (char*)argv[1];
                    char * pch = strtok (inF,",");
                    while (pch != NULL){
                      inputFiles.push_back(std::string(pch)); // printf ("%s\n",pch);
                      pch = strtok (NULL, ",");
                    }
                    reader.Analyze(inputFiles, std::string(argv[2]), selection); 
                }
            }
            else if(string(argv[2]) == "runOnOneProcess"){
                std::cerr << "please specify process to run on" << std::endl;
                if(string(argv[3]).find("selection:") == std::string::npos){
                  std::cerr << "before specifying which process to run on, please specify as well selection, use \'selection:\' before control region (ttZ3L, ttZ4L, ttZ, WZ, ZZ, ttbar, DY, Xgamma)" << std::endl;
                }
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
                    std::vector<std::string> inputFiles;
                    char * inF = (char*)argv[1];
                    char * pch = strtok (inF,",");
                    while (pch != NULL){
                      inputFiles.push_back(std::string(pch)); // printf ("%s\n",pch);
                      pch = strtok (NULL, ",");
                    }
                    reader.Analyze(inputFiles, std::string(argv[2]), selection, std::string(argv[4])); 
                }
            }
            else if(string(argv[2]) == "debug"){
                if(string(argv[3]).find("selection:") != std::string::npos){
                    std::string selection = string(argv[3]);
                    selection.erase (selection.begin(), selection.begin()+10);
                    std::cout << "output folder is set to: " << selection<< std::endl;
                    std::vector<std::string> inputFiles;
                    char * inF = (char*)argv[1];
                    char * pch = strtok (inF,",");
                    while (pch != NULL){
                      inputFiles.push_back(std::string(pch)); // printf ("%s\n",pch);
                      pch = strtok (NULL, ",");
                    }
                    reader.Analyze(inputFiles, std::string(argv[2]), selection, std::string(argv[4])); 
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
                    std::vector<std::string> inputFiles;
                    char * inF = (char*)argv[1];
                    char * pch = strtok (inF,",");
                    while (pch != NULL){
                      inputFiles.push_back(std::string(pch)); // printf ("%s\n",pch);
                      pch = strtok (NULL, ",");
                    }
                    reader.Analyze(inputFiles, std::string(argv[2]), selection, std::string(argv[4]), atol(argv[5])); 
                }
            }
        }
    }
    //rootapp->Run();
    return 0;
}
