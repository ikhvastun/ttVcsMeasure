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
#include "../interface/fillTables.h"
#include "../interface/PostFitScaler.h"

#include "tdrStyle.C"

using namespace std;
using namespace tools;

Errors LastError::lasterror = Errors::UNKNOWN;
using Output::distribs;

void treeReader::Analyze(const vector<std::string> & filesToAnalyse, const std::string option, const std::string selection, const string& sampleToDebug, long evNb){

  debug = (option == "debug" ? true : false);
//  leptonSelection = leptonSelectionAnalysis;
//  in other words, it defines at what you want to look at, e.g. 3L ttZ, WZ control region etc. 
//  This is also used to initialize only the histograms necesarry for the considered process/region/selection
  initListToPrint(selection);
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

  //  std::vector<std::string> namesOfFiles = treeReader::getNamesOfTheFiles();
  std::vector<std::string> namesOfProcesses = treeReader::getNamesOfTheProcesses();

  cout << "initiating histos...." << endl;
  initdistribs(namesOfProcesses, selection);
  cout << "finished with initiating of histos"<< endl;
  setLabelsForHistos(selection);

  // event number dump for debugging/comparing

  std::ofstream myfile;
  myfile.open("myevents.txt");

  // here you load in post fit weights, in case you have them already and want your histograms
  // to show post fit values

  PostFitScaler scaler2016, scaler2017;
  scaler2016.setInputFile("data/postFit/outputTTZ_2016_new.txt");
  scaler2016.setPostfitYields();
  scaler2017.setInputFile("data/postFit/outputTTZ_2017_new.txt");
  scaler2017.setPostfitYields();

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over all samples. Reads in from the text file. For each category, separate histograms 
  // are declared. Important: the last entry is the nonprompt data, which uses the same data root 
  // file. it's there to make sure the histograms are created.
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  
  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample("ttZ");
      //initSample("ttZ4l");
      int samCategory = processToCounterMap.at(samples[sam].getProcessName());

      Color_t color = assignColor(samples[sam].getProcessName());
      setStackColors(color, samCategory);

      //if(!(samples[sam].getFileName().find("ttHToNonbb") != std::string::npos || samples[sam].getFileName().find("ST_tWll_") != std::string::npos || samples[sam].getFileName().find("TTWJetsToLNu") != std::string::npos || samples[sam].getFileName().find("tZq_ll") != std::string::npos)) continue;
      //if(samples[sam].getProcessName() != 
      //if(samples[sam].getProcessName() != "data" && samples[sam].getProcessName() != "WZ") continue;
      //if(samples[sam].getProcessName() == "data") continue;

      if((option == "runOnOneProcess" || debug) && (samples[sam].getProcessName()) != sampleToDebug) continue;

      // adding here WZc and WZb
      // no need to run over them, just using them in the sample list file to create histograms
      // the inputs will be taken from WZ sample
      if(samples[sam].getProcessName() == "nonpromptData" || samples[sam].getProcessName() == "WZc" || samples[sam].getProcessName() == "WZb"){
          cout << "Total number of events in " << samples[sam].getProcessName() << " sample is " << distribs[0].vectorHisto[samCategory].Integral() << endl;
          continue;
      }

      std::cout<<"Entries in "<< (samples[sam].getFileName()) << " " << nEntries << std::endl;
      double progress = 0;  //for printing progress bar
		// 
		// loop over all events in a given sample
		//
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
          samCategory = processToCounterMap.at(samples[sam].getProcessName());

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
          //if(it > nEntries / 50) break;

          std::vector<unsigned> indTight, indFake, indOf2LonZ;
          //select leptons relative to the analysis
			    // start with 3 lepton selection, but it changes later
			    // the counts of tight/loose/tight4l/loose4l leptons are calculated here
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
              ind = indTight4L;
          }
          else if (lCount == 3){
              if(lCountFake != 3) continue;
              if(selection == "ZZ" || selection == "ttZ4L") continue;
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
            if(((samples[sam].getProcessName()) == "ttW" || (samples[sam].getProcessName()) == "ttH" || (samples[sam].getProcessName()) == "ttZ"     || (samples[sam].getProcessName()) == "ttX" 
                                                         || (samples[sam].getProcessName()) == "WZ"  || (samples[sam].getProcessName()) == "Xgamma"  || (samples[sam].getProcessName()) == "ZZ" 
                                                         || (samples[sam].getProcessName()) == "rare") && !allLeptonsArePrompt) continue;
          }
          if(leptonSelection == 4){
            // WZ goes to nonprompt in 4L category
            if((samples[sam].getProcessName()) == "WZ") samCategory = nonPromptSample;
            if(!allLeptonsArePrompt) samCategory = nonPromptSample;
          }

          int nLocEle = getElectronNumber(ind);
          //if(nLocEle != 0) continue;

          // lepton pt criteria
          if(leptonSelection == 4)
            if(!passPtCuts4L(ind)) continue;
          if(leptonSelection == 3)
            if(!passPtCuts3L(ind)) continue;

          // select here jets, bjets, delta from M of Z boson, HT
          std::vector<unsigned> indJets;
          std::vector<unsigned> indJetsJECUp;
          std::vector<unsigned> indJetsJECDown;
          std::vector<unsigned> indJetsJERUp;
          std::vector<unsigned> indJetsJERDown;
          std::vector<unsigned> indJetsNotB;

          std::vector<unsigned> indBJets;
          std::vector<unsigned> indBJetsJECUp;
          std::vector<unsigned> indBJetsJECDown;
          std::vector<unsigned> indBJetsJERUp;
          std::vector<unsigned> indBJetsJERDown;

          unsigned third = -9999;
          double mll = 99999;
          double mlll = 99999;
          double ptZ = 999999;
          double ptNonZ = -999999;

          nJLoc = nJets(0, true, indJets, samples[sam].is2017());
          int nJLocDown = nJets(1, true, indJetsJECDown, samples[sam].is2017());
          int nJLocUp = nJets(2, true, indJetsJECUp, samples[sam].is2017());
          int nJLocJERDown = nJets(3, true, indJetsJERDown, samples[sam].is2017());
          int nJLocJERUp = nJets(4, true, indJetsJERUp, samples[sam].is2017());
          nJLocNotB = nJetsNotB(0, true, indJetsNotB, 2, samples[sam].is2017());
          nBLoc = nBJets(0, true, true, indBJets, 1, samples[sam].is2017());
          int nBLocDown = nBJets(1, true, true, indBJetsJECDown, 1, samples[sam].is2017());
          int nBLocUp = nBJets(2, true, true, indBJetsJECUp, 1, samples[sam].is2017());
          int nBLocJERDown = nBJets(3, true, true, indBJetsJERDown, 1, samples[sam].is2017());
          int nBLocJERUp = nBJets(4, true, true, indBJetsJERUp, 1, samples[sam].is2017());

          // check how jets in WZ sample are light, c- and b-jets
          if(samples[sam].getProcessName() == "WZ" && leptonSelection == 3){
            int numberOfCJets = getNumberOfCJets(indJets);
            int numberOfBJets = getNumberOfBJets(indJets);
            if(numberOfCJets != 0) samCategory = processToCounterMap.at("WZc");
            if(numberOfBJets != 0) samCategory = processToCounterMap.at("WZb");
          }
          

          TLorentzVector Zboson, lnegative;
          double dMZ = deltaMZ(ind, third, mll, ptZ, ptNonZ, mlll, indOf2LonZ, Zboson, lnegative);
          double mll1stpair, mll2ndpair;
          double cosTSt = -999;

          if(debug) cout << "number of jets/bjets/dMZ: " << nJLoc << " " << nBLoc << " " << dMZ << endl;
          if(debug && dMZ != 999999.) cout << "index of 2 leptons that makes 1st OSSF pair: " << indOf2LonZ.at(0) << " " << indOf2LonZ.at(1) << endl;

          HTLoc = HTCalc(indJets);
          double HTLocJECUp = HTCalc(indJetsJECUp);
          double HTLocJECDown  = HTCalc(indJetsJECDown);
          double HTLocJERUp = HTCalc(indJetsJERUp);
          double HTLocJERDown  = HTCalc(indJetsJERDown);
          
          double mt1 = 9999;
          if(leptonSelection == 4){
            // used both in ttZ 4L and ZZ control region
            if(dMZ > 20) continue; 
            mll1stpair = mll;
            cosTSt = cosThetaStar(Zboson, lnegative);
            if(selection == "ttZ4L" && !passTTZ4LSelection(ind, indOf2LonZ, nJLoc)) continue;
            if(selection == "ZZ" && !passZZCRSelection(ind, indOf2LonZ, nJLoc)) continue;
            if(selection == "ttZ" && !(passTTZ4LSelection(ind, indOf2LonZ, nJLoc) || passZZCRSelection(ind, indOf2LonZ, nJLoc))) continue;
            if(selection == "ttZclean" && !passTTZCleanSelection(nJLoc, nBLoc, dMZ)) continue;
          }

          if(leptonSelection == 3){
            
            if(selection == "ttZ" && !(passTTZSelection(nJLoc, dMZ) || passWZCRSelection(nBLoc, dMZ) || passttbarCRSelection(nBLoc, dMZ, mlll))) continue;
            if(selection == "tZq" &&!passttbarCRintZqSelection(nJLoc, nBLoc, dMZ)) continue;
            if(selection == "ttZ3L" && !passTTZSelection(nJLoc, dMZ)) continue;

            // here are some thoughts: clean selection is defined with njets >= 3, nbjets >= 1 cuts, if it's defined here at selection step then no events with njets == 2 selection and having 3 jets with upward variation will not enter the selection
            // should consider njets == 2 selection and later when filling the histograms ask for 3 jets with variation to pass selection 
            // this option is used to answer the question from Giovanni, 6th of Feb 2019
            if(selection == "ttZ3Lclean" && !passTTZCleanSelection(nJLoc, nBLoc, dMZ)) continue; 
            
            // this option should be used to get the datacards to Joscha
            //if(selection == "ttZ3Lclean" && !passTTZSelection(nJLoc, dMZ)) continue;
            if(selection == "ttZclean" && !passTTZCleanSelection(nJLoc, nBLoc, dMZ)) continue;
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
          if(samples[sam].getProcessName() == "data" && samCategory == dataSample) myfile << _runNb << " " << _lumiBlock << " " << _eventNb << endl;

          // Histograms filling starts here
          // define histograms to be filled in the loop
          map<TString, double> fillVar; 
          fillVar["ptlead"] = ptCorrV[0].first;
          fillVar["sublead"] = ptCorrV[1].first;
          fillVar["trail"] = leptonSelection > 2 ? ptCorrV[2].first : 0.;
          fillVar["pt4th"] = leptonSelection > 3 ? ptCorrV[3].first : 0.;
          fillVar["etaLead"] = _lEta[ptCorrV[0].second];
          fillVar["etaSubl"] = _lEta[ptCorrV[1].second];
          fillVar["etaTrail"] = leptonSelection > 2 ? _lEta[ptCorrV[2].second] : -999.;
          fillVar["eta4th"] = leptonSelection > 3 ? _lEta[ptCorrV[3].second] : -999.;

          fillVar["ptZ"] = ptZ;
          fillVar["cosThetaStar"] = cosTSt;
          fillVar["ptNonZ"] = ptNonZ;

          fillVar["mll"] = leptonSelection != 4 ? mll:mll1stpair;
          fillVar["mll3e"] = nLocEle == 3 ? mll : -999.;
          fillVar["mll2e1mu"] = nLocEle == 2 ? mll : -999.;
          fillVar["mll1e2mu"] = nLocEle == 1 ? mll : -999.;
          fillVar["mll3mu"] = nLocEle == 0 ? mll : -999.;
          fillVar["mllnoZcut"] = leptonSelection != 4 ? mll : mll1stpair;

          fillVar["njets"] = double(nJLoc);
          fillVar["nbjets"] = double(nBLoc);
          fillVar["HT"] = HTLoc;
          fillVar["met"] = _met;
          fillVar["nPV"] = double(_nVertex);

          fillVar["SR3L"] = leptonSelection == 3 ? SRID3L(nJLoc, nBLoc, dMZ) : -999;
          fillVar["SR3L3m"] = leptonSelection == 3 && nLocEle == 0 ? SRID3L(nJLoc, nBLoc, dMZ) : -999;
          fillVar["SR3L2m1e"] = leptonSelection == 3 && nLocEle == 1 ? SRID3L(nJLoc, nBLoc, dMZ) : -999;
          fillVar["SR3L1m2e"] = leptonSelection == 3 && nLocEle == 2 ? SRID3L(nJLoc, nBLoc, dMZ) : -999;
          fillVar["SR3L3e"] = leptonSelection == 3 && nLocEle == 3 ? SRID3L(nJLoc, nBLoc, dMZ) : -999;
          fillVar["SR4L"] = leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLoc) ? SRID4L(nJLoc, nBLoc) : -999;

          fillVar["SRTTZ8SR3L"] = leptonSelection == 3 && nBLoc > 0 ? SRID8SR3L(nJLoc, nBLoc, dMZ) : -999;
          fillVar["SRallTTZ"] = SRIDTTZ(ind, indOf2LonZ, nJLoc, nBLoc, dMZ, mlll);
          fillVar["SRWZCR"] = SRIDWZCR(nJLoc, nBLoc, dMZ);
          fillVar["SRZZCR"] = SRIDZZCR(ind, indOf2LonZ, nJLoc, nBLoc);
          fillVar["SRTTCR"] = SRIDTTCR(nJLoc, nBLoc, dMZ, mlll);
          fillVar["SRttZCleanPTZ"] = leptonSelection == 3 && passTTZCleanSelection(nJLoc, nBLoc, dMZ) ? SRIDPTZ(ptZ) : -999;
          fillVar["SRttZCleanCosTheta"] = leptonSelection == 3 && passTTZCleanSelection(nJLoc, nBLoc, dMZ) ? SRIDCosTheta(cosTSt) : -999;

          fillVar["flavour3L"] = leptonSelection == 3 ? flavourCategory3L(nLocEle) : -999;
          fillVar["flavour4L"] = leptonSelection == 4 ? flavourCategory4L(nLocEle) : -999;
          fillVar["flavour4LZZ"] = leptonSelection == 4 ? flavourCategory4LZZ(nLocEle) : -999;
          fillVar["flavour3L4L"] = passTTZSRSelection(ind, indOf2LonZ, nJLoc, nBLoc, dMZ) ? flavourCategory3L4L(leptonSelection, nLocEle) : -999;

          // same as before but with JEC uncertainty up
          map<TString, double> fillVarJecUp(fillVar); 
          fillVarJecUp["njets"] = nJLocUp;
          fillVarJecUp["nbjets"] = nBLocUp;
          fillVarJecUp["HT"] = HTLocJECUp;

          fillVarJecUp["SR3L"] = leptonSelection == 3 ? SRID3L(nJLocUp, nBLocUp, dMZ) : -999;
          fillVarJecUp["SR3L3m"] = leptonSelection == 3 && nLocEle == 0 ? SRID3L(nJLocUp, nBLocUp, dMZ) : -999;
          fillVarJecUp["SR3L2m1e"] = leptonSelection == 3 && nLocEle == 1 ? SRID3L(nJLocUp, nBLocUp, dMZ) : -999;
          fillVarJecUp["SR3L1m2e"] = leptonSelection == 3 && nLocEle == 2 ? SRID3L(nJLocUp, nBLocUp, dMZ) : -999;
          fillVarJecUp["SR3L3e"] = leptonSelection == 3 && nLocEle == 3 ? SRID3L(nJLocUp, nBLocUp, dMZ) : -999;

          fillVarJecUp["SR4L"] = leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLocUp) ? SRID4L(nJLocUp, nBLocUp) : -999;
          fillVarJecUp["SRTTZ8SR3L"] = leptonSelection == 3 && nBLocUp > 0 ? SRID8SR3L(nJLocUp, nBLocUp, dMZ) : -999;
          fillVarJecUp["SRallTTZ"] = SRIDTTZ(ind, indOf2LonZ, nJLocUp, nBLocUp, dMZ, mlll);
          fillVarJecUp["SRWZCR"] = SRIDWZCR(nJLocUp, nBLocUp, dMZ);
          fillVarJecUp["SRZZCR"] = SRIDZZCR(ind, indOf2LonZ, nJLocUp, nBLocUp);
          fillVarJecUp["SRTTCR"] = SRIDTTCR(nJLocUp, nBLocUp, dMZ, mlll);
          fillVarJecUp["SRttZCleanPTZ"] = leptonSelection == 3 && passTTZCleanSelection(nJLocUp, nBLocUp, dMZ) ? SRIDPTZ(ptZ) : -999;
          fillVarJecUp["SRttZCleanCosTheta"] = leptonSelection == 3 && passTTZCleanSelection(nJLocUp, nBLocUp, dMZ) ? SRIDCosTheta(cosTSt) : -999;

          // JEC uncertainty down
          map<TString, double> fillVarJecDown(fillVar); 
          fillVarJecDown["njets"] = nJLocDown;
          fillVarJecDown["nbjets"] = nBLocDown;
          fillVarJecDown["HT"] = HTLocJECDown;

          fillVarJecDown["SR3L"] = leptonSelection == 3 ? SRID3L(nJLocDown, nBLocDown, dMZ) : -999;
          fillVarJecDown["SR3L3m"] = leptonSelection == 3 && nLocEle == 0 ? SRID3L(nJLocDown, nBLocDown, dMZ) : -999;
          fillVarJecDown["SR3L2m1e"] = leptonSelection == 3 && nLocEle == 1 ? SRID3L(nJLocDown, nBLocDown, dMZ) : -999;
          fillVarJecDown["SR3L1m2e"] = leptonSelection == 3 && nLocEle == 2 ? SRID3L(nJLocDown, nBLocDown, dMZ) : -999;
          fillVarJecDown["SR3L3e"] = leptonSelection == 3 && nLocEle == 3 ? SRID3L(nJLocDown, nBLocDown, dMZ) : -999;

          fillVarJecDown["SR4L"] = leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLocDown) ? SRID4L(nJLocDown, nBLocDown) : -999;
          fillVarJecDown["SRTTZ8SR3L"] = leptonSelection == 3 && nBLocDown > 0 ? SRID8SR3L(nJLocDown, nBLocDown, dMZ) : -999;
          fillVarJecDown["SRallTTZ"] = SRIDTTZ(ind, indOf2LonZ, nJLocDown, nBLocDown, dMZ, mlll);
          fillVarJecDown["SRWZCR"] = SRIDWZCR(nJLocDown, nBLocDown, dMZ);
          fillVarJecDown["SRZZCR"] = SRIDZZCR(ind, indOf2LonZ, nJLocDown, nBLocDown);
          fillVarJecDown["SRTTCR"] = SRIDTTCR(nJLocDown, nBLocDown, dMZ, mlll);
          fillVarJecDown["SRttZCleanPTZ"] = leptonSelection == 3 && passTTZCleanSelection(nJLocDown, nBLocDown, dMZ) ? SRIDPTZ(ptZ) : -999;
          fillVarJecDown["SRttZCleanCosTheta"] = leptonSelection == 3 && passTTZCleanSelection(nJLocDown, nBLocDown, dMZ) ? SRIDCosTheta(cosTSt) : -999;

          // JER uncertainty up
          map<TString, double> fillVarJerUp(fillVar); 
          fillVarJerUp["njets"] = nJLocJERUp;
          fillVarJerUp["nbjets"] = nBLocJERUp;
          fillVarJerUp["HT"] = HTLocJERUp;

          fillVarJerUp["SR3L"] = leptonSelection == 3 ? SRID3L(nJLocJERUp, nBLocJERUp, dMZ) : -999;
          fillVarJerUp["SR3L3m"] = leptonSelection == 3 && nLocEle == 0 ? SRID3L(nJLocJERUp, nBLocJERUp, dMZ) : -999;
          fillVarJerUp["SR3L2m1e"] = leptonSelection == 3 && nLocEle == 1 ? SRID3L(nJLocJERUp, nBLocJERUp, dMZ) : -999;
          fillVarJerUp["SR3L1m2e"] = leptonSelection == 3 && nLocEle == 2 ? SRID3L(nJLocJERUp, nBLocJERUp, dMZ) : -999;
          fillVarJerUp["SR3L3e"] = leptonSelection == 3 && nLocEle == 3 ? SRID3L(nJLocJERUp, nBLocJERUp, dMZ) : -999;

          fillVarJerUp["SR4L"] = leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLocJERUp) ? SRID4L(nJLocJERUp, nBLocJERUp) : -999;
          fillVarJerUp["SRTTZ8SR3L"] = leptonSelection == 3 && nBLocJERUp > 0 ? SRID8SR3L(nJLocJERUp, nBLocJERUp, dMZ) : -999;
          fillVarJerUp["SRallTTZ"] = SRIDTTZ(ind, indOf2LonZ, nJLocJERUp, nBLocJERUp, dMZ, mlll);
          fillVarJerUp["SRWZCR"] = SRIDWZCR(nJLocJERUp, nBLocJERUp, dMZ);
          fillVarJerUp["SRZZCR"] = SRIDZZCR(ind, indOf2LonZ, nJLocJERUp, nBLocJERUp);
          fillVarJerUp["SRTTCR"] = SRIDTTCR(nJLocJERUp, nBLocJERUp, dMZ, mlll);
          fillVarJerUp["SRttZCleanPTZ"] = leptonSelection == 3 && passTTZCleanSelection(nJLocJERUp, nBLocJERUp, dMZ) ? SRIDPTZ(ptZ) : -999;
          fillVarJerUp["SRttZCleanCosTheta"] = leptonSelection == 3 && passTTZCleanSelection(nJLocJERUp, nBLocJERUp, dMZ) ? SRIDCosTheta(cosTSt) : -999;

          // JER uncertainty down
          map<TString, double> fillVarJerDw(fillVar); 
          fillVarJerDw["njets"] = nJLocJERDown;
          fillVarJerDw["nbjets"] = nBLocJERDown;
          fillVarJerDw["HT"] = HTLocJERDown;

          fillVarJerDw["SR3L"] = leptonSelection == 3 ? SRID3L(nJLocJERDown, nBLocJERDown, dMZ) : -999;
          fillVarJerDw["SR3L3m"] = leptonSelection == 3 && nLocEle == 0 ? SRID3L(nJLocJERDown, nBLocJERDown, dMZ) : -999;
          fillVarJerDw["SR3L2m1e"] = leptonSelection == 3 && nLocEle == 1 ? SRID3L(nJLocJERDown, nBLocJERDown, dMZ) : -999;
          fillVarJerDw["SR3L1m2e"] = leptonSelection == 3 && nLocEle == 2 ? SRID3L(nJLocJERDown, nBLocJERDown, dMZ) : -999;
          fillVarJerDw["SR3L3e"] = leptonSelection == 3 && nLocEle == 3 ? SRID3L(nJLocJERDown, nBLocJERDown, dMZ) : -999;

          fillVarJerDw["SR4L"] = leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLocJERDown) ? SRID4L(nJLocJERDown, nBLocJERDown) : -999;
          fillVarJerDw["SRTTZ8SR3L"] = leptonSelection == 3 && nBLocJERDown > 0 ? SRID8SR3L(nJLocJERDown, nBLocJERDown, dMZ) : -999;
          fillVarJerDw["SRallTTZ"] = SRIDTTZ(ind, indOf2LonZ, nJLocJERDown, nBLocJERDown, dMZ, mlll);
          fillVarJerDw["SRWZCR"] = SRIDWZCR(nJLocJERDown, nBLocJERDown, dMZ);
          fillVarJerDw["SRZZCR"] = SRIDZZCR(ind, indOf2LonZ, nJLocJERDown, nBLocJERDown);
          fillVarJerDw["SRTTCR"] = SRIDTTCR(nJLocJERDown, nBLocJERDown, dMZ, mlll);
          fillVarJerDw["SRttZCleanPTZ"] = leptonSelection == 3 && passTTZCleanSelection(nJLocJERDown, nBLocJERDown, dMZ) ? SRIDPTZ(ptZ) : -999;
          fillVarJerDw["SRttZCleanCosTheta"] = leptonSelection == 3 && passTTZCleanSelection(nJLocJERDown, nBLocJERDown, dMZ) ? SRIDCosTheta(cosTSt) : -999;

                                  
          //start = std::chrono::high_resolution_clock::now();
          // SF estimation starts here
          double lepSF = 1.;
          double lepSFSystUp = 1.; double lepSFSystDown = 1.; 
          double lepSFStatUp = 1.; double lepSFStatDown = 1.;
          double lepSFRecoUp = 1.; double lepSFRecoDown = 1.;

          double puW = 1.; double puWUp = 1.; double puWDown = 1.;

          double btagL = 1.; double btagLUp = 1.; double btagLDown = 1.;
          double btagC = 1.; double btagCUp = 1.; double btagCDown = 1.;
          double btagB = 1.; double btagBUp = 1.; double btagBDown = 1.;

          double jetPrefW = 1.;

          //if((samples[sam].getProcessName()) != "data"){
          if(samCategory != dataSample){
            double lepWOnlySyst = leptonWeightOnlySyst(0); double lepWOnlySystUp = leptonWeightOnlySyst(1); double lepWOnlySystDown = leptonWeightOnlySyst(2);
            double lepWOnlyStat = leptonWeightOnlyStat(0); double lepWOnlyStatUp = leptonWeightOnlyStat(1); double lepWOnlyStatDown = leptonWeightOnlyStat(2);
            double lepWOnlyReco = leptonWeightOnlyReco(0); double lepWOnlyRecoUp = leptonWeightOnlyReco(1); double lepWOnlyRecoDown = leptonWeightOnlyReco(2);

            if(samCategory != nonPromptSample){
                lepSF = lepWOnlySyst * lepWOnlyReco; // consider only one between syst and stat, central value is the same
                lepSFSystUp = lepWOnlySystUp * lepWOnlyReco; 
                lepSFSystDown = lepWOnlySystDown * lepWOnlyReco;
                lepSFStatUp = lepWOnlyStatUp * lepWOnlyReco; 
                lepSFStatDown = lepWOnlyStatDown * lepWOnlyReco; 
                lepSFRecoUp = lepWOnlySyst * lepWOnlyRecoUp; 
                lepSFRecoDown = lepWOnlySyst * lepWOnlyRecoDown; 

                puW = puWeight(0); puWUp = puWeight(1); puWDown = puWeight(2);

                btagL = bTagWeight_udsg(0); btagLUp = bTagWeight_udsg(1); btagLDown = bTagWeight_udsg(2);
                btagC = bTagWeight_c(0); btagCUp = bTagWeight_c(1); btagCDown = bTagWeight_c(2);
                btagB = bTagWeight_b(0); btagBUp = bTagWeight_b(1); btagBDown = bTagWeight_b(2);
            }

            //jetPrefW = jetPrefiringWeight();

            // for post fit scaling
            // used previously
            // weight *= scaler.postFitScaling(samples[sam].getProcessName());
            
            if(samples[sam].is2016())
                weight *= scaler2016.postFitScaling(samples[sam].getProcessName() != "data" ? samples[sam].getProcessName() : "nonpromptData");
            else
                weight *= scaler2017.postFitScaling(samples[sam].getProcessName() != "data" ? samples[sam].getProcessName() : "nonpromptData");

          }
          if(debug){
              if(leptonSelection == 3 && nJLoc == 2 && nBLoc  == 1)
                cout << "weights are (lepSF/pu/btagL/btagBC) : " << lepSF << " " << puW << " " << btagL << " " << btagC*btagB << endl;
          }

          //finish = std::chrono::high_resolution_clock::now();
          //elapsed = finish - start;
          //std::cout << "time needed to estimate all sf and deviations: " << elapsed.count() << std::endl;

          for ( const auto & varPair : fillVar) {
            TString name = varPair.first;
            double var = varPair.second;
            int index = figNames[name].index;
            double varMaxOfVar = figNames[name].varMax;

            double varJECUp = fillVarJecUp.at(name);
            double varJECDown = fillVarJecDown.at(name);

            double varJERUp = fillVarJerUp.at(name);
            double varJERDown = fillVarJerDw.at(name);

            // check if variable is in the list to print
            if(std::find(listToPrint[selection].begin(), listToPrint[selection].end(), name) == listToPrint[selection].end()) continue;
            distribs[index].vectorHisto[samCategory].Fill(TMath::Min(var,varMaxOfVar-0.1),weight);

            if((samples[sam].getProcessName()) != "data"){

                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 0, varMaxOfVar-0.1, weight * lepSFSystUp / lepSF);
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var, 0, varMaxOfVar-0.1, weight * lepSFSystDown / lepSF);
                
                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 1, varMaxOfVar-0.1, weight * lepSFStatUp / lepSF);
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var, 1, varMaxOfVar-0.1, weight * lepSFStatDown / lepSF);
                
                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 2, varMaxOfVar-0.1, weight * lepSFRecoUp / lepSF);
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  2, varMaxOfVar-0.1, weight * lepSFRecoDown / lepSF);
                
                //distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * jetPrefW);
                //distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * jetPrefW);

                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 3, varMaxOfVar-0.1, weight*puWUp/puW);
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  3, varMaxOfVar-0.1, weight*puWDown/puW);
                
                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 4, varMaxOfVar-0.1, weight*btagLUp/btagL);
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  4, varMaxOfVar-0.1, weight*btagLDown/btagL);
                
                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 5, varMaxOfVar-0.1, weight*btagCUp*btagBUp/btagC/btagB);
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  5, varMaxOfVar-0.1, weight*btagCDown*btagBDown/btagC/btagB);
                
                distribs[index].vectorHistoUncUp[samCategory].FillUnc(varJECUp, 6, varMaxOfVar-0.1, weight);
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(varJECDown,  6, varMaxOfVar-0.1, weight);
                
                distribs[index].vectorHistoUncUp[samCategory].FillUnc(varJERUp, 7, varMaxOfVar-0.1, weight);
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(varJERDown,  7, varMaxOfVar-0.1, weight);
                
                if(samples[sam].getProcessName() == "WZ" && nBLoc > 0){ // 8 % uncertainty for WZ + bb background in high nbjets categories
                  distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 8, varMaxOfVar-0.1, weight * 1.08);
                  distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  8, varMaxOfVar-0.1, weight * 0.92);
                }
                else{
                  distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 8, varMaxOfVar-0.1, weight);
                  distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  8, varMaxOfVar-0.1, weight);
                }

                if((samples[sam].getFileName().find("ST_tWll_") != std::string::npos || samples[sam].getFileName().find("TTWJetsToLNu") != std::string::npos || samples[sam].getFileName().find("tZq_ll") != std::string::npos || samples[sam].getFileName().find("TTZToLLNuNu_M-10_") != std::string::npos) && samples[sam].is2017()){ 
                    // 10, 12 - factor 4; 6 and 8 - factor 2; 8, 9 - JERUp factor 2, 6, 7 - down factor 2, 30th Aug Daniel said in ttX chat that for FSR factor sqrt(2) should be used, indeces 3 and 5
                    // as well from https://twiki.cern.ch/twiki/bin/view/CMS/TopSystematics#Factorization_and_renormalizatio there is an instruction to use envelope uncertrtainty, i.e. largest between ISR and FSR
                    // description of the order for the PS weights: https://twiki.cern.ch/twiki/bin/view/CMS/TopModGen#Event_Generation
                    distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 9, varMaxOfVar-0.1, weight*_psWeight[8]*sumSimulatedEventWeights/sumSimulatedEventWeightsISRScaleUp); 
                    distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  9, varMaxOfVar-0.1, weight*_psWeight[6]*sumSimulatedEventWeights/sumSimulatedEventWeightsISRScaleDown); 

                    distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 10, varMaxOfVar-0.1, weight*_psWeight[5]*sumSimulatedEventWeights/sumSimulatedEventWeightsFSRScaleUp); 
                    distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  10, varMaxOfVar-0.1, weight*_psWeight[3]*sumSimulatedEventWeights/sumSimulatedEventWeightsFSRScaleDown); 
                }
                /*
                else if(samples[sam].getFileName().find("TTZToLLNuNu_M-10_") != std::string::npos && samples[sam].is2016() && index == indexSRTTZ){ // apply same weights as in 2017
                    // here let's take for the moment largest deviation from unity in a particular SR (bin)
                    for(int bin = 0; bin < figNames[name].nBins; bin++){
                      double psWUp = ttZISRUpW[bin] > ttZFSRUpW[bin] ? ttZISRUpW[bin] : ttZFSRUpW[bin];
                      double psWDown = ttZISRDownW[bin] > ttZFSRDownW[bin] ? ttZISRDownW[bin] : ttZFSRDownW[bin];

                      distribs[dist].vectorHistoUncUp[samCategory][9].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight*psWUp); 
                      distribs[dist].vectorHistoUncDown[samCategory][9].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight*psWDown); 

                      distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight*psWUp); 
                      distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight*psWDown); 
                    }
                }
                */
                else{
                    distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 9, varMaxOfVar-0.1, weight);
                    distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  9, varMaxOfVar-0.1, weight);

                    distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 10, varMaxOfVar-0.1, weight);
                    distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  10, varMaxOfVar-0.1, weight);
                }

                // renormalization and factorization uncertainties
                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 11, varMaxOfVar-0.1, weight*_lheWeight[8]*sumSimulatedEventWeights/sumSimulatedEventWeightsScaleUp);
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  11, varMaxOfVar-0.1, weight*_lheWeight[4]*sumSimulatedEventWeights/sumSimulatedEventWeightsScaleDown);

                // PDF uncertainties (will be filled later)
                if(index == indexSR3L || index == indexSR4L || index == indexSRTTZ || index == indexSRWZCR || index == indexSRZZCR || index == indexSRTTCR){
                    for(int varPDF = 0; varPDF < 100; varPDF++){
                        distribs[index].vectorHistoPDF[samCategory].var[varPDF].Fill(var, weight*_lheWeight[9+varPDF]);
                    }
                }
                
                //  all dists apart from SRs dists will have 0 uncertainty
                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 12, varMaxOfVar-0.1, weight);
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  12, varMaxOfVar-0.1, weight);


                // normalization uncertianties are stored in arrays defined in interface/readTreeSync.h
                // changed from 8 to 10, added WZc and WZb
                for(int cat = 0; cat < 10; cat++){

                    if(cat == samCategory){
                        distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 13+cat, varMaxOfVar-0.1, weight*(1+uncOnNorm[samples[sam].getProcessName()]));
                        distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  13+cat, varMaxOfVar-0.1, weight*(1-uncOnNorm[samples[sam].getProcessName()]));
                    }
                    else{
                        distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 13+cat, varMaxOfVar-0.1, weight);
                        distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  13+cat, varMaxOfVar-0.1, weight);
                    }
                }

                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var,  23, varMaxOfVar-0.1, weight * 1.025); // 2.5% for lumi
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  23, varMaxOfVar-0.1, weight * 0.975);
                
                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var,  24, varMaxOfVar-0.1, weight * 1.01); // 1% for trigger 
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  24, varMaxOfVar-0.1, weight * 0.99);

                if(samples[sam].getProcessName() == "WZ" && nJLoc > 2 && leptonSelection == 3){ // 20 % uncertainty in tails of njets, namely in njets >= 3, only in 3L category
                    distribs[index].vectorHistoUncUp[samCategory].FillUnc(var,  25, varMaxOfVar-0.1, weight * 1.2);
                    distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  25, varMaxOfVar-0.1, weight * 0.8);
                }
                else{
                    distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 25, varMaxOfVar-0.1, weight);
                    distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  25, varMaxOfVar-0.1, weight);
                }


            }
            else if(samCategory == nonPromptSample && leptonSelection != 4){
                for(int cat = 0; cat < 26; cat++){
                    if(cat == 22){
                      distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, cat, varMaxOfVar-0.1, weight*1.3);
                      distribs[index].vectorHistoUncDown[samCategory].FillUnc(var, cat, varMaxOfVar-0.1, weight*0.7);
                    }
                    else{
                      distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, cat, varMaxOfVar-0.1, weight);
                      distribs[index].vectorHistoUncDown[samCategory].FillUnc(var, cat, varMaxOfVar-0.1, weight);
                    }
                }
            }
          }
			 // break out when the event was found
          if(debug) break;
      }

      // if any histogram filled with prompt leptons input has a negative number of events
      // we set the number of events to 0
      if(samCategory != nonPromptSample && samCategory != dataSample){
        for(auto & name : listToPrint[selection]){
            int dist = figNames[name].index;
            for(int b = 1; b < distribs[dist].vectorHisto[samCategory].GetNbinsX() + 1; ++b){
               if(distribs[dist].vectorHisto[samCategory].GetBinContent(b) < 0.) distribs[dist].vectorHisto[samCategory].SetBinContent(b, 0.);
                for(int cat = 0; cat < numberOfSyst; cat++){
                    if(distribs[dist].vectorHistoUncUp[samCategory].unc[cat].GetBinContent(b) < 0.) distribs[dist].vectorHistoUncUp[samCategory].unc[cat].SetBinContent(b, 0.);
                    if(distribs[dist].vectorHistoUncDown[samCategory].unc[cat].GetBinContent(b) < 0.) distribs[dist].vectorHistoUncDown[samCategory].unc[cat].SetBinContent(b, 0.);
                }
            }
        }
      }

      cout << endl;
      samCategory = processToCounterMap.at(samples[sam].getProcessName());
      cout << "Total number of events: " << distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[samCategory].Integral() << endl;
      if(leptonSelection != 4)
        cout << "Total number of events in non prompt category: " << distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[nonPromptSample].Integral() << endl;
      cout << endl;
  }

  // after prompt subtraction from non-prompt we set all the negative yields to 0
  for(auto & name : listToPrint[selection]){
        int dist = figNames[name].index;
        for(int b = 1; b < distribs[dist].vectorHisto[nonPromptSample].GetNbinsX() + 1; ++b){
               if(distribs[dist].vectorHisto[nonPromptSample].GetBinContent(b) < 0.) distribs[dist].vectorHisto[nonPromptSample].SetBinContent(b, 0.);
                for(int cat = 0; cat < numberOfSyst; cat++){
                    if(distribs[dist].vectorHistoUncUp[nonPromptSample].unc[cat].GetBinContent(b) < 0.) distribs[dist].vectorHistoUncUp[nonPromptSample].unc[cat].SetBinContent(b, 0.);
                    if(distribs[dist].vectorHistoUncDown[nonPromptSample].unc[cat].GetBinContent(b) < 0.) distribs[dist].vectorHistoUncDown[nonPromptSample].unc[cat].SetBinContent(b, 0.);
                }
        }
  }

  // this should be done to be fully correct in PDF, takes a lot of time, an effect estimated with ttZ sample, uncertainty on acceptance is under 1%, simply assign flat uncertainty of 1 % to all signal and bkg
  // for the moment draw this uncertainty only for SR and will propagate them to datacards 
  for(auto & name : listToPrint[selection]){
    int dist = figNames[name].index; 
    if(!(dist == indexSR3L || dist == indexSR4L || dist == indexSRTTZ || dist == indexSRWZCR || dist == indexSRZZCR || dist == indexSRTTCR)) continue;

    for(unsigned proc = 0; proc < processToCounterMap.size(); ++proc){
      if(proc == dataSample) continue;
      if(proc == nonPromptSample) continue;
      for(unsigned bin = 1; bin < (unsigned) distribs[dist].vectorHistoUncUp[proc].unc.at(pdfUncIndex).GetNbinsX() + 1; ++bin){
          double pdfVarRms = 0.;
          for(unsigned pdf = 0; pdf < 100; ++pdf){
              double variedBin = distribs[dist].vectorHistoPDF[proc].var[pdf].GetBinContent(bin);
              //variedBin *= crossSectionRatio[proc][pdf]; // here should be not proc, but sample index
              // commented at the moment, because here we are looping over list of processes 
              // crossSectionRatio is stored per sample, not per process
              // therefore not clear which process correction should be taken
              double diff = (  variedBin - distribs[dist].vectorHisto[proc].GetBinContent(bin) );
              pdfVarRms += diff * diff;
          }
          pdfVarRms = sqrt( 0.01 * pdfVarRms );
          cout << "pdf rms for bin " << bin << " in process " << namesOfTheProcesses[proc] << " is equal to " << pdfVarRms << endl;
          distribs[dist].vectorHistoUncUp[proc].unc.at(pdfUncIndex).SetBinContent(bin, distribs[dist].vectorHisto[proc].GetBinContent(bin) + pdfVarRms);
          distribs[dist].vectorHistoUncDown[proc].unc.at(pdfUncIndex).SetBinContent(bin, distribs[dist].vectorHisto[proc].GetBinContent(bin) - pdfVarRms);
      }
    }
  }

  // legend to print
  TLegend* mtleg = new TLegend(0.18,0.89,0.92,0.72); 
  mtleg->SetNColumns(3);
  if(selection == "ZZ" || selection == "ttZ4L")
    mtleg->SetNColumns(4);
  if(selection == "ttZ")
    mtleg->SetNColumns(4);
  mtleg->SetFillColor(0);
  mtleg->SetFillStyle(0);
  mtleg->SetBorderSize(0);
  mtleg->SetTextFont(42);
  mtleg->SetTextSize(0.06);

  // fill the legend with entries.
  mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[dataSample],"Data","ep"); //data

  std::vector<std::string> processOrder;
  // iterator here doesn't work, map is unordered
  //for(map<std::string, int>::const_iterator it = processToCounterMap.begin();it != processToCounterMap.end(); ++it){
  for(unsigned proc = 0; proc < processToCounterMap.size(); ++proc){
    std::pair<std::string, int> it(namesOfTheProcesses[proc], proc);
    std::cout << it.first << " " << it.second << std::endl;
    processOrder.push_back(it.first);

    if(it.first == "data") continue;
    if(it.first == "ttH") continue;

    if(selection == "ZZ" || selection == "ttZ4L"){
      if(it.first == "WZ") continue;
    }

    if(it.first == "nonpromptData")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it.second],"Nonprompt","f");
    else if(it.first == "Xgamma")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it.second],"X#gamma","f");
    else if(it.first == "rare")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it.second],"Rare","f");
    else if(it.first == "ttZ")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it.second],"t#bar{t}Z","f");
    else if(it.first == "ttW")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it.second],"t#bar{t}W","f");
    else if(it.first == "ttH")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it.second],"t#bar{t}H","f");
    else if(it.first == "ttX")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it.second],"t(#bar{t})X","f");
    else if(it.first == "WZ")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it.second],"WZ+light","f");
    else if(it.first == "WZc")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it.second],"WZ+c","f");
    else if(it.first == "WZb")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it.second],"WZ+b","f");
    else
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it.second],it.first.c_str(),"f");
    
  }

  // uncertainty band on the top plot
  TH1D* histStatAndSystBand = (TH1D*)distribs[figNames[listToPrint[selection].at(0)].index].stack.GetStack()->Last()->Clone("histStatAndSystBand");
  histStatAndSystBand ->SetFillStyle(3005);
  histStatAndSystBand ->SetLineColor(kGray+2);
  histStatAndSystBand ->SetFillColor(kGray+2);
  histStatAndSystBand ->SetMarkerStyle(1);
  mtleg->AddEntry(histStatAndSystBand, "Uncertainty" ,"f"); // "Total unc." in Willem's plots

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
  gSystem->Exec("mkdir -p plotsForSave/" + folderToStorePlots + processToStore);
  double scale_num = 1.6;
  
  TCanvas* plot[nVars];
  for(int i = 0; i < nVars; i++){
      plot[i] = new TCanvas(Form("plot_%d", i),"",500,500);
  }

  std::string crToPrint = selection;

  for(int varPlot = 0; varPlot < listToPrint[crToPrint].size(); varPlot++){

    TString figName = listToPrint[crToPrint].at(varPlot);
    int distIndex = figNames[figName].index;
    bool plotInLogFlag = false;
    bool normalizedToDataFlag = false;

    plot[varPlot]->cd();
    showHist(plot[varPlot],distribs[distIndex],figNames[figName], scale_num, mtleg, plotInLogFlag, normalizedToDataFlag, showLegendOption);
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + figName + ".pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + figName + ".png");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + figName + ".root");


    plot[varPlot]->cd();
    plotInLogFlag = true;
    showHist(plot[varPlot],distribs[distIndex],figNames[figName], scale_num, mtleg, plotInLogFlag, normalizedToDataFlag, showLegendOption);
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + figName + "Log.pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + figName + "Log.png");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + figName + "Log.root");
  }
  
  if(crToPrint == "ttZ"){
    fillTablesSRTTZ(distribs[indexSRTTZ], processOrder, "SRallTTZ", showLegendOption);
    //fillTablesForFlavour(distribs[indexFlavour3L4L], processOrder, "flavour3L4L", showLegendOption);
  }
  if(crToPrint == "ttZclean"){
    fillTablesForFlavour(distribs[indexFlavour3L4L], processOrder, "flavour3L4L", showLegendOption);
  }

  // no need to fill datacards when running over 2016 and 2017 together
  if(showLegendOption == 2) return;
  if(crToPrint == "ttZ3Lclean"){
    fillDatacards(distribs[indexSRttZcleanPTZ], processOrder, "SRttZCleanPTZ", (bool)showLegendOption);
    fillDatacards(distribs[indexSRttZcleanCosTheta], processOrder, "SRttZCleanCosTheta", (bool)showLegendOption);
  }
  if(crToPrint == "ttZ"){

    fillDatacards(distribs[indexSRTTZ8SR3L], processOrder, "SRTTZ8SR3L", (bool)showLegendOption); 
    fillDatacards(distribs[indexSR3L], processOrder, "SR3L", (bool)showLegendOption); 
    fillDatacards(distribs[indexSR4L], processOrder, "SR4L", (bool)showLegendOption); 
    fillDatacards(distribs[indexSRTTZ], processOrder, "SRallTTZ", (bool)showLegendOption); 
    fillDatacards(distribs[indexSR3L3m], processOrder, "SR3L3m", (bool)showLegendOption); 
    fillDatacards(distribs[indexSR3L2m1e], processOrder, "SR3L2m1e", (bool)showLegendOption); 
    fillDatacards(distribs[indexSR3L1m2e], processOrder, "SR3L1m2e", (bool)showLegendOption); 
    fillDatacards(distribs[indexSR3L3e], processOrder, "SR3L3e", (bool)showLegendOption); 
    fillDatacards(distribs[indexSRWZCR], processOrder, "SRWZCR", (bool)showLegendOption); 
    fillDatacards(distribs[indexSRZZCR], processOrder, "SRZZCR", (bool)showLegendOption); 
    fillDatacards(distribs[indexSRTTCR], processOrder, "SRTTCR", (bool)showLegendOption); 

    fillDatacards(distribs[indexLeadPt], processOrder, "ptlead", (bool)showLegendOption); 
    fillDatacards(distribs[indexTrailPt], processOrder, "trail", (bool)showLegendOption); 
  }
  return;
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
						  // argument is "selection:..." so you need to remove the first 10 signs
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
