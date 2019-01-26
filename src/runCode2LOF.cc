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
  leptonSelection = leptonSelectionAnalysis;
  initListsToPrint(selection);
  //Set CMS plotting style
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
  initdistribs(namesOfProcesses, selection);
  cout << "finished with initiating of histos"<< endl;
  //setLabelsForHistos();

  std::ofstream myfile;
  myfile.open("myevents.txt");

  PostFitScaler scaler("data/postFit/outputTTZ.txt");

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample();
      int samCategory = processIndex.at(samples[sam].getProcessName());

      Color_t color = assignColor(samples[sam].getProcessName());
      setStackColors(color, samCategory);

      //if(!(samples[sam].getFileName().find("ST_tWll_") != std::string::npos || samples[sam].getFileName().find("TTWJetsToLNu") != std::string::npos || samples[sam].getFileName().find("tZq_ll") != std::string::npos)) continue;
      //if(samples[sam].getProcessName() != "data" && samples[sam].getProcessName() != "nonpromptData" && samples[sam].getProcessName() != "Nonprompt") continue;
      //if(samples[sam].getProcessName() != "TTto1L") continue;
      //if(samples[sam].getProcessName() != "WJets") continue;
      //if(samples[sam].getProcessName() != "TTto1L" && samples[sam].getProcessName() != "WJets") continue;

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
          // here we should ask only for single electron trigger,
          // dilepton can bias muon selection
          if(!_passTrigger_e) continue;
          //if(!_passTrigger_m) continue;
          if(debug) cout << "met filers flag: " << _passMETFilters << endl;
          if(!_passMETFilters) continue;
          
          //if(it > 10000) break;
          //if(it > nEntries / 10) break;

          // lepton selection
          // ---------------------------------------------------------------------------
          std::vector<unsigned> indTight, indFake, indLoose, indOf2LonZ;
          //select leptons relative to the analysis
          leptonSelection = 3;
          const unsigned lCount = selectLep(indTight, leptonSelection);
          const unsigned lCountFake = selectFakeLep(indFake, leptonSelection);
          const unsigned lCountLoose = selectLooseLep(indLoose, leptonSelection);

          if(lCountLoose != 2) continue;
          std::vector<unsigned> ind = indLoose;

          int nLocEle = getElectronNumber(ind);
          if(nLocEle != 1) continue;
          
          // here let's select emu pair, e should pass tight select
          int muonIndex = _lFlavor[ind[0]] == 1 ? ind[0] : ind[1]; 
          int eleIndex = _lFlavor[ind[0]] == 0 ? ind[0] : ind[1]; 

          // here explicitely we put cuts on the tag lepton
          if(!lepIsGood(eleIndex, leptonSelection) || _lPt[eleIndex] < 40) continue;
          //if(!lepIsGoodFortZq(eleIndex) || _lPt[eleIndex] < 40) continue;
          
          //if(!lepIsGood(muonIndex, leptonSelection) || _lPt[muonIndex] < 40) continue;

          if(debug) cout << "invariant mass of any fake pair is below 20 GeV: " << invMassOfAny2Lbelow20GeV(ind) << endl;
          if(invMassOfAny2Lbelow20GeV(ind)) continue;
          
          // in this block we put cuts on the charge of the leptons
          // this is used only when np is calculated from MC 
          //if(_lCharge[ind[0]] == _lCharge[ind[1]]) continue;
          // also for the check to see the contribution in ss2l final state we'll ask for 2 leptons to be the same charge
          //if(_lCharge[ind[0]] != _lCharge[ind[1]]) continue;
          /*
          if(_lCharge[ind[0]] == _lCharge[ind[1]])
              samCategory = 2;
          if(_lCharge[ind[0]] != _lCharge[ind[1]])
              samCategory = dataSample;
          */

          // ---------------------------------------------------------------------------
          // consider only prompt leptons from the MC, all nonprompt should be taken into account by DD estimation
          bool allLeptonsArePrompt = true;
          
          if((samples[sam].getProcessName()) != "data" && (samples[sam].getProcessName()) != "nonpromptData" && (samples[sam].getProcessName()) != "chargeMisIDData")
            allLeptonsArePrompt = promptLeptons(ind);

          // here we implement data driven np estimation
          if(samples[sam].getProcessName() != "data" && !allLeptonsArePrompt) continue;
          if(debug) cout << "all leptons are prompt ? " << allLeptonsArePrompt << endl;
          
          // ---------------------------------------------------------------------------
          // lepton pt criteria
          if(!passPtCuts2LOF(ind)) continue;

          // ---------------------------------------------------------------------------
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
          
          if(nJLoc < 2) continue;
          if(nBLoc != 1) continue;
          
          // ---------------------------------------------------------------------------
          // weight estimation for event
          //start = std::chrono::high_resolution_clock::now();
          double lepSF = 1.;
          double lepSFSystUp = 1.; double lepSFSystDown = 1.; 
          double lepSFStatUp = 1.; double lepSFStatDown = 1.;
          double lepSFRecoUp = 1.; double lepSFRecoDown = 1.;

          double lepWOnlySystOnlyEle = 1.; double lepWOnlySystOnlyEleUp = 1.; double lepWOnlySystOnlyEleDown = 1.;
          double lepWOnlyStatOnlyEle = 1.; double lepWOnlyStatOnlyEleUp = 1.; double lepWOnlyStatOnlyEleDown = 1.;
          double lepWOnlyRecoOnlyEle = 1.; double lepWOnlyRecoOnlyEleUp = 1.; double lepWOnlyRecoOnlyEleDown = 1.;

          double puW = 1.; double puWUp = 1.; double puWDown = 1.;

          double btagL = 1.; double btagLUp = 1.; double btagLDown = 1.;
          double btagC = 1.; double btagCUp = 1.; double btagCDown = 1.;
          double btagB = 1.; double btagBUp = 1.; double btagBDown = 1.;

          double jetPrefW = 1.;

          int muonBin = 0;
          if(_lPt[muonIndex] > 15 && _lPt[muonIndex] < 20)
              muonBin = 1;
          if(_lPt[muonIndex] > 20 && _lPt[muonIndex] < 25)
              muonBin = 2;
          if(_lPt[muonIndex] > 25 && _lPt[muonIndex] < 30)
              muonBin = 3;
          if(_lPt[muonIndex] > 30 && _lPt[muonIndex] < 35)
              muonBin = 4;
          if(_lPt[muonIndex] > 35)
              muonBin = 5;

          if((samples[sam].getProcessName()) != "data"){
            weight *= sfWeight();

            double lepWOnlySyst = leptonWeightOnlySyst(0); double lepWOnlySystUp = leptonWeightOnlySyst(1); double lepWOnlySystDown = leptonWeightOnlySyst(2);
            double lepWOnlyStat = leptonWeightOnlyStat(0); double lepWOnlyStatUp = leptonWeightOnlyStat(1); double lepWOnlyStatDown = leptonWeightOnlyStat(2);
            double lepWOnlyReco = leptonWeightOnlyReco(0); double lepWOnlyRecoUp = leptonWeightOnlyReco(1); double lepWOnlyRecoDown = leptonWeightOnlyReco(2);

            //double lepWOnlySystOnlyMuon = leptonWeightOnlySyst(0, false, true); double lepWOnlySystOnlyMuonUp = leptonWeightOnlySyst(1, false, true); double lepWOnlySystOnlyMuonDown = leptonWeightOnlySyst(2, false, true);
            //double lepWOnlyStatOnlyMuon = leptonWeightOnlyStat(0, false, true); double lepWOnlyStatOnlyMuonUp = leptonWeightOnlyStat(1, false, true); double lepWOnlyStatOnlyMuonDown = leptonWeightOnlyStat(2, false, true);
            //double lepWOnlyRecoOnlyMuon = leptonWeightOnlyReco(0, false, true); double lepWOnlyRecoOnlyMuonUp = leptonWeightOnlyReco(1, false, true); double lepWOnlyRecoOnlyMuonDown = leptonWeightOnlyReco(2, false, true);

            lepWOnlySystOnlyEle = leptonWeightOnlySyst(0, true, false); lepWOnlySystOnlyEleUp = leptonWeightOnlySyst(1, true, false); lepWOnlySystOnlyEleDown = leptonWeightOnlySyst(2, true, false);
            lepWOnlyStatOnlyEle = leptonWeightOnlyStat(0, true, false); lepWOnlyStatOnlyEleUp = leptonWeightOnlyStat(1, true, false); lepWOnlyStatOnlyEleDown = leptonWeightOnlyStat(2, true, false);
            lepWOnlyRecoOnlyEle = leptonWeightOnlyReco(0, true, false); lepWOnlyRecoOnlyEleUp = leptonWeightOnlyReco(1, true, false); lepWOnlyRecoOnlyEleDown = leptonWeightOnlyReco(2, true, false);

            lepSF = lepWOnlySyst * lepWOnlyReco; // consider only one between syst and stat, central value is the same
            lepSFSystUp = lepWOnlySystUp * lepWOnlyReco; 
            lepSFSystDown = lepWOnlySystDown * lepWOnlyReco;
            lepSFStatUp = lepWOnlyStatUp * lepWOnlyReco; 
            lepSFStatDown = lepWOnlyStatDown * lepWOnlyReco; 
            lepSFRecoUp = lepWOnlySyst * lepWOnlyRecoUp; 
            lepSFRecoDown = lepWOnlySyst * lepWOnlyRecoDown; 

            // here from data we subtract nonprompt contribution, only can be used with MC nonprompt
            if(!allLeptonsArePrompt){
                samCategory = dataSample;
                weight *= -1;
            }
            
            // we'll take ss dilepton events from MC, this event should be subtracted from SS data events, we multiply them by OS / SS ratio measured in ttbar semileptonic MC, which is 1.62733 / 1.22763 according to Willem, 18 Oct 2018, Skype chat
            // these two numbers are quite magical, real number to be used is approximately 1.15
            // P.S. later this number is measured as a function of lepton pt, here we use the measured corrections, they vary from 1.3 for first pt bin to 1 which should be applied for all pts > 35 GeV 
            if(_lCharge[ind[0]] == _lCharge[ind[1]]){
              weight *= -1 * osToss[muonBin];
              samCategory = nonPromptSample;
            } 
          }
          // here we implement data driven method of nonprompt background estimation, 
          // we'll take ss dilepton events from data and multiply them by OS / SS ratio measured in ttbar semileptonic MC, which is 1.62733 / 1.22763 according to Willem, 18 Oct 2018, Skype chat
          // P.S. later this number is measured as a function of lepton pt, here we use the measured corrections, they vary from 1.3 for first pt bin to 1 which should be applied for all pts > 35 GeV 
          if(samples[sam].getProcessName() == "data"){
            if(_lCharge[ind[0]] == _lCharge[ind[1]]){
                weight *= osToss[muonBin];
                samCategory = nonPromptSample;
            } 
          }

          if(debug) cout << "weight of event is " << weight << endl;
          // ---------------------------------------------------------------------------

          vector<double> fillVar = {
                                   _lPt[muonIndex], _lPt[eleIndex], 
                                   _lEta[muonIndex], _lEta[eleIndex], 
                                   //_lPt[eleIndex], _lPt[muonIndex], 
                                   //_lEta[eleIndex], _lEta[muonIndex], 
                                   double(nJLoc), double(nBLoc), 
                                   _met, double(_nVertex), 
                                   HTLoc,
                                   lepIsGood(muonIndex, leptonSelection) ? _lPt[muonIndex] : 0.,
                                   lepIsGood(muonIndex, leptonSelection) ? _lEta[muonIndex] : -999.,
                                   //lepIsGood(eleIndex, leptonSelection) ? _lPt[eleIndex] : 0.,
                                   //lepIsGood(eleIndex, leptonSelection) ? _lEta[eleIndex] : -999.,
                                   //lepIsGoodFortZq(muonIndex) ? _lPt[muonIndex] : 0.,
                                   //lepIsGoodFortZq(muonIndex) ? _lEta[muonIndex] : -999.,
                                   _lPt[muonIndex], _lEta[muonIndex], lepIsGood(muonIndex, leptonSelection) ? _lPt[muonIndex] : 0., lepIsGood(muonIndex, leptonSelection) ? _lEta[muonIndex] : -999.,
                                   //_lPt[muonIndex], _lEta[muonIndex], lepIsGoodFortZq(muonIndex) ? _lPt[muonIndex] : 0., lepIsGoodFortZq(muonIndex) ? _lEta[muonIndex] : -999.,
                                   _lEta[muonIndex], lepIsGood(muonIndex, leptonSelection) ? _lEta[muonIndex] : -999., _lEta[muonIndex], lepIsGood(muonIndex, leptonSelection) ? _lEta[muonIndex] : -999.,
                                   //_lEta[muonIndex], lepIsGoodFortZq(muonIndex) ? _lEta[muonIndex] : -999., _lEta[muonIndex], lepIsGoodFortZq(muonIndex) ? _lEta[muonIndex] : -999.,
                                   };

          vector<double> fillVarJecUp = fillVar; fillVarJecUp[4] = double(nJLocUp); fillVarJecUp[5] = double(nBLocUp); fillVarJecUp[8] = HTLocJECUp;
          vector<double> fillVarJecDw = fillVar; fillVarJecDw[4] = double(nJLocDown); fillVarJecDw[5] = double(nBLocDown); fillVarJecDw[8] = HTLocJECDown;
          vector<double> fillVarJerUp = fillVarJecUp;
          vector<double> fillVarJerDw = fillVarJecDw;

          vector<TString> fncName = {"ptlead", "sublead",  
                                     "etaLead", "etaSubl",  
                                     "njets", "nbjets", 
                                     "met", "nPV", 
                                     "HT",
                                     "ptMuonPassedTight",
                                     "etaMuonPassedTight",
                                     "ptLepPassedLooseForEff", "etaLepPassedLooseForEff", "ptLepPassedTightForEff", "etaLepPassedTightForEff",
                                     "etaLepPassedLooseForEffLowPt", "etaLepPassedTightForEffLowPt", "etaLepPassedLooseForEffHighPt",       "etaLepPassedTightForEffHighPt"
                                   };
                                   
          for(int counter = 0; counter < fillVar.size(); counter++){
            int dist = figNames[fncName.at(counter)].index;
            if(std::find(listToPrint[selection].begin(), listToPrint[selection].end(), fncName[counter]) == listToPrint[selection].end()) continue;
            
            // this is needed because otherwise it's impossible to draw the tgraph, number events in numerator should be always lower or equal in each bin than in denominator
            // for tZq
            /*
            if((fncName.at(counter) == "ptLepPassedTightForEff" || fncName.at(counter) == "etaLepPassedTightForEff") && !lepIsGoodFortZq(muonIndex)) continue;
            if((fncName.at(counter) == "etaLepPassedTightForEffLowPt") && (!lepIsGoodFortZq(muonIndex) || _lPt[muonIndex] > 25)) continue;
            if((fncName.at(counter) == "etaLepPassedLooseForEffLowPt") && (_lPt[muonIndex] > 25)) continue;
            if((fncName.at(counter) == "etaLepPassedTightForEffHighPt") && (!lepIsGoodFortZq(muonIndex) || _lPt[muonIndex] < 25)) continue;
            if((fncName.at(counter) == "etaLepPassedLooseForEffHighPt") && (_lPt[muonIndex] < 25)) continue;
            */
            
            if((fncName.at(counter) == "ptLepPassedTightForEff" || fncName.at(counter) == "etaLepPassedTightForEff") && !lepIsGood(muonIndex, leptonSelection)) continue;
            if((fncName.at(counter) == "etaLepPassedTightForEffLowPt") && (!lepIsGood(muonIndex, leptonSelection) || _lPt[muonIndex] > 25)) continue;
            if((fncName.at(counter) == "etaLepPassedLooseForEffLowPt") && (_lPt[muonIndex] > 25)) continue;
            if((fncName.at(counter) == "etaLepPassedTightForEffHighPt") && (!lepIsGood(muonIndex, leptonSelection) || _lPt[muonIndex] < 25)) continue;
            if((fncName.at(counter) == "etaLepPassedLooseForEffHighPt") && (_lPt[muonIndex] < 25)) continue;
            if((fncName.at(counter) == "ptLepPassedTightForEff" || fncName.at(counter) == "etaLepPassedTightForEff") && !lepIsGood(muonIndex, leptonSelection)) continue;
            //if((fncName.at(counter) == "ptMuonPassedTight" || fncName.at(counter) == "etaMuonPassedTight") && !lepIsGood(eleIndex, leptonSelection)) continue;
            
            // here we fill out the variables, the distributions for efficiency measuremt should be filled with only muon SF
            if((fncName.at(counter) == "ptLepPassedTightForEff" || fncName.at(counter) == "etaLepPassedTightForEff" || fncName.at(counter) == "ptLepPassedLooseForEff" || fncName.at(counter) == "etaLepPassedLooseForEff" || fncName.at(counter) == "etaLepPassedLooseForEffLowPt" || fncName.at(counter) == "etaLepPassedLooseForEffHighPt" || fncName.at(counter) == "etaLepPassedTightForEffLowPt" || fncName.at(counter) == "etaLepPassedTightForEffHighPt") && samples[sam].getProcessName() != "data"){
                //weight /= lepWOnlySystOnlyEle * lepWOnlyRecoOnlyEle;
                lepSF /= lepWOnlySystOnlyEle * lepWOnlyRecoOnlyEle;
                lepSFSystUp /= lepWOnlySystOnlyEleUp * lepWOnlyRecoOnlyEle;
                lepSFSystDown /= lepWOnlySystOnlyEleDown * lepWOnlyRecoOnlyEle;
                lepSFStatUp /= lepWOnlyStatOnlyEleUp * lepWOnlyRecoOnlyEle;
                lepSFStatDown /= lepWOnlyStatOnlyEleDown * lepWOnlyRecoOnlyEle;
                lepSFRecoUp /= lepWOnlyStatOnlyEle * lepWOnlyRecoOnlyEleUp;
                lepSFRecoDown /= lepWOnlyStatOnlyEle * lepWOnlyRecoOnlyEleDown;
            }
            distribs[dist].vectorHisto[samCategory].Fill(TMath::Min(fillVar.at(counter),figNames[fncName.at(counter)].varMax-0.001),weight);

            if((samples[sam].getProcessName()) != "data"){

                // here we want to fill distribution which will be later used for systematic uncertainty measurement 
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(counter), 0, figNames[fncName.at(counter)].varMax-0.1, weight * lepSFSystUp / lepSF);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(counter), 0, figNames[fncName.at(counter)].varMax-0.1, weight * lepSFSystDown / lepSF);
                
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(counter), 1, figNames[fncName.at(counter)].varMax-0.1, weight * lepSFStatUp / lepSF);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(counter), 1, figNames[fncName.at(counter)].varMax-0.1, weight * lepSFStatDown / lepSF);
                
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(counter), 2, figNames[fncName.at(counter)].varMax-0.1, weight * lepSFRecoUp / lepSF);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(counter), 2, figNames[fncName.at(counter)].varMax-0.1, weight * lepSFRecoDown / lepSF);
            }
            else if(samCategory == nonPromptSample){
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(counter), 0, figNames[fncName.at(counter)].varMax-0.1, weight * 1.3);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(counter), 0, figNames[fncName.at(counter)].varMax-0.1, weight * 0.7);
                
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(counter), 1, figNames[fncName.at(counter)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(counter), 1, figNames[fncName.at(counter)].varMax-0.1, weight);
                
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(counter), 2, figNames[fncName.at(counter)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(counter), 2, figNames[fncName.at(counter)].varMax-0.1, weight);
                
            }
          }
      }

      cout << endl;
      samCategory = processIndex.at(samples[sam].getProcessName());
      cout << "Total number of events: " << distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[samCategory].Integral() << endl;
      //cout << "Total number of events in data after nonprompt subtraction: " << distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[dataSample].Integral() << endl;
      cout << "Total number of events in nonprompt after prompt subtraction: " << distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[nonPromptSample].Integral() << endl;
      cout << endl;
  }

  // legend to print
  TLegend* mtleg = new TLegend(0.15,0.89,0.95,0.72); 
  mtleg->SetNColumns(4);
  mtleg->SetFillColor(0);
  mtleg->SetFillStyle(0);
  mtleg->SetBorderSize(0);
  mtleg->SetTextFont(42);

  mtleg->AddEntry(&distribs[0].vectorHisto[dataSample],"Data","lep"); //data
  //mtleg->AddEntry(&distribs[0].vectorHisto[dataSample],"OS ttbar 1L","lep"); //data
  //mtleg->AddEntry(&distribs[0].vectorHisto[2],"SS ttbar 1L","f"); //data

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

    if(it->second == "nonpromptData")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it->first],"Nonprompt","f");
    else if(it->second == "Xgamma")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it->first],"X#gamma","f");
    else if(it->second == "rare")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it->first],"Rare","f");
    else if(it->second == "ttZ")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it->first],"t#bar{t}Z","f");
    else if(it->second == "ttW")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it->first],"t#bar{t}W","f");
    else if(it->second == "ttH")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it->first],"t#bar{t}H","f");
    else if(it->second == "ttX")
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it->first],"t(#bar{t})X","f");
    else
      mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[it->first],it->second.c_str(),"f");
    
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

  TCanvas* plotRatio = new TCanvas("plotRatio","",500,450);
  showHistEff(plotRatio, distribs[figNames["ptLepPassedLooseForEff"].index], distribs[figNames["ptLepPassedTightForEff"].index], nonPromptSample);
  plotRatio->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/muonEffRatio.pdf");
  plotRatio->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/muonEffRatio.png");
  plotRatio->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/muonEffRatio.root");
  
  TCanvas* plotRatioEta = new TCanvas("plotRatioEta","",500,450);
  showHistEff(plotRatioEta, distribs[figNames["etaLepPassedLooseForEff"].index], distribs[figNames["etaLepPassedTightForEff"].index], nonPromptSample);
  plotRatioEta->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/muonEffRatioEta.pdf");
  plotRatioEta->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/muonEffRatioEta.png");
  plotRatioEta->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/muonEffRatioEta.root");
  
  TCanvas* plotRatioEtaLowPt = new TCanvas("plotRatioEtaLowPt","",500,450);
  showHistEff(plotRatioEtaLowPt, distribs[figNames["etaLepPassedLooseForEffLowPt"].index], distribs[figNames["etaLepPassedTightForEffLowPt"].index], nonPromptSample);
  plotRatioEtaLowPt->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/muonEffRatioEtaLowPt.pdf");
  plotRatioEtaLowPt->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/muonEffRatioEtaLowPt.png");
  plotRatioEtaLowPt->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/muonEffRatioEtaLowPt.root");
  
  TCanvas* plotRatioEtaHighPt = new TCanvas("plotRatioEtaHighPt","",500,450);
  showHistEff(plotRatioEtaHighPt, distribs[figNames["etaLepPassedLooseForEffHighPt"].index], distribs[figNames["etaLepPassedTightForEffHighPt"].index], nonPromptSample);
  plotRatioEtaHighPt->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/muonEffRatioEtaHighPt.pdf");
  plotRatioEtaHighPt->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/muonEffRatioEtaHighPt.png");
  plotRatioEtaHighPt->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/muonEffRatioEtaHighPt.root");
  
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
