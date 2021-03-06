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

#include "../interface/kinematicTools.h"
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

  leptonSelection = 2;
  if(leptonSelection == 2){

    // for reading values from xml files
    //addVariablesToBDT(is2017);

    // for adding variables to the trees
    addBranchToNNTreeVariables();
    //addBranchToBDTTreeVariables();
  }


  //PostFitScaler scaler("data/postFit/outputTTZ.txt");

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample();
      int samCategory = processToCounterMap.at(samples[sam].getProcessName());

      Color_t color = assignColor(samples[sam].getProcessName());
      setStackColors(color, samCategory);

      //if(!(samples[sam].getFileName().find("ST_tWll_") != std::string::npos || samples[sam].getFileName().find("TTWJetsToLNu") != std::string::npos || samples[sam].getFileName().find("tZq_ll") != std::string::npos)) continue;
      //if(!(samples[sam].getFileName().find("ST_tWll_") != std::string::npos)) continue;
      //if(samples[sam].getProcessName() != "data" && samples[sam].getProcessName() != "nonpromptData" && samples[sam].getProcessName() != "Nonprompt") continue;
      //if(samples[sam].getProcessName() != "data" && samples[sam].getProcessName() != "WZ") continue;
      if(samples[sam].getProcessName() == "data") continue;

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
          //if(it > nEntries / 20) break;

          std::vector<unsigned> indTight, indFake, indOf2LonZ;
          //select leptons relative to the analysis
          const unsigned lCount = selectLep(indTight, leptonSelection);
          const unsigned lCountFake = selectFakeLep(indFake, leptonSelection);
          if(lCount != 2) continue;

          /*
          std::vector<unsigned> indTight4L, indLoose4L;
          const unsigned lCount4L = selectLep(indTight4L, 4);
          const unsigned lCount4LLoose = selectFakeLep(indLoose4L, 4);

          // discard heavy flavour resonances
          if(debug) cout << "number of ttZ3L tight and fo leptons: " << lCount << " " << lCountFake << endl;
          if(debug) cout << "number of ttZ4L tight leptons: " << lCount4L << endl;
          */

          // selection of category for the event
          // 2L: possible contribution from TT, TF and FF; TTF is vetoed
          // 3L: TTT, TTF, TFF, FFF; for TTTF should consider if event pass 4L TTTT criteria
          // 4L: TTTT only combination is possible
          
          std::vector<unsigned> ind;
          ind = indTight;
          
          if(debug) cout << "invariant mass of any fake pair is below 12 GeV: " << invMassOfAny2Lbelow12GeV(ind) << endl;
          if(invMassOfAny2Lbelow12GeV(ind)) continue; 

          // consider only prompt leptons from the MC, all nonprompt should be taken into account by DD estimation
          bool allLeptonsArePrompt = true;
          
          if((samples[sam].getProcessName()) != "data" && (samples[sam].getProcessName()) != "nonpromptData" && (samples[sam].getProcessName()) != "chargeMisIDData")
            allLeptonsArePrompt = promptLeptons(ind);

          if(debug) cout << "all leptons are prompt ? " << allLeptonsArePrompt << endl;
          
          if((samples[sam].getProcessName()) == "chargeMisID" && !allLeptonsArePrompt) continue;
          if((samples[sam].getProcessName()) == "Nonprompt" && allLeptonsArePrompt) continue; // works just for MC

          if(leptonSelection == 2){
            if(((samples[sam].getProcessName()) == "ttW" || (samples[sam].getProcessName()) == "ttH" || (samples[sam].getProcessName()) == "ttZ" || (samples[sam].getProcessName()) == "ttX" 
                                                       || (samples[sam].getProcessName()) == "WZ" || (samples[sam].getProcessName()) == "Xgamma"  || (samples[sam].getProcessName()) == "ZZ" 
                                                       || (samples[sam].getProcessName()) == "rare") && !allLeptonsArePrompt) continue;
          }

          int nLocEle = getElectronNumber(ind);
          //if(nLocEle != 3) continue;

          // lepton pt criteria
          if(leptonSelection == 2)
            if(!passPtCuts2L(ind)) continue;

          if(leptonSelection == 2 && samCategory != CMIDSample)
            if(_lCharge[ind.at(0)] * _lCharge[ind.at(1)] < 0) continue;

          // select here jets, bjets, delta from M of Z boson, HT
          std::vector<unsigned> indJets;
          std::vector<unsigned> indJetsAllEta;
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
          int nJLocAllEta = nJets(0, true, indJetsAllEta, samples[sam].is2017(), 5.0);
          //if(isForwardJetPresent(samples[sam].is2017())) continue;

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
          
          double mt1 = 9999;

          double mvaVL = 0;
          double mvaVLJECUp = 0;
          double mvaVLJECDown = 0;

          //make lorentzvectors for leptons
          TLorentzVector lepV[(const unsigned int) _nEle + _nMu];
          for(unsigned l = 0; l < _nEle + _nMu; ++l) lepV[l].SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);

          TLorentzVector jetV[(const unsigned int) _nJets];
          for(unsigned j = 0; j < _nJets; ++j) jetV[j].SetPtEtaPhiE(_jetPt[j], _jetEta[j], _jetPhi[j], _jetE[j]);

          std::vector<unsigned> lepVecInd;
          for(unsigned l = 0; l < lCount; ++l) lepVecInd.push_back(ind[l]); 

          std::vector<unsigned> jetVecInd;
          for(unsigned l = 0; l < nJLocAllEta; ++l) jetVecInd.push_back(indJetsAllEta[l]); 

          std::vector<unsigned> bjetVecInd;
          for(unsigned l = 0; l < nBLoc; ++l) bjetVecInd.push_back(indBJets[l]); 

          // lepton + jet
          minMLeptonJet = kinematics::minMass(lepV, lepVecInd, jetV, jetVecInd);
          maxMLeptonJet = kinematics::maxMass(lepV, lepVecInd, jetV, jetVecInd);
          minDeltaRLeptonJet = kinematics::minDeltaR(lepV, lepVecInd, jetV, jetVecInd);
          maxDeltaRLeptonJet = kinematics::maxDeltaR(lepV, lepVecInd, jetV, jetVecInd);
          minDeltaPhiLeptonJet = kinematics::minDeltaR(lepV, lepVecInd, jetV, jetVecInd);
          maxDeltaPhiLeptonJet = kinematics::maxDeltaR(lepV, lepVecInd, jetV, jetVecInd);
          minpTLeptonJet = kinematics::minPT(lepV, lepVecInd, jetV, jetVecInd);
          maxpTLeptonJet = kinematics::maxPT(lepV, lepVecInd, jetV, jetVecInd);

          //lepton bjet
          minMLeptonbJet = kinematics::minMass(lepV, lepVecInd, jetV, bjetVecInd);
          maxMLeptonbJet = kinematics::maxMass(lepV, lepVecInd, jetV, bjetVecInd);
          minDeltaRLeptonbJet = kinematics::minDeltaR(lepV, lepVecInd, jetV, bjetVecInd);
          maxDeltaRLeptonbJet = kinematics::maxDeltaR(lepV, lepVecInd, jetV, bjetVecInd);
          minDeltaPhiLeptonbJet = kinematics::minDeltaR(lepV, lepVecInd, jetV, bjetVecInd);
          maxDeltaPhiLeptonbJet = kinematics::maxDeltaR(lepV, lepVecInd, jetV, bjetVecInd);
          minpTLeptonbJet = kinematics::minPT(lepV, lepVecInd, jetV, bjetVecInd);
          maxpTLeptonbJet = kinematics::maxPT(lepV, lepVecInd, jetV, bjetVecInd);

          //jet jet
          minMJetJet = kinematics::minMass(jetV, jetVecInd);
          maxMJetJet = kinematics::maxMass(jetV, jetVecInd);
          minDeltaRJetJet = kinematics::minDeltaR(jetV, jetVecInd);
          maxDeltaRJetJet = kinematics::maxDeltaR(jetV, jetVecInd);
          minDeltaPhiJetJet = kinematics::minDeltaPhi(jetV, jetVecInd);
          maxDeltaPhiJetJet = kinematics::maxDeltaPhi(jetV, jetVecInd);
          minpTJetJet = kinematics::minPT(jetV, jetVecInd);
          maxpTJetJet = kinematics::maxPT(jetV, jetVecInd);

          //set met index
          std::vector<unsigned> metIndex = {0};

          //make met vector 
          TLorentzVector met;
          met.SetPtEtaPhiE(_met, 0, _metPhi, _met);

          //lepton + MET 
          minDeltaPhiLeptonMET = kinematics::minDeltaPhi(lepV, lepVecInd, &met, metIndex);
          maxDeltaPhiLeptonMET = kinematics::maxDeltaPhi(lepV, lepVecInd, &met, metIndex);
          minmTLeptonMET = kinematics::minMT(lepV, lepVecInd, &met, metIndex);
          maxmTLeptonMET = kinematics::maxMT(lepV, lepVecInd, &met, metIndex);
          minpTLeptonMET = kinematics::minPT(lepV, lepVecInd, &met, metIndex);
          maxpTLeptonMET = kinematics::maxPT(lepV, lepVecInd, &met, metIndex);

          //jet + MET
          minDeltaPhiJetMET = kinematics::minDeltaPhi(jetV, jetVecInd, &met, metIndex);
          maxDeltaPhiJetMET = kinematics::maxDeltaPhi(jetV, jetVecInd, &met, metIndex);
          minmTJetMET = kinematics::minMT(jetV, jetVecInd, &met, metIndex);
          maxmTJetMET = kinematics::maxMT(jetV, jetVecInd, &met, metIndex);
          minpTJetMET = kinematics::minPT(jetV, jetVecInd, &met, metIndex);
          maxpTJetMET = kinematics::maxPT(jetV, jetVecInd, &met, metIndex);

          //bjet + MET 
          minDeltaPhiBJetMET = kinematics::minDeltaPhi(jetV, bjetVecInd, &met, metIndex);
          maxDeltaPhiBJetMET = kinematics::maxDeltaPhi(jetV, bjetVecInd, &met, metIndex);
          minmTBJetMET = kinematics::minMT(jetV, bjetVecInd, &met, metIndex);
          maxmTBJetMET = kinematics::maxMT(jetV, bjetVecInd, &met, metIndex);
          minpTBJetMET = kinematics::minPT(jetV, bjetVecInd, &met, metIndex);
          maxpTBJetMET = kinematics::maxPT(jetV, bjetVecInd, &met, metIndex);

          if(leptonSelection == 2){

            if(selection == "ttW" && !pass2Lpreselection(nJLoc, nBLoc, ind, _met, nLocEle)) continue;
            if(selection == "ttWclean" && !pass2Lcleanpreselection(nJLoc, nBLoc, ind, _met, nLocEle)) continue;

            //if(nJLoc > 6) continue;

            _jetPt1 = 0.;
            _jetEta1 = 0.;
            _jetPhi1 = 0.;
            _jetE1 = 0.;
            _jetCSV1 = 0.;

            _jetPt2 = 0.;
            _jetEta2 = 0.;
            _jetPhi2 = 0.;
            _jetE2 = 0.;
            _jetCSV2 = 0.;

            _jetPt3 = 0.;
            _jetEta3 = 0.;
            _jetPhi3 = 0.;
            _jetE3 = 0.;
            _jetCSV3 = 0.;

            _jetPt4 = 0.;
            _jetEta4 = 0.;
            _jetPhi4 = 0.;
            _jetE4 = 0.;
            _jetCSV4 = 0.;

            _jetPt5 = 0.;
            _jetEta5 = 0.;
            _jetPhi5 = 0.;
            _jetE5 = 0.;
            _jetCSV5 = 0.;

            _jetPt6 = 0.;
            _jetEta6 = 0.;
            _jetPhi6 = 0.;
            _jetE6 = 0.;
            _jetCSV6 = 0.;

            if(nJLoc > 0){
                _jetPt1 = _jetPt[indJetsAllEta[0]];
                _jetEta1 = _jetEta[indJetsAllEta[0]];
                _jetPhi1 = _jetPhi[indJetsAllEta[0]];
                _jetE1 = _jetE[indJetsAllEta[0]];
                _jetCSV1 = fabs(_jetEta[indJetsAllEta[0]]) < 2.4 ? (_closestJetDeepCsv_b[indJetsAllEta[0]] + _closestJetDeepCsv_bb[indJetsAllEta[0]] == -2 ? 0. : _closestJetDeepCsv_b[indJetsAllEta[0]] + _closestJetDeepCsv_bb[indJetsAllEta[0]]) : 0.;
                if(nJLoc > 1){
                    _jetPt2 = _jetPt[indJetsAllEta[1]];
                    _jetEta2 = _jetEta[indJetsAllEta[1]];
                    _jetPhi2 = _jetPhi[indJetsAllEta[1]];
                    _jetE2 = _jetE[indJetsAllEta[1]];
                    _jetCSV2 = fabs(_jetEta[indJetsAllEta[1]]) < 2.4 ? (_closestJetDeepCsv_b[indJetsAllEta[1]] + _closestJetDeepCsv_bb[indJetsAllEta[1]] == -2 ? 0. : _closestJetDeepCsv_b[indJetsAllEta[1]] + _closestJetDeepCsv_bb[indJetsAllEta[1]]) : 0.;
                    if(nJLoc > 2){
                        _jetPt3 = _jetPt[indJetsAllEta[2]];
                        _jetEta3 = _jetEta[indJetsAllEta[2]];
                        _jetPhi3 = _jetPhi[indJetsAllEta[2]];
                        _jetE3 = _jetE[indJetsAllEta[2]];
                        _jetCSV3 = fabs(_jetEta[indJetsAllEta[2]]) < 2.4 ? (_closestJetDeepCsv_b[indJetsAllEta[2]] + _closestJetDeepCsv_bb[indJetsAllEta[2]] == -2 ? 0. : _closestJetDeepCsv_b[indJetsAllEta[2]] + _closestJetDeepCsv_bb[indJetsAllEta[2]]) : 0.;
                        if(nJLoc > 3){
                            _jetPt4 = _jetPt[indJetsAllEta[3]];
                            _jetEta4 = _jetEta[indJetsAllEta[3]];
                            _jetPhi4 = _jetPhi[indJetsAllEta[3]];
                            _jetE4 = _jetE[indJetsAllEta[3]];
                            _jetCSV4 = fabs(_jetEta[indJetsAllEta[3]]) < 2.4 ? (_closestJetDeepCsv_b[indJetsAllEta[3]] + _closestJetDeepCsv_bb[indJetsAllEta[3]] == -2 ? 0. : _closestJetDeepCsv_b[indJetsAllEta[3]] + _closestJetDeepCsv_bb[indJetsAllEta[3]]) : 0.;
                            if(nJLoc > 4){
                                _jetPt5 = _jetPt[indJetsAllEta[4]];
                                _jetEta5 = _jetEta[indJetsAllEta[4]];
                                _jetPhi5 = _jetPhi[indJetsAllEta[4]];
                                _jetE5 = _jetE[indJetsAllEta[4]];
                                _jetCSV5 = fabs(_jetEta[indJetsAllEta[4]]) < 2.4 ? (_closestJetDeepCsv_b[indJetsAllEta[4]] + _closestJetDeepCsv_bb[indJetsAllEta[4]] == -2 ? 0. : _closestJetDeepCsv_b[indJetsAllEta[4]] + _closestJetDeepCsv_bb[indJetsAllEta[4]]) : 0.;
                                if(nJLoc > 5){
                                    _jetPt6 = _jetPt[indJetsAllEta[5]];
                                    _jetEta6 = _jetEta[indJetsAllEta[5]];
                                    _jetPhi6 = _jetPhi[indJetsAllEta[5]];
                                    _jetE6 = _jetE[indJetsAllEta[5]];
                                    _jetCSV6 = fabs(_jetEta[indJetsAllEta[5]]) < 2.4 ? (_closestJetDeepCsv_b[indJetsAllEta[5]] + _closestJetDeepCsv_bb[indJetsAllEta[5]] == -2 ? 0. : _closestJetDeepCsv_b[indJetsAllEta[5]] + _closestJetDeepCsv_bb[indJetsAllEta[5]]) : 0.;
                                }
                            }
                        }
                    }
                }
            }

            _lepPt1 = _lPt[ind[0]];
            _lepEta1 = _lEta[ind[0]];
            _lepPhi1 = _lPhi[ind[0]];
            _lepE1 = _lE[ind[0]];
            _lepCharge1 = _lCharge[ind[0]];

            _lepPt2 = _lPt[ind[1]];
            _lepEta2 = _lEta[ind[1]];
            _lepPhi2 = _lPhi[ind[1]];
            _lepE2 = _lE[ind[1]];
            _lepCharge2 = _lCharge[ind[1]];

            _metPt1 = _met;
            _metEta1 = 0;
            _metPhi1 = _metPhi;
            _metE1 = _met;
            
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

            _weightEventInTree = weight;
            if((samples[sam].getProcessName()) == "ttW")
              signalTree->Fill();
            //else if(((samples[sam].getProcessName() == "Nonprompt" || samples[sam].getProcessName() == "chargeMisID") && samCategory != nonPromptSample) || ((samples[sam].getProcessName() == "data") && samCategory == nonPromptSample))
            else
              bkgTree->Fill();

          }

          // weight estimation for event
          //auto start = std::chrono::high_resolution_clock::now();
          
          /*
          if((samples[sam].getProcessName()) != "data" && (samples[sam].getProcessName()) != "nonpromptData")
            weight *= sfWeight();
          if(samples[sam].getProcessName() == "data" && samCategory == nonPromptSample && leptonSelection == 3)
            weight *= fakeRateWeight();
          */

          //auto finish = std::chrono::high_resolution_clock::now();
          //std::chrono::duration<double> elapsed = finish - start;
          //std::cout << "time needed to estimate event weight: " << elapsed.count() << std::endl;

          if(debug) cout << "weight of event is " << weight << endl;

          int mvaValueRegion = 0;

          vector<double> fillVar = {ptCorrV[0].first, ptCorrV[1].first, leptonSelection > 2 ? ptCorrV[2].first : 0., leptonSelection > 3 ? ptCorrV[3].first : 0.,
                                   mt1, double(nJLoc), double(nBLoc), (_lCharge[ind.at(0)] == 1 ?  mvaVL : -999),
                                   // currently here we will have ttZ3L and ttZ4L categories
                                   (leptonSelection == 3 ? SRID3L(nJLoc, nBLoc, dMZ) : -999),
                                   (leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLoc) ? SRID4L(nJLoc, nBLoc) : -999),
                                   (leptonSelection != 4 ? mll:mll1stpair),ptZ,ptNonZ, (nLocEle == 3?mll:-999.), (nLocEle==2?mll:-999.), (nLocEle==1?mll:-999.), (nLocEle==0? mll: -999.),
                                   _met, minDeltaR, minDeltaRlead, mtHighest, mtLowest, leadingJetPt, trailJetPt, 0., double(_nVertex), mlll,
                                   _lEta[ptCorrV[0].second], _lEta[ptCorrV[1].second], (leptonSelection > 2 ? _lEta[ptCorrV[2].second] : -999.), (leptonSelection > 3 ? _lEta[ptCorrV[3].second] : -999.),
                                   (nLocEle == 3?mt1:-999.), (nLocEle==2?mt1:-999.), (nLocEle==1?mt1:-999.), (nLocEle==0? mt1: -999.),
                                   cosTSt, mll_ss, double(chargeOfLeptons), ll_deltaR, mt2ll_ss,
                                   (_lCharge[ind.at(0)] == -1 ?  mvaVL : -999), HTLoc,
                                   SRIDTTZ(ind, indOf2LonZ, nJLoc, nBLoc, dMZ, mlll), SRIDWZCR(nJLoc, nBLoc, dMZ), SRIDZZCR(ind, indOf2LonZ, nJLoc, nBLoc), SRIDTTCR(nJLoc, nBLoc, dMZ, mlll),
                                   (leptonSelection == 3 && nJLoc > 2 && nBLoc > 0 ? SRIDPTZ(ptZ) : -999), (leptonSelection == 3 && nJLoc > 2 && nBLoc > 0 ? SRIDCosTheta(cosTSt) : -999),
                                   (leptonSelection == 3 ? flavourCategory3L(nLocEle) : -999),
                                   (leptonSelection == 4 ? flavourCategory4L(nLocEle) : -999),
                                   (leptonSelection == 4 ? flavourCategory4LZZ(nLocEle) : -999),
                                   (leptonSelection == 3 && nLocEle == 0 ? SRID3L(nJLoc, nBLoc, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 1 ? SRID3L(nJLoc, nBLoc, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 2 ? SRID3L(nJLoc, nBLoc, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 3 ? SRID3L(nJLoc, nBLoc, dMZ) : -999),
                                   (passTTZSRSelection(ind, indOf2LonZ, nJLoc, nBLoc, dMZ) ? flavourCategory3L4L(leptonSelection, nLocEle) : -999),
                                   (leptonSelection == 3 && nBLoc > 0 ? SRID8SR3L(nJLoc, nBLoc, dMZ) : -999),
                                   minMLeptonJet, maxMLeptonJet, minDeltaRLeptonJet, maxDeltaRLeptonJet, minDeltaPhiLeptonJet, maxDeltaPhiLeptonJet, minpTLeptonJet, maxpTLeptonJet,
                                   minMLeptonbJet, maxMLeptonbJet, minDeltaRLeptonbJet, maxDeltaRLeptonbJet, minDeltaPhiLeptonbJet, maxDeltaPhiLeptonbJet, minpTLeptonbJet, maxpTLeptonbJet,
                                   minMJetJet, maxMJetJet, minDeltaRJetJet, maxDeltaRJetJet, minDeltaPhiJetJet, maxDeltaPhiJetJet, minpTJetJet, maxpTJetJet,
                                   minDeltaPhiLeptonMET, maxDeltaPhiLeptonMET, minmTLeptonMET, maxmTLeptonMET, minpTLeptonMET, maxpTLeptonMET,
                                   minDeltaPhiJetMET, maxDeltaPhiJetMET, minmTJetMET, maxmTJetMET, minpTJetMET, maxpTJetMET,
                                   minDeltaPhiBJetMET, maxDeltaPhiBJetMET, minmTBJetMET, maxmTBJetMET, minpTBJetMET, maxpTBJetMET
                                   };

          /*
          vector<double> fillVarJecUp = {ptCorrV[0].first, ptCorrV[1].first, leptonSelection > 2 ? ptCorrV[2].first : 0., leptonSelection > 3 ? ptCorrV[3].first : 0.,
                                   mt1, double(nJLocUp), double(nBLocUp), (_lCharge[ind.at(0)] == 1 ?  mvaVLJECUp : -999),
                                   (leptonSelection == 3 ? SRID3L(nJLocUp, nBLocUp, dMZ) : -999),
                                   (leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLocUp)? SRID4L(nJLocUp, nBLocUp) : -999),
                                   (leptonSelection != 4 ? mll:mll1stpair),ptZ,ptNonZ, (nLocEle == 3?mll:-999.), (nLocEle==2?mll:-999.), (nLocEle==1?mll:-999.), (nLocEle==0? mll: -999.),
                                   _met, minDeltaR, minDeltaRlead, mtHighest, mtLowest, leadingJetPt, trailJetPt, 0., double(_nVertex), mlll,
                                   _lEta[ptCorrV[0].second], _lEta[ptCorrV[1].second], (leptonSelection > 2 ? _lEta[ptCorrV[2].second] : -999.), (leptonSelection > 3 ? _lEta[ptCorrV[3].second] : -999.),
                                   (nLocEle == 3?mt1:-999.), (nLocEle==2?mt1:-999.), (nLocEle==1?mt1:-999.), (nLocEle==0? mt1: -999.),
                                   cosTSt, mll_ss, double(chargeOfLeptons), ll_deltaR, mt2ll_ss,
                                   (_lCharge[ind.at(0)] == -1 ?  mvaVLJECUp : -999), HTLocJECUp,
                                   SRIDTTZ(ind, indOf2LonZ, nJLocUp, nBLocUp, dMZ, mlll), SRIDWZCR(nJLocUp, nBLocUp, dMZ), SRIDZZCR(ind, indOf2LonZ, nJLocUp, nBLocUp), SRIDTTCR(nJLocUp, nBLocUp, dMZ, mlll),
                                   (leptonSelection == 3 && nJLocUp > 2 && nBLocUp > 0 ? SRIDPTZ(ptZ) : -999), (leptonSelection == 3 && nJLocUp > 2 && nBLocUp > 0 ? SRIDCosTheta(cosTSt) : -999),
                                   (leptonSelection == 3 ? flavourCategory3L(nLocEle) : -999),
                                   (leptonSelection == 4 ? flavourCategory4L(nLocEle) : -999),
                                   (leptonSelection == 4 ? flavourCategory4LZZ(nLocEle) : -999),
                                   (leptonSelection == 3 && nLocEle == 0 ? SRID3L(nJLocUp, nBLocUp, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 1 ? SRID3L(nJLocUp, nBLocUp, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 2 ? SRID3L(nJLocUp, nBLocUp, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 3 ? SRID3L(nJLocUp, nBLocUp, dMZ) : -999),
                                   (passTTZSRSelection(ind, indOf2LonZ, nJLocUp, nBLocUp, dMZ) ? flavourCategory3L4L(leptonSelection, nLocEle) : -999),
                                   (leptonSelection == 3 && nBLocUp > 0 ? SRID8SR3L(nJLocUp, nBLocUp, dMZ) : -999),
                                   };

          vector<double> fillVarJecDw = {ptCorrV[0].first, ptCorrV[1].first, leptonSelection > 2 ? ptCorrV[2].first : 0., leptonSelection > 3 ? ptCorrV[3].first : 0.,
                                   mt1, double(nJLocDown), double(nBLocDown), (_lCharge[ind.at(0)] == 1 ?  mvaVLJECDown : -999),
                                   (leptonSelection == 3 ? SRID3L(nJLocDown, nBLocDown, dMZ) : -999),
                                   (leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLocDown)? SRID4L(nJLocDown, nBLocDown) : -999),
                                   (leptonSelection != 4 ? mll:mll1stpair),ptZ,ptNonZ, (nLocEle == 3?mll:-999.), (nLocEle==2?mll:-999.), (nLocEle==1?mll:-999.), (nLocEle==0? mll: -999.),
                                   _met, minDeltaR, minDeltaRlead, mtHighest, mtLowest, leadingJetPt, trailJetPt, 0., double(_nVertex), mlll,
                                   _lEta[ptCorrV[0].second], _lEta[ptCorrV[1].second], (leptonSelection > 2 ? _lEta[ptCorrV[2].second] : -999.), (leptonSelection > 3 ? _lEta[ptCorrV[3].second] : -999.),
                                   (nLocEle == 3?mt1:-999.), (nLocEle==2?mt1:-999.), (nLocEle==1?mt1:-999.), (nLocEle==0? mt1: -999.),
                                   cosTSt, mll_ss, double(chargeOfLeptons), ll_deltaR, mt2ll_ss,
                                   (_lCharge[ind.at(0)] == -1 ?  mvaVLJECDown : -999), HTLocJECDown,
                                   SRIDTTZ(ind, indOf2LonZ, nJLocDown, nBLocDown, dMZ, mlll), SRIDWZCR(nJLocDown, nBLocDown, dMZ), SRIDZZCR(ind, indOf2LonZ, nJLocDown, nBLocDown), SRIDTTCR(nJLocDown, nBLocDown, dMZ, mlll),
                                   (leptonSelection == 3 && nJLocDown > 2 && nBLocDown > 0 ? SRIDPTZ(ptZ) : -999), (leptonSelection == 3 && nJLocDown > 2 && nBLocDown > 0 ? SRIDCosTheta(cosTSt) : -999),
                                   (leptonSelection == 3 ? flavourCategory3L(nLocEle) : -999),
                                   (leptonSelection == 4 ? flavourCategory4L(nLocEle) : -999),
                                   (leptonSelection == 4 ? flavourCategory4LZZ(nLocEle) : -999),
                                   (leptonSelection == 3 && nLocEle == 0 ? SRID3L(nJLocDown, nBLocDown, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 1 ? SRID3L(nJLocDown, nBLocDown, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 2 ? SRID3L(nJLocDown, nBLocDown, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 3 ? SRID3L(nJLocDown, nBLocDown, dMZ) : -999),
                                   (passTTZSRSelection(ind, indOf2LonZ, nJLocDown, nBLocDown, dMZ) ? flavourCategory3L4L(leptonSelection, nLocEle) : -999),
                                   (leptonSelection == 3 && nBLocDown > 0 ? SRID8SR3L(nJLocDown, nBLocDown, dMZ) : -999),
                                   };

          vector<double> fillVarJerUp = {ptCorrV[0].first, ptCorrV[1].first, leptonSelection > 2 ? ptCorrV[2].first : 0., leptonSelection > 3 ? ptCorrV[3].first : 0.,
                                   mt1, double(nJLocJERUp), double(nBLocJERUp), (_lCharge[ind.at(0)] == 1 ?  mvaVLJECUp : -999),
                                   (leptonSelection == 3 ? SRID3L(nJLocJERUp, nBLocJERUp, dMZ) : -999),
                                   (leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLocJERUp)? SRID4L(nJLocJERUp, nBLocJERUp) : -999),
                                   (leptonSelection != 4 ? mll:mll1stpair),ptZ,ptNonZ, (nLocEle == 3?mll:-999.), (nLocEle==2?mll:-999.), (nLocEle==1?mll:-999.), (nLocEle==0? mll: -999.),
                                   _met, minDeltaR, minDeltaRlead, mtHighest, mtLowest, leadingJetPt, trailJetPt, 0., double(_nVertex), mlll,
                                   _lEta[ptCorrV[0].second], _lEta[ptCorrV[1].second], (leptonSelection > 2 ? _lEta[ptCorrV[2].second] : -999.), (leptonSelection > 3 ? _lEta[ptCorrV[3].second] : -999.),
                                   (nLocEle == 3?mt1:-999.), (nLocEle==2?mt1:-999.), (nLocEle==1?mt1:-999.), (nLocEle==0? mt1: -999.),
                                   cosTSt, mll_ss, double(chargeOfLeptons), ll_deltaR, mt2ll_ss,
                                   (_lCharge[ind.at(0)] == -1 ?  mvaVLJECUp : -999), HTLocJERUp,
                                   SRIDTTZ(ind, indOf2LonZ, nJLocJERUp, nBLocJERUp, dMZ, mlll), SRIDWZCR(nJLocJERUp, nBLocJERUp, dMZ), SRIDZZCR(ind, indOf2LonZ, nJLocJERUp, nBLocJERUp), SRIDTTCR(nJLocJERUp, nBLocJERUp, dMZ, mlll),
                                   (leptonSelection == 3 && nJLocJERUp > 2 && nBLocJERUp > 0 ? SRIDPTZ(ptZ) : -999), (leptonSelection == 3 && nJLocJERUp > 2 && nBLocJERUp > 0 ? SRIDCosTheta(cosTSt) : -999),
                                   (leptonSelection == 3 ? flavourCategory3L(nLocEle) : -999),
                                   (leptonSelection == 4 ? flavourCategory4L(nLocEle) : -999),
                                   (leptonSelection == 4 ? flavourCategory4LZZ(nLocEle) : -999),
                                   (leptonSelection == 3 && nLocEle == 0 ? SRID3L(nJLocJERUp, nBLocJERUp, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 1 ? SRID3L(nJLocJERUp, nBLocJERUp, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 2 ? SRID3L(nJLocJERUp, nBLocJERUp, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 3 ? SRID3L(nJLocJERUp, nBLocJERUp, dMZ) : -999),
                                   (passTTZSRSelection(ind, indOf2LonZ, nJLocJERUp, nBLocJERUp, dMZ) ? flavourCategory3L4L(leptonSelection, nLocEle) : -999),
                                   (leptonSelection == 3 && nBLocJERUp > 0 ? SRID8SR3L(nJLocJERUp, nBLocJERUp, dMZ) : -999),
                                   };

          vector<double> fillVarJerDw = {ptCorrV[0].first, ptCorrV[1].first, leptonSelection > 2 ? ptCorrV[2].first : 0., leptonSelection > 3 ? ptCorrV[3].first : 0.,
                                   mt1, double(nJLocJERDown), double(nBLocJERDown), (_lCharge[ind.at(0)] == 1 ?  mvaVLJECDown : -999),
                                   (leptonSelection == 3 ? SRID3L(nJLocJERDown, nBLocJERDown, dMZ) : -999),
                                   (leptonSelection == 4 && passTTZ4LSelection(ind, indOf2LonZ, nJLocJERDown)? SRID4L(nJLocJERDown, nBLocJERDown) : -999),
                                   (leptonSelection != 4 ? mll:mll1stpair),ptZ,ptNonZ, (nLocEle == 3?mll:-999.), (nLocEle==2?mll:-999.), (nLocEle==1?mll:-999.), (nLocEle==0? mll: -999.),
                                   _met, minDeltaR, minDeltaRlead, mtHighest, mtLowest, leadingJetPt, trailJetPt, 0., double(_nVertex), mlll,
                                   _lEta[ptCorrV[0].second], _lEta[ptCorrV[1].second], (leptonSelection > 2 ? _lEta[ptCorrV[2].second] : -999.), (leptonSelection > 3 ? _lEta[ptCorrV[3].second] : -999.),
                                   (nLocEle == 3?mt1:-999.), (nLocEle==2?mt1:-999.), (nLocEle==1?mt1:-999.), (nLocEle==0? mt1: -999.),
                                   cosTSt, mll_ss, double(chargeOfLeptons), ll_deltaR, mt2ll_ss,
                                   (_lCharge[ind.at(0)] == -1 ?  mvaVLJECDown : -999), HTLocJECDown,
                                   SRIDTTZ(ind, indOf2LonZ, nJLocJERDown, nBLocJERDown, dMZ, mlll), SRIDWZCR(nJLocJERDown, nBLocJERDown, dMZ), SRIDZZCR(ind, indOf2LonZ, nJLocJERDown, nBLocJERDown), SRIDTTCR(nJLocJERDown, nBLocJERDown, dMZ, mlll),
                                   (leptonSelection == 3 && nJLocJERDown > 2 && nBLocJERDown > 0 ? SRIDPTZ(ptZ) : -999), (leptonSelection == 3 && nJLocJERDown > 2 && nBLocJERDown > 0 ? SRIDCosTheta(cosTSt) : -999),
                                   (leptonSelection == 3 ? flavourCategory3L(nLocEle) : -999),
                                   (leptonSelection == 4 ? flavourCategory4L(nLocEle) : -999),
                                   (leptonSelection == 4 ? flavourCategory4LZZ(nLocEle) : -999),
                                   (leptonSelection == 3 && nLocEle == 0 ? SRID3L(nJLocJERDown, nBLocJERDown, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 1 ? SRID3L(nJLocJERDown, nBLocJERDown, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 2 ? SRID3L(nJLocJERDown, nBLocJERDown, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 3 ? SRID3L(nJLocJERDown, nBLocJERDown, dMZ) : -999),
                                   (passTTZSRSelection(ind, indOf2LonZ, nJLocJERDown, nBLocJERDown, dMZ) ? flavourCategory3L4L(leptonSelection, nLocEle) : -999),
                                   (leptonSelection == 3 && nBLocJERDown > 0 ? SRID8SR3L(nJLocJERDown, nBLocJERDown, dMZ) : -999),
                                   };
          */

          vector<TString> fncName = {"ptlead", "sublead", "trail", "pt4th", 
                                     "mtW", "njets", "nbjets", "BDTpp", 
                                     "SR3L",
                                     "SR4L",
                                     "mll", "ptZ", "ptNonZ", "mll3e", "mll2e1mu", "mll1e2mu", "mll3mu",
                                     "met", "deltaR", "deltaRlead", "mtLeading", "mtTrailing", "leadJetPt", "trailJetPt", "SRnpCR", "nPV", "mlll",
                                     "etaLead", "etaSubl", "etaTrail", "eta4th", 
                                     "mt_3m", "mt_2m1e", "mt_1m2e", "mt_3e", 
                                     "cosThetaStar", "mll_ss", "chargeOfLeptons", "ll_deltaR", "mt2ll_ss", "BDTmm", "HT",
                                     "SRallTTZ", "SRWZCR", "SRZZCR", "SRTTCR",
                                     "SRttZCleanPTZ", "SRttZCleanCosTheta",
                                     "flavour3L", "flavour4L", "flavour4LZZ", 
                                     "SR3L3m","SR3L2m1e","SR3L1m2e","SR3L3e",
                                     "flavour3L4L",
                                     "SRTTZ8SR3L",
                                     "minMLeptonJet", "maxMLeptonJet", "minDeltaRLeptonJet", "maxDeltaRLeptonJet", "minDeltaPhiLeptonJet", "maxDeltaPhiLeptonJet", "minpTLeptonJet", "maxpTLeptonJet",
                                     "minMLeptonbJet", "maxMLeptonbJet", "minDeltaRLeptonbJet", "maxDeltaRLeptonbJet", "minDeltaPhiLeptonbJet", "maxDeltaPhiLeptonbJet", "minpTLeptonbJet", "maxpTLeptonbJet",
                                     "minMJetJet", "maxMJetJet", "minDeltaRJetJet", "maxDeltaRJetJet", "minDeltaPhiJetJet", "maxDeltaPhiJetJet", "minpTJetJet", "maxpTJetJet",
                                     "minDeltaPhiLeptonMET", "maxDeltaPhiLeptonMET", "minmTLeptonMET", "maxmTLeptonMET", "minpTLeptonMET", "maxpTLeptonMET",
                                     "minDeltaPhiJetMET", "maxDeltaPhiJetMET", "minmTJetMET", "maxmTJetMET", "minpTJetMET", "maxpTJetMET",
                                     "minDeltaPhiBJetMET", "maxDeltaPhiBJetMET", "minmTBJetMET", "maxmTBJetMET", "minpTBJetMET", "maxpTBJetMET"
                                   };
                                   
          //start = std::chrono::high_resolution_clock::now();
          /*
          double lepSF = 1.;
          double lepSFSystUp = 1.; double lepSFSystDown = 1.; 
          double lepSFStatUp = 1.; double lepSFStatDown = 1.;
          double lepSFRecoUp = 1.; double lepSFRecoDown = 1.;

          double puW = 1.; double puWUp = 1.; double puWDown = 1.;

          double btagL = 1.; double btagLUp = 1.; double btagLDown = 1.;
          double btagC = 1.; double btagCUp = 1.; double btagCDown = 1.;
          double btagB = 1.; double btagBUp = 1.; double btagBDown = 1.;

          double jetPrefW = 1.;

          if((samples[sam].getProcessName()) != "data"){
            double lepWOnlySyst = leptonWeightOnlySyst(0); double lepWOnlySystUp = leptonWeightOnlySyst(1); double lepWOnlySystDown = leptonWeightOnlySyst(2);
            double lepWOnlyStat = leptonWeightOnlyStat(0); double lepWOnlyStatUp = leptonWeightOnlyStat(1); double lepWOnlyStatDown = leptonWeightOnlyStat(2);
            double lepWOnlyReco = leptonWeightOnlyReco(0); double lepWOnlyRecoUp = leptonWeightOnlyReco(1); double lepWOnlyRecoDown = leptonWeightOnlyReco(2);

            lepSF = lepWOnlySyst * lepWOnlyStat * lepWOnlyReco; 
            lepSFSystUp = lepWOnlySystUp * lepWOnlyStat * lepWOnlyReco; 
            lepSFSystDown = lepWOnlySystDown * lepWOnlyStat * lepWOnlyReco;
            lepSFStatUp = lepWOnlySyst * lepWOnlyStatUp * lepWOnlyReco; 
            lepSFStatDown = lepWOnlySyst * lepWOnlyStatDown * lepWOnlyReco; 
            lepSFRecoUp = lepWOnlySyst * lepWOnlyStat * lepWOnlyRecoUp; 
            lepSFRecoDown = lepWOnlySyst * lepWOnlyStat * lepWOnlyRecoDown; 

            puW = puWeight(0); puWUp = puWeight(1); puWDown = puWeight(2);

            btagL = bTagWeight_udsg(0); btagLUp = bTagWeight_udsg(1); btagLDown = bTagWeight_udsg(2);
            btagC = bTagWeight_c(0); btagCUp = bTagWeight_c(1); btagCDown = bTagWeight_c(2);
            btagB = bTagWeight_b(0); btagBUp = bTagWeight_b(1); btagBDown = bTagWeight_b(2);

            //jetPrefW = jetPrefiringWeight();

            // for post fit scaling
            //weight *= scaler.postFitScaling(samples[sam].getProcessName());

          }
          if(debug){
              if(leptonSelection == 3 && nJLoc == 2 && nBLoc  == 1)
                cout << "weights are (lepSF/pu/btagL/btagBC) : " << lepSF << " " << puW << " " << btagL << " " << btagC*btagB << endl;
          }
          */

          //finish = std::chrono::high_resolution_clock::now();
          //elapsed = finish - start;
          //std::cout << "time needed to estimate all sf and deviations: " << elapsed.count() << std::endl;

          for(int dist = 0; dist < fillVar.size(); dist++){
            if(std::find(listToPrint[selection].begin(), listToPrint[selection].end(), fncName[dist]) == listToPrint[selection].end()) continue;
            //if(listToPrint[selection].find(fncName[dist]) == listToPrint[selection].end()) continue;
            distribs[dist].vectorHisto[samCategory].Fill(TMath::Min(fillVar.at(dist),figNames[fncName.at(dist)].varMax-0.1),weight);

            /*
            if((samples[sam].getProcessName()) != "data"){

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * lepSFSystUp / lepSF);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * lepSFSystDown / lepSF);
                
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 1, figNames[fncName.at(dist)].varMax-0.1, weight * lepSFStatUp / lepSF);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 1, figNames[fncName.at(dist)].varMax-0.1, weight * lepSFStatDown / lepSF);
                
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 2, figNames[fncName.at(dist)].varMax-0.1, weight * lepSFRecoUp / lepSF);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 2, figNames[fncName.at(dist)].varMax-0.1, weight * lepSFRecoDown / lepSF);
                
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 3, figNames[fncName.at(dist)].varMax-0.1, weight*puWUp/puW);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 3, figNames[fncName.at(dist)].varMax-0.1, weight*puWDown/puW);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 4, figNames[fncName.at(dist)].varMax-0.1, weight*btagLUp/btagL);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 4, figNames[fncName.at(dist)].varMax-0.1, weight*btagLDown/btagL);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 5, figNames[fncName.at(dist)].varMax-0.1, weight*btagCUp*btagBUp/btagC/btagB);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 5, figNames[fncName.at(dist)].varMax-0.1, weight*btagCDown*btagBDown/btagC/btagB);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVarJecUp.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVarJecDw.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVarJerUp.at(dist), 7, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVarJerDw.at(dist), 7, figNames[fncName.at(dist)].varMax-0.1, weight);

                if(samples[sam].getProcessName() == "WZ" && nBLoc > 0){ // 8 % uncertainty for WZ + bb background in high nbjets categories
                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 8, figNames[fncName.at(dist)].varMax-0.1, weight * 1.08);
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 8, figNames[fncName.at(dist)].varMax-0.1, weight * 0.92);
                }
                else{
                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 8, figNames[fncName.at(dist)].varMax-0.1, weight);
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 8, figNames[fncName.at(dist)].varMax-0.1, weight);
                }

                if((samples[sam].getFileName().find("ST_tWll_") != std::string::npos || samples[sam].getFileName().find("TTWJetsToLNu") != std::string::npos || samples[sam].getFileName().find("tZq_ll") != std::string::npos || samples[sam].getFileName().find("TTZToLLNuNu_M-10_") != std::string::npos) && samples[sam].is2017()){ 
                    // 10, 12 - factor 4; 6 and 8 - factor 2; 8, 9 - up factor 2, 6, 7 - down factor 2, 30th Aug Daniel said in ttX chat that for FSR factor sqrt(2) should be used, indeces 3 and 5
                    // as well from https://twiki.cern.ch/twiki/bin/view/CMS/TopSystematics#Factorization_and_renormalizatio there is an instruction to use envelope uncertrtainty, i.e. largest between ISR and FSR
                    // let's take largest deviation from unity for upward and downward fluction

                    //std::vector<double> psWeights = {_psWeight[8], _psWeight[5], _psWeight[6], _psWeight[3]};
                    //double psWUp = largestAmongAll(psWeights);
                    //double psWDown = smallestAmongAll(psWeights);

                    //if(debug) cout << "all weights, order is (ISR, FSR up) and (ISR, FSR down): " << _psWeight[8] << " " << _psWeight[5] << " " << _psWeight[6] << " " << _psWeight[3] << endl;
                    //if(debug) cout << "ps weight up and down: " << psWUp << " " << psWDown << endl;


                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight*_psWeight[8]*sumSimulatedEventWeights/sumSimulatedEventWeightsISRScaleUp); 
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight*_psWeight[6]*sumSimulatedEventWeights/sumSimulatedEventWeightsISRScaleDown); 

                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight*_psWeight[5]*sumSimulatedEventWeights/sumSimulatedEventWeightsFSRScaleUp); 
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight*_psWeight[3]*sumSimulatedEventWeights/sumSimulatedEventWeightsFSRScaleDown); 
                }
                else if(samples[sam].getFileName().find("TTZToLLNuNu_M-10_") != std::string::npos && samples[sam].is2016() && dist == indexSRTTZ){ // apply same weights as in 2017
                    // here let's take for the moment largest deviation from unity
                    double psWUp = ttZISRUpW[fillVar.at(dist)] > ttZFSRUpW[fillVar.at(dist)] ? ttZISRUpW[fillVar.at(dist)] : ttZFSRUpW[fillVar.at(dist)];
                    double psWDown = ttZISRDownW[fillVar.at(dist)] > ttZFSRDownW[fillVar.at(dist)] ? ttZISRDownW[fillVar.at(dist)] : ttZFSRDownW[fillVar.at(dist)];

                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight*psWUp); 
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight*psWDown); 

                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight*psWUp); 
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight*psWDown); 
                }
                else{
                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight);
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight);

                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight);
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight);
                }

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 11, figNames[fncName.at(dist)].varMax-0.1, weight*_lheWeight[8]*sumSimulatedEventWeights/sumSimulatedEventWeightsScaleUp);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 11,figNames[fncName.at(dist)].varMax-0.1, weight*_lheWeight[4]*sumSimulatedEventWeights/sumSimulatedEventWeightsScaleDown);

                if(dist == indexSR3L || dist == indexSR4L || dist == indexSRTTZ || dist == indexSRWZCR || dist == indexSRZZCR || dist == indexSRTTCR){
                    for(int varPDF = 0; varPDF < 100; varPDF++){
                        distribs[dist].vectorHistoPDF[samCategory].var[varPDF].Fill(fillVar.at(dist), weight*_lheWeight[9+varPDF]);
                    }
                }
                else{
                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 12, figNames[fncName.at(dist)].varMax-0.1, weight);
                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 12, figNames[fncName.at(dist)].varMax-0.1, weight);
                }

                for(int cat = 0; cat < 8; cat++){
                    if(cat == samCategory){
                        distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 13+cat, figNames[fncName.at(dist)].varMax-0.1, weight*(1+uncOnNorm[samples[sam].getProcessName()]));
                        distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 13+cat, figNames[fncName.at(dist)].varMax-0.1, weight*(1-uncOnNorm[samples[sam].getProcessName()]));
                    }
                    else{
                        distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 13+cat, figNames[fncName.at(dist)].varMax-0.1, weight);
                        distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 13+cat, figNames[fncName.at(dist)].varMax-0.1, weight);
                    }
                }

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 21, figNames[fncName.at(dist)].varMax-0.1, weight * 1.025); // 2.5% for lumi
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 21, figNames[fncName.at(dist)].varMax-0.1, weight * 0.975);
                
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 22, figNames[fncName.at(dist)].varMax-0.1, weight * 1.01); // 1% for trigger 
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 22, figNames[fncName.at(dist)].varMax-0.1, weight * 0.99);
                
            }
            else if(samCategory == nonPromptSample && leptonSelection != 4){
                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight);

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

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 11, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 11, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 12, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 12, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 13, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 13, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 14, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 14, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 15, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 15, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 16, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 16, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 17, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 17, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 18, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 18, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 19, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 19, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 20, figNames[fncName.at(dist)].varMax-0.1, weight*1.3);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 20, figNames[fncName.at(dist)].varMax-0.1, weight*0.7);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 21, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 21, figNames[fncName.at(dist)].varMax-0.1, weight);

                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 22, figNames[fncName.at(dist)].varMax-0.1, weight);
                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 22, figNames[fncName.at(dist)].varMax-0.1, weight);
            }
            */
          }
      }

      cout << endl;
      samCategory = processToCounterMap.at(samples[sam].getProcessName());
      cout << "Total number of events: " << distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[samCategory].Integral() << endl;
      //if(leptonSelection != 4)
      //  cout << "Total number of events in non prompt category: " << distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[nonPromptSample].Integral() << endl;
      cout << endl;
  }

  fileDummy->cd();
  signalTree->Write();
  bkgTree->Write();

  fileDummy->Close();


  // this should be done to be fully correct in PDF, takes a lot of time, an effect estimated with ttZ sample, uncertainty on acceptance is under 1%, simply assign flat uncertainty of 1 % to all signal and bkg
  // for the moment draw this uncertainty only for SR and will propagate them to datacards 
  /*
  for(int dist = 0; dist < figNames.size(); dist++){
    if(!(dist == indexSR3L || dist == indexSR4L || dist == indexSRTTZ || dist == indexSRWZCR || dist == indexSRZZCR || dist == indexSRTTCR)) continue;

    for(unsigned sam = 0; sam < processToCounterMap.size(); ++sam){
      if(sam == dataSample) continue;
      if(sam == nonPromptSample) continue;
      for(unsigned bin = 1; bin < (unsigned) distribs[dist].vectorHistoUncUp[sam].unc.at(pdfUncIndex).GetNbinsX() + 1; ++bin){
          double pdfVarRms = 0.;
          for(unsigned pdf = 0; pdf < 100; ++pdf){
              double variedBin = distribs[dist].vectorHistoPDF[sam].var[pdf].GetBinContent(bin);
              variedBin *= crossSectionRatio[sam][pdf];
              double diff = (  variedBin - distribs[dist].vectorHisto[sam].GetBinContent(bin) );
              pdfVarRms += diff * diff;
          }
          pdfVarRms = sqrt( 0.01 * pdfVarRms );
          //cout << "pdf rms for bin " << bin << " is equal to " << pdfVarRms << endl;
          distribs[dist].vectorHistoUncUp[sam].unc.at(pdfUncIndex).SetBinContent(bin, distribs[dist].vectorHisto[sam].GetBinContent(bin) + pdfVarRms);
          distribs[dist].vectorHistoUncDown[sam].unc.at(pdfUncIndex).SetBinContent(bin, distribs[dist].vectorHisto[sam].GetBinContent(bin) - pdfVarRms);
      }
    }
  }
  */

  // legend to print
  TLegend* mtleg = new TLegend(0.15,0.89,0.95,0.72); 
  mtleg->SetNColumns(4);
  mtleg->SetFillColor(0);
  mtleg->SetFillStyle(0);
  mtleg->SetBorderSize(0);
  mtleg->SetTextFont(42);

  //mtleg->AddEntry(&distribs[0].vectorHisto[dataSample],"Data","lep"); //data

  std::map<int, std::string> processToCounterMapReversed;
  std::vector<std::string> processOrder;
  for(map<std::string,int>::const_iterator it = processToCounterMap.begin();it != processToCounterMap.end(); ++it){
    processToCounterMapReversed.insert(std::pair<int,std::string>(it->second, it->first));
  }

  // correct order according to increase
  for(map<int,std::string>::const_iterator it = processToCounterMapReversed.begin();it != processToCounterMapReversed.end(); ++it){
    processOrder.push_back(it->second);
  }
  
  for(map<int, std::string>::const_iterator it = processToCounterMapReversed.begin();it != processToCounterMapReversed.end(); ++it){
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
    //showHist(plot[varPlot],distribs[figNames[listToPrint[crToPrint].at(varPlot)].index],figNames[listToPrint[crToPrint].at(varPlot)], scale_num, mtleg, false, false, showLegendOption);
    showSeparationHist(plot[varPlot],distribs[figNames[listToPrint[crToPrint].at(varPlot)].index],figNames[listToPrint[crToPrint].at(varPlot)], scale_num, mtleg, false, false, showLegendOption);
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + ".pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + ".png");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + ".root");
    /*
    plot[varPlot]->cd();
    //showHist(plot[varPlot],distribs[figNames[listToPrint[crToPrint].at(varPlot)].index],figNames[listToPrint[crToPrint].at(varPlot)], scale_num, mtleg, true, false, showLegendOption);
    showSeparationHist(plot[varPlot],distribs[figNames[listToPrint[crToPrint].at(varPlot)].index],figNames[listToPrint[crToPrint].at(varPlot)], scale_num, mtleg, true, false, showLegendOption);
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + "Log.pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + "Log.png");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + "Log.root");
    */
  }
  
  /*
  if(crToPrint == "ttZ3Lclean"){
    fillDatacards(distribs[indexSRttZcleanPTZ], samplesOrderNames, samplesOrder, "SRttZCleanPTZ", is2017);
    fillDatacards(distribs[indexSRttZcleanCosTheta], samplesOrderNames, samplesOrder, "SRttZCleanCosTheta", is2017);
  }
  */
  // no need to fill datacards when running over 2016 and 2017 together
  if(leptonSelection == 2) return;
  if(showLegendOption == 2) return;
  if(crToPrint == "ttZclean"){
    fillTablesForFlavour(distribs[indexFlavour3L4L], processOrder, "flavour3L4L", (bool)showLegendOption);
  }
  if(crToPrint == "ttZ"){

    fillTablesSRTTZ(distribs[indexSRTTZ], processOrder, "SRallTTZ", (bool)showLegendOption);
    fillTablesForFlavour(distribs[indexFlavour3L4L], processOrder, "flavour3L4L", (bool)showLegendOption);

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
