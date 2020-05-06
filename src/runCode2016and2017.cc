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

// a path to libKinFit.so has to be added to LD_LIBRARY_PATH
#include "/user/mniedzie/ttZ/TopKinFit/include/kinfit.h"

////////#include "TMVA/Tools.h"
////////#include "TMVA/Reader.h"
////////#include "TMVA/MethodCuts.h"

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

void treeReader::Analyze(const vector<std::string> & filesToAnalyse, const std::string option, const std::string selection, const string& sampleToDebug, long evNb ){

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

  PostFitScaler scaler2016, scaler2017, scaler2018;
  scaler2016.setInputFile("data/postFit/outputTTZ_2016_new.txt");
  scaler2016.setPostfitYields();
  scaler2017.setInputFile("data/postFit/outputTTZ_2017_new.txt");
  scaler2017.setPostfitYields();
  scaler2018.setInputFile("data/postFit/outputTTZ_2018.txt");
  scaler2018.setPostfitYields();


          KINFIT::kfit *kf = new KINFIT::kfit();
										kf->Init(TOPTOPLEPHAD); // Initialize tool for ttbar with one top decaying leptonically
										std::string pdfFileName = "/user/mniedzie/ttZ/TopKinFit/test/GenAnalysis/TopTopLepHad/pdf.root";

          kf->SetPDF("TopWMass",pdfFileName.c_str(),"TopLepWM_Fit");
          kf->SetPDF("TopMass",pdfFileName.c_str(),"TopLepRecM_Fit");
          kf->SetPDF("TopWHadMass",pdfFileName.c_str(),"TopHadWRecM_Fit");
          kf->SetPDF("TopHadMass",pdfFileName.c_str(),"TopHadRecM_Fit");
          kf->SetPDF("MetPx",pdfFileName.c_str(),"dMetPx_Gaus");
          kf->SetPDF("MetPy",pdfFileName.c_str(),"dMetPy_Gaus");
          kf->SetPDF("BJetPx",pdfFileName.c_str(),"dBJetPx_Fit");
          kf->SetPDF("BJetPy",pdfFileName.c_str(),"dBJetPy_Fit");
          kf->SetPDF("BJetPz",pdfFileName.c_str(),"dBJetPz_Fit");
          kf->SetPDF("BJetE",pdfFileName.c_str(),"dBJetE_Fit");
          kf->SetPDF("NonBJetPx",pdfFileName.c_str(),"dNonBJetPx_Fit");
          kf->SetPDF("NonBJetPy",pdfFileName.c_str(),"dNonBJetPy_Fit");
          kf->SetPDF("NonBJetPz",pdfFileName.c_str(),"dNonBJetPz_Fit");
          kf->SetPDF("NonBJetE",pdfFileName.c_str(),"dNonBJetE_Fit");
          kf->SetPDF("ElecPx",pdfFileName.c_str(),"dElecPx_Fit");
          kf->SetPDF("ElecPy",pdfFileName.c_str(),"dElecPy_Fit");
          kf->SetPDF("ElecPz",pdfFileName.c_str(),"dElecPz_Fit");
          kf->SetPDF("ElecE",pdfFileName.c_str(),"dElecE_Fit");
          kf->SetPDF("MuonPx",pdfFileName.c_str(),"dMuonPx_Fit");
          kf->SetPDF("MuonPy",pdfFileName.c_str(),"dMuonPy_Fit");
          kf->SetPDF("MuonPz",pdfFileName.c_str(),"dMuonPz_Fit");
          kf->SetPDF("MuonE",pdfFileName.c_str(),"dMuonE_Fit");
       
          kf->SetNToy(20);
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
      if(samples[sam].getProcessName() == "nonpromptData"){
          cout << "Total number of events: " << distribs[0].vectorHisto[samCategory].Integral() << endl;
          continue;
      }

      std::cout<<"Entries in "<< (samples[sam].getFileName()) << " " << nEntries << std::endl;
      double progress = 0.71;  //for printing progress bar


  //    int eeeLooseCount = 0;      int eeeeLooseCount = 0;
  //    int eeeFakeCount = 0;       int eeeeFakeCount = 0;
  //    int eeeTightCount = 0;      int eeeeTightCount = 0;
  //    int mmmLooseCount = 0;      int mmmmLooseCount = 0;
  //    int mmmFakeCount = 0;       int mmmmFakeCount = 0;
  //    int mmmTightCount = 0;      int mmmmTightCount = 0;
		// 
		// loop over all events in a given sample
		//
      for(long unsigned it = 0.71*nEntries; it < nEntries; ++it){
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

				std::cout << it << "d at 3:" << leptonMVAcutInAnalysis.at(3) << std::endl;
          GetEntry(it);
				std::cout << it << "d at 3:" << leptonMVAcutInAnalysis.at(3) << std::endl;
//										if( _eventNb == 32337502 ) debug = true;
          if(debug && (_eventNb != evNb && evNb != -999)) continue;
//										if( _eventNb == 2812206  ) { std::cout << "progress " << progress << ", event " << _eventNb << std::endl; continue; }
//										if( _eventNb == 2812197  ) { std::cout << "progress " << progress << ", event " << _eventNb << std::endl; continue; }
//										if( _eventNb == 2812200  ) { std::cout << "progress " << progress << ", event " << _eventNb << std::endl; continue; }
//										if( _eventNb == 2812191  ) { std::cout << "progress " << progress << ", event " << _eventNb << std::endl; continue; }
				std::cout << it << "d at 3:" << leptonMVAcutInAnalysis.at(3) << std::endl;
										
          if(debug) cout << "######################### New Event ###########################################" << endl;
          if(debug) cout << "event " << _eventNb << " was found" << endl;
          
          // trigger and met filters
          if(debug) cout << "trigger decision: " << _passTrigger_e << " " << _passTrigger_m << " " << _passTrigger_ee << " " << _passTrigger_em << " " << _passTrigger_mm << " " << _passTrigger_eee << " " << _passTrigger_eem << " " << _passTrigger_emm << " " << _passTrigger_mmm << endl;
          if(!(_passTrigger_e || _passTrigger_m || _passTrigger_ee || _passTrigger_em || _passTrigger_mm || _passTrigger_eee || _passTrigger_eem || _passTrigger_emm || _passTrigger_mmm)) continue;
          if(debug) cout << "met filers flag: " << _passMETFilters << endl;
          if(!_passMETFilters) continue;
          
          //if(it > 10000) break;
          //if(it > nEntries / 20) break;
          //if(it > 100) break;
//          std::cout << "ev number " << _eventNb << std::endl;

          std::vector<unsigned> indTight, indFake, indLoose, indOf2LonZ;
          //select leptons relative to the analysis
			 // start with 3 lepton selection, but it changes later
			 // the counts of tight/loose/tight4l/loose4l leptons are calculated here
          leptonSelection = 3;
          const unsigned lCount = selectLep(indTight, leptonSelection);
          const unsigned lCountFake = selectFakeLep(indFake, leptonSelection);
          const unsigned lCountLoose = selectLooseLep(indLoose);

          std::vector<unsigned> indTight4L, indLoose4L;
          const unsigned lCount4L = selectLep(indTight4L, 4);
          const unsigned lCount4LLoose = selectFakeLep(indLoose4L, 4);

         // if( selection != "ttZ4L" ){
         //   if(lCountLoose==3){
			      //     if(  _lMatchPdgId[indLoose[0]] != 22 && _lMatchPdgId[indLoose[1]] != 22 && _lMatchPdgId[indLoose[2]] != 22 ){
         //       if(_lFlavor[indLoose[0]] == 0 && _lFlavor[indLoose[1]] == 0 &&_lFlavor[indLoose[2]] == 0 ) eeeLooseCount++;
         //       if(_lFlavor[indLoose[0]] == 1 && _lFlavor[indLoose[1]] == 1 &&_lFlavor[indLoose[2]] == 1 ) mmmLooseCount++;
         //     }
         //   }
         //   if(lCountFake==3){
			      //     if(  _lMatchPdgId[indFake[0]] != 22 && _lMatchPdgId[indFake[1]] != 22 && _lMatchPdgId[indFake[2]] != 22 ){
         //       if(_lFlavor[indFake[0]] == 0 && _lFlavor[indFake[1]] == 0 &&_lFlavor[indFake[2]] == 0 ) eeeFakeCount++;
         //       if(_lFlavor[indFake[0]] == 1 && _lFlavor[indFake[1]] == 1 &&_lFlavor[indFake[2]] == 1 ) mmmFakeCount++;
         //     }
			      //  	}
         // }

         // if( selection == "ttZ4L"){
         //   if(lCount4LLoose==4){
			      //     if(  _lMatchPdgId[indLoose4L[0]] != 22 && _lMatchPdgId[indLoose4L[1]] != 22 && _lMatchPdgId[indLoose4L[2]] != 22 && _lMatchPdgId[indLoose4L[3]] != 22 ){
         //       if(_lFlavor[indLoose4L[0]] == 0 && _lFlavor[indLoose4L[1]] == 0 &&_lFlavor[indLoose4L[2]] == 0 &&_lFlavor[indLoose4L[3]] == 0) eeeeLooseCount++;
         //       if(_lFlavor[indLoose4L[0]] == 1 && _lFlavor[indLoose4L[1]] == 1 &&_lFlavor[indLoose4L[2]] == 1 &&_lFlavor[indLoose4L[3]] == 1) mmmmLooseCount++;
         //     }
         //   }
         //   
         //   if(lCountFake==4){
			      //     if(  _lMatchPdgId[indFake[0]] != 22 && _lMatchPdgId[indFake[1]] != 22 && _lMatchPdgId[indFake[2]] != 22 && _lMatchPdgId[indFake[3]] != 22 ){
         //       if(_lFlavor[indFake[0]] == 0 && _lFlavor[indFake[1]] == 0 &&_lFlavor[indFake[2]] == 0 &&_lFlavor[indFake[3]] == 0) eeeeFakeCount++;
         //       if(_lFlavor[indFake[0]] == 1 && _lFlavor[indFake[1]] == 1 &&_lFlavor[indFake[2]] == 1 &&_lFlavor[indFake[3]] == 1) mmmmFakeCount++;
         //     }
         //   }
         // }





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
          
 //  ////Marek////
 //         if( allLeptonsArePrompt== true){
 //           if( selection != "ttZ4L" && lCount==3 ){
	//		     if(  _lMatchPdgId[indTight[0]] != 22 && _lMatchPdgId[indTight[1]] != 22 && _lMatchPdgId[indTight[2]] != 22 ){
 //               if(_lFlavor[indTight[0]] == 0 && _lFlavor[indTight[1]] == 0 &&_lFlavor[indTight[2]] == 0 ) eeeTightCount++;
 //               if(_lFlavor[indTight[0]] == 1 && _lFlavor[indTight[1]] == 1 &&_lFlavor[indTight[2]] == 1 ) mmmTightCount++;
 //             }
 //           }
 //           if( selection == "ttZ4L" && lCount4L==4 ){
	//		     if( _lMatchPdgId[indTight4L[0]] != 22 && _lMatchPdgId[indTight4L[1]] != 22 && _lMatchPdgId[indTight4L[2]] != 22 && _lMatchPdgId[indTight4L[3]] != 22 ){
 //               if(_lFlavor[indTight4L[0]] == 0 && _lFlavor[indTight4L[1]] == 0 &&_lFlavor[indTight4L[2]] == 0 &&_lFlavor[indTight4L[3]] == 0) eeeeTightCount++;
 //               if(_lFlavor[indTight4L[0]] == 1 && _lFlavor[indTight4L[1]] == 1 &&_lFlavor[indTight4L[2]] == 1 &&_lFlavor[indTight4L[3]] == 1) mmmmTightCount++;
 //             }
 //           }
 //         }

	//		  //Marek//
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

          nJLoc = nJets(0, true, indJets);
          int nJLocDown = nJets(1, true, indJetsJECDown);
          int nJLocUp = nJets(2, true, indJetsJECUp);
          int nJLocJERDown = nJets(3, true, indJetsJERDown);
          int nJLocJERUp = nJets(4, true, indJetsJERUp);

          nBLoc = nBJets(0, true, true, indBJets, 1);
          int nBLocDown = nBJets(1, true, true, indBJetsJECDown, 1);
          int nBLocUp = nBJets(2, true, true, indBJetsJECUp, 1);
          int nBLocJERDown = nBJets(3, true, true, indBJetsJERDown, 1);
          int nBLocJERUp = nBJets(4, true, true, indBJetsJERUp, 1);

//          if( nBLoc != 0)  continue;

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

          std::vector<float> BJetPt,BJetEta,BJetPhi,BJetE;
          std::vector<float> NonBJetFilteredPt,NonBJetFilteredEta,NonBJetFilteredPhi,NonBJetFilteredE; 
          std::vector<float> ElectronPt,ElectronEta,ElectronPhi,ElectronE;
          std::vector<float> MuonPt,MuonEta,MuonPhi,MuonE;

          float TopLepBJetFitPt = 40.0;
          float TopLepBJetFitEta = 0.0;
          float TopLepBJetFitPhi = 0.0;
          float TopLepBJetFitE = 0.0;

          float LepTopMass = 0.0;
          float HadTopMass = 0.0;
          float HadWMass = 0.0;
          float LepTopPt = 0.0;
          float HadTopPt = 0.0;

//          cout << "number of leptons/jets/bjets/dMZ: " << lCount << " " << nJLoc << " " << nBLoc << " " << dMZ << endl;
          if( lCount>2 && dMZ < 10 && nBLoc>1 && (nJLoc-nBLoc)>1 ){
             if(debug) std::cout << "No of bjets = " << nBLoc << std::endl;
             for( int i =0; i<nBLoc; i++){
                BJetPt.push_back (_jetPt[indBJets[i]]);
                BJetEta.push_back(_jetEta[indBJets[i]]);
                BJetPhi.push_back(_jetPhi[indBJets[i]]);
                BJetE.push_back  (_jetE[indBJets[i]]);
             if(debug) std::cout <<  "indice, jet indice " << i << ", " << indBJets[i] <<
                           ", Pt  " << _jetPt[indBJets[i]] <<
                           ", Eta " << _jetEta[indBJets[i]] <<
                           ", Phi " << _jetPhi[indBJets[i]] <<
                           ", E   " << _jetE[indBJets[i]] << std::endl;
             }

             if(debug) std::cout << "No of light jets = " << nJLoc << std::endl;
             for( int i =0; i<nJLoc; i++){
                if( !((_jetDeepCsv_b[ indJets[i] ] + _jetDeepCsv_bb[ indJets[i] ]) > 0.6324) ){
                   NonBJetFilteredPt.push_back (_jetPt[indJets[i]]);
                   NonBJetFilteredEta.push_back(_jetEta[indJets[i]]);
                   NonBJetFilteredPhi.push_back(_jetPhi[indJets[i]]);
                   NonBJetFilteredE.push_back  (_jetE[indJets[i]]);
                   if(debug) std::cout <<  "indice, jet indice " << i << ", " << indJets[i] <<
                              ", Pt  " << _jetPt[indJets[i]] <<
                              ", Eta " << _jetEta[indJets[i]] <<
                              ", Phi " << _jetPhi[indJets[i]] <<
                              ", E   " << _jetE[indJets[i]] << std::endl;
                }
             }

             if(debug) std::cout << "No of light leptons = " << lCount << std::endl;
             for( int i =0; i<lCount; i++){
                      if(debug)  std::cout << i << "-th lepton flavor is " << _lFlavor[indTight[i]] << std::endl;
                if( (i != indOf2LonZ[0]) && (i != indOf2LonZ[1]) ){
                   if(_lFlavor[indTight[i]]==0){
                      if(debug) std::cout << "Lepton from top is an electron" << std::endl;
                      ElectronPt .push_back (_lPt[indTight[i]]);
                      ElectronEta.push_back(_lEta[indTight[i]]);
                      ElectronPhi.push_back(_lPhi[indTight[i]]);
                      ElectronE  .push_back  (_lE[indTight[i]]);
             if(debug) std::cout <<  ", Pt  " << _lPt[indTight[i]] <<
                                     ", Eta " << _lEta[indTight[i]] <<
                                     ", Phi " << _lPhi[indTight[i]] <<
                                     ", E   " << _lE[indTight[i]] << std::endl;
                   }else if(_lFlavor[indTight[i]]==1){
                      if(debug) std::cout << "Lepton from top is a muon" << std::endl;
                      MuonPt .push_back (_lPt[indTight[i]]);
                      MuonEta.push_back(_lEta[indTight[i]]);
                      MuonPhi.push_back(_lPhi[indTight[i]]);
                      MuonE  .push_back  (_lE[indTight[i]]);
                   }
                }
             }

             if(debug) std::cout <<  "# of bjets " << BJetPt.size() <<
                           " # of jets " << NonBJetFilteredPt.size() <<
                           " # of el " << ElectronPt.size() <<
                           " # of mu " << MuonPt.size() << 
                           " met values: " << _met*TMath::Cos(_metPhi) << ", " << _met*TMath::Sin(_metPhi) << std::endl;
             kf->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
             kf->SetNonBJet(NonBJetFilteredPt,NonBJetFilteredEta,NonBJetFilteredPhi,NonBJetFilteredE);

             // std::cout << "Setting el and mu kinematics" << std::endl;
             kf->SetElectron( ElectronPt, ElectronEta,ElectronPhi,ElectronE);
             kf->SetMuon( MuonPt, MuonEta,MuonPhi,MuonE);

             kf->SetMet(_met*TMath::Cos(_metPhi),_met*TMath::Sin(_metPhi));

             if(debug) cout << "number of jets/bjets/dMZ: " << nJLoc << " " << nBLoc << " " << dMZ << endl;
             if(debug) {
                for (int i=0; i< NonBJetFilteredPt.size(); i++){
                   if( NonBJetFilteredPt[i] != NonBJetFilteredPt[i]) 
                   std::cout << "one of jet pt is a NaN! " << std::endl;
                   std::cout << " pt is" << NonBJetFilteredPt[i]<< std::endl;
                }
             }


             // std::cout << "Running the fit " << std::endl;
             kf->Run(); // Run the tool

             int NPerm = kf->GetNPerm(); // Get number of permutations
             // std::cout << "NPerm = " << NPerm << std::endl;
             float disc = kf->GetDisc(0); // Get minimized likelihood value
             // Get index of object in input collection
             int idxEle  = kf->GetIndex( ELECTRON_TOPTOPLEPHAD, 0 );
             int idxMuon  = kf->GetIndex( MUON_TOPTOPLEPHAD, 0 );
             int idxBJetLep  = kf->GetIndex( BJETLEP_TOPTOPLEPHAD, 0); 
             int idxBJetHad  = kf->GetIndex( BJETHAD_TOPTOPLEPHAD, 0 );
             int idxJet1  = kf->GetIndex( NONBJET1_TOPTOPLEPHAD, 0 );
             int idxJet2  = kf->GetIndex( NONBJET2_TOPTOPLEPHAD, 0 );
             
             // Get reconstructed neutrino
             float NuPx = kf->GetNuPx(0,0);
             float NuPy = kf->GetNuPy(0,0);
             float NuPz = kf->GetNuPz(0,0);

             LepTopMass = kf -> GetTopMass(0,0);
             HadTopMass = kf -> GetTopMass(0,1);
             HadWMass = kf -> GetWMass(0,1);            
             LepTopPt = kf -> GetTopPt(0,0);
             HadTopPt = kf -> GetTopPt(0,1);
             // Build up b jet from leptonic top quark decay
             TopLepBJetFitPt = BJetPt[idxBJetLep];
             TopLepBJetFitEta = BJetEta[idxBJetLep];
             TopLepBJetFitPhi = BJetPhi[idxBJetLep];
             TopLepBJetFitE = BJetE[idxBJetLep];

//             std::cout << "the _jetId value for the leptonic leg b jet: " << _jetId[indBJets[idxBJetLep]] << ", was b tagged? " << bTaggedDeepCSV(indBJets[idxBJetLep], 1) << std::endl;
//             std::cout << "the leptonic top jet pt: " << TopLepBJetFitPt << std::endl;
//             std::cout << "the leptonic top mass: " << LepTopMass << std::endl;
//             std::cout << "the hadronic top mass: " << HadTopMass << std::endl;
          }


          if( TopLepBJetFitPt != TopLepBJetFitPt) 
             std::cout << "the leptonic top jet pt is a NaN! " << std::endl;


//////    TLorentzVector test;
//////    test.SetPtEtaPhiE(_lPt[0],0.0,_lPhi[0],_lE[0]);
//////    std::cout << "px, py " << test.Px() << ", " << test.Py() << "cos and sin phi times pt" << _lPt[0]*TMath::Cos(_lPhi[0]) << ", " << _lPt[0]*TMath::Sin(_lPhi[0]) << std::endl;

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
         // cout << "weight of event after est is " << weight << endl;

          int mvaValueRegion = 0;

          if(debug) cout << "lepton selection is " << leptonSelection << " total SR: " << SRIDTTZ(ind, indOf2LonZ, nJLoc, nBLoc, dMZ, mlll) << endl; 
          //if(leptonSelection == 4 && passZZCRSelection(ind, indOf2LonZ, nJLoc)) myfile << _runNb << " " << _lumiBlock << " " << _eventNb << endl;
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
          fillVar["topPt"] = TopLepBJetFitPt;

          fillVar["lepTopMass"] = LepTopMass;
          fillVar["hadTopMass"] = HadTopMass;
          fillVar["hadWMass"] = HadWMass;
          fillVar["lepTopPt"] = LepTopPt;
          fillVar["hadTopPt"] = HadTopPt;

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
//          fillVarJecUp["topPt"] = TopLepBJetFitPt*1.01;


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
//          fillVarJecDown["topPt"] = TopLepBJetFitPt*0.99;



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
//          fillVarJerUp["topPt"] = TopLepBJetFitPt*1.01;




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
//          fillVarJerDw["topPt"] = TopLepBJetFitPt*0.99;




          //start = std::chrono::high_resolution_clock::now();
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
            else if(samples[sam].is2017())
                weight *= scaler2017.postFitScaling(samples[sam].getProcessName() != "data" ? samples[sam].getProcessName() : "nonpromptData");
            else if(samples[sam].is2018())
                weight *= scaler2018.postFitScaling(samples[sam].getProcessName() != "data" ? samples[sam].getProcessName() : "nonpromptData");

          }
          if(debug){
              if(leptonSelection == 3 && nJLoc == 2 && nBLoc  == 1)
                cout << "weights are (lepSF/pu/btagL/btagBC) : " << lepSF << " " << puW << " " << btagL << " " << btagC*btagB << endl;
          }

          //finish = std::chrono::high_resolution_clock::now();
          //elapsed = finish - start;
          //std::cout << "time needed to estimate all sf and deviations: " << elapsed.count() << std::endl;

    // each histogram produced in the analysis has its version corresponding to each uncertainty (up and down)
    // below each of these uncertainty histograms is filled. The second argument of the FillUnc function corresponds to uncertainty source 

//          for(int dist = 0; dist < fillVar.size(); dist++){
//            if(std::find(listToPrint[selection].begin(), listToPrint[selection].end(), fncName[dist]) == listToPrint[selection].end()) continue;
//            //if(listToPrint[selection].find(fncName[dist]) == listToPrint[selection].end()) continue;
////            if(fncName[dist] == "topPt" ) std::cout << "filling dist : " << fncName[dist] << ", number of dist indice " << dist << ", samCategory " << samCategory <<" with value " << fillVar.at(dist) << ", top pt value : " << TopLepBJetFitPt << ", varMax thingy:  " << figNames[fncName.at(dist)].varMax-0.1  << " and the minimum " << TMath::Min(fillVar.at(dist),figNames[fncName.at(dist)].varMax-0.1) << std::endl;
//            distribs[dist].vectorHisto[samCategory].Fill(TMath::Min(fillVar.at(dist),figNames[fncName.at(dist)].varMax-0.1),weight);


          int idx1 = figNames[ "promptElMVA" ].index;
          int idx2 = figNames[ "looseElMVA"  ].index;
          int idx3 = figNames[ "promptMuMVA" ].index;
          int idx4 = figNames[ "looseMuMVA"  ].index;
          int idx5 = figNames[ "selectedByMVA"].index;

          for( int i =0; i<lCount; i++){
             if( _lFlavor[ indTight[i] ] == 0 ){
                for( int i = 0; i < distribs[ idx1 ].vectorHisto[ samCategory ].GetBin( _leptonMvatZq[i] ) ; i++ ) {
                   distribs[ idx1 ].vectorHisto[ samCategory ].Fill( distribs[idx1].vectorHisto[samCategory].GetBinCenter(i), weight );
										      }
             }else if( _lFlavor[ indTight[i] ] == 1 ){
                for( int i = 0; i < distribs[ idx3 ].vectorHisto[ samCategory ].GetBin( _leptonMvatZq[i] ) ; i++ ) {
                   distribs[ idx3 ].vectorHisto[ samCategory ].Fill( distribs[idx1].vectorHisto[samCategory].GetBinCenter(i), weight );
										      }
             }
          }

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

//                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * lepSFSystUp / lepSF);
//                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * lepSFSystDown / lepSF);
//                
//                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 1, figNames[fncName.at(dist)].varMax-0.1, weight * lepSFStatUp / lepSF);
//                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 1, figNames[fncName.at(dist)].varMax-0.1, weight * lepSFStatDown / lepSF);
//                
//                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 2, figNames[fncName.at(dist)].varMax-0.1, weight * lepSFRecoUp / lepSF);
//                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 2, figNames[fncName.at(dist)].varMax-0.1, weight * lepSFRecoDown / lepSF);

                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 0, varMaxOfVar-0.1, weight * lepSFSystUp / lepSF); 
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var, 0, varMaxOfVar-0.1, weight * lepSFSystDown / lepSF); 

                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 1, varMaxOfVar-0.1, weight * lepSFStatUp / lepSF); 
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var, 1, varMaxOfVar-0.1, weight * lepSFStatDown / lepSF); 

                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 2, varMaxOfVar-0.1, weight * lepSFRecoUp / lepSF); 
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  2, varMaxOfVar-0.1, weight * lepSFRecoDown / lepSF); 


                //distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * jetPrefW);
                //distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 0, figNames[fncName.at(dist)].varMax-0.1, weight * jetPrefW);

//                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 3, figNames[fncName.at(dist)].varMax-0.1, weight*puWUp/puW);
//                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 3, figNames[fncName.at(dist)].varMax-0.1, weight*puWDown/puW);
//
//                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 4, figNames[fncName.at(dist)].varMax-0.1, weight*btagLUp/btagL);
//                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 4, figNames[fncName.at(dist)].varMax-0.1, weight*btagLDown/btagL);
//
//                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 5, figNames[fncName.at(dist)].varMax-0.1, weight*btagCUp*btagBUp/btagC/btagB);
//                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 5, figNames[fncName.at(dist)].varMax-0.1, weight*btagCDown*btagBDown/btagC/btagB);


                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 3, varMaxOfVar-0.1, weight*puWUp/puW); 
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  3, varMaxOfVar-0.1, weight*puWDown/puW); 

                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 4, varMaxOfVar-0.1, weight*btagLUp/btagL); 
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  4, varMaxOfVar-0.1, weight*btagLDown/btagL); 

                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 5, varMaxOfVar-0.1, weight*btagCUp*btagBUp/btagC/btagB); 
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  5, varMaxOfVar-0.1, weight*btagCDown*btagBDown/btagC/btagB); 

					 // JEC and JER are applied directly to jet momenta, so rather than reweighting distributions, one files distributions with modified jet momenta
					 // thus here fillVarJecUp etc. values are used

//                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVarJecUp.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight);
//                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVarJecDw.at(dist), 6, figNames[fncName.at(dist)].varMax-0.1, weight);
//
//                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVarJerUp.at(dist), 7, figNames[fncName.at(dist)].varMax-0.1, weight);
//                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVarJerDw.at(dist), 7, figNames[fncName.at(dist)].varMax-0.1, weight);

                 distribs[index].vectorHistoUncUp[samCategory].FillUnc(varJECUp, 6, varMaxOfVar-0.1, weight); 
                 distribs[index].vectorHistoUncDown[samCategory].FillUnc(varJECDown,  6, varMaxOfVar-0.1, weight); 

                 distribs[index].vectorHistoUncUp[samCategory].FillUnc(varJERUp, 7, varMaxOfVar-0.1, weight); 
                 distribs[index].vectorHistoUncDown[samCategory].FillUnc(varJERDown,  7, varMaxOfVar-0.1, weight); 


                if(samples[sam].getProcessName() == "WZ" && nBLoc > 0){ // 8 % uncertainty for WZ + bb background in high nbjets categories
//                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 8, figNames[fncName.at(dist)].varMax-0.1, weight * 1.08);
//                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 8, figNames[fncName.at(dist)].varMax-0.1, weight * 0.92);

                   distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 8, varMaxOfVar-0.1, weight * 1.08); 
                   distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  8, varMaxOfVar-0.1, weight * 0.92); 
                }
                else{
//                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 8, figNames[fncName.at(dist)].varMax-0.1, weight);
//                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 8, figNames[fncName.at(dist)].varMax-0.1, weight);

                   distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 8, varMaxOfVar-0.1, weight); 
                   distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  8, varMaxOfVar-0.1, weight); 
                }

                if((samples[sam].getFileName().find("ST_tWll_") != std::string::npos || samples[sam].getFileName().find("TTWJetsToLNu") != std::string::npos || samples[sam].getFileName().find("tZq_ll") != std::string::npos || samples[sam].getFileName().find("TTZToLLNuNu_M-10_") != std::string::npos) && ( samples[sam].is2017() || samples[sam].is2018() )){ 
                    // 10, 12 - factor 4; 6 and 8 - factor 2; 8, 9 - JERup factor 2, 6, 7 - down factor 2, 30th Aug Daniel said in ttX chat that for FSR factor sqrt(2) should be used, indeces 3 and 5
                    // as well from https://twiki.cern.ch/twiki/bin/view/CMS/TopSystematics#Factorization_and_renormalizatio there is an instruction to use envelope uncertrtainty, i.e. largest between ISR and FSR
                    // description of the order for the PS weights: https://twiki.cern.ch/twiki/bin/view/CMS/TopModGen#Event_Generation
//                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight*_psWeight[8]*sumSimulatedEventWeights/sumSimulatedEventWeightsISRScaleUp); 
//                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight*_psWeight[6]*sumSimulatedEventWeights/sumSimulatedEventWeightsISRScaleDown); 
//
//                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight*_psWeight[5]*sumSimulatedEventWeights/sumSimulatedEventWeightsFSRScaleUp); 
//                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight*_psWeight[3]*sumSimulatedEventWeights/sumSimulatedEventWeightsFSRScaleDown); 

                     distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 9, varMaxOfVar-0.1, weight*_psWeight[8]*sumSimulatedEventWeights/sumSimulatedEventWeightsISRScaleUp); 
                     distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  9, varMaxOfVar-0.1, weight*_psWeight[6]*sumSimulatedEventWeights/sumSimulatedEventWeightsISRScaleDown); 

                     distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 10, varMaxOfVar-0.1, weight*_psWeight[5]*sumSimulatedEventWeights/sumSimulatedEventWeightsFSRScaleUp); 
                     distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  10, varMaxOfVar-0.1, weight*_psWeight[3]*sumSimulatedEventWeights/sumSimulatedEventWeightsFSRScaleDown); 

                }
//                else if(samples[sam].getFileName().find("TTZToLLNuNu_M-10_") != std::string::npos && samples[sam].is2016() && dist == indexSRTTZ){ // apply same weights as in 2017
//                    // here let's take for the moment largest deviation from unity
//        // the arrays used here contain numbers taken from 2017 dists
//                    double psWUp = ttZISRUpW[fillVar.at(dist)] > ttZFSRUpW[fillVar.at(dist)] ? ttZISRUpW[fillVar.at(dist)] : ttZFSRUpW[fillVar.at(dist)];
//                    double psWDown = ttZISRDownW[fillVar.at(dist)] > ttZFSRDownW[fillVar.at(dist)] ? ttZISRDownW[fillVar.at(dist)] : ttZFSRDownW[fillVar.at(dist)];
//
//                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight*psWUp); 
//                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight*psWDown); 
//
//                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight*psWUp); 
//                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight*psWDown); 
//                }
                else{
//                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight);
//                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 9, figNames[fncName.at(dist)].varMax-0.1, weight);
//
//                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight);
//                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 10, figNames[fncName.at(dist)].varMax-0.1, weight);


                     distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 9, varMaxOfVar-0.1, weight); 
                     distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  9, varMaxOfVar-0.1, weight); 

                     distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 10, varMaxOfVar-0.1, weight); 
                     distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  10, varMaxOfVar-0.1, weight); 
                }

      // [4] both renormalization and factorization scales down by 2
      // [8] both renormalization and factorization scales up by 2

//                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 11, figNames[fncName.at(dist)].varMax-0.1, weight*_lheWeight[8]*sumSimulatedEventWeights/sumSimulatedEventWeightsScaleUp);
//                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 11,figNames[fncName.at(dist)].varMax-0.1, weight*_lheWeight[4]*sumSimulatedEventWeights/sumSimulatedEventWeightsScaleDown);

                 // renormalization and factorization uncertainties 
                distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 11, varMaxOfVar-0.1, weight*_lheWeight[8]*sumSimulatedEventWeights/sumSimulatedEventWeightsScaleUp); 
                distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  11, varMaxOfVar-0.1, weight*_lheWeight[4]*sumSimulatedEventWeights/sumSimulatedEventWeightsScaleDown); 


                if(index == indexSR3L || index == indexSR4L || index == indexSRTTZ || index == indexSRWZCR || index == indexSRZZCR || index == indexSRTTCR){ 
//                if(dist == indexSR3L || dist == indexSR4L || dist == indexSRTTZ || dist == indexSRWZCR || dist == indexSRZZCR || dist == indexSRTTCR){
                    for(int varPDF = 0; varPDF < 100; varPDF++){
                         distribs[index].vectorHistoPDF[samCategory].var[varPDF].Fill(var, weight*_lheWeight[9+varPDF]); 
//                        distribs[dist].vectorHistoPDF[samCategory].var[varPDF].Fill(fillVar.at(dist), weight*_lheWeight[9+varPDF]);
                    }
                }
                else{
//                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 12, figNames[fncName.at(dist)].varMax-0.1, weight);
//                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 12, figNames[fncName.at(dist)].varMax-0.1, weight);

                     distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 12, varMaxOfVar-0.1, weight); 
                     distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  12, varMaxOfVar-0.1, weight); 
                }

// normalization uncertianties are stored in arrays defined in interface/readTreeSync.h
                for(int cat = 0; cat < 8; cat++){
                    if(cat == samCategory){
//                        distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 13+cat, figNames[fncName.at(dist)].varMax-0.1, weight*(1+uncOnNorm[samples[sam].getProcessName()]));
//                        distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 13+cat, figNames[fncName.at(dist)].varMax-0.1, weight*(1-uncOnNorm[samples[sam].getProcessName()]));

                         distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 13+cat, varMaxOfVar-0.1, weight*(1+uncOnNorm[samples[sam].getProcessName()])); 
                         distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  13+cat, varMaxOfVar-0.1, weight*(1-uncOnNorm[samples[sam].getProcessName()])); 
                    }
                    else{
//                        distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 13+cat, figNames[fncName.at(dist)].varMax-0.1, weight);
//                        distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 13+cat, figNames[fncName.at(dist)].varMax-0.1, weight);

                         distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 13+cat, varMaxOfVar-0.1, weight); 
                         distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  13+cat, varMaxOfVar-0.1, weight); 
                    }
                }

//                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 21, figNames[fncName.at(dist)].varMax-0.1, weight * 1.025); // 2.5% for lumi
//                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 21, figNames[fncName.at(dist)].varMax-0.1, weight * 0.975);
//                
//                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 22, figNames[fncName.at(dist)].varMax-0.1, weight * 1.01); // 1% for trigger 
//                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 22, figNames[fncName.at(dist)].varMax-0.1, weight * 0.99);

                 distribs[index].vectorHistoUncUp[samCategory].FillUnc(var,  21, varMaxOfVar-0.1, weight * 1.025); // 2.5% for lumi 
                 distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  21, varMaxOfVar-0.1, weight * 0.975); 

                 distribs[index].vectorHistoUncUp[samCategory].FillUnc(var,  22, varMaxOfVar-0.1, weight * 1.01); // 1% for trigger 
                 distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  22, varMaxOfVar-0.1, weight * 0.99); 

                if(samples[sam].getProcessName() == "WZ" && nJLoc > 2){ // 20 % uncertainty in tails of njets, namely in njets >= 3
//                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 23, figNames[fncName.at(dist)].varMax-0.1, weight * 1.2);
//                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 23, figNames[fncName.at(dist)].varMax-0.1, weight * 0.8);

                     distribs[index].vectorHistoUncUp[samCategory].FillUnc(var,  23, varMaxOfVar-0.1, weight * 1.2); 
                     distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  23, varMaxOfVar-0.1, weight * 0.8); 
                }
                else{
//                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 23, figNames[fncName.at(dist)].varMax-0.1, weight);
//                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 23, figNames[fncName.at(dist)].varMax-0.1, weight);

                     distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, 23, varMaxOfVar-0.1, weight); 
                     distribs[index].vectorHistoUncDown[samCategory].FillUnc(var,  23, varMaxOfVar-0.1, weight); 
                }


            }
            else if(samCategory == nonPromptSample && leptonSelection != 4){
                for(int cat = 0; cat < 20; cat++){
//                    distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), cat, figNames[fncName.at(dist)].varMax-0.1, weight);
//                    distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), cat, figNames[fncName.at(dist)].varMax-0.1, weight);
                     if(cat == 20){ 
                       distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, cat, varMaxOfVar-0.1, weight*1.3); 
                       distribs[index].vectorHistoUncDown[samCategory].FillUnc(var, cat, varMaxOfVar-0.1, weight*0.7); 
                     } 
                     else{ 
                       distribs[index].vectorHistoUncUp[samCategory].FillUnc(var, cat, varMaxOfVar-0.1, weight); 
                       distribs[index].vectorHistoUncDown[samCategory].FillUnc(var, cat, varMaxOfVar-0.1, weight); 
                     } 
                }

//                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 20, figNames[fncName.at(dist)].varMax-0.1, weight*1.3);
//                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 20, figNames[fncName.at(dist)].varMax-0.1, weight*0.7);
//
//                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 21, figNames[fncName.at(dist)].varMax-0.1, weight);
//                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 21, figNames[fncName.at(dist)].varMax-0.1, weight);
//
//                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 22, figNames[fncName.at(dist)].varMax-0.1, weight);
//                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 22, figNames[fncName.at(dist)].varMax-0.1, weight);
//
//                distribs[dist].vectorHistoUncUp[samCategory].FillUnc(fillVar.at(dist), 23, figNames[fncName.at(dist)].varMax-0.1, weight);
//                distribs[dist].vectorHistoUncDown[samCategory].FillUnc(fillVar.at(dist), 23, figNames[fncName.at(dist)].varMax-0.1, weight);
            }
          }
			 // break out when the event was found
          if(debug) break;
      }

		// loop through all of the histograms which zeros bins which bin contents are smaller than 0. only for prompt prediction.
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

////    std::cout << "electrons before miniIso cut     : " <<  leptons_preMiniIsoCut  << std::endl;
////    std::cout << "electrons after  miniIso cut     : " <<  leptons_postMiniIsoCut << std::endl;
////    std::cout << "electrons after  missing hits cut: " <<  leptons_postMissingHits<< std::endl;
////    std::cout << "electrons after  passEmu cut     : " <<  leptons_postPassEmu    << std::endl;
////    std::cout << "electrons after  loose cut       : " <<  leptons_loose          << std::endl;
////    std::cout << "electrons after  clean cut       : " <<  leptons_clean          << std::endl;
////    std::cout << "electrons after  MVA cut         : " <<  leptons_passMVA        << std::endl;
////
////    leptons_preMiniIsoCut  = 0;	
////    leptons_postMiniIsoCut = 0;
////    leptons_loose          = 0;
////    leptons_clean          = 0;
////    leptons_passMVA        = 0;
////    leptons_postMissingHits= 0;
////    leptons_postPassEmu    = 0;

  // ////Marek////
  //      cout << "Number of loose eee events: " << eeeLooseCount << endl; 
  //      cout << "Number of fake eee events: " << eeeFakeCount << endl;
  //      cout << "Number of tight eee events: " << eeeTightCount << endl;
  //      cout << "Number of loose mmm events: " << mmmLooseCount << endl;
  //      cout << "Number of fake mmm events: " << mmmFakeCount << endl;
  //      cout << "Number of tight mmm events: " << mmmTightCount << endl;
  //
  //      cout << "Number of loose eeee events: " << eeeeLooseCount << endl; 
  //      cout << "Number of fake eeee events: " << eeeeFakeCount << endl;
  //      cout << "Number of tight eeee events: " << eeeeTightCount << endl;
  //      cout << "Number of loose mmmm events: " << mmmmLooseCount << endl;
  //      cout << "Number of fake mmmm events: " << mmmmFakeCount << endl;
  //      cout << "Number of tight mmmm events: " << mmmmTightCount << endl;
  //			  //Marek//
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

  // legend to print
  TLegend* mtleg = new TLegend(0.18,0.89,0.92,0.72); 
  mtleg->SetNColumns(3);
  if(selection == "ZZ" || selection == "ttZ4L")
    mtleg->SetNColumns(4);
  if(selection == "ttZ")
    mtleg->SetNColumns(5);
  mtleg->SetFillColor(0);
  mtleg->SetFillStyle(0);
  mtleg->SetBorderSize(0);
  mtleg->SetTextFont(42);
  mtleg->SetTextSize(0.06);

  // fill the legend with entries.

  mtleg->AddEntry(&distribs[figNames[listToPrint[selection].at(0)].index].vectorHisto[dataSample],"Data","ep"); //data

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

  TH1D* histStatAndSystBand = (TH1D*)distribs[figNames[listToPrint[selection].at(0)].index].stack.GetStack()->Last()->Clone("histStatAndSystBand");
  histStatAndSystBand ->SetFillStyle(3005);
  histStatAndSystBand ->SetLineColor(kGray+2);
  histStatAndSystBand ->SetFillColor(kGray+2);
  histStatAndSystBand ->SetMarkerStyle(1);
  mtleg->AddEntry(histStatAndSystBand, "Uncertainty" ,"f"); // "Total unc." in Willem's plots

  // plots to make with systematics and stat uncertainty on them
  std::string processToStore = selection;
  TString folderToStorePlots;
  int showLegendOption = 0; // 0 - 2016, 1 - 2017, 2 - 2016+2017, 3 - 2018
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
    else if (filesToAnalyse.at(0).find("2018") != std::string::npos){
      folderToStorePlots = "2018/";
       showLegendOption = 3;
    }
    else{
      folderToStorePlots = "2016/";
       showLegendOption = 0;
    }
  }

  gSystem->Exec("rm plotsForSave/" + folderToStorePlots + processToStore + "/*.{pdf,png,root}");
  gSystem->Exec("rmdir plotsForSave/" + folderToStorePlots + processToStore);
  gSystem->Exec("mkdir -p plotsForSave/" + folderToStorePlots + processToStore);
// 	std::string MVAstr = std::to_string (MVAcut);
//  MVAstr.erase ( MVAstr.find_last_not_of('0') + 1, std::string::npos );
//
//  gSystem->Exec("rm output/outputMVAvalue"+ MVAstr +"/plotsForSave2/" + folderToStorePlots + processToStore + "/*.{pdf,png,root}");
//
//  gSystem->Exec("rmdir output/outputMVAvalue"+ MVAstr +"/plotsForSave2/" + folderToStorePlots + processToStore);
//  gSystem->Exec("mkdir -p output/outputMVAvalue"+ MVAstr +"/plotsForSave2/" + folderToStorePlots + processToStore);
  double scale_num = 1.6;
  
  TCanvas* plot[nVars];
  for(int i = 0; i < nVars; i++){
      plot[i] = new TCanvas(Form("plot_%d", i),"",500,500);
  }

  std::string crToPrint = selection;

// TFile f("histos.root", "recreate");
// distribs[0].vectorHisto[0].Write();
// distribs[0].vectorHisto[1].Write();
//        TFile f1("histos_before.root", "recreate");
//        distribs[figNames["topPt"].index].vectorHisto[0].Write();
//        distribs[figNames["topPt"].index].vectorHisto[1].Write();
//        distribs[figNames["topPt"].index].vectorHisto[2].Write();
//        distribs[figNames["topPt"].index].vectorHisto[3].Write();
//        distribs[figNames["topPt"].index].vectorHisto[4].Write();
//        distribs[figNames["topPt"].index].vectorHisto[5].Write();

  for(int varPlot = 0; varPlot < listToPrint[crToPrint].size(); varPlot++){
    plot[varPlot]->cd();
    showHist(plot[varPlot], distribs[figNames[listToPrint[crToPrint].at(varPlot)].index], figNames[listToPrint[crToPrint].at(varPlot)], scale_num, mtleg, false, false, showLegendOption);
//    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + ".pdf");
    std::cout << Form("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + ".pdf") << std::endl;
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + ".pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + ".png");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + ".root");
    plot[varPlot]->cd();
    showHist(plot[varPlot], distribs[figNames[listToPrint[crToPrint].at(varPlot)].index], figNames[listToPrint[crToPrint].at(varPlot)], scale_num, mtleg, true, false, showLegendOption);
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + "Log.pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + "Log.png");
    plot[varPlot]->SaveAs("plotsForSave/" + folderToStorePlots + processToStore + "/" + listToPrint[crToPrint].at(varPlot) + "Log.root");
  }

//        TFile f2("histos_after.root", "recreate");
//        distribs[figNames["topPt"].index].vectorHisto[0].Write();
//        distribs[figNames["topPt"].index].vectorHisto[1].Write();
//        distribs[figNames["topPt"].index].vectorHisto[2].Write();
//        distribs[figNames["topPt"].index].vectorHisto[3].Write();
//        distribs[figNames["topPt"].index].vectorHisto[4].Write();
//        distribs[figNames["topPt"].index].vectorHisto[5].Write();
  
  if(crToPrint == "ttZ"){
    fillTablesSRTTZ(distribs[indexSRTTZ], processOrder, "SRallTTZ", showLegendOption);
    //fillTablesForFlavour(distribs[indexFlavour3L4L], processOrder, "flavour3L4L", showLegendOption);
  }
  if(crToPrint == "ttZclean"){
    fillTablesForFlavour(distribs[indexFlavour3L4L], processOrder, "flavour3L4L", showLegendOption);
  }

  // no need to fill datacards when running over 2016 and 2017 together
		bool Analyzing2018 = false;
		if (filesToAnalyse.at(0).find("2018") != std::string::npos){
		  Analyzing2018 = true;
		}
  if(showLegendOption == 2) return;
  if(crToPrint == "ttZ3Lclean"){
    fillDatacards(distribs[indexSRttZcleanPTZ], processOrder, "SRttZCleanPTZ", (bool)showLegendOption, Analyzing2018);
    fillDatacards(distribs[indexSRttZcleanCosTheta], processOrder, "SRttZCleanCosTheta", (bool)showLegendOption, Analyzing2018);
  }
  if(crToPrint == "ttZ"){

    fillDatacards(distribs[indexSRTTZ8SR3L], processOrder, "SRTTZ8SR3L", (bool)showLegendOption, Analyzing2018); 
    fillDatacards(distribs[indexSR3L], processOrder, "SR3L", (bool)showLegendOption, Analyzing2018); 
    fillDatacards(distribs[indexSR4L], processOrder, "SR4L", (bool)showLegendOption, Analyzing2018); 
    fillDatacards(distribs[indexSRTTZ], processOrder, "SRallTTZ", (bool)showLegendOption, Analyzing2018); 
    fillDatacards(distribs[indexSR3L3m], processOrder, "SR3L3m", (bool)showLegendOption, Analyzing2018); 
    fillDatacards(distribs[indexSR3L2m1e], processOrder, "SR3L2m1e", (bool)showLegendOption, Analyzing2018); 
    fillDatacards(distribs[indexSR3L1m2e], processOrder, "SR3L1m2e", (bool)showLegendOption, Analyzing2018); 
    fillDatacards(distribs[indexSR3L3e], processOrder, "SR3L3e", (bool)showLegendOption, Analyzing2018); 
    fillDatacards(distribs[indexSRWZCR], processOrder, "SRWZCR", (bool)showLegendOption, Analyzing2018); 
    fillDatacards(distribs[indexSRZZCR], processOrder, "SRZZCR", (bool)showLegendOption, Analyzing2018); 
    fillDatacards(distribs[indexSRTTCR], processOrder, "SRTTCR", (bool)showLegendOption, Analyzing2018); 

    fillDatacards(distribs[indexLeadPt], processOrder, "ptlead", (bool)showLegendOption, Analyzing2018); 
    fillDatacards(distribs[indexTrailPt], processOrder, "trail", (bool)showLegendOption, Analyzing2018); 
  }
  return;
}

int main(int argc, const char **argv)
{
    //int rargc = 1; char// *rargv[1] = {""};
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
