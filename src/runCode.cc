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
  magicFactor = magicFactorAnalysis;
  leptonMVAcut = leptonMVAcutAnalysis;
  //Set CMS plotting style
  setTDRStyle();
  //gROOT->SetBatch(kTRUE);
  //read samples and cross sections from txt file
  //readSamples("data/samples_forTraining.txt");
  //readSamples("data/samples_forTraining_2017.txt");
  //readSamples("data/samples_MCvsData_NonpromptBoth.txt");
  //readSamples("data/samples_ForSync_2017.txt");
  //readSamples("data/samples_ZZ_sync.txt");
  //readSamples("data/samples_DYCR.txt");
  //readSamples("data/samples_DYCR_2017.txt");
  //readSamples("data/samples_ttZ.txt");
  readSamples("data/samples_ttZ_npDD.txt");
  //readSamples("data/samples_ttZ_2017.txt");
  //readSamples("data/samples_ttZ_2017_npDD.txt");
  //readSamples("data/samples_ttW.txt");
  //readSamples("data/samples_ttW_npDD.txt");
  //readSamples("data/samples_ttW_2017.txt");
  //readSamples("data/samples_ttW_2017_npDD.txt");
  setTDRStyle(); 

  std::vector<std::string> namesOfSamples = treeReader::getNamesOfTheSample();
  initdistribs(namesOfSamples);

  /*
  readerBtag[0][0].load(calib_csvv2[0], BTagEntry::FLAV_B, "comb");
  readerBtag[0][1].load(calib_csvv2[0], BTagEntry::FLAV_C, "comb");
  readerBtag[0][2].load(calib_csvv2[0], BTagEntry::FLAV_UDSG, "comb");

  readerBtag[1][0].load(calib_csvv2[1], BTagEntry::FLAV_B, "iterativefit");
  readerBtag[1][1].load(calib_csvv2[1], BTagEntry::FLAV_C, "iterativefit");
  readerBtag[1][2].load(calib_csvv2[1], BTagEntry::FLAV_UDSG, "iterativefit");
  */
   readerBtag[0][0].load(calib_csvv2[0], BTagEntry::FLAV_B, "comb");
   readerBtag[0][1].load(calib_csvv2[0], BTagEntry::FLAV_C, "comb");
   readerBtag[0][2].load(calib_csvv2[0], BTagEntry::FLAV_UDSG, "incl");

   readerBtag[1][0].load(calib_csvv2[1], BTagEntry::FLAV_B, "comb");
   readerBtag[1][1].load(calib_csvv2[1], BTagEntry::FLAV_C, "comb");
   readerBtag[1][2].load(calib_csvv2[1], BTagEntry::FLAV_UDSG, "incl");
  
  LastError::lasterror = Errors::UNKNOWN;

  vector<TH2D> fakeMaps;
  getFRmaps(fakeMaps, is2017, leptonSelection);

  //vector<vector<TH2D>> fakeMaps;
  //getFRmapsForTTZ(fakeMaps);

  if(LastError::lasterror != Errors::OK){
     cout << "FR maps not found" << endl;
     return;
  }
  
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

      Color_t color = assignColor(std::get<0>(samples[sam]));
      setStackColors(color, sam);

      //if(std::get<0>(samples[sam]) != "ZZ" && std::get<0>(samples[sam]) != "data") continue;
      //if(std::get<0>(samples[sam]) != "ZZ") continue;
      if(std::get<0>(samples[sam]) == "nonpromptData"){
          cout << "Total number of events: " << distribs[0].vectorHisto[sam].Integral() << endl;
          continue;
      }

      if(leptonSelectionAnalysis == 3 || leptonSelectionAnalysis == 4)
        if(std::get<0>(samples[sam]) == "chargeMisID") continue;

      //if(std::get<0>(samples[sam]) != "nonprompt" && std::get<0>(samples[sam]) != "ttW" && std::get<0>(samples[sam]) != "ttH" && std::get<0>(samples[sam]) != "ttZ" && std::get<0>(samples[sam]) != "chargeMisID" && std::get<0>(samples[sam]) != "WZ") continue;
      //if(std::get<0>(samples[sam]) != "Nonprompt" && std::get<0>(samples[sam]) != "ttW" && std::get<0>(samples[sam]) != "chargeMisID") continue;
      //if(!(std::get<1>(samples[sam]).find("ZZTo4L_13TeV_powheg_pythia8") != std::string::npos )) continue;
      //if(!(std::get<1>(samples[sam]).find("TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8") != std::string::npos )) continue;

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
          //if(_eventNb != 3596719) continue;
          bool printAddInfo = false;
          
          if(printAddInfo)
              cout << _passTrigger_e << " " << _passTrigger_m << " " << _passTrigger_ee << " " << _passTrigger_em  << " " << _passTrigger_mm << " " << _passTrigger_eee << " " << _passTrigger_eem << " " << _passTrigger_emm << " " << _passTrigger_mmm << endl;
          if(!(_passTrigger_e || _passTrigger_m || _passTrigger_ee || _passTrigger_em || _passTrigger_mm || _passTrigger_eee || _passTrigger_eem || _passTrigger_emm || _passTrigger_mmm)) continue;
          if(printAddInfo) cout << "pass met filters: " << _passMETFilters << endl;
          if(!_passMETFilters) continue;
          
          //if(it > 10000) break;
          //if(it > nEntries / 100) break;

          std::vector<unsigned> indTight, indFake, indOf2LonZ;
          //select leptons
          const unsigned lCount = selectLep(indTight);
          const unsigned lCountFake = selectFakeLep(indFake);

          if(printAddInfo) cout << "before 12 GeV cut" << endl;
          if(invMassOfAny2Lbelow12GeV(indFake)) continue;
          if(printAddInfo) cout << "passed 12 GeV cut" << endl;

          if(printAddInfo) cout << "number of tight and fo leptons: " << lCount << " " << lCountFake << endl;

          int samCategory = sam;

          std::vector<unsigned> ind;
          //used for main analysis
          // remove additional FO lepton from ss2l selection
          if(leptonSelectionAnalysis == 2 && lCountFake != leptonSelectionAnalysis) continue;
          if(lCount > leptonSelectionAnalysis) continue;
          if(lCountFake < leptonSelectionAnalysis) continue;
          if(lCount == leptonSelectionAnalysis) {
              samCategory = sam;
              ind = indTight;
          }
          if(lCount < leptonSelectionAnalysis){
              if(leptonSelectionAnalysis == 4) continue;
              samCategory = nonPromptSample;
              ind = indFake;
          }

          if(printAddInfo) cout << "number of tight leptons is correct" << endl;

          if(leptonSelectionAnalysis == 2)
            if(_lCharge[ind.at(0)] * _lCharge[ind.at(1)] < 0) continue;

          bool allLeptonsArePrompt = true;
          
          if(std::get<0>(samples[sam]) != "data" && std::get<0>(samples[sam]) != "nonpromptData")
            allLeptonsArePrompt = promptLeptons(ind);
          
          if(std::get<0>(samples[sam]) == "chargeMisID" && !allLeptonsArePrompt) continue;
          if(std::get<0>(samples[sam]) == "Nonprompt" && allLeptonsArePrompt) continue; // works just for MC

          if((std::get<0>(samples[sam]) == "ttW" || std::get<0>(samples[sam]) == "ttH" || std::get<0>(samples[sam]) == "ttZ" || std::get<0>(samples[sam]) == "ttX" || std::get<0>(samples[sam]) == "WZ" || std::get<0>(samples[sam]) == "Z#gamma"  || std::get<0>(samples[sam]) == "ZZ" || std::get<0>(samples[sam]) == "rare") && !allLeptonsArePrompt) continue;
          //if(!noConversionInSelection(ind)) continue;

          int nLocEle = getElectronNumber(ind);
          //if(nLocEle > 0) continue;
          //printAddInfo = true;

          /*
          unsigned electronIndex = -999;
          for(auto & l : ind){
            if(_lFlavor[l] == 0) electronIndex = l;
          }
          */

          if(leptonSelectionAnalysis == 4)
            if(!passPtCuts4L(ind)) continue;
          if(leptonSelectionAnalysis == 3)
            if(!passPtCuts3L(ind)) continue;
          if(leptonSelectionAnalysis == 2)
            if(!passPtCuts2L(ind)) continue;

          //if(!twoLeptonsInEndcap(ind)) continue;

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

          HTLoc = HTCalc(indJets);
          
          if(printAddInfo){
            cout << "number of jets/bjets/deltaMZ: " << nJLoc << " " << nBLoc << " " << dMZ << endl;
          }
          
          double mt1 = 9999;
          if(leptonSelectionAnalysis == 4){
            // ZZ control region
            if(dMZ > 15) continue; // 20 GeV originally in TOP-17-005
            mll1stpair = mll;
            std::vector<unsigned> vectorRemoved;
            vectorRemoved = ind;
            if(printAddInfo)
                cout << "index and pt of leptons that make first Z boson: ind1->" << indOf2LonZ.at(0) << " pt1->" << _lPt[indOf2LonZ.at(0)] << " ind2->" << indOf2LonZ.at(1) << " pt2->" << _lPt[indOf2LonZ.at(1)] << endl;
            vectorRemoved.erase(std::remove(vectorRemoved.begin(), vectorRemoved.end(), indOf2LonZ.at(0)), vectorRemoved.end());
            vectorRemoved.erase(std::remove(vectorRemoved.begin(), vectorRemoved.end(), indOf2LonZ.at(1)), vectorRemoved.end());
            dMZ = deltaMZ(vectorRemoved, third, mll, ptZ, ptNonZ, mlll, indOf2LonZ, Zboson, lnegative);
            if(printAddInfo) cout << "second deltaMZ is : " << dMZ << endl;
            if(dMZ > 15) continue;
            if(printAddInfo)
                cout << "index and pt of leptons that make second Z boson: ind1->" << indOf2LonZ.at(0) << " pt1->" << _lPt[indOf2LonZ.at(0)] << " ind2->" << indOf2LonZ.at(1) << " pt2->" << _lPt[indOf2LonZ.at(1)] << endl;
            mll2ndpair = mll;
          }

          if(leptonSelectionAnalysis == 3){
            
            // ttZ selection
            /*
            if(nJLoc < 2) continue;
            if(dMZ > 10) continue;
            */

            /*
            if(nJLoc < 3) continue;
            if(nBLoc < 1) continue;
            if(dMZ > 10) continue;
            */
            
            //if(false) continue;

            // WZ CR
            if(dMZ > 10) continue;
            if(nBLoc != 0) continue;
            // in origin one this cut is not used
            //if(nJLoc < 2) continue;
            //if(_met < 30) continue;
            cosTSt = cosThetaStar(Zboson, lnegative);
            //cout << "cos theta star is: " << cosTSt << endl;
            if(nJLoc < 1) continue;

            // ttbar CR
            // original one from TOP-17-005
            //if(!(dMZ == 999999 || !((dMZ < 10) || (nBLoc < 1)))) continue;
            
            // conversion CR
            //if(!(mlll > 81 && mlll < 101 && dMZ > 10)) continue; //  && nBLoc == 0

            // DY CR
            /*
            if(dMZ > 10) continue;

            TLorentzVector l0p4;
            l0p4.SetPtEtaPhiE(ptNonZ, _lEta[third], _lPhi[third], _lE[third] * ptNonZ / _lPt[third]);

            double mt1;
            mt1 = mtCalc(l0p4, _met, _metPhi);

            if(_met > 30) continue;
            if(mt1 > 30) continue;
            if(nJLoc > 1) continue;
            if(nBLoc > 0) continue;
            */

          }
          
          if(leptonSelectionAnalysis == 2){
            if(nJLoc < 2) continue;
            if(nBLoc < 1) continue;

            //if(nJLoc > 1) continue;

            TLorentzVector l0p4, l1p4;
            l0p4.SetPtEtaPhiE(ptCorrV[0].first, _lEta[ind.at(0)], _lPhi[ind.at(0)], _lE[ind.at(0)] * ptCorrV[0].first / _lPt[ind.at(0)]);
            l1p4.SetPtEtaPhiE(ptCorrV[1].first, _lEta[ind.at(1)], _lPhi[ind.at(1)], _lE[ind.at(1)] * ptCorrV[1].first / _lPt[ind.at(1)]);

            double ele_mll = (l0p4+l1p4).M();

            if(ele_mll < 12) continue;   
            if(ele_mll > 76 && ele_mll < 106 && nLocEle == 2) continue;
            if(_met < 30) continue;
          }
          
          if(printAddInfo){
            cout << "initial weight: " << weight << endl;
          }

          if(!isData)
            applyTriggerSF(weight, leptonSelection, ptCorrV[0].first);
          
          // ---------------------- here lepton SF and PU reweighing is being applied
          double dataMCSF = 1.;
          double lepSF = 1.;
          double lepSFUp = 1.;
          double lepSFDown = 1.;

          double tempValue = (double) rand() / (RAND_MAX);
          int leptonFileDicision = -99;

          if(tempValue < 20./35.9) leptonFileDicision = 0;
          else leptonFileDicision = 1;  

          //myfile << _runNb << " " << _lumiBlock << " " << _eventNb << " : ";
          myfile << _runNb << " " << _lumiBlock << " " << _eventNb << endl;
          if(std::get<0>(samples[sam]) != "nonpromptData" && std::get<0>(samples[sam]) != "data"){

            if(printAddInfo) cout << "number of true interactions are: " << _nTrueInt << endl;
            dataMCSF = (is2017 ? h_dataMC_2017->GetBinContent(h_dataMC_2017->GetXaxis()->FindBin(TMath::Min(double(_nTrueInt), 99.))) : h_dataMC_2016->GetBinContent(h_dataMC_2016->GetXaxis()->FindBin(TMath::Min(double(_nTrueInt), 49.))));

            for(unsigned int leptonInd = 0; leptonInd < leptonSelectionAnalysis; leptonInd++){

                double sfForLep = getLeptonSF(_lFlavor[ind.at(leptonInd)], _lPt[ind.at(leptonInd)], (_lFlavor[ind.at(leptonInd)] ? _lEta[ind.at(leptonInd)] : _lEtaSC[ind.at(leptonInd)]), 0, 0, leptonSelectionAnalysis, is2017);
                lepSF *= sfForLep;
                if(printAddInfo) cout << "lepton SF for lepton with pt " << _lPt[ind.at(leptonInd)] << " is " << sfForLep << endl;
                //myfile << std::setprecision(3) << sfForLep << "\t" ;

                lepSFUp *= getLeptonSF(_lFlavor[ind.at(leptonInd)], _lPt[ind.at(leptonInd)], (_lFlavor[ind.at(leptonInd)] ? _lEta[ind.at(leptonInd)] : _lEtaSC[ind.at(leptonInd)]), 1, 0, leptonSelectionAnalysis, is2017);
                lepSFDown *= getLeptonSF(_lFlavor[ind.at(leptonInd)], _lPt[ind.at(leptonInd)], (_lFlavor[ind.at(leptonInd)] ? _lEta[ind.at(leptonInd)] : _lEtaSC[ind.at(leptonInd)]), -1, 0, leptonSelectionAnalysis, is2017);
                if(printAddInfo) cout << "SF after " << leptonInd + 1 << " lepton is " << lepSF << endl; 
            }
          }
          //myfile << std::endl;
          
          weight = weight * dataMCSF;
          if(printAddInfo){
            cout << "weight after pileup reweigning is: " << weight << endl;
            cout << "lepton SF is " << lepSF << endl;
          }

          weight = weight * lepSF;

          if(printAddInfo){
            cout << "weight after lepton SF is: " << weight << endl;
          }
          
          // END ---------------------- here lepton SF and PU reweighing is being applied
          
          // ---------------------- here btag SF is being applied
          double btagSF_event = 1.;
          double btagSF_event_Up = 1.;
          double btagSF_event_Down = 1.;

          double btagSF_light_event = 1.;
          double btagSF_light_event_Up = 1.;
          double btagSF_light_event_Down = 1.;

          // see method 1a here https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
          double p_data = 1.;
          double p_mc = 1.;

          for(unsigned int i = 0; i < nJLoc; ++i){

            int j = indJets.at(i);

            double csv = _jetDeepCsv_bb[j] + _jetDeepCsv_b[j];
            if( csv < 0.0 ) csv = -0.05;
            if( csv > 1.0 ) csv = 1.0;

            double pt = _jetPt[j];
            if( pt > 1000 ) pt = 999.;

            if(std::get<0>(samples[sam]) != "nonpromptData" && std::get<0>(samples[sam]) != "data"){

              double centralValue = getBTagSF(is2017, 0, _jetHadronFlavor[j], TMath::Abs(_jetEta[j]), pt, csv);
              btagSF_event *= centralValue;
              double btagEffcalc = h_btagEff[is2017][jetFlavourNumber[_jetHadronFlavor[j]] + 3 * (leptonSelection - 2)]->GetBinContent(h_btagEff[is2017][jetFlavourNumber[_jetHadronFlavor[j]] + 3 * (leptonSelection - 2)]->FindBin(TMath::Min(_jetPt[j], 599.),fabs(_jetEta[j])));
              if(bTaggedDeepCSV(j, 1)){
               p_mc *= btagEffcalc;
               p_data *= centralValue * btagEffcalc;
              }
              else{
               p_mc *= 1 - btagEffcalc;
               p_data *= 1 - centralValue * btagEffcalc;
              }
                    
              /*
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
              */

            }

          }
          
          btagSF_event = p_data / p_mc;
          weight = weight * btagSF_event;

          if(printAddInfo){
            cout << "weight after b-tagging is: " << weight << endl;
          }
          
          // END ---------------------- here btag SF is being applied
          
          double FRloc = 1.;

          //myfile << _runNb << " " << _lumiBlock << " " << _eventNb << " " << lCount << " " << lCountFake << " "; //  << FRloc << endl; 
          if(printAddInfo){
              for(auto & i : ind)
                cout << "info about all leptons: " << _lPt[i] << " " << _lEta[i] << " " << _lFlavor[i] << " " << _lCharge[i] << " " << _leptonMvatZqTTV[i]  << endl;
          }

          //if(std::get<0>(samples[sam]) == "nonpromptData"){ 
          if(samCategory == nonPromptSample){ 

            int nFakeLepCounter = 0;

            for(int j = 0; j < ind.size(); j++){
              int i = ind.at(j);
              if(lepIsGood(i)) continue;    

              double leptFakePtCorr = magicFactor * _lPt[i] / _ptRatio[i];
              double FRloc_loc = fakeMaps.at(_lFlavor[i]).GetBinContent(fakeMaps.at(_lFlavor[i]).FindBin(TMath::Min(leptFakePtCorr,65-1.), fabs(_lEta[i])));
              //myfile << leptFakePtCorr << " " << _lEta[i] << " " << _lFlavor[i] << " ";
              if(printAddInfo) cout << "info about FO not tight lepton: " << leptFakePtCorr << " " << _lEta[i] << " " << _lFlavor[i] << " " << _leptonMvatZqTTV[i]  << endl;

              // should be applied if FR was measured as tight / FO-not-tight, here our FO definition is just FO, it can pass tight selection
              FRloc *= FRloc_loc/(1.-FRloc_loc);
              // get FR from MC maps for the moment
              //FRloc *= FRloc_loc;
              nFakeLepCounter++;

            }

            FRloc *= TMath::Power(-1, nFakeLepCounter + 1);
            // contribution of nonprompt leptons should be subtracted from nonprompt sideband category 
            if(sam != dataSample && sam != nonPromptSample)
              FRloc *= -1;
          }
          //myfile << FRloc << endl;
          
          if(printAddInfo) cout << "the FR is: " << FRloc << endl;
          
          //FRloc = 1.;
          weight = weight * FRloc;  

          if(printAddInfo) cout << "total weight is: " << weight << endl;

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
            //mvaVL = readerTTWcsttbar->EvaluateMVA( "BDTG method");
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

          if(leptonSelectionAnalysis == 3){

            // WZ CR
            TLorentzVector l0p4;
            l0p4.SetPtEtaPhiE(ptNonZ, _lEta[third], _lPhi[third], _lE[third] * ptNonZ / _lPt[third]);

            mt1 = mtCalc(l0p4, _met, _metPhi);

            distribs[4].vectorHisto[samCategory].Fill(TMath::Min(mt1,figNames["mtW"].varMax-0.1),weight);
            distribs[31 + nLocEle].vectorHisto[samCategory].Fill(TMath::Min(mt1,figNames["mtW"].varMax-0.1),weight);
            //if(mt1 < 50) continue;
          }

          distribs[7].vectorHisto[samCategory].Fill(TMath::Min(mvaVL,figNames["BDT"].varMax-0.001),weight);

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

          //myfile << _runNb << " " << _lumiBlock << " " << _eventNb << endl; //  << FRloc << endl; 

          distribs[0].vectorHisto[samCategory].Fill(TMath::Min(ptCorrV[0].first,figNames["ptlead"].varMax-0.1),weight);
          distribs[1].vectorHisto[samCategory].Fill(TMath::Min(ptCorrV[1].first,figNames["ptlead"].varMax-0.1),weight);
          if(leptonSelectionAnalysis > 2){
            distribs[2].vectorHisto[samCategory].Fill(TMath::Min(ptCorrV[2].first,figNames["trail"].varMax-0.1),weight);
            distribs[29].vectorHisto[samCategory].Fill(TMath::Min(_lEta[ptCorrV[2].second],figNames["etaSubl"].varMax-0.001),weight);
          }
          if(leptonSelectionAnalysis > 3){
            distribs[3].vectorHisto[samCategory].Fill(TMath::Min(ptCorrV[3].first,figNames["pt4th"].varMax-0.1),weight);
            distribs[30].vectorHisto[samCategory].Fill(TMath::Min(_lEta[ptCorrV[3].second],figNames["etaSubl"].varMax-0.001),weight);
          }

          distribs[27].vectorHisto[samCategory].Fill(TMath::Min(_lEta[ptCorrV[0].second],figNames["etaLead"].varMax-0.001),weight);
          distribs[28].vectorHisto[samCategory].Fill(TMath::Min(_lEta[ptCorrV[1].second],figNames["etaSubl"].varMax-0.001),weight);

          distribs[5].vectorHisto[samCategory].Fill(TMath::Min(double(nJLoc),figNames["njets"].varMax-0.1),weight);
          distribs[6].vectorHisto[samCategory].Fill(TMath::Min(double(nBLoc),figNames["nbjets"].varMax-0.1),weight);

          if(leptonSelectionAnalysis == 4){
            distribs[8].vectorHisto[samCategory].Fill(flavourCategory4L(nLocEle),weight);
          }
          if(leptonSelectionAnalysis == 3){
            distribs[8].vectorHisto[samCategory].Fill(flavourCategory3L(nLocEle),weight);
            distribs[9].vectorHisto[samCategory].Fill(SRID3L(nJLoc, nBLoc),weight);

            distribs[9].vectorHistoUp[samCategory].Fill(SRID3L(nJLoc, nBLoc),weight / lepSF * lepSFUp);
            distribs[9].vectorHistoDown[samCategory].Fill(SRID3L(nJLoc, nBLoc),weight / lepSF * lepSFDown);
          }

          if(leptonSelectionAnalysis == 2){
            distribs[8].vectorHisto[samCategory].Fill(flavourCategory2L(nLocEle, _lCharge[ind.at(0)]),weight);
            distribs[9].vectorHisto[samCategory].Fill(SRID2L(nJLoc, nBLoc, mvaValueRegion, _lCharge[ind.at(0)]),weight);
            distribs[29].vectorHisto[samCategory].Fill(SRID2L(nJLoc, nBLoc, mvaValueRegion, _lCharge[ind.at(0)]),weight);

            distribs[9].vectorHistoUp[samCategory].Fill(SRID2L(nJLoc, nBLoc, mvaValueRegion, _lCharge[ind.at(0)]),weight / lepSF * lepSFUp);
            distribs[9].vectorHistoDown[samCategory].Fill(SRID2L(nJLoc, nBLoc, mvaValueRegion, _lCharge[ind.at(0)]),weight / lepSF * lepSFDown);
          }

          if(leptonSelection != 4) 
            distribs[10].vectorHisto[samCategory].Fill(TMath::Min(mll,figNames["mll"].varMax-0.1),weight);
          else{
            distribs[10].vectorHisto[samCategory].Fill(TMath::Min(mll1stpair,figNames["mll"].varMax-0.1),weight);
            distribs[10].vectorHisto[samCategory].Fill(TMath::Min(mll2ndpair,figNames["mll"].varMax-0.1),weight);
          }
          distribs[11].vectorHisto[samCategory].Fill(TMath::Min(ptZ,figNames["ptZ"].varMax-0.1),weight);
          distribs[12].vectorHisto[samCategory].Fill(TMath::Min(ptNonZ,figNames["ptNonZ"].varMax-0.1),weight);

          distribs[16-nLocEle].vectorHisto[samCategory].Fill(TMath::Min(mll,figNames["mll"].varMax-0.1),weight);

          distribs[17].vectorHisto[samCategory].Fill(TMath::Min(_met,figNames["met"].varMax-0.1),weight);
          distribs[18].vectorHisto[samCategory].Fill(TMath::Min(minDeltaR,figNames["deltaR"].varMax-0.01),weight);
          distribs[19].vectorHisto[samCategory].Fill(TMath::Min(minDeltaRlead,figNames["deltaRlead"].varMax-0.01),weight);
          distribs[20].vectorHisto[samCategory].Fill(TMath::Min(mtHighest,figNames["mtLeading"].varMax-0.01),weight);
          distribs[21].vectorHisto[samCategory].Fill(TMath::Min(mtLowest,figNames["mtTrailing"].varMax-0.01),weight);
          distribs[22].vectorHisto[samCategory].Fill(TMath::Min(leadingJetPt,figNames["leadJetPt"].varMax-0.01),weight);
          distribs[23].vectorHisto[samCategory].Fill(TMath::Min(trailJetPt,figNames["trailJetPt"].varMax-0.01),weight);

          distribs[25].vectorHisto[samCategory].Fill(TMath::Min(double(_nVertex),figNames["nPV"].varMax-0.01),weight);
          distribs[26].vectorHisto[samCategory].Fill(TMath::Min(mlll,figNames["mlll"].varMax-0.01),weight);

          distribs[35].vectorHisto[samCategory].Fill(TMath::Min(cosTSt,figNames["cosThetaStar"].varMax-0.01),weight);
          distribs[36].vectorHisto[samCategory].Fill(TMath::Min(mll_ss,figNames["mll_ss"].varMax-0.01),weight);
          distribs[37].vectorHisto[samCategory].Fill(TMath::Min(double(chargeOfLeptons),figNames["chargeOfLeptons"].varMax-0.01),weight);
          distribs[38].vectorHisto[samCategory].Fill(TMath::Min(ll_deltaR,figNames["ll_deltaR"].varMax-0.01),weight);
          distribs[39].vectorHisto[samCategory].Fill(TMath::Min(mt2ll_ss,figNames["mt2ll_ss"].varMax-0.01),weight);

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
  
  for (int i=0; i!=nVars; ++i)  {
    
    //cout << "The sample size is " << samples.size() << endl;
    for (int sam=0; sam != samples.size(); ++sam){
      if(std::get<0>(samples[sam]) == "data") continue;
      
      // is used in MC CT 
      //if(sam == 0) continue;
      
      //cout << "the sample is added: " << std::get<0>(samples[sam]) << endl;
      distribs[i].vectorHistoTotalUnc.Add(&distribs[i].vectorHisto[sam]);
    }
        
    /*
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
    */

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
  vector<TString> listToPrint = {"ptlead", "sublead", "trail", "njets", "nbjets", "flavour", "mll", "ptZ", "ptNonZ", "mtW", "mll3e", "mll2e1mu", "mll1e2mu", "mll3mu", "met", "nPV", "mt_3m", "mt_2m1e", "mt_1m2e", "mt_3e", "cosThetaStar"};
  // list for Zgamma
  //vector<TString> listToPrint = {"ptlead", "sublead", "trail", "njets", "nbjets", "flavour", "met", "nPV", "mlll"};
  // list for ttbar CR
  //vector<TString> listToPrint = {"ptlead", "sublead", "trail", "njets", "nbjets", "flavour", "met", "nPV"};
  // list for DY nonprompt CR
  //vector<TString> listToPrint = {"ptlead", "sublead", "trail", "njets", "nbjets", "flavour", "met", "nPV", "mll", "ptZ", "ptNonZ", "mtW", "mll3e", "mll2e1mu", "mll1e2mu", "mll3mu", "mt_3m", "mt_2m1e", "mt_1m2e", "mt_3e"};
  // list for ss2l ttW
  //vector<TString> listToPrint = {"ptlead", "sublead", "njets", "nbjets", "flavour", "met", "nPV", "BDT", "deltaR", "deltaRlead", "mtLeading", "mtTrailing", "leadJetPt", "trailJetPt", "etaLead", "etaSubl", "mll_ss", "chargeOfLeptons", "ll_deltaR", "mt2ll_ss", "SR"};
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
  
  /*
  TCanvas * plotLepUnc = new TCanvas("plotLepUnc", "plotLepUnc", 500, 450);
  drawSystUnc(plotLepUnc, distribs[15], 1);
  plotLepUnc->SaveAs("plotsForSave/lepSFUnc.pdf");
  plotLepUnc->SaveAs("plotsForSave/lepSFUnc.png");
  plotLepUnc->SaveAs("plotsForSave/lepSFUnc.root");
  */

  //fillDatacards(distribs[indexSR], samplesOrderNames, samplesOrder);
}

int main(int argc, char *argv[]){

    TApplication *rootapp = new TApplication("example", &argc, argv);
    treeReader reader;
    reader.Analyze();
    rootapp->Run();

    return 0;
}

