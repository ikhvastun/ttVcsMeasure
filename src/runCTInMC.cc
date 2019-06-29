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
#include "../interface/kinematicTools.h"
#include "../interface/fillDatacards.h"

#include "tdrStyle.C"

using namespace std;
using namespace tools;

Errors LastError::lasterror = Errors::UNKNOWN;
using Output::distribs1DForCT;

void treeReader::Analyze(){

  leptonSelection = 3;
  magicFactor = 0.85;
  leptonMVAcut = 0.4;
  //Set CMS plotting style
  setTDRStyle();
  gROOT->SetBatch(kTRUE);
  //read samples and cross sections from txt file
  //readSamples("data/samples_FOtuning.txt"); // 
  //readSamples("data/samples_FOtuning_ttbar.txt"); // 
  readSamples("data/samples_FOtuning_ttbar_2017.txt"); // 
  
  std::string selection = "CTInMC";
  initListToPrint(selection);
  initdistribsForCT();

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample("ttZ");

      // only for 3L
      //if(std::get<1>(samples[sam]).find("DY") != std::string::npos ) continue;
      
      if(samples[sam].is2017() && samples[sam].getProcessName().find("TTToSemiLeptonic") != std::string::npos ) continue;
      if(samples[sam].is2017() && samples[sam].getProcessName().find("TTTo2L2Nu") != std::string::npos ) continue;

      if(samples[sam].getProcessName().find("TTJets_SingleLeptFromT") != std::string::npos ) continue;
      if(samples[sam].getProcessName().find("TTJets_DiLept") != std::string::npos ) continue;

      std::cout<<"Entries in "<< samples[sam].getProcessName() << " " << nEntries << std::endl;
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
          //if(it > nEntries / 20) break;
          
          std::vector<unsigned> ind, indFO;
          const unsigned lCountFO = selectFakeLep(indFO, leptonSelection);

          //if(!(_2017_e || _2017_m || _2017_ee || _2017_em || _2017_mm)) continue;
          //cout << "decision: " << _2017_ee << " " << _2017_mm << " " << _2017_em << " " << _2017_e << " " << _2017_m << endl;

          if(leptonSelection == 3){
            if(lCountFO < 3) continue;
            if(!passPtCuts3L(indFO)) continue;
          }
          if(leptonSelection == 2){
            if(lCountFO != 2) continue;
            if(!passPtCuts2L(indFO)) continue;
            if(_lCharge[indFO.at(0)] * _lCharge[indFO.at(1)] < 0) continue;
          }
          
          int nLocEle = getElectronNumber(indFO);
          //if(nLocEle != 2) continue;
          //if(nLocEle == 2 || nLocEle == 0) continue;
          //if(!noConversionInSelection(indFO)) continue;

          std::vector<unsigned> indJets;
          std::vector<unsigned> indBJets;

          int featureCategory = 0;
          int promptCategory = 0;
          for(auto & i : indFO){
            if (lepIsGood(i, leptonSelection)) featureCategory += 1;
            //if(_lIsPrompt[i] && _lMatchPdgId[i] != 22) promptCategory += 1;
            if(_lIsPrompt[i]) promptCategory += 1;
            //if(leptonIsPrompt(i)) promptCategory += 1;
            /*
            cout << "lepton info: " << _lPt[i] << " " << _lEta[i] << " " << _lPhi[i] << " " << _lFlavor[i] << " " << _lIsPrompt[i] << " " << _lProvenance[i] << " " << _lProvenanceCompressed[i] << endl;
            for(int j = 0; j < _gen_nL; j++){
                cout << "all gen leptons in the event: " << _gen_lPt[j] << " " << _gen_lEta[j] << " " << _gen_lPhi[j] << " " << _gen_lFlavor[j] << " " << _gen_lIsPrompt[j] << endl;
            }
            */
          }
        
          //cout << "prompt cat: " <<promptCategory << endl;
          //if(promptCategory == 2)
          //    cout << _runNb << " " << _lumiBlock << " " << _eventNb << endl; 
          
          if(leptonSelection == 3 && promptCategory > 2) continue;
          if(leptonSelection == 2 && promptCategory > 1) continue;
          //if(featureCategory < 2) continue;
          //if(featureCategory > 3) continue;
          double FRloc = 1.;

          //if(featureCategory == 3 && !eventChargeConsistent(indFO)) continue;
          //if(featureCategory < 2){ 
          bool promptInSideband = false;
          if(featureCategory < leptonSelection){ 
            weight *= fakeRateWeight();
          }

          if(promptInSideband) continue;

          //nJLoc = nJets(0, true, indJets, featureCategory < leptonSelectionAnalysis);
          //nBLoc = nBJets(0, true, true, 1, featureCategory < leptonSelectionAnalysis);
          unsigned third = -9999;
          double mll = 99999;
          double mlll = 99999;
          double ptZ = 999999;
          double ptNonZ = -999999;
          TLorentzVector Zboson, lnegative;
          std::vector<unsigned> indOf2LonZ;
          double dMZ = deltaMZ(ind, third, mll, ptZ, ptNonZ, mlll, indOf2LonZ, Zboson, lnegative);

          if(leptonSelection == 2){
            //if(nJLoc < 2) continue;
            //if(nBLoc < 1) continue;

            TLorentzVector l0p4, l1p4;
            l0p4.SetPtEtaPhiE(ptCorrV[0].first, _lEta[indFO.at(0)], _lPhi[indFO.at(0)], _lE[indFO.at(0)] * ptCorrV[0].first / _lPt[indFO.at(0)]);
            l1p4.SetPtEtaPhiE(ptCorrV[1].first, _lEta[indFO.at(1)], _lPhi[indFO.at(1)], _lE[indFO.at(1)] * ptCorrV[1].first / _lPt[indFO.at(1)]);

            double ele_mll = (l0p4+l1p4).M();

            if(ele_mll < 12) continue;
            if(ele_mll > 76 && ele_mll < 106 && nLocEle == 2) continue;
            if(_met < 30) continue;
          }

          double maxMJetJet = -999.;
          TLorentzVector recoilingJet(0,0,0,0);
          TLorentzVector taggedBJet(0,0,0,0);
          unsigned jetCount = 0;
          unsigned bJetCount = 0;
          double minDeltaRLeptonbJet = 0.;
          if(leptonSelection == 3){
             //if(dMZ > 10) continue;
             //make lorentzvectors for leptons
             TLorentzVector lepV[lCountFO];
             for(unsigned l = 0; l < lCountFO; ++l) lepV[l].SetPtEtaPhiE(ptCorrV[l].first, _lEta[indFO.at(l)], _lPhi[indFO.at(l)], _lE[indFO.at(l)] * ptCorrV[l].first / _lPt[indFO.at(l)]);

             // met value
             TLorentzVector met;
             met.SetPtEtaPhiE(_met, 0, _metPhi, _met);

             // fill njets and nbjets variables
             std::vector<unsigned> jetInd, bJetInd;
             jetCount = nJets(0, true, indJets, samples[sam].is2017());
             bJetCount = nBJets(0, true, true, indBJets, 1, samples[sam].is2017());

             //if(jetCount < 2) continue;
             // fill vector of all jets
             TLorentzVector jetV[(const unsigned) _nJets];
             for(unsigned j = 0; j < _nJets; ++j) jetV[j].SetPtEtaPhiE(_jetPt[j], _jetEta[j], _jetPhi[j], _jetE[j]);

             std::vector<unsigned> taggedJetI; //0 -> b jet from tZq, 1 -> forward recoiling jet
             TLorentzVector neutrino = findBestNeutrinoAndTop(lepV[ptCorrV[0].second], met, taggedJetI, jetInd, bJetInd, jetV);

             if(taggedJetI[0] != 99) taggedBJet = jetV[taggedJetI[0]];
             if(taggedJetI[1] != 99) recoilingJet = jetV[taggedJetI[1]];
             maxMJetJet = kinematics::maxMass(jetV, jetInd);

             std::vector<unsigned> bjetVecInd;
             for(unsigned l = 0; l < bJetCount; ++l) bjetVecInd.push_back(bJetInd[l]);

             std::vector<unsigned> lepVecInd;
             for(unsigned l = 0; l < lCountFO; ++l) lepVecInd.push_back(indFO[l]);

             minDeltaRLeptonbJet = kinematics::minDeltaR(lepV, lepVecInd, jetV, bjetVecInd);
             HTLoc = HTCalc(jetInd);

             //cout << "new values are " << fabs(recoilingJet.Eta()) << " " << std::max(maxMJetJet, 0.) << endl;
          }

          double mt1 = 9999;
          double mvaVL = 0;
          double mll1stpair = -999., mll2ndpair = -999.;
          double cosTSt = -999;

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
                                   (leptonSelection == 3 && passTTZCleanSelection(nJLoc, nBLoc, dMZ) ? SRIDPTZ(ptZ) : -999), (leptonSelection == 3 && passTTZCleanSelection(nJLoc, nBLoc, dMZ) ? SRIDCosTheta(cosTSt) : -999),
                                   (leptonSelection == 3 ? flavourCategory3L(nLocEle) : -999),
                                   (leptonSelection == 4 ? flavourCategory4L(nLocEle) : -999),
                                   (leptonSelection == 4 ? flavourCategory4LZZ(nLocEle) : -999),
                                   (leptonSelection == 3 && nLocEle == 0 ? SRID3L(nJLoc, nBLoc, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 1 ? SRID3L(nJLoc, nBLoc, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 2 ? SRID3L(nJLoc, nBLoc, dMZ) : -999),
                                   (leptonSelection == 3 && nLocEle == 3 ? SRID3L(nJLoc, nBLoc, dMZ) : -999),
                                   (passTTZSRSelection(ind, indOf2LonZ, nJLoc, nBLoc, dMZ) ? flavourCategory3L4L(leptonSelection, nLocEle) : -999),
                                   (leptonSelection == 3 && nBLoc > 0 ? SRID8SR3L(nJLoc, nBLoc, dMZ) : -999),
                                   (leptonSelection != 4 ? mll:mll1stpair),
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
                                     "flavour3L", "flavour4L", "flavour4LZZ", 
                                     "SR3L3m","SR3L2m1e","SR3L1m2e","SR3L3e",
                                     "flavour3L4L",
                                     "SRTTZ8SR3L",
                                     "mllnoZcut"
                                   };

          for(int dist = 0; dist < fillVar.size(); dist++){
            if(std::find(listToPrint[selection].begin(), listToPrint[selection].end(), fncName[dist]) == listToPrint[selection].end()) continue;
            distribs1DForCT[dist].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(TMath::Min(fillVar.at(dist),figNames[fncName.at(dist)].varMax-0.1),weight);
          }

          //distribs[0].vectorHisto[featureCategory >= 2 ? 0 : 1].Fill(TMath::Min(mtHighest,varMax[0]-0.1), weight);
          /*
          distribs[1].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(TMath::Min(_met,varMax[1]-0.1), weight);
          distribs[2].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(TMath::Min(HTLoc,varMax[2]-0.1), weight);
          distribs[3].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(TMath::Min(double(jetCount),varMax[3]-0.1), weight);
          distribs[4].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(TMath::Min(double(bJetCount),varMax[4]-0.1), weight);
          distribs[5].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(TMath::Min(ptCorrV[0].first,varMax[5]-0.1), weight);
          distribs[6].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(TMath::Min(ptCorrV[1].first,varMax[6]-0.1), weight);
          if(leptonSelection == 3)
            distribs[7].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(TMath::Min(ptCorrV[2].first,varMax[7]-0.1), weight);
          distribs[8].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(TMath::Min(_lEta[ptCorrV[0].second],varMax[8]-0.01), weight);
          distribs[9].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(TMath::Min(_lEta[ptCorrV[1].second],varMax[9]-0.01), weight);
          if(leptonSelectionAnalysis == 3)
            distribs[10].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(TMath::Min(_lEta[ptCorrV[2].second],varMax[10]-0.01), weight);
          distribs[11].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(TMath::Min(mll,varMax[11]-0.1), weight);
          //distribs[12].vectorHisto[featureCategory >= 2 ? 0 : 1].Fill(SRID(nJLoc, nBLoc, mvaValueRegion, _charges[maxPtInd]), FRloc*scale*_weight);
          if(leptonSelection == 2)
            distribs[13].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(flavourCategory2L(nLocEle, _lCharge[indFO.at(0)]), weight);
          if(leptonSelection == 3)
            distribs[13].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(flavourCategory3L(nLocEle), weight);
          if(leptonSelection == 3){
            distribs[20].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(fabs(recoilingJet.Eta()) == 0 ? -999. : fabs(recoilingJet.Eta()), weight); 
            distribs[21].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(std::max(maxMJetJet, 0.), weight); 
            distribs[22].vectorHisto[featureCategory >= leptonSelection ? 0 : 1].Fill(std::max(minDeltaRLeptonbJet, 0.), weight); 
          }
          */

      }

      std::cout << std::endl;
      //cout << "Total number of events: " << distribs[1].vectorHisto[sam].Integral() << endl;
      //std::cout << std::endl;
  }

  int indexOfFirstKinVar = figNames[listToPrint[selection].at(0)].index;

  cout << "Total number expected: " << distribs1DForCT[indexOfFirstKinVar].vectorHisto[0].Integral() << endl;
  cout << "Total number predicted: " << distribs1DForCT[indexOfFirstKinVar].vectorHisto[1].Integral() << endl;

  TLegend* mtleg = new TLegend(0.17,0.89,0.95,0.72);
  mtleg->SetNColumns(2);
  mtleg->SetFillColor(0);
  mtleg->SetFillStyle(0);
  mtleg->SetBorderSize(0);
  mtleg->SetTextFont(42);


  //mtleg->AddEntry(&distribs[0].vectorHisto[dataSample],"Data","lep"); //data
  int count = 0;

  /*
  for (std::vector<std::string>::iterator it = samplesOrderNames.begin(); it != samplesOrderNames.end(); it++) {

        cout << "count and sample name: " << count << " " << *it << " " << samplesOrder.at(count) << endl;

        mtleg->AddEntry(&distribs[0].vectorHisto[samplesOrder.at(count)],(*it).c_str(),"f");
        count++;
  }
  */
  mtleg->AddEntry(&distribs[0].vectorHisto[0],"Monte Carlo","lep");
  mtleg->AddEntry(&distribs[0].vectorHisto[1],"Tight-to-loose prediction","f");

  for (int i=0; i!=nVars; ++i)  {

    //cout << "The sample size is " << samples.size() << endl;
    for (int sam=0; sam != samples.size(); ++sam){
      if(samples[sam].getProcessName() == "data") continue;

      // is used in MC CT
      if(sam == 0) continue;

      //cout << "the sample is added: " << std::get<0>(samples[sam]) << endl;
      distribs[i].vectorHistoTotalUnc.Add(&distribs[i].vectorHisto[sam]);
    }

    for (int ibin = 1; ibin!=nBins[i]+1; ++ibin) {
      //get syst. uncertainty band:
      double err = 0.;
      for (int sam=0; sam != samples.size(); ++sam) {
        if(samples[sam].getProcessName() == "data") continue;
        // is used in MC CT
        if(sam == 0) continue;
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

  double scale_num = 1.4;

  TCanvas* plot[nVars];
  TCanvas* plot2D[nVars];

  for(int i = 0; i < nVars; i++) plot[i] = new TCanvas(Form("plot_%d", i),"",500,450);
  for(int i = 0; i < 2; i++) plot2D[i] = new TCanvas(Form("plot_2D_%d", i),"",500,450);

  vector<std::string> figNames = {"Leading lepton M_{T} [GeV]", "Missing E_{T} [GeV]", "H_{T} [GeV]", "N_{jets}", "N_{b jets}", "Leading lepton p_{T}^{corr} [GeV]", "Sub-leading p_{T}^{corr} [GeV]", "Trailing p_{T}^{corr} [GeV]", "Leading lepton #eta", "Sub-leading lepton #eta", "Trailing lepton #eta", "M_{ll} [GeV]", "SR", "Lepton flavors", "Sub-leading lepton M_{T} [GeV]", "Leading jet p_{T} [GeV]", "Sub-leading jet p_{T} [GeV]", "Min #Delta R(jet, trailing lepton)", "BDT score", "p_{T}^{corr} / p_{T}", "Recoiling jet #eta", "M_{jet-jet}^{max} [GeV]", "min #Delta R(lep, b jet)"};
  vector<TString> namesForSaveFiles = {"mtlead", "met", "HT", "njets", "nbjets", "ptCorrLead", "ptCorrSubLead", "ptCorrTrail", "etaCorrLead", "etaCorrSubLead", "etaCorrTrail", "mll", "SR", "flavour", "mtsublead", "jetptlead", "jetpttrail", "mindeltaR", "mvaVL", "ptcorToPt", "etaJetRecoil", "MjjMax", "minDeltaRLepBjet"};

  for(int varPlot = 0; varPlot < nVars; varPlot++){
    //if(varPlot == 0 || varPlot == 7 || varPlot == 10 || varPlot == 12 || varPlot > 13) continue;
    if(varPlot == 0 || varPlot == 12 || (varPlot > 13 && varPlot < 19)) continue;
    plot[varPlot]->cd();
    showHist(plot[varPlot],distribs[varPlot],"",figNames.at(varPlot),"Events", scale_num, mtleg, false, false, dataLumi); // + std::to_string(int((varMax[varPlot] - varMin[varPlot])/nBins[varPlot]))
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + ".pdf");
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + ".png");
    plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + ".root");
    //plot[varPlot]->cd();
    //showHist(plot[varPlot],distribs[varPlot],"",figNames.at(varPlot),"Events", scale_num, mtleg, true, false); //  + std::to_string(int((varMax[varPlot] - varMin[varPlot])/nBins[varPlot]))
    //plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + "Log.pdf");
    //plot[varPlot]->SaveAs("plotsForSave/" + namesForSaveFiles.at(varPlot) + "Log.png");
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

