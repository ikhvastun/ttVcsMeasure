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
using Output::distribs1DForFR;
using Output::distribs2D;

void treeReader::Analyze(){

  double leptonMVAcut = 0.4; // 0.4 for ttZ, 0.8 for 3L tZq
  double magicFactor = 0.85; // 0.85 for ttZ, 0.95 for tZq

  leptonSelection = 3;

  //Set CMS plotting style
  setTDRStyle();
  gROOT->SetBatch(kTRUE);
  //read samples and cross sections from txt file
  //readSamples("data/samples_FOtuning_ttbar.txt"); // 
  //readSamples("data/samples_FOtuning_ttbar_2017.txt"); // 
  //readSamples("data/samples_QCD.txt"); // 
  readSamples("data/samples_QCD_2017.txt"); // 
  initdistribsForFR();

  for(size_t sam = 0; sam < samples.size(); ++sam){
      initSample();

      Color_t color = assignColor(samples[sam].getProcessName());
      setStackColors(color, sam);

      if(samples[sam].is2017() && (samples[sam].getProcessName()) == "loose") continue;
      //if(!(std::get<1>(samples[sam]).find("MuEnriched") != std::string::npos )) continue;
      //if(std::get<1>(samples[sam]).find("MuEnriched") != std::string::npos ) continue;

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
          //if(it > 100000) break;
          //if(it > nEntries / 20) break;
          
          std::vector<unsigned> ind, indFO;
          const unsigned lCount = selectLep(ind, leptonSelection);
          const unsigned lCountFO = selectFakeLep(indFO, leptonSelection);

          // for QCD
          if(lCountFO != 1) continue;

          // for ttbar
          //if(lCountFO < 1) continue;

          int featureCategory = -99;
          int nLocEle = getElectronNumber(indFO);

          std::vector<unsigned> indJets;
          std::vector<unsigned> indBJets;

        
          double mll = 99999;
          double pt_Z = 999999;
          double phi_Z = 999999;
          double ptNonZ = 999999;

          nJLoc = nJets(0, true, indJets, samples[sam].is2017());
          nBLoc = nBJets(0, true, true, indBJets, 1, samples[sam].is2017());
          HTLoc = HTCalc(indJets);
          //if(HTLoc < 200) continue;

          // for QCD, maybe this is not needed in QCD MC
          if(samples[sam].getFileName().find("MuEnriched") != std::string::npos && sam != 0 && _lPt[indFO.at(0)] > 15) continue;
          if(_met > 20) continue;
          TLorentzVector l0p4;
          l0p4.SetPtEtaPhiE(_lPt[indFO.at(0)], _lEta[indFO.at(0)], _lPhi[indFO.at(0)], _lE[indFO.at(0)]);
          //l0p4.SetPtEtaPhiE(35., _lEta[indFO.at(0)], _lPhi[indFO.at(0)], _lE[indFO.at(0)]);
          double mtL = mtCalc(l0p4, _met, _metPhi);
          if(mtL > 20) continue;

          /*
          int qualityCategory = 0;
          int promptCategory = 0;
          for(auto & i : indFO){
            if (lepIsGood(i)) qualityCategory+= 1;
            if(_lIsPrompt[i] && _lMatchPdgId[i] != 22) promptCategory += 1;
          }

          if(promptCategory > 2) continue;
          if(qualityCategory< 2) continue;
          if(qualityCategory> 3) continue;
          */

          for(auto & i : indFO){
            double mvaVL = -999.;
            mvaVL =  _leptonMvatZqTTV[i];

            // for ttbar (also maybe for QCD ??? )
            if(_lIsPrompt[i]) continue;
            if(_lProvenanceCompressed[i] == 0) continue;
            if(_lProvenanceCompressed[i] == 4) continue;
            //if(_lProvenanceCompressed[i] != 1) continue;

            if(lepIsGood(i, leptonSelection)) featureCategory = 0;
            else featureCategory = 1;

            int ptAverageBin = int((mvaVL + 1) * 40);

            //hAveragePt[ptAverageBin]->Fill(double(featureCategory == 0 ? _lPt[i] : _lPt[i] / _ptRatio[i]), weight);

            int additionalFlavourIndex = 0;
            if(_lProvenanceCompressed[i] == 1)
                additionalFlavourIndex = 2;
            if(_lProvenanceCompressed[i] == 2)
                additionalFlavourIndex = 4;
            if(_lProvenanceCompressed[i] == 3)
                additionalFlavourIndex = 6;

            int additionalPositionIndex = 0;
            if(fabs(_lEta[i]) > borderOfBarrelEndcap[_lFlavor[i]])
                additionalPositionIndex = 8;

            distribs1DForFR[_lFlavor[i]].vectorHisto[featureCategory+additionalFlavourIndex+additionalPositionIndex].Fill(TMath::Min(double(featureCategory == 0 ? _lPt[i] : magicFactor * _lPt[i] / _ptRatio[i]),ptBins[nPt-1]-0.001), weight); // 
            distribs2D[_lFlavor[i]].vectorHisto[featureCategory+additionalFlavourIndex].Fill(TMath::Min(double(featureCategory == 0 ? _lPt[i] : magicFactor * _lPt[i] / _ptRatio[i]),ptBins[nPt-1]-0.001), TMath::Abs(_lEta[i]), weight);

            // for total
            distribs1DForFR[_lFlavor[i]].vectorHisto[featureCategory+additionalPositionIndex].Fill(TMath::Min(double(featureCategory == 0 ? _lPt[i] : magicFactor * _lPt[i] / _ptRatio[i]),ptBins[nPt-1]-0.001), weight);
            distribs2D[_lFlavor[i]].vectorHisto[featureCategory].Fill(TMath::Min(double(featureCategory == 0 ? _lPt[i] : magicFactor * _lPt[i] / _ptRatio[i]),ptBins[nPt-1]-0.001), TMath::Abs(_lEta[i]), weight);

            //distribsPtRatio->Fill(TMath::Min(double(magicFactor * _lPt[i] / _ptRatio[i] / _lPt[i]),varMax[16]-0.001), weight);
          }
          
      }

      std::cout << std::endl;
      //cout << "Total number of events: " << distribs[0].vectorHisto[sam].Integral() << endl;
      std::cout << std::endl;
  }

  /*
  for(int i = 0; i < 80; i++){
    cout << hAveragePt[i]->GetMean() << " ";
  }
  */
  cout << endl;

  TLegend* mtleg = new TLegend(0.77,0.89,0.95,0.62);
  //mtleg->SetNColumns(4);
  mtleg->SetFillColor(0);
  mtleg->SetFillStyle(0);
  mtleg->SetBorderSize(0);
  mtleg->SetTextFont(42);

  /*
  //mtleg->AddEntry(&distribs[0].vectorHisto[dataSample],"Data","lep"); //data
  int count = 0;

  for (std::vector<std::string>::iterator it = samplesOrderNames.begin(); it != samplesOrderNames.end(); it++) {

        cout << "count and sample name: " << count << " " << *it << " " << samplesOrder.at(count) << endl;

        mtleg->AddEntry(&distribs[0].vectorHisto[samplesOrder.at(count)],(*it).c_str(),"f");
        count++;
  }
  */

  double scale_num = 1.6;

  TCanvas* plot[2][2];
  for(int i = 0; i < 2; i++) for(int j = 0; j < 2; j++) plot[i][j] = new TCanvas(Form("plot_%d_%d", i, j),"",500,450);
  TCanvas* plot2D[8];
  for(int i = 0; i < 8; i++) plot2D[i] = new TCanvas(Form("plot_2D_%d", i),"",500,450);

  /*
  showHist2D(plot2D[0],distribs2D); 
  plot2D[0]->SaveAs("plotsForSave/correlation.pdf");
  plot2D[0]->SaveAs("plotsForSave/correlation.png");
  plot2D[0]->SaveAs("plotsForSave/correlation.root");
  */
  //showHist2D(plot2D[0],distribs[0], 0); 
  //showHist2D(plot2D[1],distribs[1], 1); 

  for(int flComp = 0; flComp < 4; flComp++){
    for(int flavour = 0; flavour < 2; flavour++){
        showHist2D(plot2D[flavour+2*flComp],distribs2D[flavour].vectorHisto[0+2*flComp], distribs2D[flavour].vectorHisto[1+2*flComp], flavour, flComp); 
    }
  }

  /*
  Color_t colorFL[4] = {kBlack, kRed, kGreen, kBlue};
  TString posString[2] = {"barrel", "endcap"};
  for(int flavor = 0; flavor < 2; flavor++){
    for(int pos = 0; pos < 2; pos++){
        //if(pos == 1) continue;
        double xmin = distribs[flavor].vectorHisto[0].GetXaxis()->GetXmin();
        double xmax = distribs[flavor].vectorHisto[0].GetXaxis()->GetXmax();

        plot[flavor][pos]->cd();
        double xPad = 0.25; // 0.25

        TPad *pad1 = new TPad("pad1","pad1",0,xPad,1,1);
        pad1->SetTopMargin(0.07);
        if(xPad != 0)
            pad1->SetBottomMargin(0.02);
        pad1->Draw();
        pad1->cd();

        for(int flComp = 0; flComp < 4; flComp++){
            distribs1DForFR[flavor].vectorHisto[2*flComp+8*pos].Divide(&distribs[flavor].vectorHisto[2*flComp+1+8*pos]);
            distribs1DForFR[flavor].vectorHisto[2*flComp+8*pos].SetLineColor(colorFL[flComp]);
            distribs1DForFR[flavor].vectorHisto[2*flComp+8*pos].SetMarkerColor(colorFL[flComp]);
            distribs1DForFR[flavor].vectorHisto[2*flComp+8*pos].SetTitle(flavorsString[flavor] + " FR (" + posString[pos] + ")");
            distribs1DForFR[flavor].vectorHisto[2*flComp+8*pos].GetXaxis()->SetTitle("p_{T}^{corr} [GeV]"); 
            distribs1DForFR[flavor].vectorHisto[2*flComp+8*pos].GetYaxis()->SetTitle("FR"); 
            distribs1DForFR[flavor].vectorHisto[2*flComp+8*pos].SetMinimum(0.);
            distribs1DForFR[flavor].vectorHisto[2*flComp+8*pos].SetMaximum(0.3);
            //distribs[flavor].vectorHisto[2*flComp+8*pos].Scale(1. / distribs[flavor].vectorHisto[2*flComp+8*pos].Integral());
            distribs1DForFR[flavor].vectorHisto[2*flComp+8*pos].Draw("same");
            if(flavor == 0 && pos == 0)
                mtleg->AddEntry(&distribs1DForFR[flavor].vectorHisto[2*flComp+8*pos], flavorComposString[flComp], "l");
        }
        */
        /*
        for(int flComp = 1; flComp < 4; flComp++){
            for(int pos = 0; pos < 2; pos++){
                distribs[flavor].vectorHisto[0].Add(&distribs[flavor].vectorHisto[2*flComp+8*pos]);
                //distribs[flavor].vectorHisto[1].Add(&distribs[flavor].vectorHisto[1+2*flComp+8*pos]);
            }
        }
        distribs[flavor].vectorHisto[0].Scale(1. / distribs[flavor].vectorHisto[0].Integral());
        distribs[flavor].vectorHisto[0].Draw("same");
        //distribs[flavor].vectorHisto[1].Scale(1. / distribs[flavor].vectorHisto[1].Integral());
        //distribs[flavor].vectorHisto[1].Draw("same");
        //
        */
        /*
        mtleg->Draw("same");

        pad1->cd();
        pad1->RedrawAxis();
        pad1->Update();

        if(xPad == 0) {
            plot[flavor][pos]->SaveAs("plotsForSave/FR_" + flavorsString[flavor] + "_" + posString[pos] + ".root");
            //distribs[flavor].vectorHisto[0].SaveAs("plotsForSave/check_" + flavorsString[flavor] + ".root");
            //distribs[flavor].vectorHisto[1].SaveAs("plotsForSave/check_" + flavorsString[flavor] + ".root");
            continue;
        }

        plot[flavor][pos]->cd();

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,xPad);

        pad2->SetBottomMargin((1.-xPad)/xPad*0.13);
        pad2->SetTopMargin(0.06);

        pad2->Draw();
        pad2->RedrawAxis();
        pad2->cd();

        TH1D *flClones[4];
        for(int flComp = 0; flComp < 4; flComp++){
           flClones[flComp] = (TH1D*)distribs[flavor].vectorHisto[2*flComp+8*pos].Clone(Form("clone_%d_%d_%d", flavor, pos, flComp));
        }
        flClones[0]->Divide(flClones[1]);
        flClones[2]->Divide(flClones[1]);
        flClones[3]->Divide(flClones[1]);

        flClones[0]->SetTitle("");
        flClones[0]->GetXaxis()->SetTitle("p_{T}^{corr} [GeV]");
        flClones[0]->GetYaxis()->SetTitle("X / b");

        flClones[0]->GetYaxis()->SetTitleOffset(1.2/((1.-xPad)/xPad));
        flClones[0]->GetYaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
        flClones[0]->GetXaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
        flClones[0]->GetYaxis()->SetLabelSize((1.-xPad)/xPad*0.05);
        flClones[0]->GetXaxis()->SetLabelSize((1.-xPad)/xPad*0.05);

        flClones[0]->Draw("axis");
        flClones[0]->SetMaximum(2.0);
        flClones[0]->SetMinimum(0.0);

        TLine *line = new TLine(xmin, 1, xmax, 1);
        line->SetLineStyle(2);

        line->Draw("same");

        flClones[0]->Draw("same");
        flClones[2]->Draw("same");
        flClones[3]->Draw("same");

        plot[flavor][pos]->SaveAs("plotsForSave/FR_" + flavorsString[flavor] + "_" + posString[pos] + ".root");
    }

  }
  */

  /*
  TCanvas * plotPtRatio = new TCanvas("plotPtRatio", "plotPtRatio", 500, 450);
  distribsPtRatio->Draw();
  plotPtRatio->SaveAs("plotsForSave/ptRatioToPt.pdf");
  plotPtRatio->SaveAs("plotsForSave/ptRatioToPt.root");
  */

  return;

}

int main(int argc, char *argv[]){

    TApplication *rootapp = new TApplication("example", &argc, argv);
    treeReader reader;
    reader.Analyze();
    rootapp->Run();

    return 0;
}

