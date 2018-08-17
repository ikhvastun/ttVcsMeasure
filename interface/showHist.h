#ifndef showHist_H
#define showHist_H

#include "TPad.h"
#include "TH1.h"
#include "TLine.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFrame.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TExec.h"
#include "TColor.h"

#include "CMS_lumi.h"
#include "Output.h"
#include "readTreeSync.h"

const int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV
const int iPos =0;

using namespace std;

using Output::distribs;
using Output::DistribsAll;

void setUpRatioFeatures(TH1D *, TGraphAsymmErrors *, histInfo & info, double);
void setUpSystUnc(DistribsAll &, TH1D *);
void calculateRatioUnc(TGraphAsymmErrors *, TH1D *, TH1D *);
void printInfoOnPlot3L();
void showHist(TVirtualPad* c1, DistribsAll & distribs, histInfo & info, double num, TLegend *leg, bool plotInLog = false, bool normalizedToData = false, const int showLegendOption = 0){ // showLegendOption 0 - 2016, 1 - 2017, 2 - 2016+2017
    double xPad = 0.25; // 0.25

    TPad *pad1 = new TPad("pad1","pad1",0,xPad,1,1);
    pad1->SetTopMargin(0.07);
    if(xPad != 0)
        pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    if(plotInLog)
        pad1->SetLogy();
    
    TH1D * dataHist = &distribs.vectorHisto[dataSample];
    //if(is2017) // keep it blinded for 2017
    //    dataHist->Reset("ICE");
    dataHist->SetMarkerSize(1);
    dataHist->SetTitle("");
    dataHist->GetXaxis()->SetTitle(info.fancyName.c_str());
    dataHist->GetYaxis()->SetTitle(("Events " + (info.isEnVar ? ("/ " + std::to_string(int((info.varMax - info.varMin) / info.nBins)) + " GeV") : "")).c_str());
    dataHist->SetMinimum(0.01);
    dataHist->SetMaximum(TMath::Max(distribs.stack.GetMaximum(), distribs.vectorHisto[dataSample].GetMaximum()) * num);
    if(plotInLog)
        dataHist->SetMaximum(TMath::Max(distribs.stack.GetMaximum(), distribs.vectorHisto[dataSample].GetMaximum()) * num * 5);
    dataHist->GetXaxis()->SetLabelOffset(0.02);
    
    dataHist->Draw("E0");
    //distribs.stack.Draw("histsame");
    //dataHist->Draw("E0same");

    leg->Draw("same");
    double lumi = 35.9;
    if(showLegendOption == 1) lumi = 41.9;
    else if (showLegendOption == 2) lumi = 77.8;
    CMS_lumi( pad1, iPeriod, iPos, lumi);

    pad1->cd();
    pad1->RedrawAxis();
    pad1->Update();

    if(xPad == 0) return;

    c1->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,xPad);

    pad2->SetBottomMargin((1.-xPad)/xPad*0.13);
    pad2->SetTopMargin(0.06);

    pad2->Draw();
    pad2->RedrawAxis();
    pad2->cd();

    TH1D *stackCopy = (TH1D*)(distribs.stack.GetStack()->Last())->Clone("stackCopy");
    TH1D *dataCopy = (TH1D*)dataHist->Clone("dataCopy");
    dataCopy->Divide(stackCopy);
    TGraphAsymmErrors * dataCopyGraph = new TGraphAsymmErrors(dataCopy);
    // calculate asymmetric uncertainties for data
    calculateRatioUnc(dataCopyGraph, dataHist, stackCopy);

    setUpRatioFeatures(stackCopy, dataCopyGraph, info, xPad);

    TH1D *histSystAndStatUnc = (TH1D*)(distribs.stack.GetStack()->Last())->Clone(Form("histSystAndStatUnc"));
    setUpSystUnc(distribs, histSystAndStatUnc);

    TLegend* mtlegRatio = new TLegend(0.17,0.39,0.85,0.58);
    mtlegRatio->SetNColumns(4);
    mtlegRatio->SetFillColor(0);
    mtlegRatio->SetFillStyle(0);
    mtlegRatio->SetBorderSize(0);
    mtlegRatio->SetTextFont(42);

    mtlegRatio->AddEntry(stackCopy, "Stat", "f");
    mtlegRatio->AddEntry(histSystAndStatUnc, "Stat+Syst", "f");

    // Draw finally the things
    stackCopy->Draw("axis");
    histSystAndStatUnc->Draw("e2same");
    stackCopy->Draw("e2same");

    double xmin = distribs.vectorHisto[dataSample].GetXaxis()->GetXmin();
    double xmax = distribs.vectorHisto[dataSample].GetXaxis()->GetXmax();
    TLine *line = new TLine(xmin, 1, xmax, 1);
    line->SetLineStyle(2);
    line->Draw("same");

    mtlegRatio->Draw("same");
    //if(!is2017)
    dataCopyGraph->Draw("p"); // dataCopyGraph = data / MC stack

    pad2->RedrawAxis();
    pad2->Update();

    c1->cd();
    pad1->cd();
    TH1D *systAndStatUnc = (TH1D*)(distribs.stack.GetStack()->Last())->Clone("systAndStatUnc");
    for(unsigned int i = 0; i < systAndStatUnc->GetNbinsX(); i++){
        systAndStatUnc->SetBinError(i+1, histSystAndStatUnc->GetBinError(i+1) * systAndStatUnc->GetBinContent(i+1));
    }
    distribs.stack.Draw("histsame");
    systAndStatUnc->SetFillStyle(3005);
    systAndStatUnc->SetFillColor(kGray+2);
    systAndStatUnc->SetMarkerStyle(1);
    systAndStatUnc->Draw("e2same");
    dataHist->Draw("E0same");

    if(info.index == indexSR3L){
       dataHist->GetXaxis()->SetTitleSize(0.07);
       dataHist->GetXaxis()->SetTitleOffset(0.8);
       if(leptonSelectionAnalysis == 3)
         printInfoOnPlot3L();
    }

    pad1->cd();
    pad1->RedrawAxis();
    pad1->Update();

}

void printInfoOnPlot3L(){

    TLine *line1 = new TLine(3.5, 0, 3.5, 1125);
    line1->SetLineStyle(2);
    line1->Draw("same");

    TLine *line2 = new TLine(5.5, 0, 5.5, 1125);
    line2->SetLineStyle(2);
    
    // need only one for the moment
    //line2->Draw("same");

    TLatex nbjetsEq0region;
    nbjetsEq0region.SetNDC();
    nbjetsEq0region.SetTextAngle(0);
    nbjetsEq0region.SetTextColor(kBlack);

    nbjetsEq0region.SetTextFont(42);
    nbjetsEq0region.SetTextAlign(31);
    nbjetsEq0region.SetTextSize(0.05);
    //nbjetsEq0region.DrawLatex(0.35, 0.56,"N_{b} = 0");

    TLatex nbjetsEq1region;
    nbjetsEq1region.SetNDC();
    nbjetsEq1region.SetTextAngle(0);
    nbjetsEq1region.SetTextColor(kBlack);

    nbjetsEq1region.SetTextFont(42);
    nbjetsEq1region.SetTextAlign(31);
    nbjetsEq1region.SetTextSize(0.05);
    nbjetsEq1region.DrawLatex(0.64, 0.56,"N_{b} = 1");

    TLatex nbjetsEq2region;
    nbjetsEq2region.SetNDC();
    nbjetsEq2region.SetTextAngle(0);
    nbjetsEq2region.SetTextColor(kBlack);

    nbjetsEq2region.SetTextFont(42);
    nbjetsEq2region.SetTextAlign(31);
    nbjetsEq2region.SetTextSize(0.05);
    nbjetsEq2region.DrawLatex(0.91, 0.56,"N_{b} > 1");
}

void setUpRatioFeatures(TH1D * stackCopy, TGraphAsymmErrors * dataCopyGraph, histInfo & info, double xPad){

    // this one will be used on the ratio plot
    stackCopy->Divide(stackCopy);

    // if there is 0 event in stack, then set uncertainty to 0 
    //for(int bin = 1; bin < stackCopy->GetNbinsX() + 1; bin++)
    //    stackCopy->SetBinError(bin, 0.);

    stackCopy->SetFillStyle(1001);
    stackCopy->SetFillColor(kCyan - 4);
    stackCopy->SetMarkerStyle(1);

    stackCopy->SetTitle("");
    stackCopy->GetXaxis()->SetTitle(info.fancyName.c_str());
    stackCopy->GetYaxis()->SetTitle("data/pred");

    stackCopy->GetYaxis()->SetTitleOffset(1.2/((1.-xPad)/xPad));
    stackCopy->GetYaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
    stackCopy->GetXaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
    stackCopy->GetYaxis()->SetLabelSize((1.-xPad)/xPad*0.05);
    if(info.index != indexSR3L && info.index != indexSR4L && info.index != indexSRTTZ) //  && info.index != indexFlavour
        stackCopy->GetXaxis()->SetLabelSize((1.-xPad)/xPad*0.05);
    else
        stackCopy->GetXaxis()->SetLabelSize(0.25);

    stackCopy->SetMaximum(2.0);
    stackCopy->SetMinimum(0.0);
    stackCopy->SetMarkerStyle(20);
    stackCopy->SetMarkerSize(0.2);

    dataCopyGraph->SetMarkerSize(0.5);

}

void calculateRatioUnc(TGraphAsymmErrors * dataGraph, TH1D * data, TH1D * stack){

    for(int i = 0; i < dataGraph->GetN(); i++){

      double dataPoint[3] = {data->GetBinContent(i+1), data->GetBinErrorUp(i+1), data->GetBinErrorLow(i+1)};
      double theMCPoint[3] = {stack->GetBinContent(i+1), stack->GetBinErrorUp(i+1), stack->GetBinErrorLow(i+1)};

      double uncRatio[2];

      // calculating the uncertainty for a / b
      // using formula: delta Unc ^ 2 = ( (a/b)'_a (delta a) ) ^ 2 + ( (a/b)'_b (delta b) ) ^ 2

      uncRatio[0] = TMath::Sqrt(TMath::Power(1 / theMCPoint[0] * dataPoint[1], 2) + TMath::Power(dataPoint[0] / TMath::Power(theMCPoint[0],2) * theMCPoint[1], 2));
      uncRatio[1] = TMath::Sqrt(TMath::Power(1 / theMCPoint[0] * dataPoint[2], 2) + TMath::Power(dataPoint[0] / TMath::Power(theMCPoint[0],2) * theMCPoint[2], 2));

      dataGraph->SetPointError(i, dataGraph->GetErrorXlow(i), dataGraph->GetErrorXhigh(i), uncRatio[1], uncRatio[0]);
    }
}

void setUpSystUnc(DistribsAll & distribs, TH1D * histSystAndStatUnc){

    // histSystAndStatUnc - a histogram with central value at 1 and with applied one of the uncertainties on top: JEC, JES and Uncl
    //TH1D *histSystAndStatUnc = (TH1D*)(distribs.stack.GetStack()->Last())->Clone(Form("histSystAndStatUnc"));
    TH1D *stackCopy = (TH1D*)(distribs.stack.GetStack()->Last())->Clone("stackCopy");
    // stack of MC with varied up and down of 3 different types of uncertainties
    TH1D *stackUncUp[numberOfSyst];
    TH1D *stackUncDown[numberOfSyst];

    for(unsigned int i = 0; i < numberOfSyst; i++){

        stackUncUp[i] = (TH1D*)distribs.vectorHisto[1].Clone(Form("stackUncUp_%d", i));
        stackUncDown[i] = (TH1D*)distribs.vectorHisto[1].Clone(Form("stackUncDown_%d", i));

        stackUncUp[i]->Reset("ICE");
        stackUncDown[i]->Reset("ICE");

        for(unsigned int j = distribs.vectorHistoUncUp.size()-1; j != 0; j--){
            stackUncUp[i]->Add((TH1D*)&distribs.vectorHistoUncUp[j].unc[i]);
            stackUncDown[i]->Add((TH1D*)&distribs.vectorHistoUncDown[j].unc[i]);
        }
    }

    for(unsigned int i = 0; i < histSystAndStatUnc->GetNbinsX(); i++){
        // content in particular bin in stack
        double stackBinContent = ((TH1D*)distribs.stack.GetStack()->Last())->GetBinContent(i+1);
        double stackBinError = ((TH1D*)distribs.stack.GetStack()->Last())->GetBinError(i+1);

        if(stackBinContent != 0.)
            stackCopy->SetBinError(i+1, stackBinError / stackBinContent);
        else
            stackCopy->SetBinError(i+1, 0.);
        double err = TMath::Power(stackCopy->GetBinError(i+1), 2);
        //cout << "stat unc is " << TMath::Sqrt(err) << endl;

        for(unsigned int j = 0; j < numberOfSyst; j++){
            // consider largest deviation between the upward and downward variations
            err += TMath::Power(TMath::Max(stackUncUp[j]->GetBinContent(i+1) - stackBinContent, stackBinContent - stackUncDown[j]->GetBinContent(i+1)) / stackBinContent, 2);
            histSystAndStatUnc->SetBinContent(i+1, 1.);
            // if uncertainty is greater than 100% consider 100% uncertainty
            if(stackBinContent != 0.){
                //cout << "unc after applying " << j << " syst: " << TMath::Sqrt(err) << endl;
                histSystAndStatUnc->SetBinError(i+1, TMath::Sqrt(err) > 1 ? 1. : TMath::Sqrt(err));
            }
            //else
            //    histSystAndStatUnc->SetBinError(i+1, 0.);
        }
    }

    histSystAndStatUnc->SetFillStyle(1001);
    histSystAndStatUnc->SetFillColor(kOrange - 4);
    histSystAndStatUnc->SetMarkerStyle(1);
    histSystAndStatUnc->Draw("same");
}

#endif  // showHist
