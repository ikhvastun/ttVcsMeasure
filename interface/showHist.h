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
#include "CMS_lumi.h"
#include "Output.h"
#include "readTreeSync.h"
#include "TExec.h"
#include "TColor.h"

const int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV
const int iPos =0;

using namespace std;

using Output::distribs;
using Output::DistribsAll;

void showHist(TVirtualPad* c1, DistribsAll & distribs, string title, string titleX, string titleY, double num, TLegend *leg, bool plotInLog = false, bool normalizedToData = false, double lumi = 35.9){
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
    // stack here corresponds to all MC
    double xmin = distribs.vectorHisto[dataSample].GetXaxis()->GetXmin();
    double xmax = distribs.vectorHisto[dataSample].GetXaxis()->GetXmax();
    //pad1->DrawFrame(xmin, -0.1, xmax, 1.1);
    
    dataHist->SetMarkerSize(1);
    dataHist->SetTitle(title.c_str());
    dataHist->GetXaxis()->SetTitle(titleX.c_str());
    dataHist->GetYaxis()->SetTitle(titleY.c_str());
    dataHist->SetMinimum(0.01);
    dataHist->SetMaximum(TMath::Max(distribs.stack.GetMaximum(), distribs.vectorHisto[dataSample].GetMaximum()) * num);
    dataHist->GetXaxis()->SetLabelOffset(0.02);
    
    dataHist->Draw("E0");
    distribs.stack.Draw("histsame");
    dataHist->Draw("E0same");
    leg->Draw("same");
    CMS_lumi( pad1, iPeriod, iPos, lumi );

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

    // here data and MC are copied to draw on ratio plot
    TH1D *dataCopy = (TH1D*)dataHist->Clone("dataCopy");
    TH1D *stackCopy = (TH1D*)(distribs.stack.GetStack()->Last())->Clone("stackCopy");
    dataCopy->Divide(stackCopy);
    // first let's calculate uncertainties for data
    TGraphAsymmErrors * dataCopyGraph = new TGraphAsymmErrors(dataCopy);

    for(int i = 0; i < dataCopyGraph->GetN(); i++){

      double dataPoint[3] = {dataHist->GetBinContent(i+1), dataHist->GetBinErrorUp(i+1), dataHist->GetBinErrorLow(i+1)};
      double theMCPoint[3] = {stackCopy->GetBinContent(i+1), stackCopy->GetBinErrorUp(i+1), stackCopy->GetBinErrorLow(i+1)};

      double uncRatio[2];

      // calculating the uncertainty for a / b
      // using formula: delta Unc ^ 2 = ( (a/b)'_a (delta a) ) ^ 2 + ( (a/b)'_b (delta b) ) ^ 2

      uncRatio[0] = TMath::Sqrt(TMath::Power(1 / theMCPoint[0] * dataPoint[1], 2) + TMath::Power(dataPoint[0] / TMath::Power(theMCPoint[0],2) * theMCPoint[1], 2));
      uncRatio[1] = TMath::Sqrt(TMath::Power(1 / theMCPoint[0] * dataPoint[2], 2) + TMath::Power(dataPoint[0] / TMath::Power(theMCPoint[0],2) * theMCPoint[2], 2));

      dataCopyGraph->SetPointError(i, dataCopyGraph->GetErrorXlow(i), dataCopyGraph->GetErrorXhigh(i), uncRatio[1], uncRatio[0]);
    }
    // this one will be used on the ratio plot
    stackCopy->Divide(stackCopy);

    stackCopy->SetFillStyle(1001);
    stackCopy->SetFillColor(kCyan - 4);
    stackCopy->SetMarkerStyle(1);

    stackCopy->SetTitle("");
    stackCopy->GetXaxis()->SetTitle(titleX.c_str());
    if(titleX == "trilep")
      stackCopy->GetXaxis()->SetTitle("");
    stackCopy->GetYaxis()->SetTitle("data/pred");

    stackCopy->GetXaxis()->SetRangeUser(xmin, xmax);

    stackCopy->GetYaxis()->SetTitleOffset(1.2/((1.-xPad)/xPad));
    stackCopy->GetYaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
    stackCopy->GetXaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
    stackCopy->GetYaxis()->SetLabelSize((1.-xPad)/xPad*0.05);
    if(titleX != "")
        stackCopy->GetXaxis()->SetLabelSize((1.-xPad)/xPad*0.05);
    else
        stackCopy->GetXaxis()->SetLabelSize(0.25);

    stackCopy->SetMaximum(2.0);
    stackCopy->SetMinimum(0.0);
    stackCopy->SetMarkerStyle(20);
    stackCopy->SetMarkerSize(0.2);

    dataCopyGraph->SetMarkerSize(0.5);

    stackCopy->Draw("axis");
    // histSystAndStatUnc - a histogram with central value at 1 and with applied one of the uncertainties on top: JEC, JES and Uncl
    TH1D *histSystAndStatUnc;
    histSystAndStatUnc =  (TH1D*)(distribs.stack.GetStack()->Last())->Clone(Form("histSystAndStatUnc"));
    // stack of MC with varied up and down of 3 different types of uncertainties
    const int numberOfSyst = 5;
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
        }
    }

    histSystAndStatUnc->SetFillStyle(1001);
    histSystAndStatUnc->SetFillColor(kOrange - 4);
    histSystAndStatUnc->SetMarkerStyle(1);

    TLegend* mtlegRatio = new TLegend(0.17,0.39,0.85,0.58);
    mtlegRatio->SetNColumns(4);
    mtlegRatio->SetFillColor(0);
    mtlegRatio->SetFillStyle(0);
    mtlegRatio->SetBorderSize(0);
    mtlegRatio->SetTextFont(42);

    mtlegRatio->AddEntry(stackCopy, "Stat", "f");
    mtlegRatio->AddEntry(histSystAndStatUnc, "Stat+Syst", "f");

    //histSystAndStatUnc->Draw("e2same");
    stackCopy->Draw("e2same");
    mtlegRatio->Draw("same");

    TLine *line = new TLine(xmin, 1, xmax, 1);
    line->SetLineStyle(2);

    line->Draw("same");

    dataCopyGraph->Draw("p"); // dataCopyGraph = data / MC stack

    pad2->RedrawAxis();
    pad2->Update();
}

/*
void drawSystUnc(TVirtualPad* c1, DistribsAll & distribs, int process){

    double xmin = distribs.vectorHisto[process].GetXaxis()->GetXmin();
    double xmax = distribs.vectorHisto[process].GetXaxis()->GetXmax();

    double xPad = 0.25; // 0.25

    TPad *pad1 = new TPad("pad1","pad1",0,xPad,1,1);
    pad1->SetTopMargin(0.07);
    if(xPad != 0)
        pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();

    distribs.vectorHisto[process].SetFillColor(0);
    distribs.vectorHisto[process].Draw("hist");

    distribs.vectorHistoUp[process].SetLineColor(kRed);
    distribs.vectorHistoUp[process].SetLineStyle(2);
    distribs.vectorHistoUp[process].Draw("histsame");

    distribs.vectorHistoDown[process].SetLineColor(kRed);
    distribs.vectorHistoDown[process].SetLineStyle(3);
    distribs.vectorHistoDown[process].Draw("histsame");

    //if(xPad == 0) return;

    c1->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,xPad);
    
    pad2->SetBottomMargin((1.-xPad)/xPad*0.13);
    pad2->SetTopMargin(0.06);

    pad2->Draw();
    pad2->RedrawAxis();
    pad2->cd();

    TH1D * var = (TH1D*)distribs.vectorHisto[process].Clone("var");
    TH1D * varUp = (TH1D*)distribs.vectorHistoUp[process].Clone("varUp");
    TH1D * varDown = (TH1D*)distribs.vectorHistoDown[process].Clone("varDown");
    
    varUp->Divide(var);
    varDown->Divide(var);

    TLine *line = new TLine(xmin, 1, xmax, 1);
    line->SetLineStyle(2);

    varUp->GetYaxis()->SetTitleOffset(1.2/((1.-xPad)/xPad));
    varUp->GetYaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
    varUp->GetXaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
    varUp->GetYaxis()->SetLabelSize((1.-xPad)/xPad*0.05);
    varUp->GetXaxis()->SetLabelSize((1.-xPad)/xPad*0.05);
    varUp->GetXaxis()->SetTitle("");
    varUp->GetXaxis()->SetLabelSize((1.-xPad)/xPad*0.055);
    varUp->GetXaxis()->SetLabelOffset(0.03);

    varUp->SetMaximum(1.2);
    varUp->SetMinimum(0.8);

    varUp->Draw("axis");
    
    line->Draw("same");

    varUp->Draw("histsame");
    varDown->Draw("histsame");

    pad2->Update();
    
}   
*/

#endif  // showHist
