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

//void showHist(TVirtualPad* c1, TH1D *hist, TH1D *hist2, THStack *stack, string title, string titleX, string titleY, double num, TLegend *leg){   
void showHist(TVirtualPad* c1, DistribsAll & distribs, string title, string titleX, string titleY, double num, TLegend *leg, bool plotInLog = false){   
 
    double xPad = 0.25; // 0.25

    TPad *pad1 = new TPad("pad1","pad1",0,xPad,1,1);
    pad1->SetTopMargin(0.07);
    if(xPad != 0)
        pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    if(plotInLog)
        pad1->SetLogy();
    
    double xmin = distribs.vectorHisto[dataSample].GetXaxis()->GetXmin();
    double xmax = distribs.vectorHisto[dataSample].GetXaxis()->GetXmax();
    //pad1->DrawFrame(xmin, -0.1, xmax, 1.1);
    
    distribs.vectorHisto[dataSample].SetMarkerSize(1);
    distribs.vectorHisto[dataSample].SetTitle(title.c_str());
    distribs.vectorHisto[dataSample].GetXaxis()->SetTitle(titleX.c_str());
    distribs.vectorHisto[dataSample].GetYaxis()->SetTitle(titleY.c_str());
    distribs.vectorHisto[dataSample].SetMinimum(0.5);
    distribs.vectorHisto[dataSample].SetMaximum(TMath::Max(distribs.stack.GetMaximum(), distribs.vectorHisto[dataSample].GetMaximum()) * num);
    distribs.vectorHisto[dataSample].GetXaxis()->SetLabelOffset(0.01);

    if(titleX == "trilep"){
      distribs.vectorHisto[dataSample].GetXaxis()->SetTitle("");
      distribs.vectorHisto[dataSample].GetXaxis()->SetLabelSize(0.1);
    }
    

    if(title == "log"){
      distribs.vectorHisto[dataSample].SetTitle("");
      distribs.vectorHisto[dataSample].SetMinimum(0.5);
      distribs.vectorHisto[dataSample].SetMaximum(1000.);
    }

    //TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.)");
    //setex2->Draw();

    distribs.vectorHisto[dataSample].Draw("E0");
    distribs.stack.Draw("histsame");

    //TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0.5)");
    //setex1->Draw();

    /*
    distribs.vectorHistoTotalUnc.SetFillStyle(3005);
    distribs.vectorHistoTotalUnc.SetFillColor(kGray+2);
    distribs.vectorHistoTotalUnc.SetMarkerStyle(1);
    distribs.vectorHistoTotalUnc.Draw("e2same");
    */

    //TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.)");
    //setex2->Draw();

    distribs.vectorHisto[dataSample].Draw("E0same");

    leg->Draw("same");
    CMS_lumi( pad1, iPeriod, iPos );

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

    TH1D * dataCopy = (TH1D*)distribs.vectorHisto[dataSample].Clone("dataCopy");
    
    TH1D *stackCopy = (TH1D*)(distribs.stack.GetStack()->Last());
    dataCopy->Divide(stackCopy);

    TGraphAsymmErrors * dataCopyGraph = new TGraphAsymmErrors(dataCopy);

    for(int i = 0; i < dataCopyGraph->GetN(); i++){

      double dataPoint[3] = {distribs.vectorHisto[dataSample].GetBinContent(i+1), distribs.vectorHisto[dataSample].GetBinErrorUp(i+1), distribs.vectorHisto[dataSample].GetBinErrorLow(i+1)};
      double theMCPoint[3] = {stackCopy->GetBinContent(i+1), stackCopy->GetBinErrorUp(i+1), stackCopy->GetBinErrorLow(i+1)};

      double uncRatio[2];

      // calculating the uncertainty for a / b
      // using formula: delta Unc ^ 2 = ( (a/b)'_a (delta a) ) ^ 2 + ( (a/b)'_b (delta b) ) ^ 2

      uncRatio[0] = TMath::Sqrt(TMath::Power(1 / theMCPoint[0] * dataPoint[1], 2) + TMath::Power(dataPoint[0] / TMath::Power(theMCPoint[0],2) * theMCPoint[1], 2));
      uncRatio[1] = TMath::Sqrt(TMath::Power(1 / theMCPoint[0] * dataPoint[2], 2) + TMath::Power(dataPoint[0] / TMath::Power(theMCPoint[0],2) * theMCPoint[2], 2));
      
      dataCopyGraph->SetPointError(i, dataCopyGraph->GetErrorXlow(i), dataCopyGraph->GetErrorXhigh(i), uncRatio[1], uncRatio[0]);
    }

    TH1D * uncHistoCopy = (TH1D*) distribs.vectorHistoTotalUnc.Clone("uncHistoCopy");
    uncHistoCopy->Divide(stackCopy);

    uncHistoCopy->SetFillStyle(3005);
    uncHistoCopy->SetFillColor(kGray+2);
    uncHistoCopy->SetMarkerStyle(1);
   

    uncHistoCopy->SetTitle("");
    uncHistoCopy->GetXaxis()->SetTitle(titleX.c_str());
    if(titleX == "trilep")
      uncHistoCopy->GetXaxis()->SetTitle("");
    uncHistoCopy->GetYaxis()->SetTitle("data/pred");

    uncHistoCopy->GetXaxis()->SetRangeUser(xmin, xmax);

    uncHistoCopy->GetYaxis()->SetTitleOffset(1.2/((1.-xPad)/xPad));
    uncHistoCopy->GetYaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
    uncHistoCopy->GetXaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
    uncHistoCopy->GetYaxis()->SetLabelSize((1.-xPad)/xPad*0.05);
    uncHistoCopy->GetXaxis()->SetLabelSize((1.-xPad)/xPad*0.05);

    uncHistoCopy->SetMaximum(1.5);
    uncHistoCopy->SetMinimum(0.5);
    uncHistoCopy->SetMarkerStyle(20);
    uncHistoCopy->SetMarkerSize(0.2);

    dataCopyGraph->SetMarkerSize(0.5);
     
    uncHistoCopy->Draw("axis");

    if(titleX=="monetplot"){
      uncHistoCopy->GetXaxis()->SetTitle("");
      uncHistoCopy->GetXaxis()->SetTitleOffset(0.6);
      pad2->SetLeftMargin(0.07);
    }

    TLine *line = new TLine(xmin, 1, xmax, 1);
    line->SetLineStyle(2);
    
    line->Draw("same");
    uncHistoCopy->Draw("e2same");
    dataCopyGraph->Draw("p"); // dataCopyGraph = data / MC stack

    pad2->Update();
}

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

#endif  // showHist
