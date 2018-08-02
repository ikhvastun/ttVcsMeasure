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

void showHist(TVirtualPad* c1, DistribsAll & distribs, string title, string titleX, string titleY, double num, TLegend *leg, bool plotLog = false, const int option = 3){ // 3 stands here for fulldataset, 0 - RunB, 1 - RunCDE, 2 - RunF
 
    double xPad = 0.25; // use 0 if ratio plot is not needed, 0.25 corresponds to 1/4 of the canvas will be filled with the ratio plot 

    TPad *pad1 = new TPad("pad1","pad1",0,xPad,1,1);
    pad1->SetTopMargin(0.07);
    if(xPad != 0)
        pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    if(plotLog)
        pad1->SetLogy();

    TH1D * dataHist;
    // stack here corresponds to all MC 
    THStack *stack = new THStack("hs","hs");
    dataHist = &distribs.vectorHistoEras[dataSample].runEras[option];
    TH1D * histoStackCopy = (TH1D*)distribs.vectorHistoEras[distribs.vectorHisto.size()-1].runEras[option].Clone("histoStackCopy");
    histoStackCopy->Reset("ICE");
    for(unsigned int j = distribs.vectorHisto.size()-1; j != 0; j--){
        histoStackCopy->Add(&distribs.vectorHistoEras[j].runEras[option]);  
    }
    double scaleForMC = distribs.vectorHistoEras[dataSample].runEras[option].Integral() / histoStackCopy->Integral();
    for(unsigned int j = distribs.vectorHisto.size()-1; j != 0; j--){      
        distribs.vectorHistoEras[j].runEras[option].Scale(scaleForMC);
        stack->Add(&distribs.vectorHistoEras[j].runEras[option]);
    }

    double xmin = dataHist->GetXaxis()->GetXmin();
    double xmax = dataHist->GetXaxis()->GetXmax();
    
    dataHist->SetMarkerSize(1);
    dataHist->SetTitle(title.c_str());
    dataHist->GetXaxis()->SetTitle(titleX.c_str());
    dataHist->GetYaxis()->SetTitle(titleY.c_str());
    dataHist->SetMinimum(0.5);
    dataHist->SetMaximum(TMath::Max(stack->GetMaximum(), dataHist->GetMaximum()) * num);
    dataHist->GetXaxis()->SetLabelOffset(0.01);

    //TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.)");
    //setex2->Draw();

    dataHist->Draw("E0");
    stack->Draw("histsame");

    //TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0.5)");
    //setex1->Draw();

    //distribs.vectorHistoTotalUnc.SetFillStyle(3005);
    //distribs.vectorHistoTotalUnc.SetFillColor(kGray+2);
    //distribs.vectorHistoTotalUnc.SetMarkerStyle(1);
    //distribs.vectorHistoTotalUnc.Draw("e2same");
    
    //TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.)");
    //setex2->Draw();

    dataHist->Draw("E0same");
    leg->Draw("same");

    // 0 - RunB, 1 - RunCDE, 2 - RunF, 3 - fullDataset
    TString lumi_13TeV = option == 0 ? "4.8 fb^{-1}" : (option == 1 ? "23.1 fb^{-1}" : (option == 2 ? "13.5 fb^{-1}" : "41.9 fb^{-1}"));
    CMS_lumi( pad1, iPeriod, iPos, lumi_13TeV );

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
    // copied in the beginning of program to normalize stack
    TH1D *stackCopy = (TH1D*)(stack->GetStack()->Last())->Clone("stackCopy");
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
    stackCopy->GetXaxis()->SetLabelSize((1.-xPad)/xPad*0.05);

    stackCopy->SetMaximum(2.0);
    stackCopy->SetMinimum(0.0);
    stackCopy->SetMarkerStyle(20);
    stackCopy->SetMarkerSize(0.2);

    dataCopyGraph->SetMarkerSize(0.5);

    stackCopy->Draw("axis");

    // histCopyUnc - a histogram with central value at 1 and with applied one of the uncertainties on top: JEC, JES and Uncl
    TH1D *histCopyUnc[3];
    // stack of MC with varied up and down of 3 different types of uncertainties
    TH1D *stackUncUp[3];
    TH1D *stackUncDown[3];

    for(unsigned int i = 0; i < 3; i++){
        
        histCopyUnc[i] =  (TH1D*)(stack->GetStack()->Last())->Clone(Form("histCopyUnc_%d", i));
        stackUncUp[i] = (TH1D*)(stack->GetStack()->Last())->Clone(Form("stackUncUp_%d", i));
        stackUncDown[i] = (TH1D*)(stack->GetStack()->Last())->Clone(Form("stackUncDown_%d", i));

        histCopyUnc[i]->Reset();
        stackUncUp[i]->Reset();
        stackUncDown[i]->Reset();

        for(unsigned int j = distribs.vectorHistoUp.size()-1; j != 0; j--){
            stackUncUp[i]->Add((TH1D*)&distribs.vectorHistoEras[j].runErasUncUp[option].unc[i]);
            stackUncDown[i]->Add((TH1D*)&distribs.vectorHistoEras[j].runErasUncDown[option].unc[i]);
        }
        stackUncUp[i]->Scale(scaleForMC);
        stackUncDown[i]->Scale(scaleForMC);
    }  

    for(unsigned int i = 0; i < histCopyUnc[0]->GetNbinsX(); i++){
        // content in particular bin in stack
        double stackBinContent = ((TH1D*)stack->GetStack()->Last())->GetBinContent(i+1);
        double stackBinError = ((TH1D*)stack->GetStack()->Last())->GetBinError(i+1);

        stackCopy->SetBinError(i+1, stackBinError / stackBinContent);
        double err = TMath::Power(stackCopy->GetBinError(i+1), 2);

        for(unsigned int j = 0; j < 3; j++){
            // consider largest deviation between the upward and downward variations 
            err += TMath::Power(TMath::Max(stackUncUp[j]->GetBinContent(i+1) - stackBinContent, stackBinContent - stackUncDown[j]->GetBinContent(i+1)) / stackBinContent, 2);
            histCopyUnc[j]->SetBinContent(i+1, 1.);
            // if uncertainty is greater than 100% consider 100% uncertainty
            histCopyUnc[j]->SetBinError(i+1, TMath::Sqrt(err) > 1 ? 1. : TMath::Sqrt(err));
        }
    }

    histCopyUnc[2]->SetFillStyle(1001);
    histCopyUnc[2]->SetFillColor(kGreen - 4);
    histCopyUnc[2]->SetMarkerStyle(1);

    histCopyUnc[1]->SetFillStyle(1001);
    histCopyUnc[1]->SetFillColor(kBlue - 4);
    histCopyUnc[1]->SetMarkerStyle(1);

    histCopyUnc[0]->SetFillStyle(1001);
    histCopyUnc[0]->SetFillColor(kOrange - 4);
    histCopyUnc[0]->SetMarkerStyle(1);

    TLegend* mtlegRatio = new TLegend(0.17,0.39,0.85,0.58); 
    mtlegRatio->SetNColumns(4);
    mtlegRatio->SetFillColor(0);
    mtlegRatio->SetFillStyle(0);
    mtlegRatio->SetBorderSize(0);
    mtlegRatio->SetTextFont(42);

    mtlegRatio->AddEntry(stackCopy, "Stat", "f");
    mtlegRatio->AddEntry(histCopyUnc[0], "JES + Stat", "f");
    mtlegRatio->AddEntry(histCopyUnc[1], "JER + JES + Stat", "f");
    mtlegRatio->AddEntry(histCopyUnc[2], "Uncl + JER + JES + Stat", "f");

    if (titleX.find("E_{T}^{miss}") != std::string::npos || titleX.find("u_") != std::string::npos){
        histCopyUnc[2]->Draw("e2same");
        histCopyUnc[1]->Draw("e2same");
        histCopyUnc[0]->Draw("e2same");
    }
    stackCopy->Draw("e2same");
    mtlegRatio->Draw("same");

    TLine *line = new TLine(xmin, 1, xmax, 1);
    line->SetLineStyle(2);
    
    line->Draw("same");
    
    dataCopyGraph->Draw("p"); // dataCopyGraph = data / MC stack

    pad2->RedrawAxis();
    pad2->Update();
}

void showDataComp(TVirtualPad* c1, DistribsAll & distribs, string title, string titleX, string titleY, double num, TLegend *leg, bool plotLog = false, int runEra = 0){
    
    double xPad = 0.25; // 0.25

    TPad *pad1 = new TPad("pad1","pad1",0,xPad,1,1);
    pad1->SetTopMargin(0.07);
    if(xPad != 0)
        pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    if(plotLog)
        pad1->SetLogy();
    
    double xmin = distribs.histDataEras[runEra][0].GetXaxis()->GetXmin();
    double xmax = distribs.histDataEras[runEra][0].GetXaxis()->GetXmax();
    //pad1->DrawFrame(xmin, -0.1, xmax, 1.1);
    
    distribs.histDataEras[runEra][0].SetMarkerSize(1);
    distribs.histDataEras[runEra][0].SetTitle(title.c_str());
    distribs.histDataEras[runEra][0].GetXaxis()->SetTitle(titleX.c_str());
    distribs.histDataEras[runEra][0].GetYaxis()->SetTitle(titleY.c_str());
    distribs.histDataEras[runEra][0].SetMinimum(0.5);
    //distribs.histDataEras[runB].SetMaximum(TMath::Max(distribs.stack.GetMaximum(), distribs.vectorHisto[dataSample].GetMaximum()) * num);
    distribs.histDataEras[runEra][0].GetXaxis()->SetLabelOffset(0.01);

    //TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.)");
    //setex2->Draw();

    distribs.histDataEras[runEra][0].Draw("E0");
    //distribs.histDataEras[runEra][1].Scale(distribs.histDataEras[runEra][0].Integral()  / distribs.histDataEras[runEra][1].Integral());
    distribs.histDataEras[runEra][1].SetMarkerColor(kRed);
    distribs.histDataEras[runEra][1].SetLineColor(kRed);
    distribs.histDataEras[runEra][1].Draw("E0same");
    /*
    distribs.histDataEras[runCDE].Scale(distribs.histDataEras[runB].Integral()  / distribs.histDataEras[runCDE].Integral());
    distribs.histDataEras[runCDE].SetMarkerColor(kRed);
    distribs.histDataEras[runCDE].SetLineColor(kRed);
    distribs.histDataEras[runCDE].Draw("E0same");
    distribs.histDataEras[runF].Scale(distribs.histDataEras[runB].Integral() / distribs.histDataEras[runF].Integral());
    distribs.histDataEras[runF].SetMarkerColor(kBlue);
    distribs.histDataEras[runF].SetLineColor(kBlue);
    distribs.histDataEras[runF].Draw("E0same");
    */
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

    /*
    TH1D * dataB = (TH1D*)distribs.histDataEras[runB].Clone("dataB");
    TH1D * dataCDE = (TH1D*)distribs.histDataEras[runCDE].Clone("dataCDE");
    TH1D * dataF = (TH1D*)distribs.histDataEras[runF].Clone("dataF");
    */
    TH1D * dataReReco = (TH1D*)distribs.histDataEras[runEra][0].Clone("dataReReco");
    TH1D * dataPromptReco = (TH1D*)distribs.histDataEras[runEra][1].Clone("dataPromptReco");
   

    dataReReco->SetTitle("");
    dataReReco->GetXaxis()->SetTitle(titleX.c_str());
    dataReReco->GetYaxis()->SetTitle("Prompt / ReReco");

    dataReReco->GetXaxis()->SetRangeUser(xmin, xmax);

    dataReReco->GetYaxis()->SetTitleOffset(1.2/((1.-xPad)/xPad));
    dataReReco->GetYaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
    dataReReco->GetXaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
    dataReReco->GetYaxis()->SetLabelSize((1.-xPad)/xPad*0.05);
    dataReReco->GetXaxis()->SetLabelSize((1.-xPad)/xPad*0.05);

    dataReReco->SetMaximum(1.5);
    dataReReco->SetMinimum(0.5);
    dataReReco->SetMarkerStyle(20);

    dataPromptReco->SetMarkerSize(0.5);
     
    dataReReco->Draw("axis");

    TLine *line = new TLine(xmin, 1, xmax, 1);
    line->SetLineStyle(2);
    
    line->Draw("same");

    dataPromptReco->Divide(dataReReco);
    dataPromptReco->Draw("psame");

    line->Draw("same");


    pad2->Update();
    
}

#endif  // showHist
