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
using Output::distribs2D;
using Output::DistribsAll;
using Output::DistribsAll2D;

void setUpRatioFeatures(TH1D *, TGraphAsymmErrors *, histInfo & info, double);
void setUpSystUnc(DistribsAll &, TH1D *);
void setUpSystUncCorr(DistribsAll &, TH1D *, DistribsAll &, TH1D *);
void calculateRatioUnc(TGraphAsymmErrors *, TH1D *, TH1D *);
void printInfoOnPlotTTZ();
void printInfoOnPlot3L();
void printInfoOnPlotNPCR();
void printInfoOnXaxisAllTTZ();
void showSeparationHist(TVirtualPad* c1, DistribsAll & distribs, histInfo & info, double num, TLegend *leg, bool plotInLog, bool normalizedToData, const int showLegendOption); // showLegendOption 0 - 2016, 1 - 2017, 2 - 2016+2017
void showHistEff(TVirtualPad* c1, DistribsAll & distribsLoose, DistribsAll & distribsTight);
void showHist(TVirtualPad* c1, DistribsAll & distribs, histInfo & info, double num, TLegend *leg, bool plotInLog = false, bool normalizedToData = false, const int showLegendOption = 0){ // showLegendOption 0 - 2016, 1 - 2017, 2 - 2016+2017
    double xPad = 0.25; // 0.25

    TPad *pad1 = new TPad("pad1","pad1",0,xPad,1,1);
    pad1->SetTopMargin(0.09);
    if(xPad != 0)
        pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    if(plotInLog)
        pad1->SetLogy();
    
    TH1D * dataHist = &distribs.vectorHisto[dataSample];
    // here is the code for blinding, number of events in data is set to 0
    // finally we are unblinded, 23 Oct 2018
    /*
    if(showLegendOption > 0 && (info.index == indexSRTTZ || info.index == indexSR3L || info.index == indexSR4L || info.index == indexSR3L3m || info.index == indexSR3L2m1e || info.index == indexSR3L1m2e || info.index == indexSR3L3e)){ // keep it blinded for 2017, showLegendOption = 1 and 2 for 2017 and comb of 2 datasets
        for(unsigned int bin = 1; bin < dataHist->GetNbinsX()+1; bin++){
            if(bin < 5) continue;
            dataHist->SetBinContent(bin, 0.);
        }
    }
    */
    dataHist->SetMarkerSize(1);
    dataHist->SetTitle("");
    dataHist->GetXaxis()->SetTitle(info.fancyName.c_str());
    dataHist->GetYaxis()->SetTitle(("Number of events " + (info.isEnVar ? ("/ " + std::to_string(int((info.varMax - info.varMin) / info.nBins)) + " GeV") : "")).c_str());
    dataHist->SetMinimum(0.01);
    dataHist->SetMaximum(TMath::Max(distribs.stack.GetMaximum(), distribs.vectorHisto[dataSample].GetMaximum()) * num);
    if(plotInLog){
        dataHist->SetMinimum(0.5);
        dataHist->SetMaximum(TMath::Max(distribs.stack.GetMaximum(), distribs.vectorHisto[dataSample].GetMaximum()) * num * 5);
    }
    dataHist->GetXaxis()->SetLabelOffset(0.02);

    TGraphAsymmErrors* dataGraph = new TGraphAsymmErrors(dataHist);
    for(int b = 1; b < dataHist->GetNbinsX() + 1; ++b){
        dataGraph->SetPointError(b - 1, 0, 0, dataHist->GetBinErrorLow(b), (dataHist->GetBinContent(b) == 0 ) ? 0 : dataHist->GetBinErrorUp(b) );
    }

    if(info.index == 5){ // Njets
        dataHist->GetXaxis()->SetNdivisions(108); // 108,505,
    }
    if(info.index == 6){ // Nbjets
        dataHist->GetXaxis()->SetNdivisions(104); // 108,505,
    }
    dataHist->Draw("axis");
    dataGraph->Draw("pe1Z same");

    leg->Draw("same");
    double lumi = 35.9;
    if(showLegendOption == 1) lumi = 41.5;
    else if (showLegendOption == 2) lumi = 77.5;
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
    //histSystAndStatUnc->Draw("same");

    TLegend* mtlegRatio = new TLegend(0.17,0.39,0.85,0.58);
    mtlegRatio->SetNColumns(4);
    mtlegRatio->SetFillColor(0);
    mtlegRatio->SetFillStyle(0);
    mtlegRatio->SetBorderSize(0);
    mtlegRatio->SetTextFont(42);

    if(info.index == 42){ // SR all TTZ
        mtlegRatio->AddEntry(stackCopy, "Stat.", "f");
        //mtlegRatio->AddEntry(histSystAndStatUnc, "Total", "f");
    }
    /*
    else{
        mtlegRatio->AddEntry(histSystAndStatUnc, "Uncertainty", "f");
    }
    */

    // Draw finally the things
    if(info.index == 5){ // Njets
        stackCopy->GetXaxis()->SetNdivisions(108); // 108,505,
    }
    if(info.index == 6){ // Njets
        stackCopy->GetXaxis()->SetNdivisions(104); // 108,505,
    }

    stackCopy->GetYaxis()->SetNdivisions(303); // 108,505,
    stackCopy->Draw("axis");
    histSystAndStatUnc->SetFillStyle(3005);
    histSystAndStatUnc->SetLineColor(kGray+2);
    histSystAndStatUnc->SetFillColor(kGray+2);
    histSystAndStatUnc->SetMarkerStyle(1);
    histSystAndStatUnc->Draw("e2same");
    if(info.index == 42) // SR all TTZ
        stackCopy->Draw("e2same");

    double xmin = distribs.vectorHisto[dataSample].GetXaxis()->GetXmin();
    double xmax = distribs.vectorHisto[dataSample].GetXaxis()->GetXmax();
    TLine *line = new TLine(xmin, 1, xmax, 1);
    line->SetLineStyle(2);
    line->Draw("same");

    mtlegRatio->Draw("same");

    dataCopyGraph->Draw("pZ"); // dataCopyGraph = data / MC stack

    if(info.index == indexSRTTZ){
        printInfoOnXaxisAllTTZ();
    }

    pad2->cd();
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

    dataGraph->Draw("pe1Z same");

    if(info.index == indexSRTTZ){
       dataHist->GetXaxis()->SetTitleSize(0.07);
       dataHist->GetXaxis()->SetTitleOffset(0.8);
       printInfoOnPlotTTZ();
    }
    if(info.index == indexSRTTCR){
       dataHist->GetXaxis()->SetTitleSize(0.07);
       dataHist->GetXaxis()->SetTitleOffset(0.8);
       printInfoOnPlotNPCR();
    }
    if(info.index == indexSR3L || info.index == indexSR3L3m || info.index == indexSR3L2m1e || info.index == indexSR3L1m2e || info.index == indexSR3L3e){
       dataHist->GetXaxis()->SetTitleSize(0.07);
       dataHist->GetXaxis()->SetTitleOffset(0.8);
       printInfoOnPlot3L();
    }

    pad1->cd();
    pad1->RedrawAxis();
    pad1->Update();

}

void showHistEff(TVirtualPad* c1, DistribsAll & distribsLoose, DistribsAll & distribsTight, const int nonPromptSample){ 
    double xPad = 0.25; // 0.25

    TPad *pad1 = new TPad("pad1","pad1",0,xPad,1,1);
    pad1->SetTopMargin(0.07);
    if(xPad != 0)
        pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    
    TH1D * dataHistLoose = (TH1D*)distribsLoose.vectorHisto[dataSample].Clone("dataHistLoose");
    TH1D * dataHistTight = (TH1D*)distribsTight.vectorHisto[dataSample].Clone("dataHistTight");
    
    TH1D *stackLoose = (TH1D*)(distribsLoose.stack.GetStack()->Last())->Clone("stackLoose");
    TH1D *stackTight = (TH1D*)(distribsTight.stack.GetStack()->Last())->Clone("stackTight");

    TH1D* npLoose = &distribsLoose.vectorHisto[nonPromptSample];
    TH1D* npTight = &distribsTight.vectorHisto[nonPromptSample];

    // subtract the nonprompt from both histos
    dataHistLoose->Add(npLoose, -1);
    dataHistTight->Add(npTight, -1);
    stackLoose->Add(npLoose, -1);
    stackTight->Add(npTight, -1);

    TGraphAsymmErrors* graph_data = new TGraphAsymmErrors( dataHistTight, dataHistLoose);
    TGraphAsymmErrors* graph_MC = new TGraphAsymmErrors( stackTight, stackLoose);

    dataHistTight->SetMarkerSize(1);
    dataHistTight->SetTitle("");
    dataHistTight->GetXaxis()->SetTitle("p_{T} [GeV]");
    dataHistTight->GetYaxis()->SetTitle("Events");
    dataHistTight->SetMinimum(0.);
    dataHistTight->SetMaximum(1.2);
    dataHistTight->GetXaxis()->SetLabelOffset(0.02);
    dataHistTight->SetFillStyle(0); 

    dataHistTight->Divide(dataHistLoose);
    dataHistTight->Draw("hist");

    stackTight->Divide(stackLoose);
    stackTight->SetFillStyle(0); 
    stackTight->Draw("histsame");

    TH1D *histSystAndStatUncLoose = (TH1D*)(distribsLoose.stack.GetStack()->Last())->Clone(Form("histSystAndStatUncLoose"));
    TH1D *histSystAndStatUncTight = (TH1D*)(distribsTight.stack.GetStack()->Last())->Clone(Form("histSystAndStatUncTight"));
    setUpSystUncCorr(distribsTight, histSystAndStatUncTight, distribsLoose, histSystAndStatUncLoose);

    for( int b = 0 ; b < histSystAndStatUncTight->GetNbinsX(); ++b){
        graph_MC->SetPointEYlow(b, histSystAndStatUncTight->GetBinError(b+1));
        graph_MC->SetPointEYhigh(b, histSystAndStatUncTight->GetBinError(b+1));
    }
    graph_data->SetFillStyle(3005);
    graph_data->SetLineWidth(2);
    graph_data->Draw("e2psame");
    graph_MC->SetFillStyle(3004);
    graph_MC->SetLineColor(kRed);
    graph_MC->SetLineWidth(2);
    graph_MC->Draw("e2psame");

    TH1D *dataHistTightUnc = (TH1D*)dataHistTight->Clone("dataHistTightUnc");
    dataHistTightUnc->SetFillStyle(3005);
    dataHistTightUnc->SetFillColor(kBlack);
    dataHistTightUnc->SetMarkerStyle(1);

    TLegend* mtlegRatio = new TLegend(0.17,0.79,0.85,0.98);
    mtlegRatio->SetNColumns(2);
    mtlegRatio->SetFillColor(0);
    mtlegRatio->SetFillStyle(0);
    mtlegRatio->SetBorderSize(0);
    mtlegRatio->SetTextFont(42);

    mtlegRatio->AddEntry(dataHistTight, "data efficiency", "lep");
    mtlegRatio->AddEntry(stackTight, "MC efficiency", "lep");
    mtlegRatio->AddEntry(graph_data, "data stat + nonprompt", "f");
    mtlegRatio->AddEntry(graph_MC, "Lepton SF unc.", "f");
    mtlegRatio->Draw("same");

    double lumi = 35.9;
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

    TH1D * dataHistTightCopy = (TH1D *) dataHistTight->Clone("dataHistTightCopy");

    dataHistTightCopy->SetTitle("");
    dataHistTightCopy->GetYaxis()->SetTitle("data/MC eff");

    dataHistTightCopy->GetYaxis()->SetTitleOffset(1.2/((1.-xPad)/xPad));
    dataHistTightCopy->GetYaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
    dataHistTightCopy->GetXaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
    dataHistTightCopy->GetYaxis()->SetLabelSize((1.-xPad)/xPad*0.05);
    dataHistTightCopy->GetXaxis()->SetLabelSize((1.-xPad)/xPad*0.05);

    dataHistTightCopy->SetMaximum(2.0);
    dataHistTightCopy->SetMinimum(0.0);
    dataHistTightCopy->SetMarkerStyle(20);
    dataHistTightCopy->SetMarkerSize(0.2);

    dataHistTightCopy->Divide(stackTight);

    dataHistTightCopy->Draw("hist");

    TGraphAsymmErrors * dataGraph = new TGraphAsymmErrors(dataHistTightCopy);
    dataGraph->SetFillStyle(3005);
    dataGraph->SetLineWidth(2);

    for(int i = 0; i < graph_data->GetN(); i++){

      double dataPoint[3] = {dataHistTight->GetBinContent(i+1), graph_data->GetErrorYhigh(i), graph_data->GetErrorYlow(i)};
      double theMCPoint[3] = {stackTight->GetBinContent(i+1), graph_MC->GetErrorYhigh(i), graph_MC->GetErrorYlow(i)};

      double uncRatio[2];

      // calculating the uncertainty for a / b
      // using formula: delta Unc ^ 2 = ( (a/b)'_a (delta a) ) ^ 2 + ( (a/b)'_b (delta b) ) ^ 2

      uncRatio[0] = TMath::Sqrt(TMath::Power(1 / theMCPoint[0] * dataPoint[1], 2) + TMath::Power(dataPoint[0] / TMath::Power(theMCPoint[0],2) * theMCPoint[1], 2));
      uncRatio[1] = TMath::Sqrt(TMath::Power(1 / theMCPoint[0] * dataPoint[2], 2) + TMath::Power(dataPoint[0] / TMath::Power(theMCPoint[0],2) * theMCPoint[2], 2));

      dataGraph->SetPointError(i, dataGraph->GetErrorXlow(i), dataGraph->GetErrorXhigh(i), uncRatio[1], uncRatio[0]);
      //dataGraph->SetPointError(i, 0., 0., uncRatio[1], uncRatio[0]);
    }

    dataGraph->Draw("e2psame");

    double xmin = dataHistTightCopy->GetXaxis()->GetXmin();
    double xmax = dataHistTightCopy->GetXaxis()->GetXmax();
    TLine *line = new TLine(xmin, 1, xmax, 1);
    line->SetLineStyle(2);
    line->Draw("same");
    
}

void showSeparationHist(TVirtualPad* c1, DistribsAll & distribs, histInfo & info, double num, TLegend *leg, bool plotInLog = false, bool normalizedToData = false, const int showLegendOption = 0){ // showLegendOption 0 - 2016, 1 - 2017, 2 - 2016+2017
    double xPad = 0.; // 0.25

    TPad *pad1 = new TPad("pad1","pad1",0,xPad,1,1);
    pad1->SetTopMargin(0.07);
    if(xPad != 0)
        pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    if(plotInLog)
        pad1->SetLogy();
    
    // here just use ttW instead of others
    //TH1D * dataHist = &distribs.vectorHisto[ttWSample];
    TH1D * dataHist = (TH1D*)distribs.vectorHisto[ttWSample].Clone("dataHist");

    TH1D *stackCopy = (TH1D*)(distribs.stack.GetStack()->Last())->Clone("stackCopy");
    dataHist->Scale(stackCopy->Integral() / dataHist->Integral());

    dataHist->SetFillStyle(0);
    dataHist->SetMarkerSize(1);
    dataHist->SetTitle("");
    dataHist->GetXaxis()->SetTitle(info.fancyName.c_str());
    dataHist->GetYaxis()->SetTitle(("Events " + (info.isEnVar ? ("/ " + std::to_string(int((info.varMax - info.varMin) / info.nBins)) + " GeV") : "")).c_str());
    dataHist->SetMinimum(0.01);
    dataHist->SetMaximum(TMath::Max(distribs.stack.GetMaximum(), dataHist->GetMaximum()) * num);
    if(plotInLog){
        dataHist->SetMinimum(0.5);
        dataHist->SetMaximum(TMath::Max(distribs.stack.GetMaximum(), dataHist->GetMaximum()) * num * 5);
    }
    dataHist->GetXaxis()->SetLabelOffset(0.02);

    TGraphAsymmErrors* dataGraph = new TGraphAsymmErrors(dataHist);
    for(int b = 1; b < dataHist->GetNbinsX() + 1; ++b){
        dataGraph->SetPointError(b - 1, 0, 0, dataHist->GetBinErrorLow(b), (dataHist->GetBinContent(b) == 0 ) ? 0 : dataHist->GetBinErrorUp(b) );
    }

    pad1->cd();
    pad1->RedrawAxis();
    pad1->Update();

    dataHist->Draw("axis");
    dataHist->Draw("histsame");
    //dataGraph->Draw("pe1 same");

    double lumi = 35.9;
    if(showLegendOption == 1) lumi = 41.5;
    else if (showLegendOption == 2) lumi = 77.5;
    CMS_lumi( pad1, iPeriod, iPos, lumi);

    pad1->cd();
    pad1->RedrawAxis();
    pad1->Update();

    c1->cd();
    pad1->cd();

    TH1D *systAndStatUnc = (TH1D*)(distribs.stack.GetStack()->Last())->Clone("systAndStatUnc");
    distribs.stack.Draw("histsame");
    systAndStatUnc->SetFillStyle(3005);
    systAndStatUnc->SetFillColor(kGray+2);
    systAndStatUnc->SetMarkerStyle(1);
    systAndStatUnc->Draw("e2same");

    //dataGraph->Draw("pe1 same");
    dataHist->Draw("histsame");
    leg->Draw("same");

    pad1->cd();
    pad1->RedrawAxis();
    pad1->Update();

}


void printInfoOnXaxisAllTTZ(){

    TLine *line1 = new TLine(3.5, -1, 3.5, 2);
    line1->SetLineStyle(2);

    TLine *line2 = new TLine(7.5, -1, 7.5, 2);
    line2->SetLineStyle(2);
    
    TLine *line3 = new TLine(11.5, -1, 11.5, 2);
    line3->SetLineStyle(2);
    
    line1->Draw("same");
    line2->Draw("same");
    line3->Draw("same");

    TLatex njetsSign;
    njetsSign.SetNDC();
    njetsSign.SetTextAngle(0);
    njetsSign.SetTextColor(kBlack);

    njetsSign.SetTextFont(42);
    njetsSign.SetTextAlign(31);
    njetsSign.SetTextSize(0.18);
    njetsSign.DrawLatex(0.84, 0.08,"N_{j}");

    TLatex njetsSign1;
    njetsSign1.SetNDC();
    njetsSign1.SetTextAngle(0);
    njetsSign1.SetTextColor(kBlack);

    njetsSign1.SetTextFont(42);
    njetsSign1.SetTextAlign(31);
    njetsSign1.SetTextSize(0.18);
    njetsSign1.DrawLatex(0.58, 0.08,"N_{j}");

    TLatex njetsSign2;
    njetsSign2.SetNDC();
    njetsSign2.SetTextAngle(0);
    njetsSign2.SetTextColor(kBlack);

    njetsSign2.SetTextFont(42);
    njetsSign2.SetTextAlign(31);
    njetsSign2.SetTextSize(0.18);
    njetsSign2.DrawLatex(0.32, 0.08,"N_{j}");

    TLatex nbjetsSign;
    nbjetsSign.SetNDC();
    nbjetsSign.SetTextAngle(0);
    nbjetsSign.SetTextColor(kBlack);

    nbjetsSign.SetTextFont(42);
    nbjetsSign.SetTextAlign(31);
    nbjetsSign.SetTextSize(0.18);
    nbjetsSign.DrawLatex(0.99, 0.08,"N_{b}");

}

void printInfoOnPlotNPCR(){

    TLine *line1 = new TLine(2.5, 0, 2.5, 550);
    line1->SetLineStyle(2);

    TLine *line2 = new TLine(5.5, 0, 5.5, 550);
    line2->SetLineStyle(2);
    
    line1->Draw("same");
    line2->Draw("same");

    TLatex nbjetsEq0region;
    nbjetsEq0region.SetNDC();
    nbjetsEq0region.SetTextAngle(0);
    nbjetsEq0region.SetTextColor(kBlack);

    nbjetsEq0region.SetTextFont(42);
    nbjetsEq0region.SetTextAlign(31);
    nbjetsEq0region.SetTextSize(0.05);
    nbjetsEq0region.DrawLatex(0.31, 0.58,"N_{b} = 0");

    TLatex nbjetsEq1region;
    nbjetsEq1region.SetNDC();
    nbjetsEq1region.SetTextAngle(0);
    nbjetsEq1region.SetTextColor(kBlack);

    nbjetsEq1region.SetTextFont(42);
    nbjetsEq1region.SetTextAlign(31);
    nbjetsEq1region.SetTextSize(0.05);
    nbjetsEq1region.DrawLatex(0.65, 0.58,"N_{b} = 1");

    TLatex nbjetsEq2region;
    nbjetsEq2region.SetNDC();
    nbjetsEq2region.SetTextAngle(0);
    nbjetsEq2region.SetTextColor(kBlack);

    nbjetsEq2region.SetTextFont(42);
    nbjetsEq2region.SetTextAlign(31);
    nbjetsEq2region.SetTextSize(0.05);
    nbjetsEq2region.DrawLatex(0.89, 0.58,"N_{b} > 1");

}

void printInfoOnPlot3L(){

    TLine *line1 = new TLine(3.5, 0, 3.5, 1125);
    line1->SetLineStyle(2);
    line1->Draw("same");

    TLine *line2 = new TLine(7.5, 0, 7.5, 1125);
    line2->SetLineStyle(2);
    
    TLine *line3 = new TLine(11.5, 0, 11.5, 1125);
    line3->SetLineStyle(2);
    
    line2->Draw("same");
    line3->Draw("same");

    TLatex threeLregion;
    threeLregion.SetNDC();
    threeLregion.SetTextAngle(0);
    threeLregion.SetTextColor(kBlack);

    threeLregion.SetTextFont(42);
    threeLregion.SetTextAlign(31);
    threeLregion.SetTextSize(0.05);
    threeLregion.DrawLatex(0.62, 0.68,"3 leptons");

    TLatex nbjetsEq0region;
    nbjetsEq0region.SetNDC();
    nbjetsEq0region.SetTextAngle(0);
    nbjetsEq0region.SetTextColor(kBlack);

    nbjetsEq0region.SetTextFont(42);
    nbjetsEq0region.SetTextAlign(31);
    nbjetsEq0region.SetTextSize(0.05);
    nbjetsEq0region.DrawLatex(0.31, 0.58,"N_{b} = 0");

    TLatex nbjetsEq1region;
    nbjetsEq1region.SetNDC();
    nbjetsEq1region.SetTextAngle(0);
    nbjetsEq1region.SetTextColor(kBlack);

    nbjetsEq1region.SetTextFont(42);
    nbjetsEq1region.SetTextAlign(31);
    nbjetsEq1region.SetTextSize(0.05);
    nbjetsEq1region.DrawLatex(0.55, 0.58,"N_{b} = 1");

    TLatex nbjetsEq2region;
    nbjetsEq2region.SetNDC();
    nbjetsEq2region.SetTextAngle(0);
    nbjetsEq2region.SetTextColor(kBlack);

    nbjetsEq2region.SetTextFont(42);
    nbjetsEq2region.SetTextAlign(31);
    nbjetsEq2region.SetTextSize(0.05);
    nbjetsEq2region.DrawLatex(0.79, 0.58,"N_{b} > 1");
}

void printInfoOnPlotTTZ(){

    TLine *line1 = new TLine(3.5, 0, 3.5, 1125);
    line1->SetLineStyle(2);
    line1->Draw("same");

    TLine *line2 = new TLine(7.5, 0, 7.5, 1125);
    line2->SetLineStyle(2);
    
    TLine *line3 = new TLine(11.5, 0, 11.5, 1125);
    line3->SetLineStyle(2);
    
    line2->Draw("same");
    line3->Draw("same");

    TLatex threeLregion;
    threeLregion.SetNDC();
    threeLregion.SetTextAngle(0);
    threeLregion.SetTextColor(kBlack);

    threeLregion.SetTextFont(42);
    threeLregion.SetTextAlign(31);
    threeLregion.SetTextSize(0.05);
    threeLregion.DrawLatex(0.62, 0.68,"3 leptons");

    TLatex nbjetsEq0region;
    nbjetsEq0region.SetNDC();
    nbjetsEq0region.SetTextAngle(0);
    nbjetsEq0region.SetTextColor(kBlack);

    nbjetsEq0region.SetTextFont(42);
    nbjetsEq0region.SetTextAlign(31);
    nbjetsEq0region.SetTextSize(0.05);
    nbjetsEq0region.DrawLatex(0.31, 0.58,"N_{b} = 0");

    TLatex nbjetsEq1region;
    nbjetsEq1region.SetNDC();
    nbjetsEq1region.SetTextAngle(0);
    nbjetsEq1region.SetTextColor(kBlack);

    nbjetsEq1region.SetTextFont(42);
    nbjetsEq1region.SetTextAlign(31);
    nbjetsEq1region.SetTextSize(0.05);
    nbjetsEq1region.DrawLatex(0.55, 0.58,"N_{b} = 1");

    TLatex nbjetsEq2region;
    nbjetsEq2region.SetNDC();
    nbjetsEq2region.SetTextAngle(0);
    nbjetsEq2region.SetTextColor(kBlack);

    nbjetsEq2region.SetTextFont(42);
    nbjetsEq2region.SetTextAlign(31);
    nbjetsEq2region.SetTextSize(0.05);
    nbjetsEq2region.DrawLatex(0.79, 0.58,"N_{b} > 1");

    TLatex fourLregion;
    fourLregion.SetNDC();
    fourLregion.SetTextAngle(0);
    fourLregion.SetTextColor(kBlack);

    fourLregion.SetTextFont(42);
    fourLregion.SetTextAlign(31);
    fourLregion.SetTextSize(0.05);
    fourLregion.DrawLatex(0.92, 0.58,"4 leptons");

    TLatex fourLregionNjets;
    fourLregionNjets.SetNDC();
    fourLregionNjets.SetTextAngle(0);
    fourLregionNjets.SetTextColor(kBlack);

    fourLregionNjets.SetTextFont(42);
    fourLregionNjets.SetTextAlign(31);
    fourLregionNjets.SetTextSize(0.05);
    fourLregionNjets.DrawLatex(0.92, 0.38,"N_{j} #geq 2");

}

void setUpRatioFeatures(TH1D * stackCopy, TGraphAsymmErrors * dataCopyGraph, histInfo & info, double xPad){

    // this one will be used on the ratio plot
    stackCopy->Divide(stackCopy);

    // if there is 0 event in stack, then set uncertainty to 0 
    //for(int bin = 1; bin < stackCopy->GetNbinsX() + 1; bin++)
    //    stackCopy->SetBinError(bin, 0.);

    stackCopy->SetFillStyle(1001);
    stackCopy->SetFillColor(kCyan - 4);
    stackCopy->SetLineColor(kCyan - 4);
    stackCopy->SetMarkerStyle(1);

    stackCopy->SetTitle("");
    stackCopy->GetXaxis()->SetTitle(info.fancyName.c_str());
    stackCopy->GetYaxis()->SetTitle("Data / Pred.");

    stackCopy->GetYaxis()->SetTitleOffset(1.2/((1.-xPad)/xPad));
    stackCopy->GetYaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
    stackCopy->GetXaxis()->SetTitleSize((1.-xPad)/xPad*0.06);
    stackCopy->GetYaxis()->SetLabelSize((1.-xPad)/xPad*0.05);
    if(info.index != indexSR3L && info.index != indexSR4L && info.index != indexSRTTZ && info.index != indexFlavour3L && info.index != indexFlavour4L && info.index != indexFlavour3L4L && info.index != indexFlavour4LZZ && info.index != indexSRTTCR)
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

      //dataGraph->SetPointError(i, dataGraph->GetErrorXlow(i), dataGraph->GetErrorXhigh(i), uncRatio[1], uncRatio[0]);
      dataGraph->SetPointError(i, 0., 0., uncRatio[1], uncRatio[0]);
    }
}

void setUpSystUnc(DistribsAll & distribs, TH1D * histSystAndStatUnc){

    /*
    histSystAndStatUnc->Reset("ICE");
    for(unsigned int bin = 0; bin < histSystAndStatUnc->GetNbinsX(); bin++){
        double err = 0.;
        cout << "#### bin number " << bin << endl;
        for(unsigned int cat = distribs.vectorHisto.size()-1; cat != 0; cat--){
            // content in particular bin in stack
            double catBinContent = distribs.vectorHisto[cat].GetBinContent(bin+1);
            if(catBinContent == 0) continue;
            double catBinError = distribs.vectorHisto[cat].GetBinError(bin+1);

            // stat uncertainty is here
            if(catBinContent != 0.) err += TMath::Power(catBinError/catBinContent,2);
            cout << "stat unc is " << TMath::Sqrt(err) << endl; 

            for(unsigned int syst = 0; syst < numberOfSyst; syst++) {
                if(syst == pdfUncIndex) continue;
                if(catBinContent != 0.) {
                    err += TMath::Power(TMath::Max(distribs.vectorHistoUncUp[cat].unc[syst].GetBinContent(bin+1) - catBinContent, catBinContent - distribs.vectorHistoUncUp[cat].unc[syst].GetBinContent(bin+1)) / catBinContent, 2);
                    cout << "unc after syst " << syst << " is " << TMath::Sqrt(err) << endl; 
                }
            }
        }
        histSystAndStatUnc->SetBinContent(bin+1, 1.);
        histSystAndStatUnc->SetBinError(bin+1, TMath::Sqrt(err) > 1 ? 1. : TMath::Sqrt(err));
    }
    */

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

    histSystAndStatUnc->SetFillStyle(3005);
    histSystAndStatUnc->SetLineColor(kGray+2);
    histSystAndStatUnc->SetFillColor(kGray+2);

    //histSystAndStatUnc->SetFillStyle(1001);
    //histSystAndStatUnc->SetLineColor(kOrange - 4);
    //histSystAndStatUnc->SetFillColor(kOrange - 4);
    histSystAndStatUnc->SetMarkerStyle(1);
    //histSystAndStatUnc->Draw("same");
}

void setUpSystUncCorr(DistribsAll & distribsTight, TH1D * histSystAndStatUncTight, DistribsAll & distribsLoose, TH1D * histSystAndStatUncLoose){

    // histSystAndStatUnc - a histogram with central value at 1 and with applied one of the uncertainties on top: JEC, JES and Uncl
    //TH1D *histSystAndStatUnc = (TH1D*)(distribs.stack.GetStack()->Last())->Clone(Form("histSystAndStatUnc"));
    TH1D *stackCopyTight = (TH1D*)(distribsTight.stack.GetStack()->Last())->Clone("stackCopyTight");
    TH1D *stackCopyLoose = (TH1D*)(distribsLoose.stack.GetStack()->Last())->Clone("stackCopyLoose");
    // stack of MC with varied up and down of 3 different types of uncertainties
    TH1D *stackUncTightUp[numberOfSyst];
    TH1D *stackUncTightDown[numberOfSyst];
    TH1D *stackUncLooseUp[numberOfSyst];
    TH1D *stackUncLooseDown[numberOfSyst];

    for(unsigned int i = 0; i < numberOfSyst; i++){

        stackUncTightUp[i] = (TH1D*)distribsTight.vectorHisto[1].Clone(Form("stackUncTightUp_%d", i));
        stackUncTightDown[i] = (TH1D*)distribsTight.vectorHisto[1].Clone(Form("stackUncTightDown_%d", i));
        stackUncLooseUp[i] = (TH1D*)distribsLoose.vectorHisto[1].Clone(Form("stackUncLooseUp_%d", i));
        stackUncLooseDown[i] = (TH1D*)distribsLoose.vectorHisto[1].Clone(Form("stackUncLooseDown_%d", i));

        stackUncTightUp[i]->Reset("ICE");
        stackUncTightDown[i]->Reset("ICE");
        stackUncLooseUp[i]->Reset("ICE");
        stackUncLooseDown[i]->Reset("ICE");

        for(unsigned int j = distribsTight.vectorHistoUncUp.size()-1; j != 0; j--){
            stackUncTightUp[i]->Add((TH1D*)&distribsTight.vectorHistoUncUp[j].unc[i]);
            stackUncTightDown[i]->Add((TH1D*)&distribsTight.vectorHistoUncDown[j].unc[i]);
            stackUncLooseUp[i]->Add((TH1D*)&distribsLoose.vectorHistoUncUp[j].unc[i]);
            stackUncLooseDown[i]->Add((TH1D*)&distribsLoose.vectorHistoUncDown[j].unc[i]);
        }
    }

    for(unsigned int i = 0; i < histSystAndStatUncTight->GetNbinsX(); i++){
        // content in particular bin in stack
        double stackBinContentTight = ((TH1D*)distribsTight.stack.GetStack()->Last())->GetBinContent(i+1);
        double stackBinErrorTight = ((TH1D*)distribsTight.stack.GetStack()->Last())->GetBinError(i+1);
        double stackBinContentLoose = ((TH1D*)distribsLoose.stack.GetStack()->Last())->GetBinContent(i+1);
        double stackBinErrorLoose = ((TH1D*)distribsLoose.stack.GetStack()->Last())->GetBinError(i+1);

        double errTight = TMath::Power(stackBinErrorTight, 2);
        double errLoose = TMath::Power(stackBinErrorLoose, 2);
        cout << "stat uncertainties are " << TMath::Sqrt(errTight) << " " << TMath::Sqrt(errLoose) << endl;

        for(unsigned int j = 0; j < numberOfSyst; j++){
            // consider largest deviation between the upward and downward variations
            errTight += TMath::Power(TMath::Max(stackUncTightUp[j]->GetBinContent(i+1) - stackBinContentTight, stackBinContentTight - stackUncTightDown[j]->GetBinContent(i+1)), 2);
            errLoose += TMath::Power(TMath::Max(stackUncLooseUp[j]->GetBinContent(i+1) - stackBinContentLoose, stackBinContentLoose - stackUncLooseDown[j]->GetBinContent(i+1)), 2);
            //cout << "unc after applying " << j << " syst tight " << TMath::Sqrt(errTight) << " and syst loose " << TMath::Sqrt(errLoose)  << endl;
        }
        // if uncertainty is greater than 100% consider 100% uncertainty
        //cout << "so all numbers in bin " << i+1 << " are " << stackBinContentTight << " " << TMath::Sqrt(errTight) << " " << stackBinContentLoose << " " << TMath::Sqrt(errLoose) << endl;
        //histSystAndStatUncTight->SetBinError(i+1, (stackBinContentTight + TMath::Sqrt(errTight))/(stackBinContentLoose + TMath::Sqrt(errLoose)) - stackBinContentTight/stackBinContentLoose);
        //histSystAndStatUncTight->SetBinError(i+1, (stackBinContentTight + TMath::Sqrt(errTight))/(stackBinContentLoose + TMath::Sqrt(errLoose)) < 1 ? (stackBinContentTight + TMath::Sqrt(errTight))/(stackBinContentLoose + TMath::Sqrt(errLoose)) - stackBinContentTight/stackBinContentLoose : 1 - stackBinContentTight/stackBinContentLoose);
        // here I simply subtract (2L + stat) - (1L + stat) uncertainty
        //histSystAndStatUncTight->SetBinError(i+1, TMath::Sqrt(errTight/stackBinContentTight - errLoose/stackBinContentLoose)*stackBinContentTight/stackBinContentLoose);

        // uncertainty will be calculated according to the formula
        // R = X / Y, dR / R = sqrt((dX/X)^2 + (dY/Y)^2)
        // this is in principle not correct for efficiency plot, 1L tight electron + 1 muon loose are the same so uncertainty should be completely divided out, the only uncertainty that has to be shown is tight to loose
        // errTight is the squared sum of uncertainty for tight selection, errLoose is the same for Loose
        //histSystAndStatUncTight->SetBinError(i+1, TMath::Sqrt(errTight/(stackBinContentTight * stackBinContentTight) + errLoose/(stackBinContentLoose * stackBinContentLoose))*stackBinContentTight/stackBinContentLoose);
        // should be a correct one 
        //histSystAndStatUncTight->SetBinError(i+1, TMath::Sqrt(errTight/(stackBinContentTight * stackBinContentTight) - errLoose/(stackBinContentLoose * stackBinContentLoose))*stackBinContentTight/stackBinContentLoose);
        // in the end we agreed to use only uncertainty on the tight
        histSystAndStatUncTight->SetBinError(i+1, TMath::Sqrt(errTight/(stackBinContentTight * stackBinContentTight))*stackBinContentTight/stackBinContentLoose);
    }

    histSystAndStatUncTight->SetFillStyle(1001);
    histSystAndStatUncTight->SetLineColor(kOrange - 4);
    histSystAndStatUncTight->SetFillColor(kOrange - 4);
    histSystAndStatUncTight->SetMarkerStyle(1);
}

void showHist2D(TVirtualPad* c1, TH2D & histo2Dpassed, TH2D & histo2Dall, int flavour = 0, int flComp = 0){

    c1->cd();
    TPad *pad1 = new TPad("pad1","pad1",0,0,1,1);
    pad1->SetRightMargin(0.07);
    TH2D * passed = (TH2D*)histo2Dpassed.Clone("passed");
    TH2D * all    = (TH2D*)histo2Dall.Clone("all");
    passed->Divide(all);
    passed->SetTitle(flavorsString[flavour] + " FR (" + flavorComposString[flComp] + ")");
    passed->GetXaxis()->SetTitle("p_{T}^{corr} [GeV]");
    passed->GetYaxis()->SetTitle("|#eta|");
    passed->GetYaxis()->SetTitleOffset(0.6);
    passed->SetMarkerSize(1.5);
    passed->Draw("etextcolz");
    passed->SaveAs("plotsForSave/" + flavorsString[flavour] + "FR_" + flavorComposString[flComp] + ".root");
    //distribs2D.vectorHisto[0].Draw("colz");

}

void showHist2D(TVirtualPad* c1, DistribsAll2D & distribs2D, int flavour = 0){

    c1->cd();
    TPad *pad1 = new TPad("pad1","pad1",0,0,1,1);
    pad1->SetRightMargin(0.07);
    TH2D * passed = (TH2D*)distribs2D.vectorHisto[0].Clone("passed");
    TH2D * all    = (TH2D*)distribs2D.vectorHisto[1].Clone("all");
    passed->Divide(all);
    passed->Draw("etextcolz");
    passed->SaveAs("plotsForSave/" + flavorsString[flavour] + "FR.root");
    //distribs2D.vectorHisto[0].Draw("colz");
    return;

}

#endif  // showHist
