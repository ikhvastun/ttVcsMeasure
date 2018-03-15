#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iomanip>

#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
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
#include "TGraphAsymmErrors.h"

//#include "mt2_bisect.h"

using namespace std;


void readTree()
{
    
    const int nLeptonsMax = 10;
    int fontToUse = 42;
    gStyle->SetOptFit(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetNdivisions(505,"XY");
    
    gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(kTRUE);
    gStyle->SetTickLength(0.03, "XYZ");
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);
    
    gStyle->SetLabelFont(fontToUse,"XYZ");
    gStyle->SetLabelSize(0.05,"XYZ");
    gStyle->SetLabelOffset(0.001,"X");
    
    gStyle->SetTitleFont(fontToUse,"");
    gStyle->SetTitleFont(fontToUse,"XYZ");
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetTitleSize(0.06,"XY");
    //gStyle->SetTitleXOffset(1.0);
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYOffset(1.2);
    
    gStyle->SetErrorX(0.);
    
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadBottomMargin(0.12);
    //gStyle->SetPadRightMargin(0.035);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadLeftMargin(0.15);
    
    gStyle->SetStatFont(fontToUse);
    gStyle->SetStatColor(10);
    gStyle->SetStatFontSize(0.08);
    gStyle->SetTitleFont(fontToUse);
    gStyle->SetTitleFontSize(0.08);
    
    gStyle->SetMarkerSize(1.);
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerColor(gStyle->GetHistLineColor());
    
    gStyle->SetPalette(1,0);
    
    gStyle->SetFuncColor(kRed);

    const double MCTruth[] = {
	    	1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,0.000140973 ,0.000240998 ,0.00071209 ,0.00130121 ,0.00245255 ,0.00502589 ,0.00919534 ,0.0146697 ,0.0204126 ,0.0267586 ,
		0.0337697 ,0.0401478 ,0.0450159 ,0.0490577 ,0.0524855 ,0.0548159 ,0.0559937 ,0.0554468 ,0.0537687 ,0.0512055 ,0.0476713 ,0.0435312 ,0.0393107 ,0.0349812 ,0.0307413 ,0.0272425 ,
		0.0237115 ,0.0208329 ,0.0182459 ,0.0160712 ,0.0142498 ,0.012804 ,0.011571 ,0.010547 ,0.00959489 ,0.00891718 ,0.00829292 ,0.0076195 ,0.0069806 ,0.0062025 ,0.00546581 ,0.00484127 ,
		0.00407168 ,0.00337681 ,0.00269893 ,0.00212473 ,0.00160208 ,0.00117884 ,0.000859662 ,0.000569085 ,0.000365431 ,0.000243565 ,0.00015688 ,9.88128e-05 ,6.53783e-05 ,3.73924e-05 ,
		2.61382e-05 ,2.0307e-05 ,1.73032e-05 ,1.435e-05 ,1.36486e-05 ,1.35555e-05 ,1.37491e-05 ,1.34255e-05 ,1.33987e-05 ,1.34061e-05 ,1.34211e-05 ,1.34177e-05 ,1.32959e-05 ,1.33287e-05
};
    

    
    TFile *file0 = TFile::Open("runAll.root");
    TH1F *hdata = (TH1F*) file0->Get("pileup");
    hdata->Sumw2();

    hdata->Scale(1./hdata->Integral());
    
    //TFile *file1 = TFile::Open("test.root");
    TH1F *htrue = new TH1F("htrue", "htrue", 80, 0, 80);
    htrue->Sumw2();

    /*
    TH1F *hdata = new TH1F("hdata", "hdata", 100, 0, 100);
    hdata->Sumw2();
    */

    TFile *hfile = new TFile("/Users/illiakhvastunov/Desktop/CERN/MCsamples/94X/ZllMET/2018MoriondMC_JECv6/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.root","read");
    //TFile *hfileData = new TFile("/Users/illiakhvastunov/Desktop/CERN/MCsamples/94X/ZllMET/DoubleEG_ReReco.root","read");

    //UChar_t _n_PV;
    Float_t _n_PV;
    Double_t        _weight;

    hfile->cd("blackJackAndHookers");
    TTree * inputTreeFake = static_cast<TTree*>( hfile->Get("blackJackAndHookers/blackJackAndHookersTree"));
	            
    inputTreeFake->SetBranchAddress("_nTrueInt", &_n_PV);
    inputTreeFake->SetBranchAddress("_weight", &_weight);

    Long64_t nEntries = inputTreeFake->GetEntries();
    TH1D* hCounter = new TH1D("hCounter", "Events counter", 1, 0, 1);
    hCounter->Read("hCounter"); 

    for (Long64_t it=0; it!=nEntries; ++it) {

	inputTreeFake->GetEntry(it);

	if(_n_PV < 10) continue;
	//htrue->Fill(_n_PV, _weight * 6025 * 41.9 * 1000 / hCounter->GetBinContent(1));
	htrue->Fill(_n_PV);
    }

/*
    hfileData->cd("blackJackAndHookers");
    TTree * inputTreeFakeData = static_cast<TTree*>( hfileData->Get("blackJackAndHookers/blackJackAndHookersTree"));
    inputTreeFakeData->SetBranchAddress("_nVertex", &_n_PV);
	            
    nEntries = inputTreeFakeData->GetEntries();

    for (Long64_t it=0; it!=nEntries; ++it) {

	inputTreeFakeData->GetEntry(it);

	if(_n_PV < 10) continue;
	hdata->Fill(_n_PV);
    }
    */

    //for(int i = 0; i < 60; i++)
    //	htrue->SetBinContent(i+1, MCTruth[i]);

    htrue->Scale(1./htrue->Integral());

    //hdata->Scale(1./hdata->Integral());
    TH1F *h3 = (TH1F*)hdata->Clone("puw");
    
    h3->Divide(htrue);

    TCanvas *c2 = new TCanvas("c2", "c2");
    TFile *file = TFile::Open("puWeights.root","RECREATE");
    h3->Draw();
    h3->Write();
    file->Close();
    

    std::cout<<"Done"<<std::endl;
    
   

}


int main(int argc, char *argv[]){

    TApplication *rootapp = new TApplication("example", &argc, argv);

    readTree();

    rootapp->Run();

    return 0;
}

