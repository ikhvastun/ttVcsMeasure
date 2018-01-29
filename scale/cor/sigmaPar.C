#include "CMS_lumi.cc"
#include "tdrStyle.C"
const int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV
const int iPos =0;

void sigmaPar(){
	gStyle->SetOptStat(0);

	setTDRStyle(); 
	TFile *_file0 = TFile::Open("sigmaParEl.root");
	TFile *_file1 = TFile::Open("sigmaParMu.root");


	
	TH1D* h1 = (TH1D*) _file0->Get("resChoice_0_1");
	TH1D* h2 = (TH1D*) _file1->Get("resChoice_0_1");


	TFile *_file2 = TFile::Open("scaleEl.root");
	TFile *_file3 = TFile::Open("scaleMu.root");

	TH1D* h3 = (TH1D*) _file2->Get("scaleChoice_0");
	TH1D* h4 = (TH1D*) _file3->Get("scaleChoice_0");



	
	h1->SetTitle("");
	h1->GetYaxis()->SetTitle("#sigma(u_{||}) [GeV]");
	h1->GetXaxis()->SetTitle("q_{T} [GeV]");
	TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
	h1->GetXaxis()->SetRangeUser(50,250);

	h1->Divide(h3);
	h2->Divide(h4);

	h1->SetMinimum(0);
	h1->SetMaximum(40);
	h1->Draw("p");
	h2->SetLineColor(kRed);
	h2->SetMarkerColor(kRed);
	h2->Draw("psame");
	TLine * line = new TLine(20, 1, 500, 1);
	line->SetLineStyle(2);
	//line->Draw("same");
	TLegend* leg = new TLegend(0.5, 0.4, 0.9, 0.2);
	leg->AddEntry(h2,"Type 1 PF E_{T}^{miss} Z#rightarrow#mu#mu", "l");
	leg->AddEntry(h1,"Type 1 PF E_{T}^{miss} Z#rightarrow ee", "l");
	leg->Draw("same");

	CMS_lumi( c1, iPeriod, iPos );
}
