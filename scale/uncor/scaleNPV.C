#include "CMS_lumi.cc"
#include "tdrStyle.C"
const int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV
const int iPos =0;

void scaleNPV(){
	gStyle->SetOptStat(0);
	setTDRStyle(); 
	/*
	TFile *_file0 = TFile::Open("17p8/npv_all/scaleMu.root");
	TFile *_file1 = TFile::Open("17p8/npv_0_18/scaleMu.root");
	TFile *_file2 = TFile::Open("17p8/npv_18_32/scaleMu.root");
	TFile *_file3 = TFile::Open("17p8/npv_32_inf/scaleMu.root");
	
	*/
	TFile *_file0 = TFile::Open("17p8/npv_all/scaleEl.root");
	TFile *_file1 = TFile::Open("17p8/npv_0_18/scaleEl.root");
	TFile *_file2 = TFile::Open("17p8/npv_18_32/scaleEl.root");
	TFile *_file3 = TFile::Open("17p8/npv_32_inf/scaleEl.root");
	
	h1 = (TH1D*) _file0->Get("scaleChoice_1");
	h2 = (TH1D*) _file1->Get("scaleChoice_1");
	h3 = (TH1D*) _file2->Get("scaleChoice_1");
	h4 = (TH1D*) _file3->Get("scaleChoice_1");
	
	h1->SetTitle("");
	h1->GetYaxis()->SetTitle("-<u_{||}>/<q_{T}>");
	h1->GetXaxis()->SetTitle("q_{T} [GeV]");
	c1 = new TCanvas("c1", "c1", 600, 600);
	h1->SetMinimum(0.85);
	h1->SetMaximum(1.05);
	h1->Draw("p");
	h2->SetLineColor(kRed);
	h2->SetMarkerColor(kRed);
	h2->Draw("psame");
	h3->SetLineColor(kGreen);
	h3->SetMarkerColor(kGreen);
	h3->Draw("psame");
	h4->SetLineColor(kYellow);
	h4->SetMarkerColor(kYellow);
	h4->Draw("psame");
	line = new TLine(20, 1, 500, 1);
	line->SetLineStyle(2);
	line->Draw("same");
	leg = new TLegend(0.5, 0.4, 0.9, 0.2);
/*
	leg->AddEntry(h1,"Raw PF E_{T}^{miss} Z#rightarrow#mu#mu", "l");
	leg->AddEntry(h2,"Raw PF E_{T}^{miss} Z#rightarrow#mu#mu NPV [0,18]", "l");
	leg->AddEntry(h3,"Raw PF E_{T}^{miss} Z#rightarrow#mu#mu NPV [18,32]", "l");
	leg->AddEntry(h4,"Raw PF E_{T}^{miss} Z#rightarrow#mu#mu NPV [32,inf]", "l");
*/
	leg->AddEntry(h1,"Raw PF E_{T}^{miss} Z#rightarrow ee", "l");
	leg->AddEntry(h2,"Raw PF E_{T}^{miss} Z#rightarrow ee NPV [0,18]", "l");
	leg->AddEntry(h3,"Raw PF E_{T}^{miss} Z#rightarrow ee NPV [18,32]", "l");
	leg->AddEntry(h4,"Raw PF E_{T}^{miss} Z#rightarrow ee NPV [32,inf]", "l");

	leg->Draw("same");

	CMS_lumi( c1, iPeriod, iPos );
}

int main(int argc, char *argv[]){

    scaleNPV();
    return 0;
}
