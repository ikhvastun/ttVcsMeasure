 
    {	
 
    	setTDRStyle();
    	gStyle->SetOptStat(0);

	    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 300);
	    TPad *pad1 = new TPad("pad1","pad1",0,0,1,1);
	    pad1->SetTopMargin(0.07);
	    pad1->SetBottomMargin(0.12);
	    pad1->SetRightMargin(0.02);
	    pad1->Draw();
	    pad1->cd();

	    TFile *_file0 = TFile::Open("shapeFile_SRallTTZ2016_outWZunc.root");

	    TString processToDraw = "ttZ";

	    hCentral = (TH1D*) _file0->Get(processToDraw);
	    
		std::vector<TString> systName = {"lepSFsyst", "lepSFstat2016", "lepSFReco"};
	    std::vector<Color_t> colorsUp = {kBlue, kBlue-8, kCyan};
	    std::vector<Color_t> colorsDown = {kRed, kRed+3, kRed+4};
	    std::vector<int> lineStyle = {1, 1, 1};
	    std::vector<TString> fancyName = {"lepSF syst", "lepSF stat", "lepSF reco"};

	    TH1D * hTotalUp; TH1D * hTotalDown;
	    std::vector<TH1D *> hUp(systName.size()); std::vector<TH1D *> hDown(systName.size());
			    
	    for(int syst = 0; syst < systName.size(); syst++){

	    	hUp[syst] = (TH1D*) _file0->Get(processToDraw + "_" + systName.at(syst) + "Up");
	   		hDown[syst] = (TH1D*) _file0->Get(processToDraw + "_" + systName.at(syst) + "Down");
	   		hUp[syst]->Divide(hCentral);
	   		hDown[syst]->Divide(hCentral);
	   		
	    	for(int bin = 1; bin < hUp[syst]->GetNbinsX()+1; bin++){
	    		cout << "up/down for " << fancyName.at(syst) << " syst are " << hUp[syst]->GetBinContent(bin) << " " << hDown[syst]->GetBinContent(bin) << endl;
	    	}
	    }

	    
		double xmin = hUp[0]->GetXaxis()->GetXmin();
    	double xmax = hUp[0]->GetXaxis()->GetXmax();
    	TLine *line = new TLine(xmin, 1, xmax, 1);


    	hUp[0]->GetXaxis()->SetLabelSize(0.09);

		hUp[0]->SetMinimum(0.9);
		hUp[0]->SetMaximum(1.1);

		hUp[0]->SetTitle("");
		hUp[0]->GetYaxis()->SetTitleOffset(0.6);
		hUp[0]->GetYaxis()->SetTitle("Relative unc.");

		TLegend* mtleg = new TLegend(0.1,0.92,0.95,0.72); 
  		mtleg->SetNColumns(7);
  		mtleg->SetFillColor(0);
  		mtleg->SetFillStyle(0);
  		mtleg->SetBorderSize(0);
  		mtleg->SetTextFont(42);

		for(int syst = 0; syst < systName.size(); syst++){
			hUp[syst]->SetLineColor(colorsUp.at(syst));
	    	hUp[syst]->SetMarkerColor(colorsUp.at(syst));
	    	hUp[syst]->SetLineWidth(3);
	    	hUp[syst]->SetLineStyle(lineStyle.at(syst));
	    	hUp[syst]->Draw("histsame");
	    	hDown[syst]->SetLineColor(colorsDown.at(syst));
	    	hDown[syst]->SetMarkerColor(colorsDown.at(syst));
	    	hDown[syst]->SetLineWidth(3);
	    	hDown[syst]->SetLineStyle(lineStyle.at(syst));
	    	hDown[syst]->Draw("histsame");

	    	mtleg->AddEntry(hUp[syst], fancyName.at(syst) + " Up", "l");
	    	mtleg->AddEntry(hDown[syst], fancyName.at(syst) + " Down", "l");
		}
		
		//mtleg->Draw("same");

	    line->SetLineStyle(2);
    	line->Draw("same");

    	TLine *lineHor1 = new TLine(3.5, 0.9, 3.5, 1.05);
    	lineHor1->SetLineStyle(2);
    	lineHor1->Draw("same");

    	TLatex regionWZ;
    	regionWZ.SetNDC();
    	regionWZ.SetTextAngle(0);
    	regionWZ.SetTextColor(kBlack);

    	regionWZ.SetTextFont(42);
    	regionWZ.SetTextAlign(31);
    	regionWZ.SetTextSize(0.09);
    	regionWZ.DrawLatex(0.29, 0.22,"3L N_{b jets} = 0");


    	TLine *lineHor2 = new TLine(7.5, 0.9, 7.5, 1.05);
    	lineHor2->SetLineStyle(2);
    	lineHor2->Draw("same");

    	TLatex regionZZ;
    	regionZZ.SetNDC();
    	regionZZ.SetTextAngle(0);
    	regionZZ.SetTextColor(kBlack);

    	regionZZ.SetTextFont(42);
    	regionZZ.SetTextAlign(31);
    	regionZZ.SetTextSize(0.09);
    	regionZZ.DrawLatex(0.52, 0.22,"3L N_{b jets} = 1");


    	TLine *lineHor3 = new TLine(11.5, 0.9, 11.5, 1.05);
    	lineHor3->SetLineStyle(2);
    	lineHor3->Draw("same");

    	TLatex regionNP;
    	regionNP.SetNDC();
    	regionNP.SetTextAngle(0);
    	regionNP.SetTextColor(kBlack);

    	regionNP.SetTextFont(42);
    	regionNP.SetTextAlign(31);
    	regionNP.SetTextSize(0.09);
    	regionNP.DrawLatex(0.79, 0.22,"3L N_{b jets} > 1");

    	TLine *lineHor4 = new TLine(24.5, 0.9, 24.5, 1.05);
    	lineHor4->SetLineStyle(2);
    	lineHor4->Draw("same");

    	TLatex regionTTZ3L;
    	regionTTZ3L.SetNDC();
    	regionTTZ3L.SetTextAngle(0);
    	regionTTZ3L.SetTextColor(kBlack);

    	regionTTZ3L.SetTextFont(42);
    	regionTTZ3L.SetTextAlign(31);
    	regionTTZ3L.SetTextSize(0.09);
    	regionTTZ3L.DrawLatex(0.93, 0.22,"4L");

    	TLatex regionTTZ4L;
    	regionTTZ4L.SetNDC();
    	regionTTZ4L.SetTextAngle(0);
    	regionTTZ4L.SetTextColor(kBlack);

    	regionTTZ4L.SetTextFont(42);
    	regionTTZ4L.SetTextAlign(31);
    	regionTTZ4L.SetTextSize(0.08);
    	regionTTZ4L.DrawLatex(0.85, 0.05,"N_{j}");

    	TLatex nbjetsSign;
    	nbjetsSign.SetNDC();
    	nbjetsSign.SetTextAngle(0);
    	nbjetsSign.SetTextColor(kBlack);

    	nbjetsSign.SetTextFont(42);
    	nbjetsSign.SetTextAlign(31);
    	nbjetsSign.SetTextSize(0.08);
    	nbjetsSign.DrawLatex(0.98, 0.05,"N_{b}");

    	TH1D * hStat;
    	hStat = (TH1D*)hCentral->Clone("hStat");
    	hStat->Reset("ICE");
		for(int bin = 1; bin < hCentral->GetNbinsX()+1; bin++){
    		TH1D * tempUp = (TH1D*) _file0->Get(processToDraw + "_" + processToDraw + "_statSRallTTZ_bin_" + std::to_string(bin) + "Up");
    		TH1D * tempDown = (TH1D*) _file0->Get(processToDraw + "_" + processToDraw + "_statSRallTTZ_bin_" + std::to_string(bin) + "Down");
			
			tempUp->Divide(hCentral);
    		tempDown->Divide(hCentral);
    		hStat->SetBinContent(bin, 1.);
    		hStat->SetBinError(bin, TMath::Max(tempUp->GetBinError(bin), tempDown->GetBinError(bin)));
    	}
    	hStat->SetLineColor(kBlack);
    	hStat->SetMarkerColor(kBlack);
    	hStat->SetMarkerSize(1);
    	hStat->Draw("psame");

    	mtleg->AddEntry(hStat, "Stat. Unc", "lep");
    	mtleg->Draw("same");

}
