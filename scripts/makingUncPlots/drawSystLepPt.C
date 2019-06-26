    {	
        /* ******
        to run it do:
        root -l
        .L ../../src/tdrStyle.C
        .x drawSyst.C
        */
 
    	setTDRStyle();
    	gStyle->SetOptStat(0);

	    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 300);
	    TPad *pad1 = new TPad("pad1","pad1",0,0,1,1);
	    pad1->SetTopMargin(0.07);
	    pad1->SetBottomMargin(0.20);
	    pad1->SetRightMargin(0.02);
        pad1->SetLeftMargin(0.08);
	    pad1->Draw();
	    pad1->cd();

        std::string name = "ptlead";
	    TFile *_file0 = TFile::Open("shapeFile_ptlead2016.root");

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
		hUp[0]->GetYaxis()->SetTitle("Relative unc.");
        hUp[0]->GetYaxis()->SetTitleOffset(0.6);
        hUp[0]->GetYaxis()->SetLabelSize(0.07);

        hUp[0]->GetXaxis()->SetTitle("Trailing lepton p_{T} [GeV]");
        hUp[0]->GetXaxis()->SetTitleOffset(1.10);
        hUp[0]->GetXaxis()->SetTitleSize(0.08);
        hUp[0]->GetXaxis()->SetLabelSize(0.08);

		TLegend* mtleg = new TLegend(0.1,0.72,0.95,0.92); 
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


    	TH1D * hStat;
    	hStat = (TH1D*)hCentral->Clone("hStat");
    	hStat->Reset("ICE");
		for(int bin = 1; bin < hCentral->GetNbinsX()+1; bin++){
    		TH1D * tempUp = (TH1D*) _file0->Get(processToDraw + "_" + processToDraw + "_stat" + name + "_bin_" + std::to_string(bin) + "Up");
    		TH1D * tempDown = (TH1D*) _file0->Get(processToDraw + "_" + processToDraw + "_stat" + name + "_bin_" + std::to_string(bin) + "Down");
			
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
