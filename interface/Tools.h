const unsigned int indexSR = 9;
const unsigned int indexFlavour = 8; 

#include "errors.h"
#include "../interface/treeReader.h"
#include "../interface/mt2_bisect.h"

bool comp (const pair<double, int> i, const pair<double, int> j) { return (i.first>j.first); }

int SRID2L (int njets, int nbjets, int mvaValueRegion, double chargesLepton) {
    int index = -1;

    int chargesLeptonIndex = (chargesLepton == 1.);
    /*
    if(mvaValueRegion == 2) chargesLeptonIndex = 0.;

    if(njets == 2)
        index = 0 + 5 * mvaValueRegion + 10 * chargesLeptonIndex;

    if(njets == 3 && nbjets == 1)
        index = 1 + 5 * mvaValueRegion + 10 * chargesLeptonIndex;

    if(njets == 3 && nbjets > 1)
        index = 2 + 5 * mvaValueRegion + 10 * chargesLeptonIndex;   

    if(njets > 3 && nbjets == 1)
        index = 3 + 5 * mvaValueRegion + 10 * chargesLeptonIndex;

    if(njets > 3 && nbjets > 1)
        index = 4 + 5 * mvaValueRegion + 10 * chargesLeptonIndex;
    
    if(mvaValueRegion == 2 && njets == 2) 
        index = 20;
    if(mvaValueRegion == 2 && njets == 3) 
        index = 21;
    if(mvaValueRegion == 2 && njets > 3) 
        index = 22;
    */
    if(njets == 2)
        index = 0 + 5 * chargesLeptonIndex;

    if(njets == 3 && nbjets == 1)
        index = 1 + 5 * chargesLeptonIndex;

    if(njets == 3 && nbjets > 1)
        index = 2 + 5 * chargesLeptonIndex;   

    if(njets > 3 && nbjets == 1)
        index = 3 + 5 * chargesLeptonIndex;

    if(njets > 3 && nbjets > 1)
        index = 4 + 5 * chargesLeptonIndex;
    
  
    return index;

}

int SRID3L (int njets, int nbjets) {
    int index = -1;
    
    const int njetsCategories = 3;

    int njetsIndex;
    if (njets == 2) njetsIndex = 0;
    else if (njets == 3) njetsIndex = 1;
    else if (njets >= 4) njetsIndex = 2;
    
    int nbjetsIndex;
    if (nbjets == 0) nbjetsIndex = 0;
    else if (nbjets == 1) nbjetsIndex = 1;
    else if (nbjets >= 2) nbjetsIndex = 2;

    index = nbjetsIndex * njetsCategories + njetsIndex;

    return index;

}

void initdistribs(std::vector<std::string> & namesOfSamples){

    for (std::map<TString, histInfo>::iterator it=figNames.begin(); it!=figNames.end(); ++it){

      histInfo hist = it->second;
      TString name = Form("varST_%d",hist.index);
      int i = hist.index;
      distribs[i].vectorHistoTotalUnc = std::move(TH1D(name,name+";",hist.nBins,hist.varMin,hist.varMax));

      for(unsigned int j = 0; j < distribs[i].vectorHisto.size(); j++){
        name = Form("var_%d_%d",i,j);
        distribs[i].vectorHisto[j] = std::move(TH1D(name,name+";",hist.nBins,hist.varMin,hist.varMax));

        /*
        name = Form("varUp_%d_%d",i,j);
        distribs[i].vectorHistoUp[j] = std::move(TH1D(name,name+";",hist.nBins,hist.varMin,hist.varMax));

        name = Form("varDown_%d_%d",i,j);
        distribs[i].vectorHistoDown[j] = std::move(TH1D(name,name+";",hist.nBins,hist.varMin,hist.varMax));
        */
        for(unsigned int k = 0; k < 6; k++){
          name = Form("varUp_%d_%d_%d",i,j,k);
          distribs[i].vectorHistoUncUp[j].unc[k] = std::move(TH1D(name,name+";",hist.nBins,hist.varMin,hist.varMax));

          name = Form("varDown_%d_%d_%d",i,j,k);
          distribs[i].vectorHistoUncDown[j].unc[k] = std::move(TH1D(name,name+";",hist.nBins,hist.varMin,hist.varMax));
        }
        
        distribs[i].vectorHisto[j].SetBinErrorOption(TH1::kPoisson);

        distribs[i].vectorHisto[j].SetMarkerStyle(20);
        distribs[i].vectorHisto[j].SetMarkerSize(0.5);
        distribs[i].vectorHisto[j].SetLineWidth(1);
        if (j < nSamples)
          distribs[i].vectorHisto[j].Sumw2();
      }
    }

    for (unsigned int i=0; i!=distribs.size();++i) {
      for(unsigned int j = namesOfSamples.size()-1; j != 0; j--){
        distribs[i].stack.Add(&distribs[i].vectorHisto[j]);
      }
    }
    
    for(auto & histo: distribs[indexFlavour].vectorHisto) {
      for(const auto & i: leptonSelectionAnalysis == 2 ? flavourLabelOptionsFor2L : (leptonSelectionAnalysis == 3 ? flavourLabelOptionsFor3L : flavourLabelOptionsFor4L)){
        histo.GetXaxis()->SetBinLabel(i.index, i.labelSR.c_str());
      }

      histo.GetXaxis()->SetLabelSize(0.1);
      histo.GetXaxis()->SetTitleSize(0.25);
      histo.GetXaxis()->SetLabelOffset(0.02);
    }

    for(const auto & i: leptonSelectionAnalysis == 2 ? flavourLabelOptionsFor2L : (leptonSelectionAnalysis == 3 ? flavourLabelOptionsFor3L : flavourLabelOptionsFor4L)){
      distribs[indexFlavour].vectorHistoTotalUnc.GetXaxis()->SetBinLabel(i.index, i.labelSR.c_str());
    }

    for(auto & histo: distribs[indexSR].vectorHisto) {
      for(const auto & i: leptonSelectionAnalysis == 2 ? theSRLabelOptionsFor2L : theSRLabelOptionsFor3L){
        histo.GetXaxis()->SetBinLabel(i.index, i.labelSR.c_str());
      }

      histo.GetXaxis()->SetTitleSize(0.15);
      histo.GetXaxis()->SetLabelOffset(0.02);
    }

    /*
    for(const auto & i: leptonSelectionAnalysis == 2 ? theSRLabelOptionsFor2L : theSRLabelOptionsFor3L){
      distribs[indexSR].vectorHistoTotalUnc.GetXaxis()->SetBinLabel(i.index, i.labelSR.c_str());
    }
    */
    for(auto & histo: distribs[indexSR].vectorHistoUncUp) 
        for(const auto & i: leptonSelectionAnalysis == 2 ? theSRLabelOptionsFor2L : theSRLabelOptionsFor3L)
            for(unsigned int k = 0; k < 6; k++)
                histo.unc[k].GetXaxis()->SetBinLabel(i.index, i.labelSR.c_str());

    for(auto & histo: distribs[indexSR].vectorHistoUncDown) 
        for(const auto & i: leptonSelectionAnalysis == 2 ? theSRLabelOptionsFor2L : theSRLabelOptionsFor3L)
            for(unsigned int k = 0; k < 6; k++)
                histo.unc[k].GetXaxis()->SetBinLabel(i.index, i.labelSR.c_str());

    for(auto & histo: distribs[indexFlavour].vectorHistoUncUp) {
        for(unsigned int k = 0; k < 6; k++){
            for(const auto & i: leptonSelectionAnalysis == 2 ? flavourLabelOptionsFor2L : (leptonSelectionAnalysis == 3 ? flavourLabelOptionsFor3L : flavourLabelOptionsFor4L))
                histo.unc[k].GetXaxis()->SetBinLabel(i.index, i.labelSR.c_str());

            histo.unc[k].GetXaxis()->SetLabelSize(0.1);
            histo.unc[k].GetXaxis()->SetTitleSize(0.25);
            histo.unc[k].GetXaxis()->SetLabelOffset(0.02);
        }
    }

    for(auto & histo: distribs[indexFlavour].vectorHistoUncDown) {
        for(unsigned int k = 0; k < 6; k++){
            for(const auto & i: leptonSelectionAnalysis == 2 ? flavourLabelOptionsFor2L : (leptonSelectionAnalysis == 3 ? flavourLabelOptionsFor3L : flavourLabelOptionsFor4L))
                histo.unc[k].GetXaxis()->SetBinLabel(i.index, i.labelSR.c_str());

            histo.unc[k].GetXaxis()->SetLabelSize(0.1);
            histo.unc[k].GetXaxis()->SetTitleSize(0.25);
            histo.unc[k].GetXaxis()->SetLabelOffset(0.02);
        }
    }

}

double flavourCategory2L(int nLocEle, int chargesTwoLepton){
  return double(1 + nLocEle + 3 * (chargesTwoLepton > 0 ? 1 : 0));
}

double flavourCategory3L(int nLocEle){
  return double(1 + nLocEle);
}

double flavourCategory4L(int nLocEle){
  return double(1 + nLocEle / 2);
}

void addBranchToBDTTreeVariables(){
    signalTree->Branch("nJLoc", &nJLoc, "nJLoc/I");
    signalTree->Branch("nBLoc", &nBLoc, "nBLoc/I");
    signalTree->Branch("HTLoc", &HTLoc, "HTLoc/D");
    signalTree->Branch("_met", &MET, "_met/D");

    signalTree->Branch("_weight", &_weightEventInTree, "_weight/D");

    signalTree->Branch("minDeltaRlead", &minDeltaRlead, "minDeltaRlead/D");
    signalTree->Branch("minDeltaR", &minDeltaR, "minDeltaR/D");
    signalTree->Branch("mt", &mtHighest, "mt/D");
    signalTree->Branch("mtlow", &mtLowest, "mtlow/D");

    signalTree->Branch("leadpt", &leadpt, "leadpt/D");
    signalTree->Branch("trailpt", &trailpt, "trailpt/D");
    signalTree->Branch("leadeta", &leadeta, "leadeta/D");
    signalTree->Branch("traileta", &traileta, "traileta/D");
    signalTree->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/D");
    signalTree->Branch("trailJetPt", &trailJetPt, "trailJetPt/D");
    signalTree->Branch("chargeOfLeptons", &chargeOfLeptons, "chargeOfLeptons/I");
    signalTree->Branch("mll_ss", &mll_ss, "mll_ss/D");
    signalTree->Branch("ll_deltaR", &ll_deltaR, "ll_deltaR/D");
    signalTree->Branch("mt2ll_ss", &mt2ll_ss, "mt2ll_ss/D");


    bkgTree->Branch("nJLoc", &nJLoc, "nJLoc/I");
    bkgTree->Branch("nBLoc", &nBLoc, "nBLoc/I");
    bkgTree->Branch("HTLoc", &HTLoc, "HTLoc/D");
    bkgTree->Branch("_met", &MET, "_met/D");

    bkgTree->Branch("_weight", &_weightEventInTree, "_weight/D");

    bkgTree->Branch("minDeltaRlead", &minDeltaRlead, "minDeltaRlead/D");
    bkgTree->Branch("minDeltaR", &minDeltaR, "minDeltaR/D");
    bkgTree->Branch("mt", &mtHighest, "mt/D");
    bkgTree->Branch("mtlow", &mtLowest, "mtlow/D");

    bkgTree->Branch("leadpt", &leadpt, "leadpt/D");
    bkgTree->Branch("trailpt", &trailpt, "trailpt/D");
    bkgTree->Branch("leadeta", &leadeta, "leadeta/D");
    bkgTree->Branch("traileta", &traileta, "traileta/D");
    bkgTree->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/D");  
    bkgTree->Branch("trailJetPt", &trailJetPt, "trailJetPt/D");
    bkgTree->Branch("chargeOfLeptons", &chargeOfLeptons, "chargeOfLeptons/I");
    bkgTree->Branch("mll_ss", &mll_ss, "mll_ss/D");
    bkgTree->Branch("ll_deltaR", &ll_deltaR, "ll_deltaR/D");
    bkgTree->Branch("mt2ll_ss", &mt2ll_ss, "mt2ll_ss/D");
}

void addVariablesToBDT(const bool is2017 = false){

    readerTTWcsttbar->AddVariable( "nJLoc", &usernJLoc );
    readerTTWcsttbar->AddVariable( "nBLoc", &usernBLoc );
    readerTTWcsttbar->AddVariable( "HTLoc", &userHTLoc ); 
    readerTTWcsttbar->AddVariable( "_met", &user_met );
    readerTTWcsttbar->AddVariable( "minDeltaRlead", &userminDeltaRlead );
    readerTTWcsttbar->AddVariable( "minDeltaR", &userminDeltaR );
    readerTTWcsttbar->AddVariable( "ll_deltaR", &userll_deltaR );
    readerTTWcsttbar->AddVariable( "mt", &usermt );
    readerTTWcsttbar->AddVariable( "mtlow", &usermtlow );
    readerTTWcsttbar->AddVariable( "leadpt", &userleadpt );
    readerTTWcsttbar->AddVariable( "trailpt", &usertrailpt );
    readerTTWcsttbar->AddVariable( "leadingJetPt", &userleadingjetpt );
    readerTTWcsttbar->AddVariable( "trailJetPt", &usertrailjetpt );  
    readerTTWcsttbar->AddVariable( "leadeta", &userleadeta );
    readerTTWcsttbar->AddVariable( "traileta", &usertraileta );
    //readerTTWcsttbar->AddVariable( "chargeOfLeptons", &userchargeOfLeptons);
    readerTTWcsttbar->AddVariable( "mll_ss", &usermll_ss );
    readerTTWcsttbar->AddVariable( "mt2ll_ss", &usermt2ll_ss );

    // the one used for leptonMVA
    // obtained from TMVA training
    TString dir    = "MVAtrainings/" + (TString)(is2017 ? "2017MC" : "2016MC") + "/ttVvsNPchargeMisIDplusDDnonprompt/dataset/weights/"; 

    // used for cut based
    //TString dir    = "/Users/illiakhvastunov/Desktop/CERN/ss2l_2016_fulldataset/analysis/ttWvsttbar_MC_newJEC_fixed/weights/";
    
    TString prefix = "TMVAClassification";
      
    TString methodName = TString("BDTG") + TString(" method");
    TString weightfile = dir + prefix + TString("_") + TString("BDTG") + TString(".weights.xml");
    //TString methodName = TString("BDT::BDT");
    //TString weightfile = TString("MVAtrainings/2017MC/ttVvsNPchargeMisID/trainingInSKlearn/bdt.weights.xml");
    readerTTWcsttbar->BookMVA( methodName, weightfile ); 
}

void fillBDTvariables(vector<Float_t> & varForBDT){

    usernJLoc = varForBDT.at(0);
    usernBLoc = varForBDT.at(1);
    userHTLoc = varForBDT.at(2);
    user_met = varForBDT.at(3);
    userminDeltaRlead = varForBDT.at(4);
    userminDeltaR = varForBDT.at(5);
    userll_deltaR = varForBDT.at(6);
    usermt = varForBDT.at(7);
    usermtlow = varForBDT.at(8);            
    userleadpt = varForBDT.at(9);
    usertrailpt = varForBDT.at(10);
    userleadingjetpt = varForBDT.at(11);
    usertrailjetpt = varForBDT.at(12);
    userleadeta = fabs(varForBDT.at(13));
    usertraileta = fabs(varForBDT.at(14));
    userchargeOfLeptons = varForBDT.at(15);
    usermll_ss = varForBDT.at(16);
    usermt2ll_ss = varForBDT.at(17);
}


void setStackColors(Color_t & color, int sam){

    for(unsigned int i = 0; i < distribs.size(); i++){
        distribs[i].vectorHisto[sam].SetLineColor(color);
        distribs[i].vectorHisto[sam].SetFillColor(color);
        distribs[i].vectorHisto[sam].SetMarkerColor(color);
    }
}

double mt2ll(const TLorentzVector& l1, const TLorentzVector& l2, const TLorentzVector& metVec){
    return  asymm_mt2_lester_bisect::get_mT2(l1.M(), l1.Px(), l1.Py(), l2.M(), l2.Px(), l2.Py(), metVec.Px(), metVec.Py(), 0, 0);
}

double getTriggerSF(int leptonSelection, const double leptPt, bool is2017 = false){
    if(leptonSelection == 3 && leptPt < 120 && is2017)
        //weight *= 0.985; // got it from Daniel's trigger measurement, agreed with him on skype, 5 Jul 2018
        return 0.985;
    if(leptonSelection == 3 && leptPt < 80 && !is2017)
        return 0.966000020504;
    return 1.;
}
