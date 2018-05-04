const unsigned int indexSR = 12;
const unsigned int indexFlavour = 13; 

#include "errors.h"
#include "../interface/treeReader.h"

bool comp (const pair<double, int> i, const pair<double, int> j) { return (i.first>j.first); }

int SRID2L (int njets, int nbjets, int mvaValueRegion, double chargesLepton) {
    int index = -1;

    int chargesLeptonIndex = (chargesLepton == 1.);
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

void getFRmaps(vector<TH2D> & fakeMaps, bool is2017 = false){

    TFile *fakerate = NULL;

    for (int i=0; i!=nFlavors; ++i){

      //fakerate = TFile::Open("data/FRmaps/" + flavorsString[i] + "FR_ttbar_" + (TString)(is2017 ? "2017" : "2016") + "_noChargeAgreementCut.root","READ");
      fakerate = TFile::Open("data/FRmaps/" + flavorsString[i] + "FR.root","READ");
      TH2D * tempPtr = (TH2D*) (fakerate->Get("passed"));

      if(tempPtr != NULL){
        fakeMaps.push_back(*tempPtr);
      }
      else{
        LastError::lasterror = Errors::FRmapsReturnNULL;
        return;
      }

    }
      
    LastError::lasterror = Errors::OK;
    
}


void initdistribs(std::vector<std::string> & namesOfSamples){

    for(unsigned int i = 0; i < distribs.size(); i++){
      TString name = Form("varST_%d",i);
      //distribs[i].colsStack = std::move(THStack(name,varN[i]));
      distribs[i].stack.SetName(name);
      distribs[i].stack.SetTitle(varN[i]);

      distribs[i].vectorHistoTotalUnc = std::move(TH1D(name,name+";",nBins[i],varMin[i],varMax[i]));

      for(unsigned int j = 0; j < distribs[i].vectorHisto.size(); j++){
        name = Form("var_%d_%d",i,j);
        distribs[i].vectorHisto[j] = std::move(TH1D(name,name+";",nBins[i],varMin[i],varMax[i]));

        name = Form("varUp_%d_%d",i,j);
        distribs[i].vectorHistoUp[j] = std::move(TH1D(name,name+";",nBins[i],varMin[i],varMax[i]));

        name = Form("varDown_%d_%d",i,j);
        distribs[i].vectorHistoDown[j] = std::move(TH1D(name,name+";",nBins[i],varMin[i],varMax[i]));
        
        distribs[i].vectorHisto[j].SetBinErrorOption(TH1::kPoisson);
        /*
        distribs[i].vectorHisto[j].SetLineColor(colsStack[j]);
        if (j < nSamples-1)
          distribs[i].vectorHisto[j].SetFillColor(colsStack[j]);
        distribs[i].vectorHisto[j].SetMarkerColor(colsStack[j]);
        */
        distribs[i].vectorHisto[j].SetMarkerStyle(20);
        distribs[i].vectorHisto[j].SetMarkerSize(0.5);
        distribs[i].vectorHisto[j].SetLineWidth(1);
        if (j < nSamples)
          distribs[i].vectorHisto[j].Sumw2();
      }
    }

    // 0 - fakes, 1-2 ttZ, 3-WZ, 4--21 -rares, 22 - data

    
    for (unsigned int i=0; i!=distribs.size();++i) {
      /*
      for(unsigned int j = 1; j != distribsOrder.size(); j++){
        distribs[i].stack.Add(&distribs[i].vectorHisto[distribsOrder[j]]);
      }
      */
      for(unsigned int j = namesOfSamples.size()-1; j != 0; j--){
      //for(unsigned int j = distribs.size()-1; j != 0; j--){
        distribs[i].stack.Add(&distribs[i].vectorHisto[j]);
      }
    }

    TString name = Form("var2D");
    distribs2D.vectorHisto[0] = std::move(TH2D(name,name+";",nBins2D[0],varMin2D[0],varMax2D[0],nBins2D[1],varMin2D[1],varMax2D[1]));

    name = Form("FRele");
    distribs[0].vectorHisto2D[0] = std::move(TH2D(name,           name+";",nPt-1, ptBins, nEta-1, etaBins[0]));
    distribs[0].vectorHisto2D[1] = std::move(TH2D(name + "passed",name+";",nPt-1, ptBins, nEta-1, etaBins[0]));
    name = Form("FRmu");
    distribs[1].vectorHisto2D[0] = std::move(TH2D(name,           name+";",nPt-1, ptBins, nEta-1, etaBins[1]));
    distribs[1].vectorHisto2D[1] = std::move(TH2D(name + "passed",name+";",nPt-1, ptBins, nEta-1, etaBins[1]));

    for(auto & histo: distribs[indexFlavour].vectorHisto) {
        for(const auto & i: leptonSelectionAnalysis == 2 ? flavourLabelOptionsFor2L : flavourLabelOptionsFor3L){
            histo.GetXaxis()->SetBinLabel(i.index, i.labelSR.c_str());
        }

        histo.GetXaxis()->SetLabelSize(0.1);
        histo.GetXaxis()->SetTitleSize(0.25);
        histo.GetXaxis()->SetLabelOffset(0.01);
    }

    for(const auto & i: leptonSelectionAnalysis == 2 ? flavourLabelOptionsFor2L : flavourLabelOptionsFor3L){
        distribs[indexFlavour].vectorHistoTotalUnc.GetXaxis()->SetBinLabel(i.index, i.labelSR.c_str());
    }
    distribs[indexFlavour].vectorHistoTotalUnc.GetXaxis()->SetLabelSize(0.1);
    distribs[indexFlavour].vectorHistoTotalUnc.GetXaxis()->SetTitleSize(0.25);
    distribs[indexFlavour].vectorHistoTotalUnc.GetXaxis()->SetLabelOffset(0.02);

    for(auto & histo: distribs[indexSR].vectorHisto) {
      
      for(const auto & i: leptonSelectionAnalysis == 2 ? theSRLabelOptionsFor2L : theSRLabelOptionsFor3L){
        histo.GetXaxis()->SetBinLabel(i.index, i.labelSR.c_str());
      }

      histo.GetXaxis()->SetTitleSize(0.15);
      histo.GetXaxis()->SetLabelOffset(0.02);
    }

}


int flavourCategory2L(int nLocEle, int chargesTwoLepton){
  
  return 1+nLocEle + 3 * (chargesTwoLepton > 0 ? 1 : 0);

}


int flavourCategory3L(int nLocEle){
  
  return 1+nLocEle;

}

float getLeptonSF(int flavour, Float_t pt, Float_t eta, float var, int eraDecision){

    float lepSF = 1.;

    if(flavour == 0){
      TH2F * hist = lepSFMapsElectron[0];
      int etabin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(eta))); // careful, different convention
      int ptbin  = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(pt)));
      lepSF *= hist->GetBinContent(etabin,ptbin) + var * hist->GetBinError(etabin,ptbin) ;

      
      for(int sfFile = 4; sfFile < 5; sfFile++){
        TH2F *hist = lepSFMapsElectron[sfFile];
        ptbin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
        etabin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(TMath::Abs(eta))));
        lepSF *= hist->GetBinContent(ptbin,etabin) + var * hist->GetBinError(ptbin, etabin) ;
      }
      
    }
                
    if(flavour == 1){
                  
      lepSF *= lepSFMaps1DMuon[eraDecision]->Eval(eta);

      //lepSF *= lepSFMaps1DMuon[eraDecision]->Eval(eta) + var * TMath::Max(lepSFMaps1DMuon[eraDecision]->GetErrorY(eta), lepSFMaps1DMuon[eraDecision]->GetErrorY(eta));

      
      for(int sfFile = 0; sfFile < 4; sfFile++){
        TH2F *hist = lepSFMapsMuon[sfFile];
        int ptbin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
        int etabin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(TMath::Abs(eta))));
        lepSF *= hist->GetBinContent(ptbin,etabin) + var * hist->GetBinError(ptbin, etabin);
      } 
        
    }

    return lepSF;
}


float getBTagSF(int btagFileDicision, float var, int jf, float eta, float pt, float csv){

   bool isBFlav = false;
   bool isCFlav = false;
   bool isLFlav = false;

   BTagEntry::JetFlavor jfInput = BTagEntry::FLAV_UDSG;
   int jetFlavorVariation = 0;
   if( abs(jf) == 5 ){
    jfInput = BTagEntry::FLAV_B;
    jetFlavorVariation = 0;
   }
   else if( abs(jf) == 4 ){
    jfInput = BTagEntry::FLAV_C;
    jetFlavorVariation = 1;
   } 
   else{
    jfInput = BTagEntry::FLAV_UDSG;
    jetFlavorVariation = 2;
   }

   std::string varStr;

   if(var == 0) varStr = "central";
   else if(var == -1) varStr = "down_hf";
   else if(var == -2) varStr = "down_lf";
   else if(var == 1) varStr = "up_hf";
   else if(var == 2) varStr = "up_lf";

   float sf = readerBtag[btagFileDicision][jetFlavorVariation].eval_auto_bounds(varStr, jfInput, fabs(eta), pt, csv);
   //float sf = readerBtag[btagFileDicision][0].eval_auto_bounds(varStr, jfInput, fabs(eta), pt, csv);
                                              
   return sf;      
}



double mtCalc(TLorentzVector Vect, double MET, double MET_Phi){

    double MT=sqrt(2* Vect.Pt() * MET * ( 1 - (TMath::Cos(Vect.Phi() - MET_Phi )) ) );

    return MT;
}

void setStackColors(Color_t & color, int sam){

    for(unsigned int i = 0; i < distribs.size(); i++){
        distribs[i].vectorHisto[sam].SetLineColor(color);
        distribs[i].vectorHisto[sam].SetFillColor(color);
        distribs[i].vectorHisto[sam].SetMarkerColor(color);
    }
}



/*
void addBranchToBDTTreeVariables(){

    signalTree->Branch("pt", &LepGood_pt, "pt/D");
    signalTree->Branch("eta", &LepGood_eta, "eta/D");
    //signalTree->Branch("lepSelTrackMult", &LepGood_jetNDauChargedMVASel, "trackMult/D");
    signalTree->Branch("miniIsoCharged", &LepGood_miniRelIsoCharged, "miniIsoCharged/D");
    signalTree->Branch("miniIsoNeutral", &LepGood_miniRelIsoNeutral, "miniIsoNeutral/D");
    signalTree->Branch("ptrel", &LepGood_jetPtRelv2, "ptrel/D");
    signalTree->Branch("ptratio", &LepGood_jetPtRatio, "ptratio/D");
    signalTree->Branch("jetBtagCSV", &LepGood_jetBTagCSV, "btagCSV/D");
    signalTree->Branch("sip3d", &LepGood_sip3d, "sip3d/D");
    signalTree->Branch("dxy", &LepGood_dxy, "dxy/D");
    signalTree->Branch("dz", &LepGood_dz, "dz/D");
    signalTree->Branch("segmComp", &LepGood_segmentCompatibility, "segmComp/D");
    signalTree->Branch("_weight", &_weightEventInTree, "_weight/D");

    bkgTree->Branch("pt", &LepGood_pt, "pt/D");
    bkgTree->Branch("eta", &LepGood_eta, "eta/D");
    //bkgTree->Branch("lepSelTrackMult", &LepGood_jetNDauChargedMVASel, "trackMult/D");
    bkgTree->Branch("miniIsoCharged", &LepGood_miniRelIsoCharged, "miniIsoCharged/D");
    bkgTree->Branch("miniIsoNeutral", &LepGood_miniRelIsoNeutral, "miniIsoNeutral/D");
    bkgTree->Branch("ptrel", &LepGood_jetPtRelv2, "ptrel/D");
    bkgTree->Branch("ptratio", &LepGood_jetPtRatio, "ptratio/D");
    bkgTree->Branch("jetBtagCSV", &LepGood_jetBTagCSV, "btagCSV/D");
    bkgTree->Branch("sip3d", &LepGood_sip3d, "sip3d/D");
    bkgTree->Branch("dxy", &LepGood_dxy, "dxy/D");
    bkgTree->Branch("dz", &LepGood_dz, "dz/D");
    bkgTree->Branch("segmComp", &LepGood_segmentCompatibility, "segmComp/D");
    bkgTree->Branch("_weight", &_weightEventInTree, "_weight/D");

}
*/

void addVariablesToBDT(){

    readerLeptonMVAele->AddVariable( "pt", &user_pt); 
    readerLeptonMVAele->AddVariable( "eta", &user_eta); 
    readerLeptonMVAele->AddVariable( "trackMult", &user_trackMult); 
    readerLeptonMVAele->AddVariable( "miniIsoCharged", &user_miniIsoCharged); 
    readerLeptonMVAele->AddVariable( "miniIsoNeutral", &user_miniIsoNeutral); 
    readerLeptonMVAele->AddVariable( "ptrel", &user_ptrel); 
    readerLeptonMVAele->AddVariable( "min(ptratio,1.5)", &user_ptratio); 
    readerLeptonMVAele->AddVariable( "relIso0p3", &user_relIso);
    readerLeptonMVAele->AddVariable( "max(jetbtagCSV,0)", &user_jetBtagCSV); 
    readerLeptonMVAele->AddVariable( "sip3d", &user_sip3d); 
    //readerLeptonMVAele->AddVariable( "log(abs(dxy))", &user_dxy); 
    //readerLeptonMVAele->AddVariable( "log(abs(dz))", &user_dz); 
    readerLeptonMVAele->AddVariable( "eleMVA", &user_eleMVA); 

    readerLeptonMVAmu->AddVariable( "pt", &user_pt); 
    readerLeptonMVAmu->AddVariable( "eta", &user_eta); 
    readerLeptonMVAmu->AddVariable( "trackMult", &user_trackMult); 
    readerLeptonMVAmu->AddVariable( "miniIsoCharged", &user_miniIsoCharged); 
    readerLeptonMVAmu->AddVariable( "miniIsoNeutral", &user_miniIsoNeutral); 
    readerLeptonMVAmu->AddVariable( "ptrel", &user_ptrel); 
    readerLeptonMVAmu->AddVariable( "min(ptratio,1.5)", &user_ptratio); 
    readerLeptonMVAmu->AddVariable( "relIso0p3", &user_relIso);
    readerLeptonMVAmu->AddVariable( "max(jetbtagCSV,0)", &user_jetBtagCSV); 
    readerLeptonMVAmu->AddVariable( "sip3d", &user_sip3d); 
    //readerLeptonMVAmu->AddVariable( "log(abs(dxy))", &user_dxy); 
    //readerLeptonMVAmu->AddVariable( "log(abs(dz))", &user_dz); 
    readerLeptonMVAmu->AddVariable( "segmComp", &user_segmComp); 

    //TString dirEle    = "../checkLeptonMVAvar/MVAtrainings/2017MC/ele_withoutIsoReq_removePromptFromTaus/dataset/weights/"; 
    //TString dirMu    = "../checkLeptonMVAvar/MVAtrainings/2017MC/muon_removePromptFromTaus/dataset/weights/"; 

    //TString dirEle    = "../checkLeptonMVAvar/MVAtrainings/2016MC/ele_dRMatchingToGen_removePromptFromTaus_removedIP_relIso_tryBOTHptratioANDrelISO/dataset/weights/";  // _removePromptFromTaus
    //TString dirMu    = "../checkLeptonMVAvar/MVAtrainings/2016MC/muon_dRMatchingToGen_removePromptFromTaus_removedIP_relIso_tryBOTHptratioANDrelISO/dataset/weights/"; 

    TString dirEle    = "../checkLeptonMVAvar/MVAtrainings/2016MC/ele_relISO0p3_ptGR10/dataset/weights/";
    TString dirMu    = "../checkLeptonMVAvar/MVAtrainings/2016MC/muon_relISO0p3_ptGR10/dataset/weights/";

    TString prefix = "TMVAClassification";
      
    TString methodName = TString("BDTG") + TString(" method");
    TString weightfile = dirEle + prefix + TString("_") + TString("BDTG") + TString(".weights.xml");
    readerLeptonMVAele->BookMVA( methodName, weightfile ); 

    weightfile = dirMu + prefix + TString("_") + TString("BDTG") + TString(".weights.xml");
    readerLeptonMVAmu->BookMVA( methodName, weightfile ); 

} 

void fillBDTvariables(vector<Float_t> & varForBDT, int flavor){

    user_pt =  varForBDT.at(0);
    user_eta = varForBDT.at(1);
    user_trackMult = varForBDT.at(2);
    user_miniIsoCharged = varForBDT.at(3);
    user_miniIsoNeutral = varForBDT.at(4);
    user_ptrel = varForBDT.at(5);
    user_ptratio = varForBDT.at(6);
    user_relIso = varForBDT.at(7);
    user_jetBtagCSV = varForBDT.at(8);
    user_sip3d = varForBDT.at(9);
    user_dxy = varForBDT.at(10);
    user_dz = varForBDT.at(11);
    if(flavor == 0)
       user_eleMVA = varForBDT.at(12);
    if(flavor == 1)
       user_segmComp = varForBDT.at(12);
}

