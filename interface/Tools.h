const unsigned int indexSR = 15;
const unsigned int indexFlavour = 14; 

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

void getFRmaps(vector<TH2D> & fakeMaps){

    //TFile *fakerate = TFile::Open("/Users/illiakhvastunov/Desktop/CERN/ss2l_2016_fulldataset/analysis_withMVATTH/FRmaps/FR_data_ttH_mva.root","READ");
    TFile *fakerate = NULL;
    

    for (int i=0; i!=nFlavors; ++i){

      fakerate = TFile::Open("FRmapsQCDmine/fakerate_" + flavorsString[i] + "_QCD.root","READ");
      TH2D * tempPtr = (TH2D*) (fakerate->Get("fakerate_" + flavorsString[i]));

      //TH2D * tempPtr = (TH2D*) (fakerate->Get("FR_mva090_" + flavorsString[i] + "_data_comb" + additionalString[i]));
      //TH2D * tempPtr = (TH2D*) (fakerate->Get("FR_mva090_" + flavorsString[i] + "_QCD")); //  + additionalString[i]
      
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

}


int flavourCategory2L(int nLocEle, int chargesTwoLepton){
  
  return 1 + nLocEle + 3 * (chargesTwoLepton > 0 ? 1 : 0);

}


int flavourCategory3L(int nLocEle){
  
  return 1 + nLocEle;

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

void addBranchToBDTTreeVariables(){
    signalTree->Branch("nJLoc", &nJLoc, "nJLoc/I");
    signalTree->Branch("nBLoc", &nBLoc, "nBLoc/I");
    signalTree->Branch("HTLoc", &HTLoc, "HTLoc/D");
    signalTree->Branch("_met", &MET, "_met/D");

    signalTree->Branch("_weight", &_weightEventInTree, "_weight/D");

    signalTree->Branch("minDeltaR", &minDeltaR, "minDeltaR/D");
    signalTree->Branch("mt", &mtHighest, "mt/D");
    signalTree->Branch("mtlow", &mtLowest, "mtlow/D");

    signalTree->Branch("leadpt", &leadpt, "leadpt/D");
    signalTree->Branch("trailpt", &trailpt, "trailpt/D");
    signalTree->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/D");
    signalTree->Branch("trailJetPt", &trailJetPt, "trailJetPt/D");


    bkgTree->Branch("nJLoc", &nJLoc, "nJLoc/I");
    bkgTree->Branch("nBLoc", &nBLoc, "nBLoc/I");
    bkgTree->Branch("HTLoc", &HTLoc, "HTLoc/D");
    bkgTree->Branch("_met", &MET, "_met/D");

    bkgTree->Branch("_weight", &_weightEventInTree, "_weight/D");

    bkgTree->Branch("minDeltaR", &minDeltaR, "minDeltaR/D");
    bkgTree->Branch("mt", &mtHighest, "mt/D");
    bkgTree->Branch("mtlow", &mtLowest, "mtlow/D");

    bkgTree->Branch("leadpt", &leadpt, "leadpt/D");
    bkgTree->Branch("trailpt", &trailpt, "trailpt/D");
    bkgTree->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/D");  
    bkgTree->Branch("trailJetPt", &trailJetPt, "trailJetPt/D");
}

void addVariablesToBDT(){

    reader->AddVariable( "HTLoc", &userHTLoc ); 
    reader->AddVariable( "nJLoc", &usernJLoc );
    reader->AddVariable( "nBLoc", &usernBLoc );
    reader->AddVariable( "_met", &user_met );
    reader->AddVariable( "minDeltaR", &userminDeltaR );
    reader->AddVariable( "mt", &usermt );
    reader->AddVariable( "mtlow", &usermtlow );
    reader->AddVariable( "leadpt", &userleadpt );
    reader->AddVariable( "trailpt", &usertrailpt );
    reader->AddVariable( "leadingJetPt", &userleadingjetpt );
    reader->AddVariable( "trailJetPt", &usertrailjetpt );  

    TString dir    = "ttWvsttbarMC_LeptonMVA_0p9/weights/";
    TString prefix = "TMVAClassification";
      
    TString methodName = TString("BDTG") + TString(" method");
    TString weightfile = dir + prefix + TString("_") + TString("BDTG") + TString(".weights.xml");
    reader->BookMVA( methodName, weightfile ); 
}

void fillBDTvariables(vector<Float_t> & varForBDT){

    usernJLoc = varForBDT.at(0);
    usernBLoc = varForBDT.at(1);
    userHTLoc = varForBDT.at(2);
    user_met = varForBDT.at(3);
    userminDeltaR = varForBDT.at(4);
    userleadpt = varForBDT.at(5);
    usertrailpt = varForBDT.at(6);
    usermt = varForBDT.at(7);
    usermtlow = varForBDT.at(8);            
    userleadingjetpt = varForBDT.at(9);
    usertrailjetpt = varForBDT.at(10);
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
std::map<std::string, unsigned> pair_to_map(std::vector<std::string> & a, std::vector<unsigned> b)
{
    std::map<std::string, unsigned> my_map;
    for (unsigned i = 0; i < a.size(); ++i)
    {
        my_map[a] = b;
    }
    return my_map;
}
*/
