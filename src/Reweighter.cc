#include "../interface/Reweighter.h"

//include c++ library classes

//include ROOT classes
#include "TFile.h"
#include "TROOT.h"

//include other parts of code 
#include "../interface/analysisTools.h"

Reweighter::Reweighter(const std::vector<Sample>& samples, const bool sampleIs2016, const int lepSel) : is2016(sampleIs2016), leptonSelection(lepSel) {
    initializeAllWeights(samples);
}

void Reweighter::initializeAllWeights(const std::vector<Sample>& samples){

    //initialize pu weights
    initializePuWeights(samples);

    //initialize b-tag weights
    initializeBTagWeights();

    //initialize electron weights 
    initializeElectronWeights();

    //initialize muon weights 
    initializeMuonWeights();

    //initialize fake-rate
    initializeFakeRate();

    //initialize charge mis ID rate
    initializeCMIDRate();

    //initialize prefiring probabilities
    initializePrefiringProbabilities();
}

void Reweighter::initializePuWeights(const std::vector< Sample >& sampleList){

    static const std::string minBiasVariation[3] = {"central", "down", "up"};
    for(auto& sample : sampleList){
        //no pu weights for data 
        if( sample.isData() ) continue;

        //open root file corresponding to sample
        std::string str = sample.getFileName();
        TFile* puFile = TFile::Open( (const TString&) "data/pileUpReweighing/puWeights_" + sample.getFileName());

        //extract pu weights 
        for(unsigned var = 0; var < 3; ++var){
            std::string histName = "puw_Run";
            histName += (sample.is2016() ? "2016" : "2017");
            histName += "Inclusive_" + minBiasVariation[var];
            puWeights[sample.getUniqueName()].push_back( std::shared_ptr<TH1D> ( (TH1D*)puFile->Get( (const TString&) histName ) ) );

            //make sure histogram does not get deleted when closing file
            puWeights[sample.getUniqueName()].back()->SetDirectory(gROOT);
        }
        puFile->Close();
    }
}

void Reweighter::initializeBTagWeights(){

    //assuming medium WP of deepCSV tagger 
    std::string sfFileName;
    if(is2016){
        sfFileName = "DeepCSV_Moriond17_B_H.csv";
    } else {
        sfFileName = "DeepCSV_94XSF_V3_B_F.csv";
    }
    bTagCalib = std::shared_ptr<BTagCalibration>( new BTagCalibration("deepCsv", "data/btagSF/" + sfFileName) );

    //WARNING: b-tagging efficiencies currently assume medium deepCSV tagger!
    bTagCalibReader = std::shared_ptr<BTagCalibrationReader>( new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"}) );
    bTagCalibReader->load(*bTagCalib, BTagEntry::FLAV_B, "comb");
    bTagCalibReader->load(*bTagCalib, BTagEntry::FLAV_C, "comb");
    bTagCalibReader->load(*bTagCalib, BTagEntry::FLAV_UDSG, "incl");    

    //initialize b-tag efficiencies
    std::string effFileName;
    if(is2016){
        if(leptonSelection == 2)
            effFileName = "bTagEff_deepCSV_cleaned_ttW_2016.root";
        if(leptonSelection == 3)
            effFileName = "bTagEff_deepCSV_cleaned_ttZ3l_2016.root";
            //effFileName = "bTagEff_deepCSV_cleaned_tZq_2016.root";
        if(leptonSelection == 4)
            effFileName = "bTagEff_deepCSV_cleaned_ttZ4l_2016.root";
    } else {
        if(leptonSelection == 2)
            effFileName = "bTagEff_deepCSV_cleaned_ttW_2017.root";
        if(leptonSelection == 3)
            effFileName = "bTagEff_deepCSV_cleaned_ttZ3l_2017.root";
            //effFileName = "bTagEff_deepCSV_cleaned_tZq_2017.root";
        if(leptonSelection == 4)
            effFileName = "bTagEff_deepCSV_cleaned_ttZ4l_2017.root";
    }

    TFile* bTagFile = TFile::Open( (const TString&) "data/btagSF/" + effFileName);
    const std::string quarkFlavors[3] = {"udsg", "charm", "beauty"};
    for(unsigned flav = 0; flav < 3; ++flav){
        
        //WARNING: b-tagging efficiencies currently assume medium deepCSV tagger!
        bTagEffHist[flav] = std::shared_ptr<TH2D>( (TH2D*) bTagFile->Get( (const TString&) "bTagEff_medium" + quarkFlavors[flav]) );
        bTagEffHist[flav]->SetDirectory(gROOT);
    }
    bTagFile->Close();
}

void Reweighter::initializeElectronWeights(){   

    //read electron reco SF weights
    if( is2016 ){
        TFile* electronRecoFile_pT0to20 = TFile::Open("data/leptonSF/EGM2D_BtoH_low_RecoSF_Legacy2016.root");
        electronRecoSF_pT0to20 = std::shared_ptr<TH2D>( (TH2D*) electronRecoFile_pT0to20->Get("EGamma_SF2D") );
        electronRecoSF_pT0to20->SetDirectory(gROOT);
        electronRecoFile_pT0to20->Close();
        TFile* electronRecoFile_pT20toInf = TFile::Open("data/leptonSF/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root");
        electronRecoSF_pT20toInf = std::shared_ptr<TH2D>( (TH2D*) electronRecoFile_pT20toInf->Get("EGamma_SF2D") );
        electronRecoSF_pT20toInf->SetDirectory(gROOT);
        electronRecoFile_pT20toInf->Close();
    } else {
        //low pT SF
        TFile* electronRecoFile_pT0to20 = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt_2017.root");
        electronRecoSF_pT0to20 = std::shared_ptr<TH2D>( (TH2D*) electronRecoFile_pT0to20->Get("EGamma_SF2D") );
        electronRecoSF_pT0to20->SetDirectory(gROOT);
        electronRecoFile_pT0to20->Close();

        //high pT SF
        TFile* electronRecoFile_pT20toInf = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_2017.root");
        electronRecoSF_pT20toInf = std::shared_ptr<TH2D>( (TH2D*) electronRecoFile_pT20toInf->Get("EGamma_SF2D") );
        electronRecoSF_pT20toInf->SetDirectory(gROOT);
        electronRecoFile_pT20toInf->Close();
    }

    //read electron ID SF weights
    TFile* electronIdFile; //  = TFile::Open("data/leptonSF/scaleFactors_electrons_2016.root");
    if(is2016){
        electronIdFile = TFile::Open("data/leptonSF/scaleFactors_electrons_2016_upd.root");
    } else {
        electronIdFile = TFile::Open("data/leptonSF/scaleFactors_electrons_2017_upd.root");
    } 
    /*
    electronLooseToRecoSF = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("EleToTTVLoose") );
    electronLooseToRecoSF->SetDirectory(gROOT);
    if(leptonSelection == 2)
        electronTightToLooseSF = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("TTVLooseToTTVLeptonMvattW") );
    if(leptonSelection == 3)
        electronTightToLooseSF = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("TTVLooseToTTVLeptonMvattZ3l") );
        //electronTightToLooseSF = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("TTVLooseToTTVLeptonMvatZq") );
    if(leptonSelection == 4)
        electronTightToLooseSF = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("TTVLooseToTTVLeptonMvattZ4l") );
    electronTightToLooseSF->SetDirectory(gROOT);
    if(leptonSelection == 2){
        electronChargeConsToTightSF = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("TTVLeptonMvattWToTightCharge") );
        electronChargeConsToTightSF->SetDirectory(gROOT);
    }
    */
    if(leptonSelection == 2){
        electronTightToRecoSF_syst = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("TTVLooseToTTVLeptonMvattW_syst") );
        electronTightToRecoSF_stat = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("TTVLooseToTTVLeptonMvattW_stat") );
    }
    if(leptonSelection == 3){
        //electronTightToRecoSF_syst = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("EleToTTVLeptonMvattZ3l_sys") );
        //electronTightToRecoSF_stat = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("EleToTTVLeptonMvattZ3l_stat") );
        electronTightToRecoSF_syst = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("EleToTTVLeptonMvatZq_sys") );
        electronTightToRecoSF_stat = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("EleToTTVLeptonMvatZq_stat") );
    }
    if(leptonSelection == 4){
        electronTightToRecoSF_syst = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("EleToTTVLeptonMvattZ4l_sys") );
        electronTightToRecoSF_stat = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("EleToTTVLeptonMvattZ4l_stat") );
    }
    electronTightToRecoSF_syst->SetDirectory(gROOT);
    electronTightToRecoSF_stat->SetDirectory(gROOT);
    electronIdFile->Close();
}


void Reweighter::initializeMuonWeights(){

    //read muon reco SF weights
    // no reco SF, email from Sergio and Sandro on 30th of July, 2018, they removed as well recommendations from twiki
    if(is2016){
        //WARNING: not clear how ownership works for TGraph, can not set directory
        //make sure the TGraph is not DELETED when file is closed!
        TFile* muonRecoFile = TFile::Open("data/leptonSF/Tracking_EfficienciesAndSF_BCDEFGH.root");
        muonRecoSF = std::shared_ptr<TGraph>( (TGraph*) muonRecoFile->Get("ratio_eff_eta3_dr030e030_corr") );
        muonRecoFile->Close();
    } else {

    }   

    //read muon ID SF weights
    TFile* muonIdFile;
    if(is2016){
        muonIdFile = TFile::Open("data/leptonSF/scaleFactors_muons_2016_upd.root");
    } else {
        muonIdFile = TFile::Open("data/leptonSF/scaleFactors_muons_2017_upd.root"); // use _upd when Tom will update on muon SFs, updated on 14th of September, 2018
    }
    
    muonLooseToRecoSF = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("MuonToTTVLoose") );
    muonLooseToRecoSF->SetDirectory(gROOT);
    /*
    if(leptonSelection == 2)
        muonTightToLooseSF = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("TTVLooseToTTVLeptonMvattW") );
    if(leptonSelection == 3)
        muonTightToLooseSF = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("TTVLooseToTTVLeptonMvattZ3l") );
        //muonTightToLooseSF = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("TTVLooseToTTVLeptonMvatZq") );
    if(leptonSelection == 4)
        muonTightToLooseSF = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("TTVLooseToTTVLeptonMvattZ4l") );
    muonTightToLooseSF->SetDirectory(gROOT);
    if(leptonSelection == 2){
        muonChargeConsToTightSF = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("TTVLeptonMvattWTotkSigmaPtOverPtCut") );
        muonChargeConsToTightSF->SetDirectory(gROOT);
    }
    */
    if(leptonSelection == 2){
        muonTightToRecoSF_syst = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("MuonToTTVLeptonMvattW_sys") );
        muonTightToRecoSF_stat = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("MuonToTTVLeptonMvattW_stat") );
    }
    if(leptonSelection == 3){
        //muonTightToRecoSF_syst = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("MuonToTTVLeptonMvattZ3l_sys") );
        //muonTightToRecoSF_stat = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("MuonToTTVLeptonMvattZ3l_stat") );
        muonTightToRecoSF_syst = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("MuonToTTVLeptonMvatZq_sys") );
        muonTightToRecoSF_stat = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("MuonToTTVLeptonMvatZq_stat") );
    }
    if(leptonSelection == 4){
        muonTightToRecoSF_syst = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("MuonToTTVLeptonMvattZ4l_sys") );
        muonTightToRecoSF_stat = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("MuonToTTVLeptonMvattZ4l_stat") );
    }
    muonTightToRecoSF_syst->SetDirectory(gROOT);
    muonTightToRecoSF_stat->SetDirectory(gROOT);
    muonIdFile->Close();    
}

void Reweighter::initializeFakeRate(){

    //WARNING : To be updated with new fake-rate in 2016/2017 splitting
    TFile* frFileEl = TFile::Open("data/FRmaps/elFR_data_" + TString(leptonSelection == 2 ? "2L_" : "") + (is2016 ? "2016" : "2017") + ".root");
    TFile* frFileMu = TFile::Open("data/FRmaps/muFR_data_" + TString(leptonSelection == 2 ? "2L_" : "") + (is2016 ? "2016" : "2017") + ".root");
    // tZq
    //TFile* frFileEl = TFile::Open("data/FRmaps/elFR_data_" + TString(leptonSelection == 2 ? "2L_" : "3L_tZq_") + (is2016 ? "2016" : "2017") + ".root");
    //TFile* frFileMu = TFile::Open("data/FRmaps/muFR_data_" + TString(leptonSelection == 2 ? "2L_" : "3L_tZq_") + (is2016 ? "2016" : "2017") + ".root");
    //const std::string frUnc[3] = {"", "_down", "_up"};
    const std::string frUnc[3] = {"", "", ""};
    for(unsigned unc = 0; unc < 3; ++unc){
        frMapEle[unc] = std::shared_ptr<TH2D>( (TH2D*) frFileEl->Get((const TString&) "passed" + frUnc[unc]) );
        frMapEle[unc]->SetDirectory(gROOT);
        frMapMu[unc] = std::shared_ptr<TH2D>( (TH2D*) frFileMu->Get((const TString&) "passed" + frUnc[unc]) );
        frMapMu[unc]->SetDirectory(gROOT);
    }
    frFileEl->Close();
    frFileMu->Close();
}

void Reweighter::initializeCMIDRate(){

    //WARNING : To be updated with new fake-rate in 2016/2017 splitting
    TFile* frFileEl = TFile::Open("data/CMID/chargeMisID_" + (TString)(is2016 ? "2016" : "2017") + ".root");

    const std::string frUnc[3] = {"", "", ""};
    for(unsigned unc = 0; unc < 3; ++unc){
        CMIDMapEle[unc] = std::shared_ptr<TH2D>( (TH2D*) frFileEl->Get((const TString&) "passed" + frUnc[unc]) );
        CMIDMapEle[unc]->SetDirectory(gROOT);
    }
    frFileEl->Close();
}

void Reweighter::initializePrefiringProbabilities(){
    std::string year = ( is2016 ? "2016" : "2017" );
    TFile* prefiringFile = TFile::Open( (const TString&) "data/L1prefiring/L1prefiring_eff_jet_" + year + ".root");
    prefiringMap = std::shared_ptr<TH2D>( (TH2D*) prefiringFile->Get( is2016? "prefireEfficiencyMap" : "L1prefiring_jet_2017BtoF" ) );
    prefiringMap->SetDirectory(gROOT);
    prefiringFile->Close();
}

Reweighter::~Reweighter(){}

double Reweighter::puWeight(const double nTrueInt, const Sample& sample, const unsigned unc) const{
    if(unc < 3){

        //find weights for given sample
        const auto& weightVectorIter = puWeights.find(sample.getUniqueName() );

        //check if pu weights are available for this sample
        if( weightVectorIter == puWeights.cend() ){
            std::cerr << "Error: no pu weights found for sample : " << sample.getUniqueName() << " returning weight 0  " << std::endl;
            return 0.;
        }

        double maxBin = (sample.is2016() ?  49.5 : 99.5 );
        TH1D* weights = ( (*weightVectorIter).second)[unc].get();
        if(std::isinf(nTrueInt)) return 1.;
        return weights->GetBinContent(weights->FindBin(std::min(std::max(nTrueInt, 0.), maxBin) ) );
    }
    else {
        std::cerr << "Error: invalid pu uncertainty requested: returning weight 0" << std::endl;
        return 0.;
    }
}

double Reweighter::bTagWeight(const unsigned jetFlavor, const double jetPt, const double jetEta, const double jetCSV, const unsigned unc) const{
    static const BTagEntry::JetFlavor flavorEntries[3] = {BTagEntry::FLAV_UDSG, BTagEntry::FLAV_C, BTagEntry::FLAV_B};
    static const std::string uncName[3] = {"central", "down", "up"};
    if(unc < 3){
        return bTagCalibReader->eval_auto_bounds(uncName[unc], flavorEntries[flavorInd(jetFlavor)], jetEta, jetPt, jetCSV);
    } else {
        std::cerr << "Error: invalid b-tag SF uncertainty requested: returning weight 1" << std::endl;
        return 1.;
    }
}

double Reweighter::bTagEff(const unsigned jetFlavor, const double jetPt, const double jetEta) const{
    double croppedPt = std::min( std::max(jetPt, 25.), 599.);
    double croppedAbsEta = std::min( fabs(jetEta), 2.39 ); 
    return bTagEffHist[flavorInd(jetFlavor)]->GetBinContent(bTagEffHist[flavorInd(jetFlavor)]->FindBin(croppedPt, croppedAbsEta) );
}

double Reweighter::muonRecoWeight(const double eta, const unsigned unc) const{
    //!!!! To be split for 2016 and 2017 data !!!!
    // no reco SF, email from Sergio and Sandro on 30th of July, 2018, they removed as well recommendations from twiki
    if( is2016 ){
        return muonRecoSF->Eval(std::max(-2.4,std::min(eta, 2.4) ) );
    } else {
        //WARNING temporary, implement this later
        return 1.;
    }
}

double Reweighter::electronRecoWeight(const double superClusterEta, const double pt, const unsigned unc) const{
    double croppedSuperClusterEta = std::max(-2.49, std::min(superClusterEta, 2.49) );
    int var = unc == 2 ? -1 : int(unc);
    if( is2016 ){
        if( pt <= 20 ){
            double croppedPt = std::max(10.01, pt);
            return electronRecoSF_pT0to20->GetBinContent( electronRecoSF_pT0to20->FindBin( croppedSuperClusterEta, croppedPt ) ) + (var != 0 ? var * electronRecoSF_pT0to20->GetBinError( electronRecoSF_pT0to20->FindBin( croppedSuperClusterEta, croppedPt ) ) : 0.);
        } else {
            double croppedPt = std::min(pt, 499.); 
            return electronRecoSF_pT20toInf->GetBinContent( electronRecoSF_pT20toInf->FindBin( croppedSuperClusterEta, croppedPt ) ) + (var != 0 ? var * electronRecoSF_pT20toInf->GetBinError( electronRecoSF_pT20toInf->FindBin( croppedSuperClusterEta, croppedPt ) ) : 0.);
        }
    } else {
        if( pt <= 20 ){
            double croppedPt = std::max(10.01, pt);
            return electronRecoSF_pT0to20->GetBinContent( electronRecoSF_pT0to20->FindBin( croppedSuperClusterEta, croppedPt ) ) + (var != 0 ? var * electronRecoSF_pT0to20->GetBinError( electronRecoSF_pT0to20->FindBin( croppedSuperClusterEta, croppedPt ) ) : 0.);
        } else {
            double croppedPt = std::min(pt, 499.); 
            return electronRecoSF_pT20toInf->GetBinContent( electronRecoSF_pT20toInf->FindBin( croppedSuperClusterEta, croppedPt ) ) + (var != 0 ? var * electronRecoSF_pT20toInf->GetBinError( electronRecoSF_pT20toInf->FindBin( croppedSuperClusterEta, croppedPt ) ) : 0.);
        }
    }
}

double Reweighter::muonLooseIdWeight(const double pt, const double eta, const unsigned unc) const{
    double croppedPt = std::min(pt, 199.);
    double croppedEta = std::min( fabs(eta), 2.39 ); 
    return muonLooseToRecoSF->GetBinContent( muonLooseToRecoSF->FindBin( croppedPt, croppedEta) ) + (unc != 0 ? unc * muonLooseToRecoSF->GetBinError( muonLooseToRecoSF->FindBin( croppedPt, croppedEta) ) : 0.);
}

double Reweighter::electronLooseIdWeight(const double pt, const double eta, const unsigned unc) const{
    double croppedPt = std::min(pt, 199.);
    double croppedEta;
    if( is2016 ){   //asymmetric scale factors are currently only available for 2016 data!
        croppedEta = std::min( std::max( -2.49, eta ), 2.49 ); 
    } else {
        croppedEta = std::min( fabs(eta), 2.49); 
    }
    return electronLooseToRecoSF->GetBinContent( electronLooseToRecoSF->FindBin( croppedPt, croppedEta) ) + (unc != 0 ? unc * electronLooseToRecoSF->GetBinError( electronLooseToRecoSF->FindBin( croppedPt, croppedEta) ) : 0.);
}

// old implementation, split into tight to loose and loose to reco, switched to 1-go implementation, 29 Aug 2018
/*
double Reweighter::muonTightIdWeight(const double pt, const double eta, const unsigned unc) const{
    double croppedPt = std::min(pt, 199.);
    double croppedEta = std::min( fabs(eta), 2.39);
    double varOne = (unc != 0 ? muonLooseToRecoSF->GetBinError( muonLooseToRecoSF->FindBin( croppedPt, croppedEta) ) : 0.);
    double varTwo = (unc != 0 ? muonTightToLooseSF->GetBinError( muonTightToLooseSF->FindBin( croppedPt, croppedEta) ) : 0.);
    //if(unc == 2) std::cout << "muon with pt " << pt << " has var 1 and 2 : " << varOne << " " << varTwo << std::endl;
    double sfloose = muonLooseToRecoSF->GetBinContent( muonLooseToRecoSF->FindBin( croppedPt, croppedEta) );
    double sftight = muonTightToLooseSF->GetBinContent( muonTightToLooseSF->FindBin( croppedPt, croppedEta) );
    //if(unc == -1) std::cout << "muon with pt " << pt << " has sf(l,t) " << sfloose << " and " << sftight << "; unc is " << (unc != 0 ? unc * uncorUncCalc(sfloose, sftight, varOne, varTwo) : 0.) << std::endl;
    int var = unc == 2 ? -1 : int(unc);
    return sfloose*sftight + (var != 0 ? var * uncorUncCalc(sfloose, sftight, varOne, varTwo) : 0.);
}
*/
double Reweighter::muonTightIdWeightOnlySyst(const double pt, const double eta, const unsigned unc) const{
    double croppedPt = std::min(pt, 199.);
    double croppedEta = std::min( fabs(eta), 2.39);

    double sftight = muonTightToRecoSF_syst->GetBinContent( muonTightToRecoSF_syst->FindBin( croppedPt, croppedEta) );
    int var = unc == 2 ? -1 : int(unc);
    return sftight + (var != 0 ? var * muonTightToRecoSF_syst->GetBinError( muonTightToRecoSF_syst->FindBin( croppedPt, croppedEta) ) : 0.);
}

double Reweighter::muonTightIdWeightOnlyStat(const double pt, const double eta, const unsigned unc) const{
    double croppedPt = std::min(pt, 199.);
    double croppedEta = std::min( fabs(eta), 2.39);

    double sftight = muonTightToRecoSF_stat->GetBinContent( muonTightToRecoSF_stat->FindBin( croppedPt, croppedEta) );
    int var = unc == 2 ? -1 : int(unc);
    return sftight + (var != 0 ? var * muonTightToRecoSF_stat->GetBinError( muonTightToRecoSF_stat->FindBin( croppedPt, croppedEta) ) : 0.);
}

double Reweighter::muonChargeConsWeight(const double pt, const double eta, const unsigned unc) const{
    double croppedPt = std::min(pt, 199.);
    double croppedEta = std::min( fabs(eta), 2.39);
    double sf = muonChargeConsToTightSF->GetBinContent( muonChargeConsToTightSF->FindBin( croppedPt, croppedEta) ) + unc * muonChargeConsToTightSF->GetBinError( muonChargeConsToTightSF->FindBin( croppedPt, croppedEta) );
    //std::cout << "lepton with pt " << pt << " has sf " << sf << std::endl;
    return sf;
}

/* // old implementation, split into tight to loose and loose to reco, switched to 1-go implementation, 29 Aug 2018
double Reweighter::electronTightIdWeight(const double pt, const double eta, const unsigned unc) const{
    double croppedPt = std::min(pt, 199.);
    double croppedEta;
    if( is2016 ){   //asymmetric scale factors are currently only available for 2016 data!
        croppedEta = std::min( std::max( -2.49, eta ), 2.49 ); 
    } else {
        croppedEta = std::min( fabs(eta), 2.49); 
    }
    double varOne = unc != 0 ? electronLooseToRecoSF->GetBinError( electronLooseToRecoSF->FindBin( croppedPt,  croppedEta) ) : 0.;
    double varTwo = unc != 0 ? electronTightToLooseSF->GetBinError( electronTightToLooseSF->FindBin( croppedPt, croppedEta) ) : 0.;
    double sfloose = electronLooseToRecoSF->GetBinContent( electronLooseToRecoSF->FindBin( croppedPt, croppedEta) );
    double sftight = electronTightToLooseSF->GetBinContent( electronTightToLooseSF->FindBin( croppedPt, croppedEta) );
    int var = unc == 2 ? -1 : int(unc);
    return sfloose*sftight + (var != 0 ? var * uncorUncCalc(sfloose, sftight, varOne, varTwo) : 0.);
}
*/

double Reweighter::electronTightIdWeightOnlySyst(const double pt, const double eta, const unsigned unc) const{
    double croppedPt = std::min(pt, 199.);
    double croppedEta;
    if( is2016 ){   //asymmetric scale factors are currently only available for 2016 data!
        croppedEta = std::min( std::max( -2.49, eta ), 2.49 ); 
    } else {
        croppedEta = std::min( fabs(eta), 2.49); 
    }

    double sftight = electronTightToRecoSF_syst->GetBinContent( electronTightToRecoSF_syst->FindBin( croppedPt, croppedEta) );

    int var = unc == 2 ? -1 : int(unc);
    return sftight + (var != 0 ? var * electronTightToRecoSF_syst->GetBinError( electronTightToRecoSF_syst->FindBin( croppedPt, croppedEta)) : 0.);
}

double Reweighter::electronTightIdWeightOnlyStat(const double pt, const double eta, const unsigned unc) const{
    double croppedPt = std::min(pt, 199.);
    double croppedEta;
    if( is2016 ){   //asymmetric scale factors are currently only available for 2016 data!
        croppedEta = std::min( std::max( -2.49, eta ), 2.49 ); 
    } else {
        croppedEta = std::min( fabs(eta), 2.49); 
    }

    double sftight = electronTightToRecoSF_stat->GetBinContent( electronTightToRecoSF_stat->FindBin( croppedPt, croppedEta) );

    int var = unc == 2 ? -1 : int(unc);
    return sftight + (var != 0 ? var * electronTightToRecoSF_stat->GetBinError( electronTightToRecoSF_stat->FindBin( croppedPt, croppedEta)) : 0.);
}

double Reweighter::electronChargeConsWeight(const double pt, const double eta, const unsigned unc) const{
    double croppedPt = std::min(pt, 199.);
    double croppedEta;
    if( is2016 ){   //asymmetric scale factors are currently only available for 2016 data!
        croppedEta = std::min( std::max( -2.49, eta ), 2.49 ); 
    } else {
        croppedEta = std::min( fabs(eta), 2.49); 
    }
    double sf = electronChargeConsToTightSF->GetBinContent( electronChargeConsToTightSF->FindBin( croppedPt, croppedEta) ) + (unc != 0 ? unc * electronChargeConsToTightSF->GetBinError( electronChargeConsToTightSF->FindBin( croppedPt, croppedEta) ) : 0.);
    return sf;
}

double Reweighter::muonFakeRate(const double pt, const double eta, const unsigned unc) const{
    //!!!! To be split for 2016 and 2017 data !!!! 
    if(unc < 3){
        double etaFixed = eta;
        // empty bin in tZq FR maps
        //if(pt < 15 && fabs(eta) > 2.1) etaFixed = 2.09;
        return frMapMu[unc]->GetBinContent(frMapMu[unc]->FindBin(std::min(pt, (leptonSelection == 2 ? 64. : 44.)), std::min(fabs(etaFixed), 2.4) ) );
    } else {
        std::cerr << "Error: invalid muon fake-rate uncertainty requested: returning fake-rate 99" << std::endl;
        return 99;
    }
}

double Reweighter::electronFakeRate(const double pt, const double eta, const unsigned unc) const{
    //!!!! To be split for 2016 and 2017 data !!!!
    if(unc < 3){ 
        return frMapEle[unc]->GetBinContent(frMapEle[unc]->FindBin(std::min(pt, (leptonSelection == 2 ? 64. : 44.)), std::min(fabs(eta), 2.5) ) );
    } else {
        std::cerr << "Error: invalid electron fake-rate uncertainty requested: returning fake-rate 99" << std::endl;
        return 99;
    }
}

double Reweighter::electronCMIDRate(const double pt, const double eta, const unsigned unc) const{
    //!!!! To be split for 2016 and 2017 data !!!!
    if(unc < 3){ 
        return CMIDMapEle[unc]->GetBinContent(CMIDMapEle[unc]->FindBin(std::min(pt, 299.), std::min(fabs(eta), 2.5) ) );
    } else {
        std::cerr << "Error: invalid electron charge mis ID rate uncertainty requested: returning charge mis ID rate 99" << std::endl;
        return 99;
    }
}

double Reweighter::getTriggerSF(const double pt) const{
    if(leptonSelection == 3 && pt < 120 && is2016) // cross checked from plots in AN, figure 5, 29th Aug 2018
        //weight *= 0.985; // got it from Daniel's trigger measurement, agreed with him on skype, 5 Jul 2018
        return 0.985;
    if(leptonSelection == 3 && pt < 80 && !is2016) // 
        return 0.966000020504;
    return 1.;
}

double Reweighter::uncorUncCalc(const double & centralOne, const double & centralTwo, const double & uncOne, const double & uncTwo) const{
    // http://lectureonline.cl.msu.edu/~mmp/labs/error/e2.htm, point 5
    return sqrt(centralTwo * centralTwo * uncOne * uncOne + centralOne * centralOne * uncTwo * uncTwo);
}

double Reweighter::jetPrefiringProbability(const double pt, const double eta) const{

    //consider probabilities outside of map to be zero
    //this is implicitly implemented by not requiring the pt and eta values to fall in the map

    //abseta binning for 2016, eta binning for 2017
    double croppedEta = eta;
    if( is2016 ){
        croppedEta = fabs(eta);
    }
    return prefiringMap->GetBinContent( prefiringMap->FindBin(croppedEta, pt ) );
}
