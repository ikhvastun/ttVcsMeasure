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
        sfFileName = "DeepCSV_94XSF_V2_B_F.csv";
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
        effFileName = "bTagEff_deepCSV_cleaned_ttZ3l_2016.root";
    } else {
        effFileName = "bTagEff_deepCSV_cleaned_ttZ3l_2017.root";
    }

    TFile* bTagFile = TFile::Open( (const TString&) "data/btagSF/" + effFileName);
    const std::string quarkFlavors[3] = {"udsg", "charm", "beauty"};
    for(unsigned flav = 0; flav < 3; ++flav){
        
        //WARNING: b-tagging efficiencies currently assume medium deepCSV tagger!
        bTagEffHist[flav] = std::shared_ptr<TH1D>( (TH1D*) bTagFile->Get( (const TString&) "bTagEff_medium" + quarkFlavors[flav]) );
        bTagEffHist[flav]->SetDirectory(gROOT);
    }
    bTagFile->Close();
}

void Reweighter::initializeElectronWeights(){   

    //read electron reco SF weights
    if( is2016 ){
        TFile* electronRecoFile = TFile::Open("data/leptonSF/egammaEffi.txt_EGM2D.root");
        electronRecoSF = std::shared_ptr<TH2D>( (TH2D*) electronRecoFile->Get("EGamma_SF2D") );
        electronRecoSF->SetDirectory(gROOT);
        electronRecoFile->Close();
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
        electronIdFile = TFile::Open("data/leptonSF/scaleFactors_electrons_2016.root");
    } else {
        electronIdFile = TFile::Open("data/leptonSF/scaleFactors_electrons_2017.root");
    } 
    electronLooseToRecoSF = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("EleToTTVLoose") );
    electronLooseToRecoSF->SetDirectory(gROOT);
    electronTightToLooseSF = std::shared_ptr<TH2D>( (TH2D*) electronIdFile->Get("TTVLooseToTTVLeptonMvattZ3l") );
    electronTightToLooseSF->SetDirectory(gROOT);
    electronIdFile->Close();
}


void Reweighter::initializeMuonWeights(){

    //read muon reco SF weights
    if(is2016){

        //WARNING: not clear how ownership works for TGraph, can not set directory
        //make sure the TGraph is not DELETED when file is closed!
        TFile* muonRecoFile = TFile::Open("data/leptonSF/Tracking_EfficienciesAndSF_BCDEFGH.root");
        muonRecoSF = std::shared_ptr<TGraph>( (TGraph*) muonRecoFile->Get("ratio_eff_eta3_dr030e030_corr") );
        muonRecoFile->Close();
    } else {

    }   

    //read muon ID SF weights
    TFile* muonIdFile; //  = TFile::Open("data/leptonSF/scaleFactors_electrons_2016.root");
    if(is2016){
        muonIdFile = TFile::Open("data/leptonSF/scaleFactors_muons_2016.root");
    } else {
        muonIdFile = TFile::Open("data/leptonSF/scaleFactors_muons_2017.root");
    }
    
    muonLooseToRecoSF = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("MuonToTTVLoose") );
    muonLooseToRecoSF->SetDirectory(gROOT);
    muonTightToLooseSF = std::shared_ptr<TH2D>( (TH2D*) muonIdFile->Get("TTVLooseToTTVLeptonMvattZ3l") );
    muonTightToLooseSF->SetDirectory(gROOT);
    muonIdFile->Close();    
}

void Reweighter::initializeFakeRate(){

    //WARNING : To be updated with new fake-rate in 2016/2017 splitting
    TFile* frFileEl = TFile::Open("data/FRmaps/elFR_data_" + TString(leptonSelection == 2 ? "2L_" : "") + (is2016 ? "2016" : "2017") + (leptonSelection == 2 ? "" : "_test") + ".root");
    TFile* frFileMu = TFile::Open("data/FRmaps/muFR_data_" + TString(leptonSelection == 2 ? "2L_" : "") + (is2016 ? "2016" : "2017") + (leptonSelection == 2 ? "" : "_test") + ".root");
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
    if( is2016 ){
        return muonRecoSF->Eval(std::max(-2.4,std::min(eta, 2.4) ) );
    } else {
        //WARNING temporary, implement this later
        return 1.;
    }
}

double Reweighter::electronRecoWeight(const double superClusterEta, const double pt, const unsigned unc) const{
    double croppedSuperClusterEta = std::max(-2.49, std::min(superClusterEta, 2.49) );
    if( is2016 ){
        double croppedPt = std::max(40., std::min(pt, 499.) );
        return electronRecoSF->GetBinContent( electronRecoSF->FindBin( croppedSuperClusterEta , croppedPt ) ) + var.at(unc) * electronRecoSF->GetBinError( electronRecoSF->FindBin( croppedSuperClusterEta , croppedPt ) );
    } else {
        if( pt <= 20 ){
            double croppedPt = std::max(10.01, pt);
            return electronRecoSF_pT0to20->GetBinContent( electronRecoSF_pT0to20->FindBin( croppedSuperClusterEta, croppedPt ) ) + var.at(unc) * electronRecoSF_pT0to20->GetBinError( electronRecoSF_pT0to20->FindBin( croppedSuperClusterEta, croppedPt ) );
        } else {
            double croppedPt = std::min(pt, 499.); 
            return electronRecoSF_pT20toInf->GetBinContent( electronRecoSF_pT20toInf->FindBin( croppedSuperClusterEta, croppedPt ) ) + var.at(unc) * electronRecoSF_pT20toInf->GetBinError( electronRecoSF_pT20toInf->FindBin( croppedSuperClusterEta, croppedPt ) );
        }
    }
}

double Reweighter::muonLooseIdWeight(const double pt, const double eta, const unsigned unc) const{
    double croppedPt = std::min(pt, 199.);
    double croppedEta = std::min( fabs(eta), 2.39 ); 
    return muonLooseToRecoSF->GetBinContent( muonLooseToRecoSF->FindBin( croppedPt, croppedEta) ) + var.at(unc) * muonLooseToRecoSF->GetBinError( muonLooseToRecoSF->FindBin( croppedPt, croppedEta) );
}

double Reweighter::electronLooseIdWeight(const double pt, const double eta, const unsigned unc) const{
    double croppedPt = std::min(pt, 199.);
    double croppedEta;
    if( is2016 ){   //asymmetric scale factors are currently only available for 2016 data!
        croppedEta = std::min( std::max( -2.49, eta ), 2.49 ); 
    } else {
        croppedEta = std::min( fabs(eta), 2.49); 
    }
    return electronLooseToRecoSF->GetBinContent( electronLooseToRecoSF->FindBin( croppedPt, croppedEta) ) + var.at(unc) * electronLooseToRecoSF->GetBinError( electronLooseToRecoSF->FindBin( croppedPt, croppedEta) );
}

double Reweighter::muonTightIdWeight(const double pt, const double eta, const unsigned unc) const{
    double croppedPt = std::min(pt, 199.);
    double croppedEta = std::min( fabs(eta), 2.39);
    double sf = muonLooseToRecoSF->GetBinContent( muonLooseToRecoSF->FindBin( croppedPt, croppedEta) );
    sf*= muonTightToLooseSF->GetBinContent( muonTightToLooseSF->FindBin( croppedPt, croppedEta) ) + var.at(unc) * muonTightToLooseSF->GetBinError( muonTightToLooseSF->FindBin( croppedPt, croppedEta) );
    //std::cout << "lepton with pt " << pt << " has sf " << sf << std::endl;
    return sf;
}

double Reweighter::electronTightIdWeight(const double pt, const double eta, const unsigned unc) const{
    double croppedPt = std::min(pt, 199.);
    double croppedEta;
    if( is2016 ){   //asymmetric scale factors are currently only available for 2016 data!
        croppedEta = std::min( std::max( -2.49, eta ), 2.49 ); 
    } else {
        croppedEta = std::min( fabs(eta), 2.49); 
    }
    double sf = electronLooseToRecoSF->GetBinContent( electronLooseToRecoSF->FindBin( croppedPt, croppedEta) );
    sf*= electronTightToLooseSF->GetBinContent( electronTightToLooseSF->FindBin( croppedPt, croppedEta) ) + var.at(unc) * electronTightToLooseSF->GetBinError( electronTightToLooseSF->FindBin( croppedPt, croppedEta) );
    return sf;
}

double Reweighter::muonFakeRate(const double pt, const double eta, const unsigned unc) const{
    //!!!! To be split for 2016 and 2017 data !!!! 
    if(unc < 3){
        return frMapMu[unc]->GetBinContent(frMapMu[unc]->FindBin(std::min(pt, 99.), std::min(fabs(eta), 2.4) ) );
    } else {
        std::cerr << "Error: invalid muon fake-rate uncertainty requested: returning fake-rate 99" << std::endl;
        return 99;
    }
}

double Reweighter::electronFakeRate(const double pt, const double eta, const unsigned unc) const{
    //!!!! To be split for 2016 and 2017 data !!!!
    if(unc < 3){ 
        return frMapEle[unc]->GetBinContent(frMapEle[unc]->FindBin(std::min(pt, 99.), std::min(fabs(eta), 2.5) ) );
    } else {
        std::cerr << "Error: invalid electron fake-rate uncertainty requested: returning fake-rate 99" << std::endl;
        return 99;
    }
}
