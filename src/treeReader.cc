#include <iostream>
#include <fstream>

#include "../interface/treeReader.h"
#include "../interface/analysisTools.h"

using namespace std;

treeReader::treeReader(TTree *tree) : fChain(nullptr) 
{
    if (tree != nullptr){
        initTree(tree);
    }
}

void treeReader::readSamples(const std::string& list, std::vector<Sample>& sampleVector){
    //sampleVector.clear();    //clear current sample list
    //read sample info (names and xSec) from txt file
    std::ifstream file(list);
    int sampleCounter = 0;
    int processCounter = 0;
    do {
        sampleVector.push_back(Sample(file));
        namesOfTheFiles.push_back(sampleVector.back().getProcessName());
        if(std::find(namesOfTheProcesses.begin(), namesOfTheProcesses.end(), sampleVector.back().getProcessName()) == namesOfTheProcesses.end() && sampleVector.back().getProcessName() != "") {
            namesOfTheProcesses.push_back(sampleVector.back().getProcessName());
            processIndex.insert(std::pair<std::string,int>(sampleVector.back().getProcessName(), processCounter));
            if(sampleVector.back().getProcessName() == "nonpromptData")
                nonPromptSample = processCounter;
            processCounter++;
        } 

        //if(sampleVector.back().getProcessName() == "chargeMisIDData")
        //    CMIDSample = sampleCounter;
        sampleCounter++;
    } while(!file.eof());
    sampleVector.pop_back();
    file.close();       //close file after usage
    //display samples that have been read 
    //is2017folder = list.find("2017") != std::string::npos;
    //dataLumi = is2017 ? 41.9 : 35.9;

}

void treeReader::readSamples(const std::string& list){
    readSamples(list, this->samples);
}

void treeReader::initSample(const Sample& samp){ 

    //update current sample
    currentSample = samp;
    sampleFile = samp.getFile("/user/ikhvastu/Work/ntuples_ttV_" + std::string(samp.is2017() ? "2017/" : "2016/")); //  + (TString)(is2017 ? "" : "newReReco/") 
    //sampleFile = samp.getFile("/Users/illiakhvastunov/Desktop/CERN/MCsamples/94X/ntuples_ttV_" + std::string(samp.is2017() ? "2017/" : "2016/")); //  + (TString)(is2017 ? "" : "newReReco/") 
    sampleFile->cd("blackJackAndHookers");
    fChain = (TTree*) sampleFile->Get("blackJackAndHookers/blackJackAndHookersTree");
    initTree(fChain, samp.isData());
    nEntries = fChain->GetEntries();
    
    isData = (samples[currentSampleIndex].getProcessName()) == "data";
    isDataNonprompt = (samples[currentSampleIndex].getProcessName()) == "nonpromptData";
    isChargeMisIDSample = (samples[currentSampleIndex].getProcessName()) == "chargeMisIDData";

    if(!samp.isData()){

        //read sum of simulated event weights
        TH1D* hCounter = new TH1D("hCounter", "Events counter", 1, 0, 1);
        hCounter->Read("hCounter"); 
        sumSimulatedEventWeights = hCounter->GetBinContent(1);
        delete hCounter;

        TH1D* lheCounter = new TH1D("lheCounter", "Events counter", 110, 0, 110);
        lheCounter->Read("lheCounter"); 

        // for some samples these lhe weights are 0 (not stored in root file), simply take initial number of simulated events
        sumSimulatedEventWeightsScaleUp = lheCounter->GetBinContent(9) != 0 ? lheCounter->GetBinContent(9) : sumSimulatedEventWeights;
        sumSimulatedEventWeightsScaleDown = lheCounter->GetBinContent(5) != 0 ? lheCounter->GetBinContent(5) : sumSimulatedEventWeights;
        
        for(unsigned lhe = 9; lhe < 110; ++lhe){
            double variedSumOfWeights = lheCounter->GetBinContent(lhe + 1) != 0 ? lheCounter->GetBinContent(lhe + 1) : sumSimulatedEventWeights;
            crossSectionRatio[currentSampleIndex][lhe-9] = sumSimulatedEventWeights /( variedSumOfWeights );
        }

        delete lheCounter;
        //event weights set with lumi depending on sample's era 
        double dataLumi;
        if( is2016() ){
            dataLumi = lumi2016;
        } else {
            dataLumi = lumi2017;
        } 
        scale = samp.getXSec()*dataLumi*1000/sumSimulatedEventWeights;       //xSec*lumi divided by total sum of simulated event weights
    }

    if(std::find(samplesOrderNames.begin(), samplesOrderNames.end(), (samples[currentSampleIndex].getProcessName())) == samplesOrderNames.end()){
        std::cout << "New collection: " << (samples[currentSampleIndex].getProcessName()) << "; with number: " << currentSampleIndex << std::endl;
        samplesOrder.push_back(currentSampleIndex);
        samplesOrderNames.push_back((samples[currentSampleIndex].getProcessName()));
    }
    if(currentSampleIndex == samples.size()-1){
        samplesOrder.push_back(currentSampleIndex+1);
    }

    //++currentSampleIndex;    //increment the current sample for the next iteration
}

void treeReader::initSample(){ //initialize the next sample in the list 
    initSample(samples[++currentSampleIndex]);
}

void treeReader::GetEntry(long unsigned entry)
{
    if (!fChain) return;
    fChain->GetEntry(entry);
    //Set up correct weights
    if(!isData && !isDataNonprompt && !isChargeMisIDSample) {
        weight = _weight*scale; //MC
        /*
        if(isChargeMisIDSample && is2017)
            weight *= 1.2; // apply 20% correction to observed charge mis ID in 2017
        */
    }
    else weight = 1;                               //data
}

void treeReader::initTree(TTree *tree, const bool isData)
{
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("_runNb", &_runNb, &b__runNb);
    fChain->SetBranchAddress("_lumiBlock", &_lumiBlock, &b__lumiBlock);
    fChain->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
    fChain->SetBranchAddress("_nVertex", &_nVertex, &b__nVertex);    

    fChain->SetBranchAddress("_passTrigger_e", &_passTrigger_e, &b__passTrigger_e);
    fChain->SetBranchAddress("_passTrigger_m", &_passTrigger_m, &b__passTrigger_m);
    fChain->SetBranchAddress("_passTrigger_ee", &_passTrigger_ee, &b__passTrigger_ee);
    fChain->SetBranchAddress("_passTrigger_em", &_passTrigger_em, &b__passTrigger_em);
    fChain->SetBranchAddress("_passTrigger_mm", &_passTrigger_mm, &b__passTrigger_mm);
    fChain->SetBranchAddress("_passTrigger_eee", &_passTrigger_eee, &b__passTrigger_eee);
    fChain->SetBranchAddress("_passTrigger_eem", &_passTrigger_eem, &b__passTrigger_eem);
    fChain->SetBranchAddress("_passTrigger_emm", &_passTrigger_emm, &b__passTrigger_emm);
    fChain->SetBranchAddress("_passTrigger_mmm", &_passTrigger_mmm, &b__passTrigger_mmm);

    fChain->SetBranchAddress("_passMETFilters", &_passMETFilters, &b__passMETFilters);
    fChain->SetBranchAddress("_nL", &_nL, &b__nL);
    fChain->SetBranchAddress("_nMu", &_nMu, &b__nMu);
    fChain->SetBranchAddress("_nEle", &_nEle, &b__nEle);
    fChain->SetBranchAddress("_nLight", &_nLight, &b__nLight);
    //fChain->SetBranchAddress("_nTau", &_nTau, &b__nTau);
    fChain->SetBranchAddress("_lPt", _lPt, &b__lPt);
    fChain->SetBranchAddress("_lEta", _lEta, &b__lEta);
    fChain->SetBranchAddress("_lEtaSC", _lEtaSC, &b__lEtaSC);
    fChain->SetBranchAddress("_lPhi", _lPhi, &b__lPhi);
    fChain->SetBranchAddress("_lE", _lE, &b__lE);
    fChain->SetBranchAddress("_lFlavor", _lFlavor, &b__lFlavor);
    fChain->SetBranchAddress("_lCharge", _lCharge, &b__lCharge);
    fChain->SetBranchAddress("_dxy", _dxy, &b__dxy);
    fChain->SetBranchAddress("_dz", _dz, &b__dz);
    fChain->SetBranchAddress("_3dIP", _3dIP, &b__3dIP);
    fChain->SetBranchAddress("_3dIPSig", _3dIPSig, &b__3dIPSig);
    fChain->SetBranchAddress("_lElectronMva", _lElectronMva, &b__lElectronMva);
    fChain->SetBranchAddress("_lElectronMvaHZZ", _lElectronMvaHZZ, &b__lElectronMvaHZZ);
    fChain->SetBranchAddress("_lElectronMvaFall17NoIso", _lElectronMvaFall17NoIso, &b__lElectronMvaFall17NoIso);
    fChain->SetBranchAddress("_lElectronPassEmu", _lElectronPassEmu, &b__lElectronPassEmu);
    fChain->SetBranchAddress("_lElectronPassConvVeto", _lElectronPassConvVeto, &b__lElectronPassConvVeto);
    fChain->SetBranchAddress("_lElectronChargeConst", _lElectronChargeConst, &b__lElectronChargeConst);
    fChain->SetBranchAddress("_lElectronMissingHits", _lElectronMissingHits, &b__lElectronMissingHits);
    fChain->SetBranchAddress("_leptonMvaSUSY", _leptonMvaSUSY, &b__leptonMvaSUSY);
    fChain->SetBranchAddress("_leptonMvaTTH", _leptonMvaTTH, &b__leptonMvaTTH);
    fChain->SetBranchAddress("_leptonMvatZqTTV", _leptonMvatZqTTV, &b__leptonMvatZqTTV);

    fChain->SetBranchAddress("_lPOGLoose", _lPOGLoose, &b__lPOGLoose);
    fChain->SetBranchAddress("_lPOGMedium", _lPOGMedium, &b__lPOGMedium);
    fChain->SetBranchAddress("_lPOGTight", _lPOGTight, &b__lPOGTight);
    
    //fChain->SetBranchAddress("_tauMuonVeto", _tauMuonVeto, &b__tauMuonVeto);
    //fChain->SetBranchAddress("_tauEleVeto", _tauEleVeto, &b__tauEleVeto);
    //fChain->SetBranchAddress("_decayModeFindingNew", _decayModeFindingNew, &b__decayModeFindingNew);
    //fChain->SetBranchAddress("_tauVLooseMvaNew", _tauVLooseMvaNew, &b__tauVLooseMvaNew);
    //fChain->SetBranchAddress("_tauLooseMvaNew", _tauLooseMvaNew, &b__tauLooseMvaNew);
    //fChain->SetBranchAddress("_tauMediumMvaNew", _tauMediumMvaNew, &b__tauMediumMvaNew);
    //fChain->SetBranchAddress("_tauTightMvaNew", _tauTightMvaNew, &b__tauTightMvaNew);
    //fChain->SetBranchAddress("_tauVTightMvaNew", _tauVTightMvaNew, &b__tauVTightMvaNew);
    //fChain->SetBranchAddress("_tauVTightMvaOld", _tauVTightMvaOld, &b__tauVTightMvaOld);
    fChain->SetBranchAddress("_relIso", _relIso, &b__relIso);
    fChain->SetBranchAddress("_relIso0p4Mu", _relIso0p4Mu, &b__relIso0p4Mu);
    fChain->SetBranchAddress("_miniIso", _miniIso, &b__miniIso);
    fChain->SetBranchAddress("_miniIsoCharged", _miniIsoCharged, &b__miniIsoCharged);
    fChain->SetBranchAddress("_ptRel", _ptRel, &b__ptRel);
    fChain->SetBranchAddress("_ptRatio", _ptRatio, &b__ptRatio);
    fChain->SetBranchAddress("_closestJetCsvV2", _closestJetCsvV2, &b__closestJetCsvV2);
    fChain->SetBranchAddress("_closestJetDeepCsv_b", _closestJetDeepCsv_b, &b__closestJetDeepCsv_b);
    fChain->SetBranchAddress("_closestJetDeepCsv_bb", _closestJetDeepCsv_bb, &b__closestJetDeepCsv_bb);
    fChain->SetBranchAddress("_selectedTrackMult", _selectedTrackMult, &b__selectedTrackMult);
    fChain->SetBranchAddress("_lMuonSegComp", _lMuonSegComp, &b__lMuonSegComp);
    fChain->SetBranchAddress("_lMuonTrackPt", _lMuonTrackPt, &b__lMuonTrackPt);
    fChain->SetBranchAddress("_lMuonTrackPtErr", _lMuonTrackPtErr, &b__lMuonTrackPtErr);
    fChain->SetBranchAddress("_nJets", &_nJets, &b__nJets);
    fChain->SetBranchAddress("_jetPt", _jetPt, &b__jetPt);
    fChain->SetBranchAddress("_jetPt_JECUp", _jetPt_JECUp, &b__jetPt_JECUp);
    fChain->SetBranchAddress("_jetPt_JECDown", _jetPt_JECDown, &b__jetPt_JECDown);
    fChain->SetBranchAddress("_jetPt_JERUp", _jetPt_JERUp, &b__jetPt_JERUp);
    fChain->SetBranchAddress("_jetPt_JERDown", _jetPt_JERDown, &b__jetPt_JERDown);
    fChain->SetBranchAddress("_jetEta", _jetEta, &b__jetEta);
    fChain->SetBranchAddress("_jetPhi", _jetPhi, &b__jetPhi);
    fChain->SetBranchAddress("_jetE", _jetE, &b__jetE);
    fChain->SetBranchAddress("_jetCsvV2", _jetCsvV2, &b__jetCsvV2);
    fChain->SetBranchAddress("_jetDeepCsv_udsg", _jetDeepCsv_udsg, &b__jetDeepCsv_udsg);
    fChain->SetBranchAddress("_jetDeepCsv_b", _jetDeepCsv_b, &b__jetDeepCsv_b);
    fChain->SetBranchAddress("_jetDeepCsv_c", _jetDeepCsv_c, &b__jetDeepCsv_c);
    fChain->SetBranchAddress("_jetDeepCsv_bb", _jetDeepCsv_bb, &b__jetDeepCsv_bb);
    fChain->SetBranchAddress("_jetHadronFlavor", _jetHadronFlavor, &b__jetHadronFlavor);
    //fChain->SetBranchAddress("_jetId", _jetId, &b__jetId);
    fChain->SetBranchAddress("_jetIsTight", _jetIsTight, &b__jetIsTight);
    fChain->SetBranchAddress("_met", &_met, &b__met);
    fChain->SetBranchAddress("_metPhi", &_metPhi, &b__metPhi);
    /*
    fChain->SetBranchAddress("_metJECDown", &_metJECDown, &b__metJECDown);
    fChain->SetBranchAddress("_metJECUp", &_metJECUp, &b__metJECUp);
    fChain->SetBranchAddress("_metUnclDown", &_metUnclDown, &b__metUnclDown);
    fChain->SetBranchAddress("_metUnclUp", &_metUnclUp, &b__metUnclUp);
    fChain->SetBranchAddress("_metPhiJECDown", &_metPhiJECDown, &b__metPhiJECDown);
    fChain->SetBranchAddress("_metPhiJECUp", &_metPhiJECUp, &b__metPhiJECUp);
    fChain->SetBranchAddress("_metPhiUnclDown", &_metPhiUnclDown, &b__metPhiUnclDown);
    fChain->SetBranchAddress("_metPhiUnclUp", &_metPhiUnclUp, &b__metPhiUnclUp);
    */

    if(!isData){
        fChain->SetBranchAddress("_weight", &_weight, &b__weight);
        //fChain->SetBranchAddress("_nLheWeights", &_nLheWeights, &b__nLheWeights);
        fChain->SetBranchAddress("_lheWeight", _lheWeight, &b__lheWeight);
        fChain->SetBranchAddress("_nTrueInt", &_nTrueInt, &b__nTrueInt);
        fChain->SetBranchAddress("_gen_met", &_gen_met, &b__gen_met);
        fChain->SetBranchAddress("_gen_metPhi", &_gen_metPhi, &b__gen_metPhi);
        fChain->SetBranchAddress("_gen_nL", &_gen_nL, &b__gen_nL);
        fChain->SetBranchAddress("_gen_lPt", _gen_lPt, &b__gen_lPt);
        fChain->SetBranchAddress("_gen_lEta", _gen_lEta, &b__gen_lEta);
        fChain->SetBranchAddress("_gen_lPhi", _gen_lPhi, &b__gen_lPhi);
        fChain->SetBranchAddress("_gen_lE", _gen_lE, &b__gen_lE);
        fChain->SetBranchAddress("_gen_lFlavor", _gen_lFlavor, &b__gen_lFlavor);
        fChain->SetBranchAddress("_gen_lCharge", _gen_lCharge, &b__gen_lCharge);
        fChain->SetBranchAddress("_gen_lMomPdg", _gen_lMomPdg, &b__gen_lMomPdg);
        fChain->SetBranchAddress("_gen_lIsPrompt", _gen_lIsPrompt, &b__gen_lIsPrompt);
        fChain->SetBranchAddress("_lIsPrompt", _lIsPrompt, &b__lIsPrompt);
        fChain->SetBranchAddress("_lMatchPdgId", _lMatchPdgId, &b__lMatchPdgId);
    }

    
}

void treeReader::setOutputTree(TTree* outputTree, const bool isData){
    outputTree->Branch("_runNb",                        &_runNb,                        "_runNb/l");
    outputTree->Branch("_lumiBlock",                    &_lumiBlock,                    "_lumiBlock/l");
    outputTree->Branch("_eventNb",                      &_eventNb,                      "_eventNb/l");
    outputTree->Branch("_nVertex",                      &_nVertex,                      "_nVertex/b");
    outputTree->Branch("_met",                          &_met,                          "_met/D");
    outputTree->Branch("_metJECDown",                   &_metJECDown,                   "_metJECDown/D");
    outputTree->Branch("_metJECUp",                     &_metJECUp,                     "_metJECUp/D");
    outputTree->Branch("_metUnclDown",                  &_metUnclDown,                  "_metUnclDown/D");
    outputTree->Branch("_metUnclUp",                    &_metUnclUp,                    "_metUnclUp/D");
    outputTree->Branch("_metPhi",                       &_metPhi,                       "_metPhi/D");
    outputTree->Branch("_metPhiJECDown",                &_metPhiJECDown,                "_metPhiJECDown/D");
    outputTree->Branch("_metPhiJECUp",                  &_metPhiJECUp,                  "_metPhiJECUp/D");
    outputTree->Branch("_metPhiUnclDown",               &_metPhiUnclDown,               "_metPhiUnclDown/D");
    outputTree->Branch("_metPhiUnclUp",                 &_metPhiUnclUp,                 "_metPhiUnclUp/D");
    if(isData){                 //Temporarily only store 2017 triggers for data, to be updated when 2017 MC is available
        outputTree->Branch("_passTrigger_e", &_passTrigger_e, "_passTrigger_e/O");
        outputTree->Branch("_HLT_Ele35_WPTight_Gsf", &_HLT_Ele35_WPTight_Gsf, "_HLT_Ele35_WPTight_Gsf/O");
        outputTree->Branch("_HLT_Ele35_WPTight_Gsf_prescale", &_HLT_Ele35_WPTight_Gsf_prescale, "_HLT_Ele35_WPTight_Gsf_prescale/O");
        outputTree->Branch("_HLT_Ele40_WPTight_Gsf", &_HLT_Ele40_WPTight_Gsf, "_HLT_Ele40_WPTight_Gsf/O");
        outputTree->Branch("_HLT_Ele40_WPTight_Gsf_prescale", &_HLT_Ele40_WPTight_Gsf_prescale, "_HLT_Ele40_WPTight_Gsf_prescale/O");
        outputTree->Branch("_passTrigger_ee", &_passTrigger_ee, "_passTrigger_ee/O");
        outputTree->Branch("_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350", &_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350, "_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350/O");
        outputTree->Branch("_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_prescale", &_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_prescale, "_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_prescale/O");
        outputTree->Branch("_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, "_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL/O");
        outputTree->Branch("_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_prescale", &_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_prescale, "_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_prescale/O");
        outputTree->Branch("_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, "_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ/O");
        outputTree->Branch("_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale", &_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale, "_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale/O");
        outputTree->Branch("_passTrigger_eee", &_passTrigger_eee, "_passTrigger_eee/O");
        outputTree->Branch("_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", &_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL, "_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL/O");
        outputTree->Branch("_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale", &_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale, "_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale/O");
        outputTree->Branch("_passTrigger_em", &_passTrigger_em, "_passTrigger_em/O");
        outputTree->Branch("_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ", &_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ, "_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ/O");
        outputTree->Branch("_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_prescale", &_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_prescale, "_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_prescale/O");
        outputTree->Branch("_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, "_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ/O");
        outputTree->Branch("_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale", &_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale, "_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale/O");
        outputTree->Branch("_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, "_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ/O");
        outputTree->Branch("_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale", &_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale, "_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale/O");
        outputTree->Branch("_passTrigger_m", &_passTrigger_m, "_passTrigger_m/O");
        outputTree->Branch("_HLT_IsoMu27", &_HLT_IsoMu27, "_HLT_IsoMu27/O");
        outputTree->Branch("_HLT_IsoMu27_prescale", &_HLT_IsoMu27_prescale, "_HLT_IsoMu27_prescale/O");
        outputTree->Branch("_HLT_IsoMu30", &_HLT_IsoMu30, "_HLT_IsoMu30/O");
        outputTree->Branch("_HLT_IsoMu30_prescale", &_HLT_IsoMu30_prescale, "_HLT_IsoMu30_prescale/O");
        outputTree->Branch("_passTrigger_mee", &_passTrigger_mee, "_passTrigger_mee/O");
        outputTree->Branch("_HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &_HLT_Mu8_DiEle12_CaloIdL_TrackIdL, "_HLT_Mu8_DiEle12_CaloIdL_TrackIdL/O");
        outputTree->Branch("_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_prescale", &_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_prescale, "_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_prescale/O");
        outputTree->Branch("_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ", &_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ, "_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ/O");
        outputTree->Branch("_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_prescale", &_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_prescale, "_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_prescale/O");
        outputTree->Branch("_passTrigger_mm", &_passTrigger_mm, "_passTrigger_mm/O");
        outputTree->Branch("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, "_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL/O");
        outputTree->Branch("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_prescale", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_prescale, "_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_prescale/O");
        outputTree->Branch("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, "_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ/O");
        outputTree->Branch("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_prescale", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_prescale, "_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_prescale/O");
        outputTree->Branch("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_Mass8", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_Mass8, "_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_Mass8/O");
        outputTree->Branch("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_Mass8_prescale", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_Mass8_prescale, "_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_Mass8_prescale/O");
        outputTree->Branch("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8, "_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8/O");
        outputTree->Branch("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_prescale", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_prescale, "_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_prescale/O");
        outputTree->Branch("_HLT_DoubleMu4_Mass8_DZ_PFHT350", &_HLT_DoubleMu4_Mass8_DZ_PFHT350, "_HLT_DoubleMu4_Mass8_DZ_PFHT350/O");
        outputTree->Branch("_HLT_DoubleMu4_Mass8_DZ_PFHT350_prescale", &_HLT_DoubleMu4_Mass8_DZ_PFHT350_prescale, "_HLT_DoubleMu4_Mass8_DZ_PFHT350_prescale/O");
        outputTree->Branch("_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", &_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8, "_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8/O");
        outputTree->Branch("_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_prescale", &_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_prescale, "_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_prescale/O");
        outputTree->Branch("_passTrigger_mme", &_passTrigger_mme, "_passTrigger_mme/O");
        outputTree->Branch("_HLT_DiMu9_Ele9_CaloIdL_TrackIdL", &_HLT_DiMu9_Ele9_CaloIdL_TrackIdL, "_HLT_DiMu9_Ele9_CaloIdL_TrackIdL/O");
        outputTree->Branch("_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale", &_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale, "_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale/O");
        outputTree->Branch("_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ", &_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ, "_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ/O");
        outputTree->Branch("_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_prescale", &_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_prescale, "_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_prescale/O");
        outputTree->Branch("_passTrigger_mmm", &_passTrigger_mmm, "_passTrigger_mmm/O");
        outputTree->Branch("_HLT_TripleMu_10_5_5_DZ", &_HLT_TripleMu_10_5_5_DZ, "_HLT_TripleMu_10_5_5_DZ/O");
        outputTree->Branch("_HLT_TripleMu_10_5_5_DZ_prescale", &_HLT_TripleMu_10_5_5_DZ_prescale, "_HLT_TripleMu_10_5_5_DZ_prescale/O");
        outputTree->Branch("_HLT_TripleMu_5_3_3_Mass3p8to60_DZ", &_HLT_TripleMu_5_3_3_Mass3p8to60_DZ, "_HLT_TripleMu_5_3_3_Mass3p8to60_DZ/O");
        outputTree->Branch("_HLT_TripleMu_5_3_3_Mass3p8to60_DZ_prescale", &_HLT_TripleMu_5_3_3_Mass3p8to60_DZ_prescale, "_HLT_TripleMu_5_3_3_Mass3p8to60_DZ_prescale/O");
        outputTree->Branch("_TripleMu_12_10_5", &_TripleMu_12_10_5, "_TripleMu_12_10_5/O");
    }
    outputTree->Branch("_passMETFilters", &_passMETFilters, "_passMETFilters/O");
    outputTree->Branch("_nL",                           &_nL,                           "_nL/b");
    outputTree->Branch("_nMu",                          &_nMu,                          "_nMu/b");
    outputTree->Branch("_nEle",                         &_nEle,                         "_nEle/b");
    outputTree->Branch("_nLight",                       &_nLight,                       "_nLight/b");
    outputTree->Branch("_nTau",                         &_nTau,                         "_nTau/b");
    outputTree->Branch("_lPt",                          &_lPt,                          "_lPt[_nL]/D");
    outputTree->Branch("_lEta",                         &_lEta,                         "_lEta[_nL]/D");
    outputTree->Branch("_lEtaSC",                       &_lEtaSC,                       "_lEtaSC[_nLight]/D");
    outputTree->Branch("_lPhi",                         &_lPhi,                         "_lPhi[_nL]/D");
    outputTree->Branch("_lE",                           &_lE,                           "_lE[_nL]/D");
    outputTree->Branch("_lFlavor",                      &_lFlavor,                      "_lFlavor[_nL]/i");
    outputTree->Branch("_lCharge",                      &_lCharge,                      "_lCharge[_nL]/I");
    outputTree->Branch("_dxy",                          &_dxy,                          "_dxy[_nL]/D");
    outputTree->Branch("_dz",                           &_dz,                           "_dz[_nL]/D");
    outputTree->Branch("_3dIP",                         &_3dIP,                         "_3dIP[_nL]/D");
    outputTree->Branch("_3dIPSig",                      &_3dIPSig,                      "_3dIPSig[_nL]/D");
    outputTree->Branch("_lElectronMva",                 &_lElectronMva,                 "_lElectronMva[_nLight]/F");
    outputTree->Branch("_lElectronPassEmu",             &_lElectronPassEmu,             "_lElectronPassEmu[_nLight]/O");
    outputTree->Branch("_leptonMvaSUSY",                &_leptonMvaSUSY,                "_leptonMvaSUSY[_nLight]/D");
    outputTree->Branch("_leptonMvaTTH",                 &_leptonMvaTTH,                 "_leptonMvaTTH[_nLight]/D");
    outputTree->Branch("_lHNLoose",                     &_lHNLoose,                     "_lHNLoose[_nLight]/O");
    outputTree->Branch("_lHNFO",                        &_lHNFO,                        "_lHNFO[_nLight]/O");
    outputTree->Branch("_lHNTight",                     &_lHNTight,                     "_lHNTight[_nLight]/O");
    outputTree->Branch("_lEwkLoose",                    &_lEwkLoose,                    "_lEwkLoose[_nL]/O");
    outputTree->Branch("_lEwkFO",                       &_lEwkFO,                       "_lEwkFO[_nL]/O");
    outputTree->Branch("_lEwkTight",                    &_lEwkTight,                    "_lEwkTight[_nL]/O");
    outputTree->Branch("_lPOGVeto",                     &_lPOGVeto,                     "_lPOGVeto[_nL]/O");
    outputTree->Branch("_lPOGLoose",                    &_lPOGLoose,                    "_lPOGLoose[_nL]/O");
    outputTree->Branch("_lPOGMedium",                   &_lPOGMedium,                   "_lPOGMedium[_nL]/O");
    outputTree->Branch("_lPOGTight",                    &_lPOGTight,                    "_lPOGTight[_nL]/O");
    outputTree->Branch("_tauMuonVeto",                  &_tauMuonVeto,                  "_tauMuonVeto[_nL]/O");
    outputTree->Branch("_tauEleVeto",                   &_tauEleVeto,                   "_tauEleVeto[_nL]/O");
    outputTree->Branch("_decayModeFindingNew",          &_decayModeFindingNew,          "_decayModeFindingNew[_nL]/O");
    outputTree->Branch("_tauVLooseMvaNew",              &_tauVLooseMvaNew,              "_tauVLooseMvaNew[_nL]/O");
    outputTree->Branch("_tauLooseMvaNew",               &_tauLooseMvaNew,               "_tauLooseMvaNew[_nL]/O");
    outputTree->Branch("_tauMediumMvaNew",              &_tauMediumMvaNew,              "_tauMediumMvaNew[_nL]/O");
    outputTree->Branch("_tauTightMvaNew",               &_tauTightMvaNew,               "_tauTightMvaNew[_nL]/O");
    outputTree->Branch("_tauVTightMvaNew",              &_tauVTightMvaNew,              "_tauVTightMvaNew[_nL]/O");
    outputTree->Branch("_tauVTightMvaOld",              &_tauVTightMvaOld,              "_tauVTightMvaOld[_nL]/O");
    outputTree->Branch("_relIso",                       &_relIso,                       "_relIso[_nLight]/D");
    outputTree->Branch("_miniIso",                      &_miniIso,                      "_miniIso[_nLight]/D");
    outputTree->Branch("_miniIsoCharged",               &_miniIsoCharged,               "_miniIsoCharged[_nLight]/D");
    outputTree->Branch("_ptRel",                        &_ptRel,                        "_ptRel[_nLight]/D");
    outputTree->Branch("_ptRatio",                      &_ptRatio,                      "_ptRatio[_nLight]/D");
    outputTree->Branch("_closestJetCsvV2",              &_closestJetCsvV2,              "_closestJetCsvV2[_nLight]/D");
    outputTree->Branch("_closestJetDeepCsv_b",          &_closestJetDeepCsv_b,          "_closestJetDeepCsv_b[_nLight]/D");
    outputTree->Branch("_closestJetDeepCsv_bb",         &_closestJetDeepCsv_bb,         "_closestJetDeepCsv_bb[_nLight]/D");
    outputTree->Branch("_selectedTrackMult",            &_selectedTrackMult,            "_selectedTrackMult[_nLight]/i");
    outputTree->Branch("_lMuonSegComp",                  &_lMuonSegComp,                  "_lMuonSegComp[_nMu]/D");
    outputTree->Branch("_nJets",                     &_nJets,                    "_nJets/b");
    outputTree->Branch("_jetPt",                     &_jetPt,                    "_jetPt[_nJets]/D");
    outputTree->Branch("_jetPt_JECUp",               &_jetPt_JECUp,              "_jetPt_JECUp[_nJets]/D");
    outputTree->Branch("_jetPt_JECDown",             &_jetPt_JECDown,            "_jetPt_JECDown[_nJets]/D");
    outputTree->Branch("_jetPt_JERUp",               &_jetPt_JERUp,              "_jetPt_JERUp[_nJets]/D");
    outputTree->Branch("_jetPt_JERDown",             &_jetPt_JERDown,            "_jetPt_JERDown[_nJets]/D");
    outputTree->Branch("_jetEta",                    &_jetEta,                   "_jetEta[_nJets]/D");
    outputTree->Branch("_jetPhi",                    &_jetPhi,                   "_jetPhi[_nJets]/D");
    outputTree->Branch("_jetE",                      &_jetE,                     "_jetE[_nJets]/D");
    outputTree->Branch("_jetCsvV2",                  &_jetCsvV2,                 "_jetCsvV2[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_udsg",           &_jetDeepCsv_udsg,          "_jetDeepCsv_udsg[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_b",              &_jetDeepCsv_b,             "_jetDeepCsv_b[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_c",              &_jetDeepCsv_c,             "_jetDeepCsv_c[_nJets]/D");
    outputTree->Branch("_jetDeepCsv_bb",             &_jetDeepCsv_bb,            "_jetDeepCsv_bb[_nJets]/D");
    //outputTree->Branch("_jetDeepCsv_cc",             &_jetDeepCsv_cc,            "_jetDeepCsv_cc[_nJets]/D");
    outputTree->Branch("_jetHadronFlavor",           &_jetHadronFlavor,          "_jetHadronFlavor[_nJets]/i");
    //outputTree->Branch("_jetId",                     &_jetId,                    "_jetId[_nJets]/i");

    outputTree->Branch("_HLT_Ele27_WPTight_Gsf",        &_HLT_Ele27_WPTight_Gsf,        "_HLT_Ele27_WPTight_Gsf/O");
    outputTree->Branch("_HLT_IsoMu24",        &_HLT_IsoMu24,        "_HLT_IsoMu24/O");
    outputTree->Branch("_HLT_IsoTkMu24",        &_HLT_IsoTkMu24,        "_HLT_IsoTkMu24/O");
    


    if(!isData){
        outputTree->Branch("_nLheWeights",               &_nLheWeights,               "_nLheWeights/b");
        outputTree->Branch("_lheWeight",                 &_lheWeight,                 "_lheWeight[_nLheWeights]/D");
        outputTree->Branch("_weight",                    &_weight,                    "_weight/D");
        outputTree->Branch("_lIsPrompt",                 &_lIsPrompt,                 "_lIsPrompt[_nL]/O");
        outputTree->Branch("_lMatchPdgId",               &_lMatchPdgId,               "_lMatchPdgId[_nL]/I");
        outputTree->Branch("_nTrueInt",                  &_nTrueInt,                  "_nTrueInt/F");
        outputTree->Branch("_gen_met",                   &_gen_met,                   "_gen_met/D");
        outputTree->Branch("_gen_metPhi",                &_gen_metPhi,                "_gen_metPhi/D");
        outputTree->Branch("_gen_nL",                    &_gen_nL,                    "_gen_nL/b");
        outputTree->Branch("_gen_lPt",                   &_gen_lPt,                   "_gen_lPt[_gen_nL]/D");
        outputTree->Branch("_gen_lEta",                  &_gen_lEta,                  "_gen_lEta[_gen_nL]/D");
        outputTree->Branch("_gen_lPhi",                  &_gen_lPhi,                  "_gen_lPhi[_gen_nL]/D");
        outputTree->Branch("_gen_lE",                    &_gen_lE,                    "_gen_lE[_gen_nL]/D");
        outputTree->Branch("_gen_lFlavor",               &_gen_lFlavor,               "_gen_lFlavor[_gen_nL]/i");
        outputTree->Branch("_gen_lCharge",               &_gen_lCharge,               "_gen_lCharge[_gen_nL]/I");
        outputTree->Branch("_gen_lMomPdg",               &_gen_lMomPdg,               "_gen_lMomPdg[_gen_nL]/I");
        outputTree->Branch("_gen_lIsPrompt",             &_gen_lIsPrompt,             "_gen_lIsPrompt[_gen_nL]/O");

    }
}

