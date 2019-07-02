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
	 // also creates map between process and some index/counter, which is used e.g. to identify process
    std::ifstream file(list);
    int sampleCounter = 0;
    int processCounter = 0;
    do {
        sampleVector.push_back(Sample(file));
        namesOfTheFiles.push_back(sampleVector.back().getProcessName());
        if(std::find(namesOfTheProcesses.begin(), namesOfTheProcesses.end(), sampleVector.back().getProcessName()) == namesOfTheProcesses.end() && sampleVector.back().getProcessName() != "") {
            namesOfTheProcesses.push_back(sampleVector.back().getProcessName());
            processToCounterMap.insert(std::pair<std::string,int>(sampleVector.back().getProcessName(), processCounter));
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

void treeReader::initSample(const Sample& samp, std::string option){ 

	 ///////////////////////////////////////////////////////////////////
    // update current sample
	 // hard coded path with input root files. !!!! IMPORTANT !!!!
	 ///////////////////////////////////////////////////////////////////
    currentSample = samp;
    if(option == "ttZ")
        sampleFile = samp.getFile("/user/ikhvastu/Work/ntuples_ttV_" + std::string(samp.is2017() ? "2017/" : "2016/")); //  + (TString)(is2017 ? "" : "newReReco/") 
    else if(option == "ttZ4l")
//        sampleFile = samp.getFile("/user/ikhvastu/Work/ntuples_ttz_4l"); //" + istd::string(samp.is2017() ? "2017/" : "2016/"));
        sampleFile = samp.getFile("/eos/user/m/mniedzie/ttZ4l/ntuples_ttV_2016/"); //" + std::string(samp.is2017() ? "2017/" : "2016/"));
    else if(option == "FR")
        sampleFile = samp.getFile("/user/ikhvastu/Work/ntuples_FR/" + std::string(samp.getFileName().find("TT") != std::string::npos ? "ttbar" : "QCD") + "/" + std::string(samp.is2017() ?    "2017" : "2016") + "MC/");
    //sampleFile = samp.getFile("/pnfs/iihe/cms/store/user/wverbeke/ntuples_ewkino/"); //  + (TString)(is2017 ? "" : "newReReco/") 
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
        
        TH1D* psCounter = new TH1D("psCounter", "Events counter", 14, 0, 14);
        psCounter->Read("psCounter"); 

        sumSimulatedEventWeightsISRScaleUp = psCounter->GetBinContent(8) != 0 ? psCounter->GetBinContent(8) : sumSimulatedEventWeights;
        sumSimulatedEventWeightsISRScaleDown = psCounter->GetBinContent(6) != 0 ? psCounter->GetBinContent(6) : sumSimulatedEventWeights;
        
        sumSimulatedEventWeightsFSRScaleUp = psCounter->GetBinContent(5) != 0 ? psCounter->GetBinContent(5) : sumSimulatedEventWeights;
        sumSimulatedEventWeightsFSRScaleDown = psCounter->GetBinContent(3) != 0 ? psCounter->GetBinContent(3) : sumSimulatedEventWeights;
        
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

void treeReader::initSample(std::string option){ //initialize the next sample in the list 
    initSample(samples[++currentSampleIndex], option);
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
    //fChain->SetBranchAddress("_leptonMvatZqTTV16", _leptonMvatZqTTV, &b__leptonMvatZqTTV);

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
    fChain->SetBranchAddress("_jetSmearedPt", _jetSmearedPt, &b__jetSmearedPt);
    fChain->SetBranchAddress("_jetSmearedPt_JECDown", _jetSmearedPt_JECDown, &b__jetSmearedPt_JECDown);
    fChain->SetBranchAddress("_jetSmearedPt_JECUp", _jetSmearedPt_JECUp, &b__jetSmearedPt_JECUp);
    fChain->SetBranchAddress("_jetSmearedPt_JERDown", _jetSmearedPt_JERDown, &b__jetSmearedPt_JERDown);
    fChain->SetBranchAddress("_jetSmearedPt_JERUp", _jetSmearedPt_JERUp, &b__jetSmearedPt_JERUp);
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
        fChain->SetBranchAddress("_nLheWeights", &_nLheWeights, &b__nLheWeights);
        fChain->SetBranchAddress("_lheWeight", _lheWeight, &b__lheWeight);
        fChain->SetBranchAddress("_nPsWeights", &_nPsWeights, &b__nPsWeights);
        fChain->SetBranchAddress("_psWeight", _psWeight, &b__psWeight);
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
        fChain->SetBranchAddress("_lProvenance", _lProvenance, &b__lProvenance);
        fChain->SetBranchAddress("_lProvenanceCompressed", _lProvenanceCompressed, &b__lProvenanceCompressed);
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
    //outputTree->Branch("_metJetResDown",                &_metJetResDown,                "_metJetResDown/D");
    //outputTree->Branch("_metJetResUp",                  &_metJetResUp,                  "_metJetResUp/D");
    outputTree->Branch("_metUnclDown",                  &_metUnclDown,                  "_metUnclDown/D");
    outputTree->Branch("_metUnclUp",                    &_metUnclUp,                    "_metUnclUp/D");

    outputTree->Branch("_metPhi",                       &_metPhi,                       "_metPhi/D");

    outputTree->Branch("_metPhiJECDown",                &_metPhiJECDown,                "_metPhiJECDown/D");
    outputTree->Branch("_metPhiJECUp",                  &_metPhiJECUp,                  "_metPhiJECUp/D");
    //outputTree->Branch("_metPhiJetResDown",             &_metPhiJetResDown,             "_metPhiJetResDown/D");
    //outputTree->Branch("_metPhiJetResUp",               &_metPhiJetResUp,               "_metPhiJetResUp/D");
    outputTree->Branch("_metPhiUnclDown",               &_metPhiUnclDown,               "_metPhiUnclDown/D");
    outputTree->Branch("_metPhiUnclUp",                 &_metPhiUnclUp,                 "_metPhiUnclUp/D");

    outputTree->Branch("_passTrigger_e", &_passTrigger_e, "_passTrigger_e/O");
    outputTree->Branch("_passTrigger_m", &_passTrigger_m, "_passTrigger_m/O");
    outputTree->Branch("_passTrigger_ee", &_passTrigger_ee, "_passTrigger_ee/O");
    outputTree->Branch("_passTrigger_em", &_passTrigger_em, "_passTrigger_em/O");
    outputTree->Branch("_passTrigger_mm", &_passTrigger_mm, "_passTrigger_mm/O");
    outputTree->Branch("_passTrigger_eee", &_passTrigger_eee, "_passTrigger_eee/O");
    outputTree->Branch("_passTrigger_eem", &_passTrigger_eem, "_passTrigger_eem/O");
    outputTree->Branch("_passTrigger_emm", &_passTrigger_emm, "_passTrigger_emm/O");
    outputTree->Branch("_passTrigger_mmm", &_passTrigger_mmm, "_passTrigger_mmm/O");

    outputTree->Branch("_passMETFilters",               &_passMETFilters,               "_passMETFilters/O");

    outputTree->Branch("_nL",                           &_nL,                           "_nL/b");
    outputTree->Branch("_nMu",                          &_nMu,                          "_nMu/b");
    outputTree->Branch("_nEle",                         &_nEle,                         "_nEle/b");
    outputTree->Branch("_nLight",                       &_nLight,                       "_nLight/b");

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
    outputTree->Branch("_lElectronMvaHZZ",              &_lElectronMvaHZZ,              "_lElectronMvaHZZ[_nLight]/F");
    //outputTree->Branch("_lElectronMvaFall17Iso",        &_lElectronMvaFall17Iso,        "_lElectronMvaFall17Iso[_nLight]/F");
    outputTree->Branch("_lElectronMvaFall17NoIso",      &_lElectronMvaFall17NoIso,      "_lElectronMvaFall17NoIso[_nLight]/F");
    outputTree->Branch("_lElectronPassEmu",             &_lElectronPassEmu,             "_lElectronPassEmu[_nLight]/O");
    outputTree->Branch("_lElectronPassConvVeto",        &_lElectronPassConvVeto,        "_lElectronPassConvVeto[_nLight]/O");
    outputTree->Branch("_lElectronChargeConst",         &_lElectronChargeConst,         "_lElectronChargeConst[_nLight]/O");
    outputTree->Branch("_lElectronMissingHits",         &_lElectronMissingHits,         "_lElectronMissingHits[_nLight]/i");
    outputTree->Branch("_leptonMvaSUSY",                &_leptonMvaSUSY,                "_leptonMvaSUSY[_nLight]/D");
    outputTree->Branch("_leptonMvaTTH",                 &_leptonMvaTTH,                 "_leptonMvaTTH[_nLight]/D");
    outputTree->Branch("_leptonMvatZqTTV",              &_leptonMvatZqTTV,              "_leptonMvatZqTTV[_nLight]/D");
 
     // the only string needed for ZllMET
    outputTree->Branch("_lPOGTight",                    &_lPOGTight,                    "_lPOGTight[_nL]/O");

    outputTree->Branch("_lPOGLoose",                    &_lPOGLoose,                    "_lPOGLoose[_nL]/O");
    outputTree->Branch("_lPOGMedium",                   &_lPOGMedium,                   "_lPOGMedium[_nL]/O");
    //outputTree->Branch("_lPOGLooseWOIso",               &_lPOGLooseWOIso,               "_lPOGLooseWOIso[_nL]/O");
    //outputTree->Branch("_lPOGMediumWOIso",              &_lPOGMediumWOIso,              "_lPOGMediumWOIso[_nL]/O");
    //outputTree->Branch("_lPOGTightWOIso",               &_lPOGTightWOIso,               "_lPOGTightWOIso[_nL]/O");

    outputTree->Branch("_relIso",                       &_relIso,                       "_relIso[_nLight]/D");
    outputTree->Branch("_relIso0p4Mu",                  &_relIso0p4Mu,                  "_relIso0p4Mu[_nMu]/D");

    //outputTree->Branch("_relIso0p4",                  &_relIso0p4,                  "_relIso0p4[_nLight]/D");
    //outputTree->Branch("_relIso0p6",                  &_relIso0p6,                  "_relIso0p6[_nLight]/D");
    //outputTree->Branch("_relIso0p8",                  &_relIso0p8,                  "_relIso0p8[_nLight]/D");
    //outputTree->Branch("_relIso1p0",                  &_relIso1p0,                  "_relIso1p0[_nLight]/D");
    outputTree->Branch("_miniIso",                      &_miniIso,                      "_miniIso[_nLight]/D");
    outputTree->Branch("_miniIsoCharged",               &_miniIsoCharged,               "_miniIsoCharged[_nLight]/D");
    outputTree->Branch("_ptRel",                        &_ptRel,                        "_ptRel[_nLight]/D");
    outputTree->Branch("_ptRatio",                      &_ptRatio,                      "_ptRatio[_nLight]/D");
    outputTree->Branch("_closestJetCsvV2",              &_closestJetCsvV2,              "_closestJetCsvV2[_nLight]/D");
    outputTree->Branch("_closestJetDeepCsv_b",          &_closestJetDeepCsv_b,          "_closestJetDeepCsv_b[_nLight]/D");
    outputTree->Branch("_closestJetDeepCsv_bb",         &_closestJetDeepCsv_bb,         "_closestJetDeepCsv_bb[_nLight]/D");
    outputTree->Branch("_selectedTrackMult",            &_selectedTrackMult,            "_selectedTrackMult[_nLight]/i");
    outputTree->Branch("_lMuonSegComp",                 &_lMuonSegComp,                 "_lMuonSegComp[_nMu]/D");
    outputTree->Branch("_lMuonTrackPt",                 &_lMuonTrackPt,                 "_lMuonTrackPt[_nMu]/D");
    outputTree->Branch("_lMuonTrackPtErr",              &_lMuonTrackPtErr,              "_lMuonTrackPtErr[_nMu]/D");

    outputTree->Branch("_nJets",                     &_nJets,                    "_nJets/b");
    outputTree->Branch("_jetPt",                     &_jetPt,                    "_jetPt[_nJets]/D");
    outputTree->Branch("_jetSmearedPt",              &_jetSmearedPt,             "_jetSmearedPt[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JECDown",      &_jetSmearedPt_JECDown,     "_jetSmearedPt_JECDown[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JECUp",        &_jetSmearedPt_JECUp,       "_jetSmearedPt_JECUp[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JERDown",      &_jetSmearedPt_JERDown,     "_jetSmearedPt_JERDown[_nJets]/D");
    outputTree->Branch("_jetSmearedPt_JERUp",        &_jetSmearedPt_JERUp,       "_jetSmearedPt_JERUp[_nJets]/D");

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
    outputTree->Branch("_jetIsTight",                &_jetIsTight,               "_jetIsTight[_nJets]/O");
    //outputTree->Branch("_jetIsTightLepVeto",         &_jetIsTightLepVeto,        "_jetIsTightLepVeto[_nJets]/O");

    if(!isData){
        outputTree->Branch("_nLheWeights",               &_nLheWeights,               "_nLheWeights/b");
        outputTree->Branch("_lheWeight",                 &_lheWeight,                 "_lheWeight[_nLheWeights]/D");
        outputTree->Branch("_nPsWeights",                &_nPsWeights,                "_nPsWeights/b");
        outputTree->Branch("_psWeight",                  &_psWeight,                  "_psWeight[_nPsWeights]/D");

        outputTree->Branch("_lIsPrompt",                 &_lIsPrompt,                 "_lIsPrompt[_nL]/O");
        outputTree->Branch("_lMatchPdgId",               &_lMatchPdgId,               "_lMatchPdgId[_nL]/I");
        //outputTree->Branch("_lMomPdgId",                 &_lMomPdgId,                 "_lMomPdgId[_nL]/I");

        // only these 2 lines are needed for ZMET
        outputTree->Branch("_weight",                    &_weight,                    "_weight/D");
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
        //outputTree->Branch("_gen_partonPt",              &_gen_partonPt,              "_gen_partonPt[_gen_nL]/D");

        //outputTree->Branch("_lProvenance",               &_lProvenance,               "_lProvenance[_nL]/i");
        //outputTree->Branch("_lProvenanceCompressed",     &_lProvenanceCompressed,     "_lProvenanceCompressed[_nL]/i");

    }

}

