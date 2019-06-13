//include ROOT classes
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TROOT.h"

//include C++ library classes
#include <sstream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <algorithm>

//include other parts of the code
#include "../interface/treeReader.h"
#include "../interface/analysisTools.h"


void treeReader::skimTree(const std::string& fileName, std::string outputDirectory, const bool isData){//std::string outputFileName){
    //Read tree	
    std::shared_ptr<TFile> sampleFile = std::make_shared<TFile>( (const TString& ) fileName,"read");	
    sampleFile->cd("blackJackAndHookers");
    //Determine hcounter and lheCounter for MC cross section scaling and uncertainties
    TH1D* hCounter;
    TH1D* lheCounter;
    TH1D* nTrueInteractions;
    TH1D* psCounter;
    if(!isData){
        hCounter = (TH1D*) sampleFile->Get("blackJackAndHookers/hCounter");
        lheCounter = (TH1D*) sampleFile->Get("blackJackAndHookers/lheCounter");
        nTrueInteractions = (TH1D*) sampleFile->Get("blackJackAndHookers/nTrueInteractions");
        psCounter = (TH1D*) sampleFile->Get("blackJackAndHookers/psCounter");
    }
    //Get Tree
    TTree* sampleTree = (TTree*) (sampleFile->Get("blackJackAndHookers/blackJackAndHookersTree"));
    initTree(sampleTree, isData);
    outputDirectory = (outputDirectory == "") ? "~/Work/ntuples_temp/" : outputDirectory;
    std::string outputFileName = fileName;
    /*
    auto it = outputFileName.find_last_of("/");
    outputFileName.erase(0, it + 1);
    */
    outputFileName.erase(std::remove(outputFileName.begin(), outputFileName.end(), '/'), outputFileName.end());
    outputFileName.insert(0, outputDirectory);
    auto it = outputFileName.find(".root");
    outputFileName.erase(it, outputFileName.size());
    outputFileName.append("_dilepSkim.root");
    std::cout << "output file : " << outputFileName << std::endl;
    
    std::shared_ptr<TFile> outputFile = std::make_shared<TFile>((const TString&) outputFileName ,"RECREATE");
    std::shared_ptr<TTree> outputTree = std::make_shared<TTree>("blackJackAndHookersTree","blackJackAndHookersTree");
    outputFile->mkdir("blackJackAndHookers");
    outputFile->cd("blackJackAndHookers"); 
    setOutputTree(outputTree.get(), isData);
    //
    // here is the event list to be checked
    //std::set<ULong64_t> eventsList = {1492195608, 146025307, 508933991, 321197298, 277740241, 892590768, 26652538, 315966641, 320474674, 59566750, 66424217, 144766203, 252985782, 65013805, 133493915,       102651335, 455363823, 111273860, 373187240, 239729165, 239078987, 195340703, 1856389, 461581059, 506460369, 656450946, 184079568, 963204027, 963441575, 336829296, 115879231, 41506348, 640566084,          78165014, 1255495340, 64807712, 2666356, 107674345, 529788719, 293132512, 39548847, 267438002, 55582650, 111863351, 1061848933, 217235889, 562481829, 704561987, 480406039, 744805798, 202616779,           327142299, 167976892, 391402062, 327599224, 503087687, 98671267, 198248356, 90666260, 406200233, 272203871, 163786055, 253703343, 142520294, 382842138, 405922138, 111982906, 210867535, 135955294,         252681643, 1150234279, 1103452740, 1262721932, 36893531, 312773489, 253682178, 450349644, 37685662, 38382410, 400846856, 552906843, 82761276, 664083895, 325380020, 585216792, 102898869, 18194142,         637834750, 573681949, 318118288, 529995837, 144910725, 14517387, 226434855, 192484572, 1617293625, 196222566, 75686784, 1147847462, 108502231, 1534207907, 1306527199, 131111757, 48742638, 131076926,      385195332, 694375811, 384522728, 894287433, 370558190, 21779924, 1099098547, 898202951, 1216158608, 253240568, 963202886, 367361713, 587481963, 115793233, 322860364, 211654778, 1206592648, 1655550800,    1350511364, 1695619142, 1754680495, 267327098, 1110306499, 5545114, 27071751, 211503876, 3377288405, 55721552, 618991682, 1117123544, 142703839, 517934154, 834180370, 1371104318, 1269000161, 2763422568,  1083197265, 1444157466, 1421307868, 1027796765, 1861225520, 1372234311, 172620241, 1040650747, 342353721, 234086389, 236104091, 76206227, 204990413, 39938667, 191786789, 325353201, 183963204, 158559514,  96976446, 157532648, 73280873, 97907529, 197833912, 697733424, 418136939, 305060877, 1273299716, 778786389, 131142593, 842492785, 332241563, 760618473, 244050843, 675289609, 694369509, 1053432589,        1266575101, 598168272, 638247408, 14814758, 1151944380, 1301114924, 101525287, 274827271, 1041988330, 321245229, 390563015, 508574825, 344972376, 522431310, 221056243, 40355959, 409853777, 200561756,     74081286, 89251178, 119661967, 207168464, 58509301, 251183962, 28591491, 30036062, 1124290562, 289965738, 1090947268, 699472740, 75939013, 411070835, 536170736, 144462225, 528847863, 260371307,           202531851, 1106760047, 82655785, 773526112, 248826088, 83999606, 297716316, 18276366, 1085186745, 372661440, 729677885, 1318455870};

    addVariablesToBDT();
    double progress = 0; 	//For printing progress bar
    long nEntries = sampleTree->GetEntries();
    for (long it=0; it <nEntries; ++it) {
        if(it%100 == 0 && it != 0){
            progress += (double) (100./ (double) nEntries);
            tools::printProgress(progress);
        } else if(it == nEntries -1){
            progress = 1.;
            tools::printProgress(progress);
        }
        sampleTree->GetEntry(it);
        std::vector<unsigned> ind, indJets;
        unsigned lCount = selectLepGoodForLeptonMVA(ind);

        if(lCount < 1) continue;
        if(_nJets < 2) continue;

        /*
        bool allPrompt = true;
        for(auto & i : ind){
            if(!_lIsPrompt[i])
                allPrompt = false;
        }
        if(allPrompt) continue;
        */
        
        /*
        bool allChargeConsistent = true;
        for(auto & i : ind){
            if(leptonChargeIsInconsistent(i))
                allChargeConsistent = false;
        }
        if(allChargeConsistent) continue;
        */
        
        // this is used for FR measurement in QCD
        /*
        if(lCount < 1) continue;
        unsigned jetCount = nJets(indJets, 0, true);
        if(!deltaRCalcForWS(indJets, ind.at(0))) continue;
        */

        // additional requirement for FR in data
        //if(!(_HLT_Mu3_PFJet40 || _HLT_Mu8 || _HLT_Mu17 || _HLT_Mu27 || _HLT_Ele8_CaloIdM_TrackIdM_PFJet30 || _HLT_Ele17_CaloIdM_TrackIdM_PFJet30 || _HLT_Ele23_CaloIdM_TrackIdM_PFJet30)) continue;

        // same sign selection
        /*
        if(lCount < 2) continue;
        if(lCount == 2){
            if(_lCharge[ind.at(0)] * _lCharge[ind.at(1)] < 0)
                continue;
        }
        */

        // is used for ttbar emu
        /*
        if(lCount < 2) continue;
        int nLocEle = getElectronNumber(ind);
        int nLocMu = lCount - nLocEle;

        if(nLocEle < 1) continue;
        if(nLocMu < 1) continue;
        */

        /* DY-> ee
        if(lCount < 2) continue;
        int nLocEle = getElectronNumber(ind);
        if(nLocEle < 2) continue;
        int positiveCharge = 0;
        int negativeCharge = 0;
        for(auto & i : ind){
            if(_lFlavor[i] == 0){
                if(_lCharge[i] == 1)
                    positiveCharge+=1;
                if(_lCharge[i] == -1)
                    negativeCharge+=1;
            }
        }

        if(positiveCharge < 2 && negativeCharge < 2) continue;
        */

        //Zll MET task and also DY -> 2L on-Z selection
        /*
        double deltaMZ = 99999.;
        //std::cout << "Number of leptons: " << lCount << std::endl;
        if(lCount < 2) continue;
        */
        
        // we are interested only in ee events, muon are not needed here
        /*
        int nLocEle = getElectronNumber(ind);
        if(nLocEle < 2) continue;
        */

        /*
        for (int l0 = 0; l0<lCount; ++l0) {
            TLorentzVector l0p4;
            l0p4.SetPtEtaPhiE(_lPt[ind.at(l0)],_lEta[ind.at(l0)],_lPhi[ind.at(l0)],_lE[ind.at(l0)]);
            for(int l1 = l0 + 1; l1 < lCount; ++l1){
                // needed for DY prompt CR and Zll MET task, not needed for charge mis ID measurement, where basically we should meausre ee / e+e-, so both same charge and opposite charge events are nedeed
                if (_lCharge[ind.at(l0)] != _lCharge[ind.at(l1)]) {
                    TLorentzVector l1p4;
                    l1p4.SetPtEtaPhiE(_lPt[ind.at(l1)],_lEta[ind.at(l1)],_lPhi[ind.at(l1)],_lE[ind.at(l1)]);
                    l1p4+=l0p4;
                    double mdiL = l1p4.M();
                    if (_lFlavor[ind.at(l0)] == _lFlavor[ind.at(l1)] ) {
                        if (fabs(mdiL - 91.2) < deltaMZ) {
                            deltaMZ = fabs(mdiL - 91.2);
                        }
                    }
                }
            }
        }

        if(deltaMZ > 10) continue;
        */

        // what is needed to run over list of events
        /*
        const bool is_in = eventsList.find(_eventNb) != eventsList.end();
        if(!is_in) continue;
        */

        for(int i = 0; i < _nLight; i++){
            std::vector<Float_t> varForBDT = { (Float_t)_lPt[i], (Float_t)fabs(_lEta[i]), (Float_t)_selectedTrackMult[i], (Float_t)_miniIsoCharged[i], (Float_t)(_miniIso[i] - _miniIsoCharged[i]), (Float_t)_ptRel[i], (Float_t)_ptRatio[i], (Float_t)_relIso[i],(Float_t)(std::isnan(_closestJetDeepCsv_b[i] + _closestJetDeepCsv_bb[i]) ? 0. : std::max(_closestJetDeepCsv_b[i] + _closestJetDeepCsv_bb[i], 0.)),(Float_t)_3dIPSig[i], (Float_t)log( fabs(_dxy[i])), (Float_t)log( fabs(_dz[i])), 
            _lFlavor[i] == 0 ? (Float_t)_lElectronMva[i] : (Float_t)_lMuonSegComp[i]};
            //_lFlavor[i] == 0 ? (Float_t)_lElectronMvaFall17NoIso[i] : (Float_t)_lMuonSegComp[i]};
            fillBDTvariables(varForBDT, _lFlavor[i]);
            //if(_closestJetDeepCsv_bb[i] + _closestJetDeepCsv_b[i] != _closestJetDeepCsv_bb[i] + _closestJetDeepCsv_b[i]) continue;
            double mvaVL =  _lFlavor[i] == 0 ? readerLeptonMVAele->EvaluateMVA( "BDTG method") : readerLeptonMVAmu->EvaluateMVA( "BDTG method");
            _leptonMvatZqTTV[i] = mvaVL;
            //std::cout << "mva value is " << mvaVL << std::endl;
        }
        //std::vector<unsigned> indTightTTW;
        //unsigned lCountTightTTW = selectLepTightTTW(indTightTTW);
        //if(lCountTightTTW < 2) continue;

        outputTree->Fill();
    }   
    std::cout << std::endl;
    if(!isData){
        hCounter->Write();
        lheCounter->Write();
        nTrueInteractions->Write();
        //psCounter->Write();
    }
    outputTree->Write("",  BIT(2));
    outputFile->Close(); 
}

int main(int argc, char* argv[]){
    treeReader reader;
    bool isData = false;
    if(argc != 0){
        std::vector<std::string> datasets = {"SingleElectron", "SingleMuon", "DoubleEG", "DoubleMuon", "MuonEG"}; 
        for(auto it = datasets.cbegin(); it != datasets.cend(); ++it){
            std::string name(argv[1]);
            auto pos = name.find(*it);
            if(pos < name.size()) isData = true;
        }
    }
    switch(argc){
        case 2:{
                   reader.skimTree(argv[1], "", isData);
                   return 0;
               }
        case 3:{
                   reader.skimTree(argv[1], argv[2], isData);
                   return 0;
               }
        default:{
                    std::cerr << "Error: Wrong number of options given!" << std::endl;
                    return 1;
                }
    }
}



