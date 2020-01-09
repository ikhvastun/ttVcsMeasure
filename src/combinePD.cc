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

//include other parts of the code
#include "../interface/treeReader.h"
#include "../interface/analysisTools.h"


void treeReader::combinePD(const std::vector<std::string>& datasets, std::string outputDirectory){
    std::set<std::tuple<long unsigned, long unsigned, long unsigned> > usedEvents;
    //Set output file and tree
    const std::string outputFileName = "~/Work/ntuples_ttz_dilep_2018/data_2018.root";
    TFile* outputFile = new TFile((const TString&) outputFileName ,"RECREATE");
    outputFile->mkdir("blackJackAndHookers");
    outputFile->cd("blackJackAndHookers"); 
    TTree* outputTree = new TTree("blackJackAndHookersTree","blackJackAndHookersTree");
    setOutputTree(outputTree, true);
    for(std::vector<std::string>::const_iterator it = datasets.cbegin(); it != datasets.cend(); ++it){
        std::cout << *it << std::endl;
        //Read tree	
        TFile* sampleFile = new TFile( (const TString&)"~/Work/ntuples_ttz_dilep_2018/" + *it,"read");//ntuples_ttV, ntuples_ttV_2017	
        //Determine hcounter for cross section scaling
        sampleFile->cd("blackJackAndHookers");	
        TTree* sampleTree = (TTree*) (sampleFile->Get("blackJackAndHookers/blackJackAndHookersTree"));
        initTree(sampleTree, true);

        double progress = 0; 	//For printing progress bar
        long nEntries = sampleTree->GetEntries();
        for (long it=0; it <nEntries; ++it){
            sampleTree->GetEntry(it);
            if(it%100 == 0 && it != 0){
                progress += (double) (100./ (double) nEntries);
                tools::printProgress(progress);
            } else if(it == nEntries -1){
                progress = 1.;
                tools::printProgress(progress);
            }
            if(usedEvents.find(std::make_tuple(_runNb, _lumiBlock, _eventNb) ) == usedEvents.end()){
                usedEvents.insert(std::make_tuple(_runNb, _lumiBlock, _eventNb) );
            } else{
                continue;
            }
            outputTree->Fill();
        }
        //sampleFile->Close();
        std::cout << std::endl;
    }
    outputFile->cd("blackJackAndHookers"); 
    outputTree->Write("",  BIT(2));
    outputFile->Close();
}

int main(int argc, char* argv[]){	
//    std::vector<std::string> datasets = {"SingleElectron.root", "SingleMuon.root", "DoubleEG.root", "DoubleMuon.root", "MuonEG.root"}; 
				std::vector<std::string> datasets = {"DoubleMuon_MiniAOD2018.root", "EGamma_MiniAOD2018.root", "MuonEG_MiniAOD2018.root", "SingleMuon_MiniAOD2018.root"};
    //std::vector<std::string> datasets = {"SingleElectron.root", "DoubleEG.root", "MuonEG.root"}; 
    //std::vector<std::string> datasets = {"data_SEDE.root", "MuonEG.root"}; 
    treeReader reader;
    switch(argc){
        case 1:{
                   reader.combinePD(datasets);
                   return 0;
               }
        case 2:{
                   reader.combinePD(datasets, argv[1]);
                   return 0;
               }
        default:{
                    std::cerr << "Error: Wrong number of options given!" << std::endl;
                    return 1;
                }
    }
}

