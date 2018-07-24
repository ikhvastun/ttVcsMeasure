#ifndef filldatacards_H
#define filldatacards_H

#include "Output.h"
#include "readTreeSync.h"
#include "TString.h"

#include <iostream>
#include <cstdlib>

using Output::distribs;
using Output::DistribsAll;

using namespace std;

void fillString(ofstream &, vector<TString> &);
void fillString(ofstream &, vector<std::string> &);
void fillString(ofstream & file, vector<double> & str);
void fillString(ofstream & file, vector<int> & str);
vector<double> formEmptyString(int);
vector<int> formEmptyStringInt(int);
vector<int> formUnityStringInt(int);
void fillExperUnc(ofstream &, vector<std::string> & );

std::vector<double> experUnc      = {1.025}; 
std::vector<TString> experUncName = {"lumi"}; //"JES", "btagl", "btagb", "PDF", "Q2"};
std::vector<TString> ttVprocesses = {"ttW", "ttZ", "ttH", "ttX"};

const int SRNumber = leptonSelectionAnalysis == 2 ? theSRLabelOptionsFor2L.size() : theSRLabelOptionsFor3L.size();

// just an example, nameOfProcessForDatacards: data, ttZ, ttW, etc.
// distribsOrderForYields = 0, 1, 3, ...., 
// 0 - data
// so 1 + 2 = ttZ
// 3 = ttW
void fillDatacards(DistribsAll & distribs, vector<std::string> & nameOfProcessesForDatacard, vector<unsigned> & distribsOrderForYields){

  const int nCategories = nameOfProcessesForDatacard.size(); // nameOfProcessesForDatacard - all MC + 1 for data
  std::vector<double> intYield;

  // first let's calculate total yield in each category
  for (int i = 0; i != distribsOrderForYields.size()-1; ++i) {

     double yield = 0;
     for (int k = 0; k != SRNumber; ++k) {

        for(unsigned int j = distribsOrderForYields.at(i); j != distribsOrderForYields.at(i+1); j++){
          yield += distribs.vectorHisto[j].GetBinContent(k+1);
        }
     }
     intYield.push_back(yield);
  }

  cout << "number of SR and Categories: " << SRNumber << " " << nCategories << endl;
    
  nameOfProcessesForDatacard.erase(nameOfProcessesForDatacard.begin()); // here we erase 1 element - data
  const int numberOfBKG = nameOfProcessesForDatacard.size() - 1; // -1 for signal

   gSystem->Exec("rm datacards/shapes/shapeFile_ttZ3L.root"); // delete previous tex file
   gSystem->Exec("rm datacards/datacard_ttZ3L.txt"); // delete previous tex file
   ofstream fileout;
   fileout .open ( "datacards/datacard_ttZ3L.txt", ios_base::app); // create a new tex file
   fileout << fixed << showpoint << setprecision(2);
   fileout << "imax 1 number of channels " <<  endl;
   fileout << "jmax " << numberOfBKG << " number of backgrounds " <<  endl;
   fileout << "kmax " << (leptonSelectionAnalysis == 2 ? 24 : 81) <<  " number of nuisance parameters (sources of systematical uncertainties) " <<  endl;
   fileout << "----------- " <<  endl;
   TFile *file = TFile::Open("datacards/shapes/shapeFile_ttZ3L.root", "RECREATE");
   fileout << "shapes * * shapes/shapeFile_ttZ3L.root  $PROCESS $PROCESS_$SYSTEMATIC" << endl;
   fileout << "-----------  " <<  endl;
   fileout << "bin  1" <<  endl;
   fileout << "observation  " << intYield[0] << endl;
   fileout << "-----------  " <<  endl;

   vector<TString> datastrVector;
   vector<TString> processNumberVector;
   vector<double> rateVector;

   for(int i = 1; i < nCategories; i++){ // here -1 we don't consider data
      datastrVector.push_back(std::to_string(1));
      processNumberVector.push_back(std::to_string(i - 1));
      rateVector.push_back(intYield[i] < 0.01 ? 0.01 : intYield[i]);
   }     
     
   fileout << "bin     " ;
   fillString(fileout, datastrVector);

   fileout << "process     " ;
   fillString(fileout, nameOfProcessesForDatacard);

   fileout << "process     " ;
   fillString(fileout, processNumberVector);

   fileout << "rate        " ;
   fillString(fileout, rateVector);

   // first fill stat unc
   for(int cat = 0; cat < nCategories; cat++){ 
       TH1D *hist, *histStUp, *histStDown;
       if(cat == 0){
          hist = (TH1D*)distribs.vectorHisto[cat].Clone("hist");
          hist->SetName("data_obs");
          hist->Write();
          continue;
       }
       for(int sam = distribsOrderForYields[cat]; sam < distribsOrderForYields[cat+1]; sam++){
          if(sam == distribsOrderForYields[cat]){
              hist = (TH1D*)distribs.vectorHisto[sam].Clone("hist");
              histStUp = (TH1D*)distribs.vectorHisto[sam].Clone("histUp");
              histStDown  = (TH1D*)distribs.vectorHisto[sam].Clone("histDown");
          }
          else{
            hist->Add(&distribs.vectorHisto[sam]); 
            histStUp->Add(&distribs.vectorHisto[sam]); 
            histStDown->Add(&distribs.vectorHisto[sam]); 
          }
       }
       hist->SetName(nameOfProcessesForDatacard[cat-1].c_str());
       hist->Write();
       for(int sr = 0; sr < SRNumber; sr++){
       
         histStUp->SetBinContent(sr+1, hist->GetBinContent(sr+1) + hist->GetBinError(sr+1));
         histStDown->SetBinContent(sr+1, hist->GetBinContent(sr+1) - hist->GetBinError(sr+1));

         histStUp->SetName((nameOfProcessesForDatacard[cat-1] + "_" + nameOfProcessesForDatacard[cat-1] + "_stat_bin_" + std::to_string(sr+1) + "Up").c_str());
         histStDown->SetName((nameOfProcessesForDatacard[cat-1] + "_" + nameOfProcessesForDatacard[cat-1] + "_stat_bin_" + std::to_string(sr+1) + "Down").c_str());
         histStUp->Write();
         histStDown->Write();

         histStUp->SetBinContent(sr+1, hist->GetBinContent(sr+1));
         histStDown->SetBinContent(sr+1, hist->GetBinContent(sr+1));

         fileout << nameOfProcessesForDatacard[cat-1] + "_stat_bin_" + std::to_string(sr+1) + " shape     " ;
         vector<int> newString = formEmptyStringInt(numberOfBKG + 1); // + 1 for signal
         newString[cat-1] = 1;
         fillString(fileout, newString);
       }
   }

   fileout << "----------- " <<  endl;
   fillExperUnc(fileout, nameOfProcessesForDatacard);
   
   // here fill all syst shape uncertainties 
   std::vector<std::string> systShapeNames = {"lepSF", "pileup", "scale", "bTag_udsg", "bTag_bc"};
   for(int syst = 0; syst < systShapeNames.size(); syst++){
       for(int cat = 1; cat < nCategories; cat++){ 
          TH1D *histStUp, *histStDown;
          for(int sam = distribsOrderForYields[cat]; sam < distribsOrderForYields[cat+1]; sam++){
             if(sam == distribsOrderForYields[cat]){
                histStUp = (TH1D*)distribs.vectorHistoUncUp[sam].unc[syst].Clone("histUp");
                histStDown  = (TH1D*)distribs.vectorHistoUncDown[sam].unc[syst].Clone("histDown");
             }
             else{
                histStUp->Add(&distribs.vectorHistoUncUp[sam].unc[syst]); 
                histStDown->Add(&distribs.vectorHistoUncDown[sam].unc[syst]); 
             }
          }
          histStUp->SetName((nameOfProcessesForDatacard[cat-1] + "_" + systShapeNames.at(syst) + "Up").c_str());
          histStDown->SetName((nameOfProcessesForDatacard[cat-1] + "_" + systShapeNames.at(syst) + "Down").c_str());

          histStUp->Write();
          histStDown->Write();
      }
      fileout <<  systShapeNames.at(syst) << " shape     " ;
      vector<int> newString = formUnityStringInt(numberOfBKG + 1); // + 1 for signal
      newString[2] = 999; // don't consider uncertainty on nonprompt
      fillString(fileout, newString);
   }

   fileout << ("fake     lnN ");
   vector<double> newString = formEmptyString(numberOfBKG + 1);
   for(int i = 0; i < nameOfProcessesForDatacard.size(); i++){
      if(nameOfProcessesForDatacard.at(i) == "nonpromptData" || nameOfProcessesForDatacard.at(i) == "nonprompt")
        newString[i] = 1.3;
   }
   fillString(fileout, newString);

   if(leptonSelectionAnalysis == 2){
      fileout << ("charge     lnN ");
      newString = formEmptyString(numberOfBKG + 1);
      for(int i = 0; i < nameOfProcessesForDatacard.size(); i++){
          if(nameOfProcessesForDatacard.at(i) == "chargeMisID")
            newString[i] = 1.2;
      }
      fillString(fileout, newString);
   }

   // this should be covered by scale and pdf 
   /*
        fileout << ("ttX     lnN ");
        newString = formEmptyString(numberOfBKG + 1);
        for(int i = 1; i < nameOfProcessesForDatacard.size(); i++){
          if(std::find(ttVprocesses.begin(), ttVprocesses.end(), nameOfProcessesForDatacard.at(i)) != ttVprocesses.end() )
            newString[i] = 1.1;
        }
        fillString(fileout, newString);
        */

   fileout << ("WZ     lnN ");
   newString = formEmptyString(numberOfBKG + 1);
   for(int i = 0; i < nameOfProcessesForDatacard.size(); i++){
      if(nameOfProcessesForDatacard.at(i) == "WZ")
        newString[i] = 1.1;
   }
   fillString(fileout, newString);

   fileout << ("rare    lnN ");
   newString = formEmptyString(numberOfBKG + 1);
   for(int i = 0; i < nameOfProcessesForDatacard.size(); i++){
       if(nameOfProcessesForDatacard.at(i) == "rare")
         newString[i] = 1.5;
   }
   fillString(fileout, newString);

   file->Close();
    
    std::cout << "datacard DONE" << std::endl;

}


void fillString(ofstream & file, vector<TString> & str){
  file << fixed << showpoint << setprecision(2);
  for(int i = 0; i < str.size(); i++)
    file <<  str.at(i) << '\t';
  file << endl;
}

void fillString(ofstream & file, vector<std::string> & str){
  file << fixed << showpoint << setprecision(2);
  for(int i = 0; i < str.size(); i++)
    file <<  str.at(i) << '\t';
  file << endl;
}

void fillString(ofstream & file, vector<double> & str){
  file << fixed << showpoint << setprecision(2);
  for(int i = 0; i < str.size(); i++){
    if(str.at(i) != 1.0)
      file <<  str.at(i) << '\t';
    else
      file <<  "-" << '\t';
  }
  file << endl;
}

void fillString(ofstream & file, vector<int> & str){
  file << fixed << showpoint << setprecision(2);
  for(int i = 0; i < str.size(); i++){
    if(str.at(i) != 999)
      file <<  str.at(i) << '\t';
    else
      file <<  "-" << '\t';
  }
  file << endl;
}

vector<double> formEmptyString(int numberOfProcesses){
  vector<double> formString;
  formString.clear();
  for(int i = 0; i < numberOfProcesses; i++)
    formString.push_back(1.0);
  return formString;

}

vector<int> formEmptyStringInt(int numberOfProcesses){
  vector<int> formString;
  formString.clear();
  for(int i = 0; i < numberOfProcesses; i++)
    formString.push_back(999);
  return formString;

}

vector<int> formUnityStringInt(int numberOfProcesses){
  vector<int> formString;
  formString.clear();
  for(int i = 0; i < numberOfProcesses; i++)
    formString.push_back(1);
  return formString;

}


void fillExperUnc(ofstream & file, vector<std::string> & nameOfProcessesForDatacard){
  for(int i = 0; i < experUnc.size(); i++){
    file << experUncName.at(i) << "      lnN  ";
    for(int statInd = 0; statInd < nameOfProcessesForDatacard.size(); statInd++){
      if(nameOfProcessesForDatacard.at(statInd) == "nonpromptData" || nameOfProcessesForDatacard.at(statInd) == "nonprompt")
        file << "-" << '\t' ;
      else if (nameOfProcessesForDatacard.at(statInd) == "chargeMisID")
        file << "-" << '\t' ;
      //else if (nameOfProcessesForDatacard.at(statInd) == "WZ" && experUncName.at(i) == "LeptonId")
      //  file << "-" << '\t' ;
      else
        file << experUnc.at(i) << '\t' ;
    }
    file << endl;
  }


}


#endif  // fillDatacards
