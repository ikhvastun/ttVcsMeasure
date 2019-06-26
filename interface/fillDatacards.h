#ifndef filldatacards_H
#define filldatacards_H

#include "Output.h"
#include "readTreeSync.h"
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"

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
void fillExperUnc(ofstream &, vector<std::string> &, std::vector<TString> &, std::vector<double> &);
double largestAmongAll(const std::vector<double> & weights);
double smallestAmongAll(const std::vector<double> & weights);

void fillDatacards(DistribsAll & distribs, vector<std::string> nameOfProcessesForDatacard, const TString name, bool is2017 = false){

  std::vector<double> experUnc      = {1.025, 1.01};  // for trigger agreed to reduce it to 1%
  std::vector<TString> experUncName = {"lumi" + (std::string)(is2017 ? "2017" : "2016"), "trigger" + (std::string)(is2017 ? "2017" : "2016")}; //"JES", "btagl", "btagb", "PDF", "Q2"};
  std::vector<TString> ttVprocesses = {"ttW", "ttZ", "ttH", "ttX"};

  const int SRNumber = figNames[name].nBins;
  const int nCategories = nameOfProcessesForDatacard.size(); // nameOfProcessesForDatacard - all MC + 1 for data
  std::vector<double> intYield;

  // first let's calculate total yield in each category
  for (int i = 0; i != nameOfProcessesForDatacard.size(); ++i) {

     double yield = 0;
     for (int k = 0; k != SRNumber; ++k) {
        yield += distribs.vectorHisto[i].GetBinContent(k+1);
     }
     intYield.push_back(yield);
  }

  cout << "number of SR and Categories: " << SRNumber << " " << nCategories << endl;
    
  nameOfProcessesForDatacard.erase(nameOfProcessesForDatacard.begin()); // here we erase 1 element - data
  const int numberOfBKG = nameOfProcessesForDatacard.size() - 1; // -1 for signal

   gSystem->Exec("rm datacards/shapes/shapeFile_" + name + (TString)(is2017 ? "2017" : "2016") + ".root"); // delete previous tex file
   gSystem->Exec("rm datacards/datacard_" + name + (TString)(is2017 ? "2017" : "2016") + ".txt"); // delete previous tex file
   gSystem->Exec("mkdir -p datacards/"); // delete previous tex file
   gSystem->Exec("mkdir -p datacards/shapes"); // delete previous tex file
   ofstream fileout;
   fileout .open ( "datacards/datacard_" + name + (TString)(is2017 ? "2017" : "2016") + ".txt", ios_base::app); // create a new tex file
   fileout << fixed << showpoint << setprecision(2);
   fileout << "imax 1 number of channels " <<  endl;
   fileout << "jmax " << numberOfBKG << " number of backgrounds " <<  endl;
   //fileout << "kmax " << (lepSel == 2 ? 195 : (lepSel == 34 ? 129 : (lepSel == 3 ? 96 : 33))) <<  " number of nuisance parameters (sources of systematical uncertainties) " <<  endl;
   fileout << "kmax 130 number of nuisance parameters (sources of systematical uncertainties) " <<  endl;
   fileout << "----------- " <<  endl;
   TFile *file = TFile::Open("datacards/shapes/shapeFile_" + name + (TString)(is2017 ? "2017" : "2016") + ".root", "RECREATE");
   fileout << "shapes * * shapes/shapeFile_" + name + (TString)(is2017 ? "2017" : "2016") + ".root  $PROCESS $PROCESS_$SYSTEMATIC" << endl;
   fileout << "-----------  " <<  endl;
   fileout << "bin  bin1" <<  endl;
   fileout << "observation  " << intYield[0] << endl;
   fileout << "-----------  " <<  endl;

   vector<TString> datastrVector;
   vector<TString> processNumberVector;
   vector<double> rateVector;

   for(int i = 1; i < nCategories; i++){ // here -1 we don't consider data
      datastrVector.push_back("bin" + std::to_string(1));
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
       hist = (TH1D*)distribs.vectorHisto[cat].Clone("hist");
       histStUp = (TH1D*)distribs.vectorHisto[cat].Clone("histUp");
       histStDown  = (TH1D*)distribs.vectorHisto[cat].Clone("histDown");

       hist->SetName(nameOfProcessesForDatacard[cat-1].c_str());
       hist->Write();

       for(int sr = 0; sr < SRNumber; sr++){
       
         histStUp->SetBinContent(sr+1, hist->GetBinContent(sr+1) + hist->GetBinError(sr+1));
         histStDown->SetBinContent(sr+1, hist->GetBinContent(sr+1) - hist->GetBinError(sr+1));

         histStUp->SetName((nameOfProcessesForDatacard[cat-1] + "_" + nameOfProcessesForDatacard[cat-1] + "_stat" + std::string(name) + "_bin_" + std::to_string(sr+1) + "Up").c_str());
         histStDown->SetName((nameOfProcessesForDatacard[cat-1] + "_" + nameOfProcessesForDatacard[cat-1] + "_stat" + std::string(name) + "_bin_" + std::to_string(sr+1) + "Down").c_str());
         histStUp->Write();
         histStDown->Write();

         histStUp->SetBinContent(sr+1, hist->GetBinContent(sr+1));
         histStDown->SetBinContent(sr+1, hist->GetBinContent(sr+1));

         fileout << nameOfProcessesForDatacard[cat-1] + "_stat" + name + "_bin_" + std::to_string(sr+1) + " shape     " ;
         vector<int> newString = formEmptyStringInt(numberOfBKG + 1); // + 1 for signal
         newString[cat-1] = 1;
         fillString(fileout, newString);
       }
   }

   fileout << "----------- " <<  endl;
   fillExperUnc(fileout, nameOfProcessesForDatacard, experUncName, experUnc);
   
   // here fill all syst shape uncertainties 
   std::vector<std::string> systShapeNames = {"lepSFsyst", "lepSFstat" + (std::string)(is2017 ? "2017" : "2016"), "lepSFReco", "pileup", "bTag_udsg" + (std::string)(is2017 ? "2017" : "2016"), "bTag_bc" + (std::string)(is2017 ? "2017" : "2016"), "jec", "jer", "WZbb"}; // , "ISRandFSR"
   for(int syst = 0; syst < systShapeNames.size(); syst++){
       for(int cat = 1; cat < nCategories; cat++){ 
          TH1D *histStUp, *histStDown;
          histStUp = (TH1D*)distribs.vectorHistoUncUp[cat].unc[syst].Clone("histUp");
          histStDown  = (TH1D*)distribs.vectorHistoUncDown[cat].unc[syst].Clone("histDown");

          histStUp->SetName((nameOfProcessesForDatacard[cat-1] + "_" + systShapeNames.at(syst) + "Up").c_str());
          histStDown->SetName((nameOfProcessesForDatacard[cat-1] + "_" + systShapeNames.at(syst) + "Down").c_str());

          histStUp->Write();
          histStDown->Write();
      }
      fileout <<  systShapeNames.at(syst) << " shape     " ;
      vector<int> newString;
      if(systShapeNames[syst] == "WZbb"){
        newString = formEmptyStringInt(numberOfBKG + 1);
        newString[3] = 1;
      }
      /*
      else if(systShapeNames[syst] == "ISRandFSR"){
        newString = formEmptyStringInt(numberOfBKG + 1);
        newString[0] = 1.;
      }
      */
      else{
        newString = formUnityStringInt(numberOfBKG + 1); // + 1 for signal
        newString[7] = 999; // don't consider uncertainty on nonprompt
      }
      fillString(fileout, newString);
   }

   std::vector<std::string> systISRFSRNames = {"ISR", "FSR"};
   for(int cat = 1; cat < nCategories; cat++){ 
        if(cat != 1) continue;
        TH1D *histISRandFSRUp, *histISRandFSRDown, *histISRUp, *histISRDown, *histFSRUp, *histFSRDown;

        histISRandFSRUp = (TH1D*)distribs.vectorHistoUncUp[cat].unc[systShapeNames.size()].Clone("histISRandFSRUp");
        histISRandFSRDown  = (TH1D*)distribs.vectorHistoUncDown[cat].unc[systShapeNames.size()].Clone("histISRandFSRDown");

        histISRUp = (TH1D*)distribs.vectorHistoUncUp[cat].unc[systShapeNames.size()].Clone("histISRUp");
        histISRDown  = (TH1D*)distribs.vectorHistoUncDown[cat].unc[systShapeNames.size()].Clone("histISRDown");
        histFSRUp = (TH1D*)distribs.vectorHistoUncUp[cat].unc[systShapeNames.size()+1].Clone("histFSRUp");
        histFSRDown  = (TH1D*)distribs.vectorHistoUncDown[cat].unc[systShapeNames.size()+1].Clone("histFSRDown");

        for(int sr = 0; sr < SRNumber; sr++){
       
            std::vector<double> valueDev;
            valueDev.clear();
            valueDev.push_back(histISRUp->GetBinContent(sr+1));
            valueDev.push_back(histISRDown->GetBinContent(sr+1));
            valueDev.push_back(histFSRUp->GetBinContent(sr+1));
            valueDev.push_back(histFSRDown->GetBinContent(sr+1));

            histISRandFSRUp->SetBinContent(sr+1, largestAmongAll(valueDev));
            histISRandFSRDown->SetBinContent(sr+1, smallestAmongAll(valueDev));
        }

        histISRandFSRUp->SetName((nameOfProcessesForDatacard[cat-1] + "_ISRandFSRUp").c_str());
        histISRandFSRDown->SetName((nameOfProcessesForDatacard[cat-1] + "_ISRandFSRDown").c_str());

        histISRandFSRUp->Write();
        histISRandFSRDown->Write();

        fileout <<  "ISRandFSR shape     " ;
        vector<int> newString;
        newString = formEmptyStringInt(numberOfBKG + 1);
        newString[0] = 1.;
        fillString(fileout, newString);
   }

   // here fill all syst shape uncertainties that have acceptance uncertainty 
   std::vector<std::string> systShapeAcceptNames = {"scaleAcc", "pdfAcc"};
   for(int syst = systShapeNames.size() + systISRFSRNames.size(); syst < systShapeAcceptNames.size() + systShapeNames.size() + systISRFSRNames.size(); syst++){
       for(int cat = 1; cat < nCategories; cat++){ 
          TH1D *hist, *histStUp, *histStDown, *histAccUp, *histAccDown;
          hist = (TH1D*)distribs.vectorHisto[cat].Clone("hist");

          histStUp = (TH1D*)distribs.vectorHistoUncUp[cat].unc[syst].Clone("histUp");
          histStDown  = (TH1D*)distribs.vectorHistoUncDown[cat].unc[syst].Clone("histDown");
                
          histAccUp = (TH1D*)distribs.vectorHistoUncUp[cat].unc[syst].Clone("histAccUp");
          histAccDown = (TH1D*)distribs.vectorHistoUncDown[cat].unc[syst].Clone("histAccDown");

          histAccUp->Reset("ICE");
          histAccDown->Reset("ICE");
          
          for(int sr = 0; sr < SRNumber; sr++){
             histAccUp->SetBinContent(sr+1, (histStUp->GetBinContent(sr+1) / histStUp->Integral()) / (hist->GetBinContent(sr+1) / hist->Integral()) * hist->GetBinContent(sr+1));
             histAccDown->SetBinContent(sr+1, (histStDown->GetBinContent(sr+1) / histStDown->Integral()) / (hist->GetBinContent(sr+1) / hist->Integral()) * hist->GetBinContent(sr+1));
          }


          histAccUp->SetName((nameOfProcessesForDatacard[cat-1] + "_" + systShapeAcceptNames.at(syst - systShapeNames.size() - systISRFSRNames.size()) + "Up").c_str());
          histAccDown->SetName((nameOfProcessesForDatacard[cat-1] + "_" + systShapeAcceptNames.at(syst - systShapeNames.size() - systISRFSRNames.size()) + "Down").c_str());

          histAccUp->Write();
          histAccDown->Write();
      }
      fileout <<  systShapeAcceptNames.at(syst - systShapeNames.size() - systISRFSRNames.size()) << " shape     " ;
      vector<int> newString = formUnityStringInt(numberOfBKG + 1); // + 1 for signal
      newString[7] = 999; // don't consider uncertainty on nonprompt
      //if(lepSel == 2)
      //    newString[3] = 999;
      fillString(fileout, newString);
   }

   fileout << ("nonprompt     lnN ");
   vector<double> newString = formEmptyString(numberOfBKG + 1);
   for(int i = 0; i < nameOfProcessesForDatacard.size(); i++){
      if(nameOfProcessesForDatacard.at(i) == "nonpromptData" || nameOfProcessesForDatacard.at(i) == "nonprompt")
        newString[i] = 1.3;
   }
   fillString(fileout, newString);

   /*
   if(lepSel == 2){
      fileout << ("charge     lnN ");
      newString = formEmptyString(numberOfBKG + 1);
      for(int i = 0; i < nameOfProcessesForDatacard.size(); i++){
          if(nameOfProcessesForDatacard.at(i) == "chargeMisID" || nameOfProcessesForDatacard.at(i) == "chargeMisIDData")
            newString[i] = 1.2;
      }
      fillString(fileout, newString);
   }
   */

   // this should be covered by scale and pdf 
   fileout << ("ttX     lnN ");
   newString = formEmptyString(numberOfBKG + 1);
   for(int i = 1; i < nameOfProcessesForDatacard.size(); i++){
       if(std::find(ttVprocesses.begin(), ttVprocesses.end(), nameOfProcessesForDatacard.at(i)) != ttVprocesses.end() )
       newString[i] = 1.11;
   }
   fillString(fileout, newString);

   fileout << ("WZ     lnN ");
   newString = formEmptyString(numberOfBKG + 1);
   for(int i = 0; i < nameOfProcessesForDatacard.size(); i++){
      if(nameOfProcessesForDatacard.at(i) == "WZ")
        newString[i] = 1.1;
   }
   fillString(fileout, newString);

   fileout << ("ZZ     lnN ");
   newString = formEmptyString(numberOfBKG + 1);
   for(int i = 0; i < nameOfProcessesForDatacard.size(); i++){
      if(nameOfProcessesForDatacard.at(i) == "ZZ")
        newString[i] = 1.1;
   }
   fillString(fileout, newString);

   //if(lepSel != 2){
     fileout << ("Xgamma     lnN ");
     newString = formEmptyString(numberOfBKG + 1);
     for(int i = 0; i < nameOfProcessesForDatacard.size(); i++){
       if(nameOfProcessesForDatacard.at(i) == "Xgamma")
         newString[i] = 1.2;
     }
     fillString(fileout, newString);
   //}

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


void fillExperUnc(ofstream & file, vector<std::string> & nameOfProcessesForDatacard, std::vector<TString> & experUncName, std::vector<double> & experUnc){
  for(int i = 0; i < experUnc.size(); i++){
    file << experUncName.at(i) << "      lnN  ";
    for(int statInd = 0; statInd < nameOfProcessesForDatacard.size(); statInd++){
      if(nameOfProcessesForDatacard.at(statInd) == "nonpromptData" || nameOfProcessesForDatacard.at(statInd) == "nonprompt")
        file << "-" << '\t' ;
      else if (nameOfProcessesForDatacard.at(statInd) == "chargeMisID" || nameOfProcessesForDatacard.at(statInd) == "chargeMisIDData")
        file << "-" << '\t' ;
      //else if (nameOfProcessesForDatacard.at(statInd) == "WZ" && experUncName.at(i) == "LeptonId")
      //  file << "-" << '\t' ;
      else
        file << experUnc.at(i) << '\t' ;
    }
    file << endl;
  }


}

double largestAmongAll(const std::vector<double> & weights){

    double max_value = weights.front() ;
    for( std::size_t i = 1 ; i < weights.size() ; ++i ) 
        if( weights[i] > max_value ) max_value = weights[i] ;
    return max_value;

}

double smallestAmongAll(const std::vector<double> & weights){

    double min_value = weights.front() ;
    for( std::size_t i = 1 ; i < weights.size() ; ++i ) 
        if( weights[i] < min_value ) min_value = weights[i] ;
    return min_value;

}

#endif  // fillDatacards
