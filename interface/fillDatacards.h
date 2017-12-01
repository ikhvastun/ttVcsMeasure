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
vector<double> formEmptyString(int);
void fillExperUnc(ofstream &, vector<std::string> & );

std::vector<double> experUnc      = {1.025,  1.01, 1.03, (leptonSelectionAnalysis == 2 ? 1.03 : 1.04), 1.01}; // 1.01,  1.01,    1.01,    1.01,  1.01};
std::vector<TString> experUncName = {"lumi", "PU", "trigger", "LeptonId", "JER", }; //"JES", "btagl", "btagb", "PDF", "Q2"};
std::vector<TString> ttVprocesses = {"ttW", "ttZ", "ttH", "ttX"};

const int SRNumber = leptonSelectionAnalysis == 2 ? theSRLabelOptionsFor2L.size() : theSRLabelOptionsFor3L.size();

//std::vector<int> distribsOrderForYields;
//std::vector<TString> nameOfProcessesForDatacard;
std::vector<int> distribsOrderForDatacards; 

void fillDatacards(DistribsAll & distribs, vector<std::string> & nameOfProcessesForDatacard, vector<unsigned> & distribsOrderForYields){

  /*
  There are 3 vectors here:
    distribsOrderForDatacards -> order that is used for datacard, from 1 to 9, where 1 stands for ttW for 2L, 9 for rare

  */

  const int nCategories = nameOfProcessesForDatacard.size(); // nameOfProcessesForDatacard - all MC + 1 for data

  vector<double> wz_jec, wz_bl, wz_bb, ttx_jec, ttx_bl, ttx_bb, tth_jec, tth_bl, tth_bb, ttw_jec, ttw_bl, ttw_bb, ttz_jec, ttz_bl, ttz_bb, WZbkg;

  if(leptonSelectionAnalysis == 3){
    wz_jec = {1.03, 1.05, 1.06, 1.02, 1.05, 1.08, 1.00, 1.06, 1.04};
    wz_bl = {1.02, 1.02, 1.04, 1.04, 1.05, 1.06, 1.03, 1.04, 1.06};
    wz_bb = {1.00, 1.00, 1.00, 1.01, 1.01, 1.01, 1.03, 1.03, 1.02};

    ttx_jec = {1.00, 1.03, 1.05, 1.00, 1.03, 1.05, 1.00, 1.02, 1.05};
    ttx_bl = {1.01, 1.02, 1.03, 1.01, 1.02, 1.03, 1.00, 1.02, 1.03};
    ttx_bb = {1.00, 1.00, 1.00, 1.01, 1.01, 1.01, 1.02, 1.02, 1.03};

    tth_jec = {1.00, 0.97, 1.05, 0.97, 0.98, 1.03, 0.96, 1.03, 1.03};
    tth_bl = {1.01, 1.01, 1.03, 1.01, 1.01, 1.03, 1.00, 1.02, 1.03};
    tth_bb = {1.00, 1.00, 1.00, 1.01, 1.01, 1.01, 1.02, 1.03, 1.03};

    ttw_jec = {1.00, 1.00, 1.03, 0.99, 1.01, 1.04, 0.99, 1.01, 1.03};
    ttw_bl = {1.01, 1.01, 1.03, 1.01, 1.01, 1.03, 1.00, 1.01, 1.03};
    ttw_bb = {1.00, 1.00, 1.00, 1.01, 1.01, 1.01, 1.02, 1.02, 1.03};

    ttz_jec = {0.96, 0.98, 1.02, 0.95, 0.99, 1.02, 0.96, 0.98, 1.03};
    ttz_bl = {1.01, 1.01, 1.03, 1.01, 1.01, 1.03, 1.01, 1.01, 1.03};
    ttz_bb = {1.00, 1.00, 1.00, 1.01, 1.01, 1.01, 1.02, 1.03, 1.03};

    WZbkg = {1.10, 1.10, 1.20, 1.10, 1.10, 1.20, 1.10, 1.10, 1.20};
  }

  if(leptonSelectionAnalysis == 2){
    ttw_jec = {1.05, 1.03, 1.01, 1.02, 1.03, 1.03, 1.01, 1.01, 1.01, 1.03, 1.05, 1.03, 1.01, 1.02, 1.03, 1.03, 1.01, 1.01, 1.01, 1.03, 1.1, 1.1, 1.1};
    ttw_bl = {1.01, 1.01, 1.01, 1.02, 1.02, 1.01, 1.01, 1.01, 1.03, 1.02, 1.01, 1.01, 1.01, 1.02, 1.02, 1.01, 1.01, 1.01, 1.03, 1.02, 1.1, 1.1, 1.1};
    ttw_bb = {1.02, 1.03, 1.05, 1.07, 1.03, 1.04, 1.03, 1.06, 1.07, 1.05, 1.02, 1.03, 1.05, 1.07, 1.03, 1.04, 1.03, 1.06, 1.07, 1.05, 1.1, 1.1, 1.1};

    ttz_jec = {1.07, 1.01, 1.02, 1.03, 1.01, 1.09, 1.05, 1.01, 1.03, 1.04, 1.07, 1.01, 1.02, 1.03, 1.01, 1.09, 1.05, 1.01, 1.03, 1.04, 1.1, 1.1, 1.1};
    ttz_bl = {1.01, 1.01, 1.01, 1.02, 1.02, 1.01, 1.01, 1.01, 1.02, 1.01, 1.01, 1.01, 1.01, 1.02, 1.02, 1.01, 1.01, 1.01, 1.02, 1.01, 1.1, 1.1, 1.1};
    ttz_bb = {1.01, 1.01, 1.02, 1.01, 1.02, 1.02, 1.01, 1.02, 1.01, 1.02, 1.01, 1.01, 1.02, 1.01, 1.02, 1.02, 1.01, 1.02, 1.01, 1.02, 1.1, 1.1, 1.1};

    tth_jec = {1.08, 1.02, 1.02, 1.01, 1.03, 1.05, 1.04, 1.01, 1.03, 1.03, 1.08, 1.02, 1.02, 1.01, 1.03, 1.05, 1.04, 1.01, 1.03, 1.03, 1.1, 1.1, 1.1};
    tth_bl = {1.01, 1.01, 1.01, 1.02, 1.02, 1.01, 1.01, 1.01, 1.02, 1.02, 1.01, 1.01, 1.01, 1.02, 1.02, 1.01, 1.01, 1.01, 1.02, 1.02, 1.1, 1.1, 1.1};
    tth_bb = {1.01, 1.01, 1.02, 1.01, 1.02, 1.01, 1.01, 1.02, 1.01, 1.02, 1.01, 1.01, 1.02, 1.01, 1.02, 1.01, 1.01, 1.02, 1.01, 1.02, 1.1, 1.1, 1.1};

    ttx_jec = {1.01, 1.05, 1.05, 1.09, 1.02, 1.03, 1.01, 1.02, 1.01, 1.03, 1.01, 1.05, 1.05, 1.09, 1.02, 1.03, 1.01, 1.02, 1.01, 1.03, 1.1, 1.1, 1.1};
    ttx_bl = {1.01, 1.01, 1.03, 1.03, 1.03, 1.01, 1.01, 1.02, 1.03, 1.03, 1.01, 1.01, 1.03, 1.03, 1.03, 1.01, 1.01, 1.02, 1.03, 1.03, 1.1, 1.1, 1.1};
    ttx_bb = {1.01, 1.01, 1.02, 1.01, 1.02, 1.02, 1.01, 1.02, 1.01, 1.02, 1.01, 1.01, 1.02, 1.01, 1.02, 1.02, 1.01, 1.02, 1.01, 1.02, 1.1, 1.1, 1.1};

    WZbkg = {1.10, 1.10, 1.10, 1.20, 1.20, 1.10, 1.10, 1.10, 1.20, 1.20, 1.10, 1.10, 1.10, 1.20, 1.20, 1.10, 1.10, 1.10, 1.20, 1.20, 1.10, 1.10, 1.20};
  }
 
  int indSR[SRNumber];

  for(int i = 0; i < SRNumber; i++){
    indSR[i] = i+1;
  }

  std::vector<std::vector<std::pair<double, double> >> yields;

  /*
  for(auto & i : distribsOrderForYields){
    cout << i << " ";
  }
  cout << endl;
  */

    for (int i = 0; i != distribsOrderForYields.size()-1; ++i) {

        if(i != 0){
          //cout << "what we push to distribsOrderForDatacards: " << i << endl;
          distribsOrderForDatacards.push_back(i); // = { 1    , 2     , 3       , 4    , 5    , 6    , 7   , 8  , 9};
        }

        vector<pair<double, double> > yieldForEachSR;

        for (int k = 0; k != SRNumber; ++k) {
          double yield = 0;
          double unc = 0;

          for(unsigned int j = distribsOrderForYields.at(i); j != distribsOrderForYields.at(i+1); j++){
            yield += distribs.vectorHisto[j].GetBinContent(indSR[k]);
            unc += TMath::Power(distribs.vectorHisto[j].GetBinError(indSR[k]), 2);

          }
          yieldForEachSR.push_back((make_pair(yield, unc))); // uncertainty on stats    
        }
        
        yields.push_back(yieldForEachSR);
    }

    cout << "number of SR and Categories: " << SRNumber << " " << nCategories << endl;

    // datacards
    double statUnc[SRNumber][nCategories];
    TH1F* hSRyield[nCategories]; 
    for(int yieldNumber = 0; yieldNumber < nCategories; yieldNumber++){

        TString name = Form("hSRyield_%d",yieldNumber);
        
        hSRyield[yieldNumber] = new TH1F(name,name,400,0,400);
        hSRyield[yieldNumber]->Sumw2();
        for (int k=0; k!=SRNumber; ++k) {
            
            hSRyield[yieldNumber]->SetBinContent(k+1,yields.at(yieldNumber).at(k).first);
            hSRyield[yieldNumber]->SetBinError(k+1,TMath::Sqrt(yields.at(yieldNumber).at(k).second));
            
            double uncy = hSRyield[yieldNumber]->GetBinContent(k+1) != 0 ? TMath::Abs(hSRyield[yieldNumber]->GetBinError(k+1) / hSRyield[yieldNumber]->GetBinContent(k+1)) : 1.00;
            //if(yieldNumber != nCategories)
            statUnc[k][yieldNumber] = 1 + uncy;

        }
    }
    
    nameOfProcessesForDatacard.erase(nameOfProcessesForDatacard.begin());
    const int numberOfBKG = nameOfProcessesForDatacard.size() - 1; // -1 for signal and -1 for data
        
    for(int id = 0; id < SRNumber; id++){

        TString datastr;
        if(leptonSelectionAnalysis == 2)
          datastr  = "A" + std::to_string(id); 
        if(leptonSelectionAnalysis == 3)
          datastr  = "B" + std::to_string(id); 

        gSystem->Exec("rm datacards/" + datastr + ".txt"); // delete previous tex file
        ofstream fileout;
        fileout .open ( "datacards/" + datastr + ".txt", ios_base::app); // create a new tex file
        fileout << fixed << showpoint << setprecision(2);
        fileout << "#  " << datastr  << endl;
        fileout << "imax 1 number of channels " <<  endl;
        fileout << "jmax " << numberOfBKG << " number of backgrounds " <<  endl;
        fileout << "kmax " << (leptonSelectionAnalysis == 2 ? 24 : 22) <<  " number of nuisance parameters (sources of systematical uncertainties) " <<  endl;
        fileout << "----------- " <<  endl;
        fileout << "shapes * * FAKE" << endl;
        fileout << "-----------  " <<  endl;
        fileout << "bin          "  << datastr <<  endl;
        fileout << "observation  " << int(hSRyield[0]->GetBinContent(id+1)) << endl;
        fileout << "-----------  " <<  endl;
        fileout << "bin          ";

        vector<TString> datastrVector;
        vector<TString> processNumberVector;
        vector<double> rateVector;

        for(int i = 0; i < nCategories-1; i++){ // here -1 we don't consider data
          datastrVector.push_back(datastr);
          processNumberVector.push_back(std::to_string(i));
          rateVector.push_back((hSRyield[distribsOrderForDatacards.at(i)]->GetBinContent(id+1) > 0 ? hSRyield[distribsOrderForDatacards.at(i)]->GetBinContent(id+1) : 0.0001));
        }
        
        fillString(fileout, datastrVector);

        fileout << "process     " ;
        fillString(fileout, nameOfProcessesForDatacard);
        
        fileout << "process     " ;
        fillString(fileout, processNumberVector);

        fileout << "rate        " ;
        fillString(fileout, rateVector);

        fileout << "----------- " <<  endl;

        
        
        for(int statInd = 0; statInd < numberOfBKG + 1; statInd++){
          fileout << ("st" + nameOfProcessesForDatacard.at(statInd) + datastr + "    lnN  ");
          vector<double> newString = formEmptyString(numberOfBKG + 1);
          newString[statInd] = statUnc[id][distribsOrderForDatacards.at(statInd)];
          fillString(fileout, newString);
        }

        fileout << "----------- " <<  endl;
        fillExperUnc(fileout, nameOfProcessesForDatacard);

        fileout << "JES      lnN  ";
        vector<double> newString = formEmptyString(numberOfBKG + 1);
        for(int i = 0; i <  nameOfProcessesForDatacard.size(); i++){
          if(nameOfProcessesForDatacard.at(i) == "ttW")
            newString[i] = ttw_jec[id];
          if(nameOfProcessesForDatacard.at(i) == "ttZ")
            newString[i] = ttz_jec[id];
          if(nameOfProcessesForDatacard.at(i) == "ttH")
            newString[i] = tth_jec[id];
          if(nameOfProcessesForDatacard.at(i) == "ttX")
            newString[i] = ttx_jec[id];
        }
        fillString(fileout, newString);

        fileout << "btagl      lnN  ";
        newString.clear();
        newString = formEmptyString(numberOfBKG + 1);
        for(int i = 0; i <  nameOfProcessesForDatacard.size(); i++){
          if(nameOfProcessesForDatacard.at(i) == "ttW")
            newString[i] = ttw_bl[id];
          if(nameOfProcessesForDatacard.at(i) == "ttZ")
            newString[i] = ttz_bl[id];
          if(nameOfProcessesForDatacard.at(i) == "ttH")
            newString[i] = tth_bl[id];
          if(nameOfProcessesForDatacard.at(i) == "ttX")
            newString[i] = ttx_bl[id];
        }
        fillString(fileout, newString);

        fileout << "btagb      lnN  ";
        newString.clear();
        newString = formEmptyString(numberOfBKG + 1);
        for(int i = 0; i <  nameOfProcessesForDatacard.size(); i++){
          if(nameOfProcessesForDatacard.at(i) == "ttW")
            newString[i] = ttw_bb[id];
          if(nameOfProcessesForDatacard.at(i) == "ttZ")
            newString[i] = ttz_bb[id];
          if(nameOfProcessesForDatacard.at(i) == "ttH")
            newString[i] = tth_bb[id];
          if(nameOfProcessesForDatacard.at(i) == "ttX")
            newString[i] = ttx_bb[id];
        }
        fillString(fileout, newString);

        fileout << endl;
        
        fileout << "PDF      lnN  ";
        newString = formEmptyString(numberOfBKG + 1);
        for(int i = 0; i < nameOfProcessesForDatacard.size(); i++){
          if(std::find(ttVprocesses.begin(), ttVprocesses.end(), nameOfProcessesForDatacard.at(i)) != ttVprocesses.end() )
            newString[i] = 1.01;
        }
        fillString(fileout, newString);

        fileout << "Q2      lnN  ";
        newString = formEmptyString(numberOfBKG + 1);
        for(int i = 0; i < nameOfProcessesForDatacard.size(); i++){
          if(std::find(ttVprocesses.begin(), ttVprocesses.end(), nameOfProcessesForDatacard.at(i)) != ttVprocesses.end() )
            newString[i] = 1.01;
        }
        fillString(fileout, newString);
        fileout << endl;

        fileout << ("fake     lnN ");
        newString = formEmptyString(numberOfBKG + 1);
        for(int i = 0; i < nameOfProcessesForDatacard.size(); i++){
          if(nameOfProcessesForDatacard.at(i) == "nonpromptData")
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

        fileout << ("ttX     lnN ");
        newString = formEmptyString(numberOfBKG + 1);
        for(int i = 1; i < nameOfProcessesForDatacard.size(); i++){
          if(std::find(ttVprocesses.begin(), ttVprocesses.end(), nameOfProcessesForDatacard.at(i)) != ttVprocesses.end() )
            newString[i] = 1.1;
        }
        fillString(fileout, newString);

        fileout << ("WZ     lnN ");
        newString = formEmptyString(numberOfBKG + 1);
        for(int i = 0; i < nameOfProcessesForDatacard.size(); i++){
          if(nameOfProcessesForDatacard.at(i) == "WZ")
            newString[i] = WZbkg[id];
        }
        fillString(fileout, newString);

        fileout << ("rare    lnN ");
        newString = formEmptyString(numberOfBKG + 1);
        for(int i = 0; i < nameOfProcessesForDatacard.size(); i++){
          if(nameOfProcessesForDatacard.at(i) == "rare")
            newString[i] = 1.5;
        }
        fillString(fileout, newString);

        fileout << endl;
        
        
    }
    
    
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

vector<double> formEmptyString(int numberOfProcesses){
  vector<double> formString;
  formString.clear();
  for(int i = 0; i < numberOfProcesses; i++)
    formString.push_back(1.0);
  return formString;

}


void fillExperUnc(ofstream & file, vector<std::string> & nameOfProcessesForDatacard){
  for(int i = 0; i < experUnc.size(); i++){
    file << experUncName.at(i) << "      lnN  ";
    for(int statInd = 0; statInd < nameOfProcessesForDatacard.size(); statInd++){
      if(nameOfProcessesForDatacard.at(statInd) == "nonpromptData")
        file << "-" << '\t' ;
      else if (nameOfProcessesForDatacard.at(statInd) == "chargeMisID")
        file << "-" << '\t' ;
      else if (nameOfProcessesForDatacard.at(statInd) == "WZ" && experUncName.at(i) == "LeptonId")
        file << "-" << '\t' ;
      else
        file << experUnc.at(i) << '\t' ;
    }
    file << endl;
  }


}


#endif  // fillDatacards