#ifndef filltables_H
#define filltables_H

#include "Output.h"
#include "readTreeSync.h"
#include "TString.h"

#include <iostream>
#include <cstdlib>

using Output::distribs;
using Output::DistribsAll;

using namespace std;

void setUpSystUncForCategory(DistribsAll & distribs, TH1D * histSystAndStatUnc, const int cat);
void setUpSystUncForIntegral(DistribsAll & distribs, TH1D * histSystAndStatUnc, int);
double setUpSystUncForIntegralInCategory(DistribsAll & distribs, int cat);
double setUpSystUncForIntegralTotal(DistribsAll & distribs, int nCategories);

void fillTablesSRTTZ(DistribsAll & distribs, vector<std::string> nameOfProcessesForDatacard, const TString name, int period = 0){
  const int SRNumber = figNames[name].nBins;
  const int nCategories = nameOfProcessesForDatacard.size(); // nameOfProcessesForDatacard - all MC + 1 for data

  TH1D *stackCopy = (TH1D*)(distribs.stack.GetStack()->Last())->Clone("stackCopy");
  TH1D *histSystAndStatUnc = (TH1D*)(distribs.stack.GetStack()->Last())->Clone(Form("histSystAndStatUnc"));
  setUpSystUncForIntegral(distribs, histSystAndStatUnc, nCategories);

  std::vector<std::string> srCategoriesLabel = {
      // here let's split into 4 tables, 3L nb = 0, nb = 1, nb > 1, 4L nj >= 2
      "3 leptons, $\\rm N_{b} = 0$",
      "3 leptons, $\\rm N_{b} = 1$",
      "3 leptons, $\\rm N_{b} \\geq 2$",
      "4 leptons, $\\rm N_{j} \\geq 2$",
  };

  std::vector<std::string> srSubCatLabelWZ3L = {"\\rm $N_{j} = 1$", "\\rm $N_{j} = 2$", "\\rm $N_{j} = 3$", "\\rm $N_{j} > 3$"};
  std::vector<std::string> srSubCatLabelTTZ3L = {"\\rm $N_{j} = 2$", "\\rm $N_{j} = 3$", "\\rm $N_{j} = 4$", "\\rm $N_{j} > 4$"};
  std::vector<std::string> srSubCatLabelTTZ4L = {"\\rm $N_{b} = 0$", "\\rm $N_{b} > 0$"};

  int srCatNumber[5] = {0, 4, 8, 12, 14};

  for(int srCat = 0; srCat < srCategoriesLabel.size(); srCat++){

    ofstream tableBkg;
    tableBkg.open("tables/tableSR" + std::to_string(srCat) + "_" + (TString)(period == 1 ? "2017" : (period == 0 ? "2016" : "comb")) + ".tex");
    tableBkg<<"\\begin{table}\n";
    tableBkg<<"\\begin{adjustbox}{width=1\\textwidth}\n";
    tableBkg<<"\\begin{tabular}{|c|c|c|c|c|}\\hline\n";
    tableBkg << std::fixed << setprecision(1) << "\n";

    tableBkg << " & \\multicolumn{4}{c|}{" << srCategoriesLabel.at(srCat) << "}" ;
    tableBkg << "\\\\ \\hline\n";

    tableBkg << "Process";

    for(const auto & label: srCat == 0 ? srSubCatLabelWZ3L : (srCat == 3 ? srSubCatLabelTTZ4L : srSubCatLabelTTZ3L))
        tableBkg << " & " << label;
    tableBkg << "\\\\ \\hline\n";

    for(int cat = 1; cat < nCategories; cat++){ // cat == 0 is data
        tableBkg << nameOfProcessesForDatacard.at(cat) << " ";
        TH1D *histSystAndStatUncInd = (TH1D*)(distribs.stack.GetStack()->Last())->Clone(Form("histSystAndStatUncInd"));
        histSystAndStatUncInd->Reset("ICE");
        setUpSystUncForCategory(distribs, histSystAndStatUncInd, cat);
        for(int sr = srCatNumber[srCat]; sr < srCatNumber[srCat+1]; sr++){
            tableBkg << " & $ " << distribs.vectorHisto[cat].GetBinContent(sr+1) << " \\pm "  << histSystAndStatUncInd->GetBinError(sr+1) << " $  ";
        }
        tableBkg << "\\\\ \\hline\n";
    }

    tableBkg << "Total";
    for(int sr = srCatNumber[srCat]; sr < srCatNumber[srCat+1]; sr++)
      tableBkg << " & $ " << stackCopy->GetBinContent(sr+1) << " \\pm "  << stackCopy->GetBinContent(sr+1) * histSystAndStatUnc->GetBinError(sr+1) << " $  ";
    tableBkg << "\\\\ \\hline\n";

    tableBkg << "Observed";
    for(int sr = srCatNumber[srCat]; sr < srCatNumber[srCat+1]; sr++)
        tableBkg << " & $ " << int(distribs.vectorHisto[0].GetBinContent(sr+1)) << " $  ";

    tableBkg << " \\\\ \\hline\n";

    tableBkg <<"\\end{tabular}\n";
    tableBkg <<"\\end{adjustbox}\n";
    tableBkg <<"\\end{table}\n";
    tableBkg.close();
  }

  std::cout << "tables DONE" << std::endl;

}

void fillTablesForFlavour(DistribsAll & distribs, vector<std::string> nameOfProcessesForDatacard, const TString name, int period = 0){

  const int SRNumber = figNames[name].nBins;
  const int nCategories = nameOfProcessesForDatacard.size(); // nameOfProcessesForDatacard - all MC + 1 for data

  TH1D *stackCopy = (TH1D*)(distribs.stack.GetStack()->Last())->Clone("stackCopy");
  TH1D *histSystAndStatUnc = (TH1D*)(distribs.stack.GetStack()->Last())->Clone(Form("histSystAndStatUnc"));
  setUpSystUncForIntegral(distribs, histSystAndStatUnc, nCategories);

  ofstream tableBkg;
  tableBkg.open("tables/tableFlavour" + (TString)(period == 1 ? "2017" : (period == 0 ? "2016" : "comb")) + ".tex");
  tableBkg<<"\\begin{table}\n";
  tableBkg<<"\\begin{adjustbox}{width=1\\textwidth}\n";
  tableBkg<<"\\begin{tabular}{|c|c|c|c|c||c|}\\hline\n";
  tableBkg << std::fixed << setprecision(1) << "\n";

  tableBkg << "Process";

  std::vector<std::string> flavourLabel = {

      "$(\\mu)\\mu\\mu\\mu$",
      "$(\\mu)\\mu\\mu$ e",
      "($\\mu$/e)$\\mu$ ee",
      "(e)eee",

  };

  for(const auto & label: flavourLabel)
    tableBkg << " & " << label;
  tableBkg << "\\\\ \\hline\n";

  for(int cat = 1; cat < nCategories; cat++){ // cat == 0 is data
     tableBkg << nameOfProcessesForDatacard.at(cat) << " ";
     TH1D *histSystAndStatUncInd = (TH1D*)(distribs.stack.GetStack()->Last())->Clone(Form("histSystAndStatUncInd"));
     histSystAndStatUncInd->Reset("ICE");
     setUpSystUncForCategory(distribs, histSystAndStatUncInd, cat);
     for(int sr = 0; sr < SRNumber; sr++)
        tableBkg << " & $ " << distribs.vectorHisto[cat].GetBinContent(sr+1) << " \\pm "  << histSystAndStatUncInd->GetBinError(sr+1) << " $  ";
     tableBkg << " & $ " << distribs.vectorHisto[cat].Integral() << " \\pm "  << setUpSystUncForIntegralInCategory(distribs, cat) << " $  ";
     tableBkg << "\\\\ \\hline\n";
  }

  tableBkg << "Total";
  for(int sr = 0; sr < SRNumber; sr++)
     tableBkg << " & $ " << stackCopy->GetBinContent(sr+1) << " \\pm "  << stackCopy->GetBinContent(sr+1) * histSystAndStatUnc->GetBinError(sr+1) << " $  ";
  tableBkg << " & $ " << stackCopy->Integral() << " \\pm "  << setUpSystUncForIntegralTotal(distribs, nCategories) << " $  ";
  tableBkg << "\\\\ \\hline\n";

  tableBkg << "Observed";
  for(int sr = 0; sr < SRNumber; sr++)
     tableBkg << " & $ " << int(distribs.vectorHisto[0].GetBinContent(sr+1)) << " $  ";
  tableBkg << " & $ " << int(distribs.vectorHisto[0].GetBinContent(1) + distribs.vectorHisto[0].GetBinContent(2) + distribs.vectorHisto[0].GetBinContent(3) + distribs.vectorHisto[0].GetBinContent(4)) << " $  ";

  tableBkg << " \\\\ \\hline\n";

  tableBkg <<"\\end{tabular}\n";
  tableBkg <<"\\end{adjustbox}\n";
  tableBkg <<"\\end{table}\n";
  tableBkg.close();

   std::cout << "tables DONE" << std::endl;
}

void setUpSystUncForCategory(DistribsAll & distribs, TH1D * histSystAndStatUnc, const int cat){

    TH1D *stackUncUp[numberOfSyst];
    TH1D *stackUncDown[numberOfSyst];

    for(unsigned int i = 0; i < numberOfSyst; i++){

        stackUncUp[i] = (TH1D*)distribs.vectorHisto[cat].Clone(Form("stackUncUp_%d", i));
        stackUncDown[i] = (TH1D*)distribs.vectorHisto[cat].Clone(Form("stackUncDown_%d", i));

        stackUncUp[i]->Reset("ICE");
        stackUncDown[i]->Reset("ICE");

        stackUncUp[i]->Add((TH1D*)&distribs.vectorHistoUncUp[cat].unc[i]);
        stackUncDown[i]->Add((TH1D*)&distribs.vectorHistoUncDown[cat].unc[i]);
    }

    for(unsigned int i = 0; i < histSystAndStatUnc->GetNbinsX(); i++){
        // content in particular bin in stack
        double stackBinContent = distribs.vectorHisto[cat].GetBinContent(i+1);
        double stackBinError = distribs.vectorHisto[cat].GetBinError(i+1);

        double err = 0.;
        if(stackBinContent != 0.)
            err += TMath::Power(distribs.vectorHisto[cat].GetBinError(i+1), 2);

        for(unsigned int j = 0; j < numberOfSyst; j++){
            // consider largest deviation between the upward and downward variations
            err += TMath::Power(TMath::Max(fabs(stackUncUp[j]->GetBinContent(i+1) - stackBinContent), fabs(stackBinContent - stackUncDown[j]->GetBinContent(i+1))) , 2);
        }
        // if uncertainty is greater than 100% consider 100% uncertainty
        if(stackBinContent != 0.){
            histSystAndStatUnc->SetBinContent(i+1, TMath::Sqrt(err));
            histSystAndStatUnc->SetBinError(i+1, TMath::Sqrt(err));
        }
    }
}

void setUpSystUncForIntegral(DistribsAll & distribs, TH1D * histSystAndStatUnc, int nCategories){

    // histSystAndStatUnc - a histogram with central value at 1 and with applied one of the uncertainties on top: JEC, JES and Uncl
    //TH1D *histSystAndStatUnc = (TH1D*)(distribs.stack.GetStack()->Last())->Clone(Form("histSystAndStatUnc"));
    TH1D *stackCopy = (TH1D*)(distribs.stack.GetStack()->Last())->Clone("stackCopy");
    // stack of MC with varied up and down of 3 different types of uncertainties
    TH1D *stackUncUp[numberOfSyst];
    TH1D *stackUncDown[numberOfSyst];

    for(unsigned int i = 0; i < numberOfSyst; i++){

        stackUncUp[i] = (TH1D*)distribs.vectorHisto[1].Clone(Form("stackUncUp_%d", i));
        stackUncDown[i] = (TH1D*)distribs.vectorHisto[1].Clone(Form("stackUncDown_%d", i));

        stackUncUp[i]->Reset("ICE");
        stackUncDown[i]->Reset("ICE");

        for(unsigned int j = nCategories; j != 0; j--){
            stackUncUp[i]->Add((TH1D*)&distribs.vectorHistoUncUp[j].unc[i]);
            stackUncDown[i]->Add((TH1D*)&distribs.vectorHistoUncDown[j].unc[i]);
        }
    }

    for(unsigned int i = 0; i < histSystAndStatUnc->GetNbinsX(); i++){
        // content in particular bin in stack
        double stackBinContent = ((TH1D*)distribs.stack.GetStack()->Last())->GetBinContent(i+1);
        double stackBinError = ((TH1D*)distribs.stack.GetStack()->Last())->GetBinError(i+1);

        if(stackBinContent != 0.)
            stackCopy->SetBinError(i+1, stackBinError / stackBinContent);
        else
            stackCopy->SetBinError(i+1, 0.);
        double err = TMath::Power(stackCopy->GetBinError(i+1), 2);
        //cout << "stat unc is " << TMath::Sqrt(err) << endl;

        for(unsigned int j = 0; j < numberOfSyst; j++){
            // consider largest deviation between the upward and downward variations
            err += TMath::Power(TMath::Max(stackUncUp[j]->GetBinContent(i+1) - stackBinContent, stackBinContent - stackUncDown[j]->GetBinContent(i+1)) / stackBinContent, 2);
            histSystAndStatUnc->SetBinContent(i+1, 1.);
            // if uncertainty is greater than 100% consider 100% uncertainty
            if(stackBinContent != 0.){
                //cout << "unc after applying " << j << " syst: " << TMath::Sqrt(err) << endl;
                histSystAndStatUnc->SetBinError(i+1, TMath::Sqrt(err) > 1 ? 1. : TMath::Sqrt(err));
            }
            //else
            //    histSystAndStatUnc->SetBinError(i+1, 0.);
        }
    }

}

double setUpSystUncForIntegralTotal(DistribsAll & distribs, int nCategories){

    // stack of MC with varied up and down of 3 different types of uncertainties
    TH1D *stackUncUp[numberOfSyst];
    TH1D *stackUncDown[numberOfSyst];

    for(unsigned int i = 0; i < numberOfSyst; i++){

        stackUncUp[i] = (TH1D*)distribs.vectorHisto[1].Clone(Form("stackUncUp_%d", i));
        stackUncDown[i] = (TH1D*)distribs.vectorHisto[1].Clone(Form("stackUncDown_%d", i));

        stackUncUp[i]->Reset("ICE");
        stackUncDown[i]->Reset("ICE");

        for(unsigned int j = nCategories; j != 0; j--){
            stackUncUp[i]->Add((TH1D*)&distribs.vectorHistoUncUp[j].unc[i]);
            stackUncDown[i]->Add((TH1D*)&distribs.vectorHistoUncDown[j].unc[i]);
        }
    }

    // content in particular bin in stack
    double stackBinContent = ((TH1D*)distribs.stack.GetStack()->Last())->Integral();

    double err = 0.;
    if(stackBinContent != 0.) {
        ((TH1D*)distribs.stack.GetStack()->Last())->IntegralAndError(1, 4, err);
        err = TMath::Power(err, 2);
    }

    for(unsigned int j = 0; j < numberOfSyst; j++){
       // consider largest deviation between the upward and downward variations
       err += TMath::Power(TMath::Max(stackUncUp[j]->Integral() - stackBinContent, stackBinContent - stackUncDown[j]->Integral()), 2);
    }
    return TMath::Sqrt(err);

}

double setUpSystUncForIntegralInCategory(DistribsAll & distribs, int cat){

    TH1D *stackUncUp[numberOfSyst];
    TH1D *stackUncDown[numberOfSyst];

    double err = 0.;
    double stackBinContent = distribs.vectorHisto[cat].Integral();
    if(stackBinContent != 0.){
       distribs.vectorHisto[cat].IntegralAndError(1, 4, err);
       err = TMath::Power(err, 2);
    }

    for(unsigned int i = 0; i < numberOfSyst; i++){
        stackUncUp[i] = (TH1D*)distribs.vectorHistoUncUp[cat].unc[i].Clone(Form("stackUncUp_%d", i));
        stackUncDown[i] = (TH1D*)distribs.vectorHistoUncDown[cat].unc[i].Clone(Form("stackUncDown_%d", i));
        err += TMath::Power(TMath::Max(fabs(stackUncUp[i]->Integral() - stackBinContent), fabs(stackBinContent - stackUncDown[i]->Integral())) , 2);
    }
    return TMath::Sqrt(err);
}

#endif  // filltables
