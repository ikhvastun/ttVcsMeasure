#include "errors.h"
#include "treeReader.h"
#include "mt2_bisect.h"
#include "readTreeSync.h"

#include "Output.h"
using Output::distribs1DForFR;
using Output::distribs2D;
using Output::distribs;
using Output::distribs1DForCT;
using namespace std;

bool comp (const pair<double, int> i, const pair<double, int> j) { return (i.first>j.first); }

double SRID2L (double mvaVL, double chargesLepton) {
    double index = -1;
    int chargesLeptonIndex = (chargesLepton == 1.);
    if(mvaVL < -1) return -1;
    else return floor((mvaVL + 1) / 0.2) + 10 * chargesLeptonIndex;
}

void initdistribs(std::vector<std::string> & namesOfProcesses, const std::string & selection){

    for (std::map<TString, histInfo>::iterator it=figNames.begin(); it!=figNames.end(); ++it){
      if(std::find(listToPrint[selection].begin(), listToPrint[selection].end(), it->first) == listToPrint[selection].end()) continue;

      histInfo hist = it->second;
      TString name = Form("varST_%d",hist.index);
      int i = hist.index;

		// we have distribs object for each process i. 
		// for each such process we have all histograms required for the selection
      for(unsigned int j = 0; j < distribs[i].vectorHisto.size(); j++){
        name = Form("var_%d_%d",i,j);
        distribs[i].vectorHisto[j] = std::move(TH1D(name,name+";",hist.nBins,hist.varMin,hist.varMax));

			// we make histog for each uncertainty up and down. 
			// we acces the ditributions later in code by calling:
			// distribs[i].vectorHistoUncUp[j].unc[k] where i refers to variable/hist, j to process and k to uncertianty
        for(unsigned int k = 0; k < numberOfSyst; k++){
          name = Form("varUp_%d_%d_%d",i,j,k);
          distribs[i].vectorHistoUncUp[j].unc[k] = std::move(TH1D(name,name+";",hist.nBins,hist.varMin,hist.varMax));

          name = Form("varDown_%d_%d_%d",i,j,k);
          distribs[i].vectorHistoUncDown[j].unc[k] = std::move(TH1D(name,name+";",hist.nBins,hist.varMin,hist.varMax));
        }

        // pdf uncertainties, taes quite a lot for time to initialize them, for convenience atm are commented 
        if(i == indexSR3L || i == indexSR4L || i == indexSRTTZ || i == indexSRTTCR || i == indexSRWZCR || i == indexSRZZCR){
            for(unsigned int pdf = 0; pdf < 100; pdf++){
                name = Form("pdf_%d_%d_%d",i,j,pdf);
                distribs[i].vectorHistoPDF[j].var[pdf] = std::move(TH1D(name,name+";",hist.nBins,hist.varMin,hist.varMax));
            }
        }
        
		  // we set for all processes poissonian errors.

        distribs[i].vectorHisto[j].SetBinErrorOption(TH1::kPoisson);

        distribs[i].vectorHisto[j].SetMarkerStyle(20);
        distribs[i].vectorHisto[j].SetMarkerSize(0.5);
        distribs[i].vectorHisto[j].SetLineWidth(1);
        if (j < nProcesses)
          distribs[i].vectorHisto[j].Sumw2();
      }
    }

    for (unsigned int i=0; i!=distribs.size();++i) {
      for(unsigned int j = namesOfProcesses.size()-1; j != 0; j--){
//        if(selection == "ttW" && j == ttWSample) continue;
        distribs[i].stack.Add(&distribs[i].vectorHisto[j]);
      }
    }
}

void initdistribsForCT(){

    for(unsigned int i = 0; i < distribs1DForCT.size(); i++){
      TString name = Form("varST_%d",i);

      for(unsigned int j = 0; j < distribs1DForCT[i].vectorHisto.size(); j++){
        name = Form("var_%d_%d",i,j);
        distribs1DForFR[i].vectorHisto[j] = std::move(TH1D(name,name+";",nPt-1, ptBins));

        distribs1DForFR[i].vectorHisto[j].SetBinErrorOption(TH1::kPoisson);

        distribs1DForFR[i].vectorHisto[j].SetMarkerStyle(20);
        distribs1DForFR[i].vectorHisto[j].SetMarkerSize(0.5);
        distribs1DForFR[i].vectorHisto[j].SetLineWidth(1);
        distribs1DForFR[i].vectorHisto[j].Sumw2();
      }
    }
}

void initdistribsForFR(){

    for(unsigned int i = 0; i < distribs1DForFR.size(); i++){
      TString name = Form("varST_%d",i);

      for(unsigned int j = 0; j < distribs1DForFR[i].vectorHisto.size(); j++){
        name = Form("var_%d_%d",i,j);
        distribs1DForFR[i].vectorHisto[j] = std::move(TH1D(name,name+";",nPt-1, ptBins));

        distribs1DForFR[i].vectorHisto[j].SetBinErrorOption(TH1::kPoisson);

        distribs1DForFR[i].vectorHisto[j].SetMarkerStyle(20);
        distribs1DForFR[i].vectorHisto[j].SetMarkerSize(0.5);
        distribs1DForFR[i].vectorHisto[j].SetLineWidth(1);
        distribs1DForFR[i].vectorHisto[j].Sumw2();
      }
    }

    for(unsigned int i = 0; i < distribs2D.size(); i++){
      for(unsigned int j = 0; j < distribs2D[i].vectorHisto.size(); j++){

        TString name = Form("FRflavour_%d_%d", i, j);
        distribs2D[i].vectorHisto[j] = std::move(TH2D(name + (TString)(j%2 == 1 ? "passed" : ""),name+";",nPt-1, ptBins, nEta-1, etaBins[i%2]));
        distribs2D[i].vectorHisto[j].Sumw2();

      }
    }
}

    
void setLabelsForHistos(const std::string & selection){

    std::vector<std::vector<BinLabelOptions>> labelVector = {flavourLabelOptionsFor3L, flavourLabelOptionsFor4L, flavourLabelOptionsFor4LZZ, theSRLabelOptionsForWZCR, theSRLabelOptionsForZZCR, theSRLabelOptionsForTTCR, theSRLabelOptionsFor3L, theSRLabelOptionsFor4L, theSRLabelOptionsForTTZ, theSRLabelOptionsFor3L, theSRLabelOptionsFor3L, theSRLabelOptionsFor3L, theSRLabelOptionsFor3L, flavourLabelOptionsFor3L4L, theSRLabelOptionsForTTZ8SR3L, theSRLabelOptionsForttZCleanPTZ, theSRLabelOptionsForttZCleanCosTheta};
    std::vector<unsigned int> indices = {indexFlavour3L, indexFlavour4L, indexFlavour4LZZ, indexSRWZCR, indexSRZZCR, indexSRTTCR, indexSR3L, indexSR4L, indexSRTTZ, indexSR3L3m, indexSR3L2m1e, indexSR3L1m2e, indexSR3L3e, indexFlavour3L4L, indexSRTTZ8SR3L, indexSRttZcleanPTZ, indexSRttZcleanCosTheta};
    std::vector<TString> namesOfIndices = {"flavour3L", "flavour4L", "flavour4LZZ", "SRWZCR", "SRZZCR", "SRTTCR", "SR3L", "SR4L", "SRallTTZ", "SR3L3m", "SR3L2m1e", "SR3L1m2e", "SR3L3e", "flavour3L4L", "SRTTZ8SR3L", "SRttZCleanPTZ", "SRttZCleanCosTheta"};
    for(int ind = 0; ind < labelVector.size(); ind++){
        if(std::find(listToPrint[selection].begin(), listToPrint[selection].end(), namesOfIndices.at(ind)) == listToPrint[selection].end()) continue;
        for(auto & histo: distribs[indices.at(ind)].vectorHisto) {
            for(const auto & i: labelVector.at(ind)){
                histo.GetXaxis()->SetBinLabel(i.index, i.labelSR.c_str());
            }

            histo.GetXaxis()->SetLabelSize(0.1);
            histo.GetXaxis()->SetTitleSize(0.25);
            histo.GetXaxis()->SetLabelOffset(0.02);
        }

        for(auto & histo: distribs[indices.at(ind)].vectorHistoUncUp) 
            for(const auto & i: labelVector.at(ind))
                for(unsigned int k = 0; k < numberOfSyst; k++)
                    histo.unc[k].GetXaxis()->SetBinLabel(i.index, i.labelSR.c_str());

        for(auto & histo: distribs[indices.at(ind)].vectorHistoUncDown) 
            for(const auto & i: labelVector.at(ind))
                for(unsigned int k = 0; k < numberOfSyst; k++)
                    histo.unc[k].GetXaxis()->SetBinLabel(i.index, i.labelSR.c_str());
    }

}

double flavourCategory2L(int nLocEle, int chargesTwoLepton){
  return double(1 + nLocEle + 3 * (chargesTwoLepton > 0 ? 1 : 0));
}

double flavourCategory3L4L(int lepSel, int nLocEle){
  if(lepSel == 3)
    return double(1 + nLocEle);
  else if(lepSel == 4){
    int elCat = nLocEle < 2 ? nLocEle : (nLocEle == 4 ? 3 : 2);
    return 1 + elCat;
  }
  return -999.;
}

double flavourCategory3L(int nLocEle){
  return double(1 + nLocEle);
}

double flavourCategory4L(int nLocEle){
  return double(1 + nLocEle);
}

double flavourCategory4LZZ(int nLocEle){
  return double(1 + nLocEle / 2);
}

void addBranchToNNTreeVariables(){

    signalTree->Branch("_jetPt1", &_jetPt1, "_jetPt1/D");
    signalTree->Branch("_jetEta1", &_jetEta1, "_jetEta1/D");
    signalTree->Branch("_jetPhi1", &_jetPhi1, "_jetPhi1/D");
    signalTree->Branch("_jetE1", &_jetE1, "_jetE1/D");
    signalTree->Branch("_jetCSV1", &_jetCSV1, "_jetCSV1/D");

    signalTree->Branch("_jetPt2", &_jetPt2, "_jetPt2/D");
    signalTree->Branch("_jetEta2", &_jetEta2, "_jetEta2/D");
    signalTree->Branch("_jetPhi2", &_jetPhi2, "_jetPhi2/D");
    signalTree->Branch("_jetE2", &_jetE2, "_jetE2/D");
    signalTree->Branch("_jetCSV2", &_jetCSV2, "_jetCSV2/D");

    signalTree->Branch("_jetPt3", &_jetPt3, "_jetPt3/D");
    signalTree->Branch("_jetEta3", &_jetEta3, "_jetEta3/D");
    signalTree->Branch("_jetPhi3", &_jetPhi3, "_jetPhi3/D");
    signalTree->Branch("_jetE3", &_jetE3, "_jetE3/D");
    signalTree->Branch("_jetCSV3", &_jetCSV3, "_jetCSV3/D");

    signalTree->Branch("_jetPt4", &_jetPt4, "_jetPt4/D");
    signalTree->Branch("_jetEta4", &_jetEta4, "_jetEta4/D");
    signalTree->Branch("_jetPhi4", &_jetPhi4, "_jetPhi4/D");
    signalTree->Branch("_jetE4", &_jetE4, "_jetE4/D");
    signalTree->Branch("_jetCSV4", &_jetCSV4, "_jetCSV4/D");

    signalTree->Branch("_jetPt5", &_jetPt5, "_jetPt5/D");
    signalTree->Branch("_jetEta5", &_jetEta5, "_jetEta5/D");
    signalTree->Branch("_jetPhi5", &_jetPhi5, "_jetPhi5/D");
    signalTree->Branch("_jetE5", &_jetE5, "_jetE5/D");
    signalTree->Branch("_jetCSV5", &_jetCSV5, "_jetCSV5/D");

    signalTree->Branch("_jetPt6", &_jetPt6, "_jetPt6/D");
    signalTree->Branch("_jetEta6", &_jetEta6, "_jetEta6/D");
    signalTree->Branch("_jetPhi6", &_jetPhi6, "_jetPhi6/D");
    signalTree->Branch("_jetE6", &_jetE6, "_jetE6/D");
    signalTree->Branch("_jetCSV6", &_jetCSV6, "_jetCSV6/D");

    signalTree->Branch("_lepPt1", &_lepPt1, "_lepPt1/D");
    signalTree->Branch("_lepEta1", &_lepEta1, "_lepEta1/D");
    signalTree->Branch("_lepPhi1", &_lepPhi1, "_lepPhi1/D");
    signalTree->Branch("_lepE1", &_lepE1, "_lepE1/D");
    signalTree->Branch("_lepCharge1", &_lepCharge1, "_lepCharge1/D");

    signalTree->Branch("_lepPt2", &_lepPt2, "_lepPt2/D");
    signalTree->Branch("_lepEta2", &_lepEta2, "_lepEta2/D");
    signalTree->Branch("_lepPhi2", &_lepPhi2, "_lepPhi2/D");
    signalTree->Branch("_lepE2", &_lepE2, "_lepE2/D");
    signalTree->Branch("_lepCharge2", &_lepCharge2, "_lepCharge2/D");

    signalTree->Branch("_metPt1", &_metPt1, "_metPt1/D");
    signalTree->Branch("_metEta1", &_metEta1, "_metEta1/D");
    signalTree->Branch("_metPhi1", &_metPhi1, "_metPhi1/D");
    signalTree->Branch("_metE1", &_metE1, "_metE1/D");

    signalTree->Branch("nJLoc", &nJLoc, "nJLoc/I");
    signalTree->Branch("nBLoc", &nBLoc, "nBLoc/I");
    signalTree->Branch("HTLoc", &HTLoc, "HTLoc/D");
    signalTree->Branch("_met", &MET, "_met/D");

    signalTree->Branch("minDeltaR", &minDeltaR, "minDeltaR/D");
    signalTree->Branch("mt", &mtHighest, "mt/D");

    signalTree->Branch("leadpt", &leadpt, "leadpt/D");
    signalTree->Branch("trailpt", &trailpt, "trailpt/D");
    signalTree->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/D");
    signalTree->Branch("trailJetPt", &trailJetPt, "trailJetPt/D");
    signalTree->Branch("chargeOfLeptons", &chargeOfLeptons, "chargeOfLeptons/I");
    signalTree->Branch("mll_ss", &mll_ss, "mll_ss/D");
    signalTree->Branch("ll_deltaR", &ll_deltaR, "ll_deltaR/D");
    signalTree->Branch("mt2ll_ss", &mt2ll_ss, "mt2ll_ss/D");

    signalTree->Branch("maxMLeptonJet", &maxMLeptonJet, "maxMLeptonJet/D");
    signalTree->Branch("maxpTLeptonJet", &maxpTLeptonJet, "maxpTLeptonJet/D");

    signalTree->Branch("maxMLeptonbJet", &maxMLeptonbJet, "maxMLeptonbJet/D");
    signalTree->Branch("maxpTLeptonbJet", &maxpTLeptonbJet, "maxpTLeptonbJet/D");

    signalTree->Branch("minMJetJet", &minMJetJet, "minMJetJet/D");
    signalTree->Branch("maxMJetJet", &maxMJetJet, "maxMJetJet/D");
    signalTree->Branch("minDeltaRJetJet", &minDeltaRJetJet, "minDeltaRJetJet/D");
    signalTree->Branch("maxDeltaRJetJet", &maxDeltaRJetJet, "maxDeltaRJetJet/D");
    signalTree->Branch("minDeltaPhiJetJet", &minDeltaPhiJetJet, "minDeltaPhiJetJet/D");
    signalTree->Branch("maxDeltaPhiJetJet", &maxDeltaPhiJetJet, "maxDeltaPhiJetJet/D");
    signalTree->Branch("maxpTJetJet", &maxpTJetJet, "maxpTJetJet/D");

    signalTree->Branch("maxmTLeptonMET", &maxmTLeptonMET, "maxmTLeptonMET/D");
    signalTree->Branch("minpTLeptonMET", &minpTLeptonMET, "minpTLeptonMET/D");
    signalTree->Branch("maxpTLeptonMET", &maxpTLeptonMET, "maxpTLeptonMET/D");

    signalTree->Branch("maxmTJetMET", &maxmTJetMET, "maxmTJetMET/D");
    signalTree->Branch("maxpTJetMET", &maxpTJetMET, "maxpTJetMET/D");

    signalTree->Branch("maxmTBJetMET", &maxmTBJetMET, "maxmTBJetMET/D");
    signalTree->Branch("maxpTBJetMET", &maxpTBJetMET, "maxpTBJetMET/D");

    signalTree->Branch("_weight", &_weightEventInTree, "_weight/D");

    bkgTree->Branch("_jetPt1", &_jetPt1, "_jetPt1/D");
    bkgTree->Branch("_jetEta1", &_jetEta1, "_jetEta1/D");
    bkgTree->Branch("_jetPhi1", &_jetPhi1, "_jetPhi1/D");
    bkgTree->Branch("_jetE1", &_jetE1, "_jetE1/D");
    bkgTree->Branch("_jetCSV1", &_jetCSV1, "_jetCSV1/D");

    bkgTree->Branch("_jetPt2", &_jetPt2, "_jetPt2/D");
    bkgTree->Branch("_jetEta2", &_jetEta2, "_jetEta2/D");
    bkgTree->Branch("_jetPhi2", &_jetPhi2, "_jetPhi2/D");
    bkgTree->Branch("_jetE2", &_jetE2, "_jetE2/D");
    bkgTree->Branch("_jetCSV2", &_jetCSV2, "_jetCSV2/D");

    bkgTree->Branch("_jetPt3", &_jetPt3, "_jetPt3/D");
    bkgTree->Branch("_jetEta3", &_jetEta3, "_jetEta3/D");
    bkgTree->Branch("_jetPhi3", &_jetPhi3, "_jetPhi3/D");
    bkgTree->Branch("_jetE3", &_jetE3, "_jetE3/D");
    bkgTree->Branch("_jetCSV3", &_jetCSV3, "_jetCSV3/D");

    bkgTree->Branch("_jetPt4", &_jetPt4, "_jetPt4/D");
    bkgTree->Branch("_jetEta4", &_jetEta4, "_jetEta4/D");
    bkgTree->Branch("_jetPhi4", &_jetPhi4, "_jetPhi4/D");
    bkgTree->Branch("_jetE4", &_jetE4, "_jetE4/D");
    bkgTree->Branch("_jetCSV4", &_jetCSV4, "_jetCSV4/D");

    bkgTree->Branch("_jetPt5", &_jetPt5, "_jetPt5/D");
    bkgTree->Branch("_jetEta5", &_jetEta5, "_jetEta5/D");
    bkgTree->Branch("_jetPhi5", &_jetPhi5, "_jetPhi5/D");
    bkgTree->Branch("_jetE5", &_jetE5, "_jetE5/D");
    bkgTree->Branch("_jetCSV5", &_jetCSV5, "_jetCSV5/D");

    bkgTree->Branch("_jetPt6", &_jetPt6, "_jetPt6/D");
    bkgTree->Branch("_jetEta6", &_jetEta6, "_jetEta6/D");
    bkgTree->Branch("_jetPhi6", &_jetPhi6, "_jetPhi6/D");
    bkgTree->Branch("_jetE6", &_jetE6, "_jetE6/D");
    bkgTree->Branch("_jetCSV6", &_jetCSV6, "_jetCSV6/D");

    bkgTree->Branch("_lepPt1", &_lepPt1, "_lepPt1/D");
    bkgTree->Branch("_lepEta1", &_lepEta1, "_lepEta1/D");
    bkgTree->Branch("_lepPhi1", &_lepPhi1, "_lepPhi1/D");
    bkgTree->Branch("_lepE1", &_lepE1, "_lepE1/D");
    bkgTree->Branch("_lepCharge1", &_lepCharge1, "_lepCharge1/D");

    bkgTree->Branch("_lepPt2", &_lepPt2, "_lepPt2/D");
    bkgTree->Branch("_lepEta2", &_lepEta2, "_lepEta2/D");
    bkgTree->Branch("_lepPhi2", &_lepPhi2, "_lepPhi2/D");
    bkgTree->Branch("_lepE2", &_lepE2, "_lepE2/D");
    bkgTree->Branch("_lepCharge2", &_lepCharge2, "_lepCharge2/D");

    bkgTree->Branch("_metPt1", &_metPt1, "_metPt1/D");
    bkgTree->Branch("_metEta1", &_metEta1, "_metEta1/D");
    bkgTree->Branch("_metPhi1", &_metPhi1, "_metPhi1/D");
    bkgTree->Branch("_metE1", &_metE1, "_metE1/D");

    bkgTree->Branch("nJLoc", &nJLoc, "nJLoc/I");
    bkgTree->Branch("nBLoc", &nBLoc, "nBLoc/I");
    bkgTree->Branch("HTLoc", &HTLoc, "HTLoc/D");
    bkgTree->Branch("_met", &MET, "_met/D");

    bkgTree->Branch("minDeltaR", &minDeltaR, "minDeltaR/D");
    bkgTree->Branch("mt", &mtHighest, "mt/D");

    bkgTree->Branch("leadpt", &leadpt, "leadpt/D");
    bkgTree->Branch("trailpt", &trailpt, "trailpt/D");
    bkgTree->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/D");
    bkgTree->Branch("trailJetPt", &trailJetPt, "trailJetPt/D");
    bkgTree->Branch("chargeOfLeptons", &chargeOfLeptons, "chargeOfLeptons/I");
    bkgTree->Branch("mll_ss", &mll_ss, "mll_ss/D");
    bkgTree->Branch("ll_deltaR", &ll_deltaR, "ll_deltaR/D");
    bkgTree->Branch("mt2ll_ss", &mt2ll_ss, "mt2ll_ss/D");

    bkgTree->Branch("maxMLeptonJet", &maxMLeptonJet, "maxMLeptonJet/D");
    bkgTree->Branch("maxpTLeptonJet", &maxpTLeptonJet, "maxpTLeptonJet/D");

    bkgTree->Branch("maxMLeptonbJet", &maxMLeptonbJet, "maxMLeptonbJet/D");
    bkgTree->Branch("maxpTLeptonbJet", &maxpTLeptonbJet, "maxpTLeptonbJet/D");

    bkgTree->Branch("minMJetJet", &minMJetJet, "minMJetJet/D");
    bkgTree->Branch("maxMJetJet", &maxMJetJet, "maxMJetJet/D");
    bkgTree->Branch("minDeltaRJetJet", &minDeltaRJetJet, "minDeltaRJetJet/D");
    bkgTree->Branch("maxDeltaRJetJet", &maxDeltaRJetJet, "maxDeltaRJetJet/D");
    bkgTree->Branch("minDeltaPhiJetJet", &minDeltaPhiJetJet, "minDeltaPhiJetJet/D");
    bkgTree->Branch("maxDeltaPhiJetJet", &maxDeltaPhiJetJet, "maxDeltaPhiJetJet/D");
    bkgTree->Branch("maxpTJetJet", &maxpTJetJet, "maxpTJetJet/D");

    bkgTree->Branch("maxmTLeptonMET", &maxmTLeptonMET, "maxmTLeptonMET/D");
    bkgTree->Branch("minpTLeptonMET", &minpTLeptonMET, "minpTLeptonMET/D");
    bkgTree->Branch("maxpTLeptonMET", &maxpTLeptonMET, "maxpTLeptonMET/D");

    bkgTree->Branch("maxmTJetMET", &maxmTJetMET, "maxmTJetMET/D");
    bkgTree->Branch("maxpTJetMET", &maxpTJetMET, "maxpTJetMET/D");

    bkgTree->Branch("maxmTBJetMET", &maxmTBJetMET, "maxmTBJetMET/D");
    bkgTree->Branch("maxpTBJetMET", &maxpTBJetMET, "maxpTBJetMET/D");

    bkgTree->Branch("_weight", &_weightEventInTree, "_weight/D");

}

void addBranchToBDTTreeVariables(){
    signalTree->Branch("nJLoc", &nJLoc, "nJLoc/I");
    signalTree->Branch("nBLoc", &nBLoc, "nBLoc/I");
    signalTree->Branch("HTLoc", &HTLoc, "HTLoc/D");
    signalTree->Branch("_met", &MET, "_met/D");

    signalTree->Branch("_weight", &_weightEventInTree, "_weight/D");

    signalTree->Branch("minDeltaRlead", &minDeltaRlead, "minDeltaRlead/D");
    signalTree->Branch("minDeltaR", &minDeltaR, "minDeltaR/D");
    signalTree->Branch("mt", &mtHighest, "mt/D");
    signalTree->Branch("mtlow", &mtLowest, "mtlow/D");

    signalTree->Branch("leadpt", &leadpt, "leadpt/D");
    signalTree->Branch("trailpt", &trailpt, "trailpt/D");
    signalTree->Branch("leadeta", &leadeta, "leadeta/D");
    signalTree->Branch("traileta", &traileta, "traileta/D");
    signalTree->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/D");
    signalTree->Branch("trailJetPt", &trailJetPt, "trailJetPt/D");
    signalTree->Branch("chargeOfLeptons", &chargeOfLeptons, "chargeOfLeptons/I");
    signalTree->Branch("mll_ss", &mll_ss, "mll_ss/D");
    signalTree->Branch("ll_deltaR", &ll_deltaR, "ll_deltaR/D");
    signalTree->Branch("mt2ll_ss", &mt2ll_ss, "mt2ll_ss/D");

    signalTree->Branch("maxMLeptonJet", &maxMLeptonJet, "maxMLeptonJet/D");
    signalTree->Branch("maxpTLeptonJet", &maxpTLeptonJet, "maxpTLeptonJet/D");

    signalTree->Branch("maxMLeptonbJet", &maxMLeptonbJet, "maxMLeptonbJet/D");
    signalTree->Branch("maxpTLeptonbJet", &maxpTLeptonbJet, "maxpTLeptonbJet/D");

    signalTree->Branch("minMJetJet", &minMJetJet, "minMJetJet/D");
    signalTree->Branch("maxMJetJet", &maxMJetJet, "maxMJetJet/D");
    signalTree->Branch("minDeltaRJetJet", &minDeltaRJetJet, "minDeltaRJetJet/D");
    signalTree->Branch("maxDeltaRJetJet", &maxDeltaRJetJet, "maxDeltaRJetJet/D");
    signalTree->Branch("minDeltaPhiJetJet", &minDeltaPhiJetJet, "minDeltaPhiJetJet/D");
    signalTree->Branch("maxDeltaPhiJetJet", &maxDeltaPhiJetJet, "maxDeltaPhiJetJet/D");
    signalTree->Branch("maxpTJetJet", &maxpTJetJet, "maxpTJetJet/D");

    signalTree->Branch("minDeltaPhiLeptonMET", &minDeltaPhiLeptonMET, "minDeltaPhiLeptonMET/D");
    signalTree->Branch("maxDeltaPhiLeptonMET", &maxDeltaPhiLeptonMET, "maxDeltaPhiLeptonMET/D");
    signalTree->Branch("minmTLeptonMET", &minmTLeptonMET, "minmTLeptonMET/D");
    signalTree->Branch("maxmTLeptonMET", &maxmTLeptonMET, "maxmTLeptonMET/D");
    signalTree->Branch("minpTLeptonMET", &minpTLeptonMET, "minpTLeptonMET/D");
    signalTree->Branch("maxpTLeptonMET", &maxpTLeptonMET, "maxpTLeptonMET/D");

    signalTree->Branch("minDeltaPhiJetMET", &minDeltaPhiJetMET, "minDeltaPhiJetMET/D");
    signalTree->Branch("maxDeltaPhiJetMET", &maxDeltaPhiJetMET, "maxDeltaPhiJetMET/D");
    signalTree->Branch("minmTJetMET", &minmTJetMET, "minmTJetMET/D");
    signalTree->Branch("maxmTJetMET", &maxmTJetMET, "maxmTJetMET/D");
    signalTree->Branch("minpTJetMET", &minpTJetMET, "minpTJetMET/D");
    signalTree->Branch("maxpTJetMET", &maxpTJetMET, "maxpTJetMET/D");

    signalTree->Branch("minDeltaPhiBJetMET", &minDeltaPhiBJetMET, "minDeltaPhiBJetMET/D");
    signalTree->Branch("maxDeltaPhiBJetMET", &maxDeltaPhiBJetMET, "maxDeltaPhiBJetMET/D");
    signalTree->Branch("minmTBJetMET", &minmTBJetMET, "minmTBJetMET/D");
    signalTree->Branch("maxmTBJetMET", &maxmTBJetMET, "maxmTBJetMET/D");
    signalTree->Branch("minpTBJetMET", &minpTBJetMET, "minpTBJetMET/D");
    signalTree->Branch("maxpTBJetMET", &maxpTBJetMET, "maxpTBJetMET/D");

    bkgTree->Branch("nJLoc", &nJLoc, "nJLoc/I");
    bkgTree->Branch("nBLoc", &nBLoc, "nBLoc/I");
    bkgTree->Branch("HTLoc", &HTLoc, "HTLoc/D");
    bkgTree->Branch("_met", &MET, "_met/D");

    bkgTree->Branch("_weight", &_weightEventInTree, "_weight/D");

    bkgTree->Branch("minDeltaRlead", &minDeltaRlead, "minDeltaRlead/D");
    bkgTree->Branch("minDeltaR", &minDeltaR, "minDeltaR/D");
    bkgTree->Branch("mt", &mtHighest, "mt/D");
    bkgTree->Branch("mtlow", &mtLowest, "mtlow/D");

    bkgTree->Branch("leadpt", &leadpt, "leadpt/D");
    bkgTree->Branch("trailpt", &trailpt, "trailpt/D");
    bkgTree->Branch("leadeta", &leadeta, "leadeta/D");
    bkgTree->Branch("traileta", &traileta, "traileta/D");
    bkgTree->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/D");  
    bkgTree->Branch("trailJetPt", &trailJetPt, "trailJetPt/D");
    bkgTree->Branch("chargeOfLeptons", &chargeOfLeptons, "chargeOfLeptons/I");
    bkgTree->Branch("mll_ss", &mll_ss, "mll_ss/D");
    bkgTree->Branch("ll_deltaR", &ll_deltaR, "ll_deltaR/D");
    bkgTree->Branch("mt2ll_ss", &mt2ll_ss, "mt2ll_ss/D");

    bkgTree->Branch("maxMLeptonJet", &maxMLeptonJet, "maxMLeptonJet/D");
    bkgTree->Branch("maxpTLeptonJet", &maxpTLeptonJet, "maxpTLeptonJet/D");

    bkgTree->Branch("maxMLeptonbJet", &maxMLeptonbJet, "maxMLeptonbJet/D");
    bkgTree->Branch("maxpTLeptonbJet", &maxpTLeptonbJet, "maxpTLeptonbJet/D");

    bkgTree->Branch("minMJetJet", &minMJetJet, "minMJetJet/D");
    bkgTree->Branch("maxMJetJet", &maxMJetJet, "maxMJetJet/D");
    bkgTree->Branch("minDeltaRJetJet", &minDeltaRJetJet, "minDeltaRJetJet/D");
    bkgTree->Branch("maxDeltaRJetJet", &maxDeltaRJetJet, "maxDeltaRJetJet/D");
    bkgTree->Branch("minDeltaPhiJetJet", &minDeltaPhiJetJet, "minDeltaPhiJetJet/D");
    bkgTree->Branch("maxDeltaPhiJetJet", &maxDeltaPhiJetJet, "maxDeltaPhiJetJet/D");
    bkgTree->Branch("maxpTJetJet", &maxpTJetJet, "maxpTJetJet/D");

    bkgTree->Branch("minDeltaPhiLeptonMET", &minDeltaPhiLeptonMET, "minDeltaPhiLeptonMET/D");
    bkgTree->Branch("maxDeltaPhiLeptonMET", &maxDeltaPhiLeptonMET, "maxDeltaPhiLeptonMET/D");
    bkgTree->Branch("minmTLeptonMET", &minmTLeptonMET, "minmTLeptonMET/D");
    bkgTree->Branch("maxmTLeptonMET", &maxmTLeptonMET, "maxmTLeptonMET/D");
    bkgTree->Branch("minpTLeptonMET", &minpTLeptonMET, "minpTLeptonMET/D");
    bkgTree->Branch("maxpTLeptonMET", &maxpTLeptonMET, "maxpTLeptonMET/D");

    bkgTree->Branch("minDeltaPhiJetMET", &minDeltaPhiJetMET, "minDeltaPhiJetMET/D");
    bkgTree->Branch("maxDeltaPhiJetMET", &maxDeltaPhiJetMET, "maxDeltaPhiJetMET/D");
    bkgTree->Branch("minmTJetMET", &minmTJetMET, "minmTJetMET/D");
    bkgTree->Branch("maxmTJetMET", &maxmTJetMET, "maxmTJetMET/D");
    bkgTree->Branch("minpTJetMET", &minpTJetMET, "minpTJetMET/D");
    bkgTree->Branch("maxpTJetMET", &maxpTJetMET, "maxpTJetMET/D");

    bkgTree->Branch("minDeltaPhiBJetMET", &minDeltaPhiBJetMET, "minDeltaPhiBJetMET/D");
    bkgTree->Branch("maxDeltaPhiBJetMET", &maxDeltaPhiBJetMET, "maxDeltaPhiBJetMET/D");
    bkgTree->Branch("minmTBJetMET", &minmTBJetMET, "minmTBJetMET/D");
    bkgTree->Branch("maxmTBJetMET", &maxmTBJetMET, "maxmTBJetMET/D");
    bkgTree->Branch("minpTBJetMET", &minpTBJetMET, "minpTBJetMET/D");
    bkgTree->Branch("maxpTBJetMET", &maxpTBJetMET, "maxpTBJetMET/D");

}

/*
void addVariablesToBDT(const bool is2017 = false){

    readerTTWcsttbar->AddVariable( "nJLoc", &usernJLoc );
    readerTTWcsttbar->AddVariable( "nBLoc", &usernBLoc );
    readerTTWcsttbar->AddVariable( "HTLoc", &userHTLoc ); 
    readerTTWcsttbar->AddVariable( "_met", &user_met );
    readerTTWcsttbar->AddVariable( "minDeltaRlead", &userminDeltaRlead );
    readerTTWcsttbar->AddVariable( "minDeltaR", &userminDeltaR );
    readerTTWcsttbar->AddVariable( "ll_deltaR", &userll_deltaR );
    readerTTWcsttbar->AddVariable( "mt", &usermt );
    readerTTWcsttbar->AddVariable( "mtlow", &usermtlow );
    readerTTWcsttbar->AddVariable( "leadpt", &userleadpt );
    readerTTWcsttbar->AddVariable( "trailpt", &usertrailpt );
    readerTTWcsttbar->AddVariable( "leadingJetPt", &userleadingjetpt );
    readerTTWcsttbar->AddVariable( "trailJetPt", &usertrailjetpt );  
    readerTTWcsttbar->AddVariable( "leadeta", &userleadeta );
    readerTTWcsttbar->AddVariable( "traileta", &usertraileta );
    //readerTTWcsttbar->AddVariable( "chargeOfLeptons", &userchargeOfLeptons);
    readerTTWcsttbar->AddVariable( "mll_ss", &usermll_ss );
    readerTTWcsttbar->AddVariable( "mt2ll_ss", &usermt2ll_ss );

    // the one used for leptonMVA
    // obtained from TMVA training
    TString dir    = "MVAtrainings/" + (TString)(is2017 ? "2017MC" : "2016MC") + "/ttVvsNPchargeMisIDplusDDnonprompt/dataset/weights/"; 

    // used for cut based
    //TString dir    = "/Users/illiakhvastunov/Desktop/CERN/ss2l_2016_fulldataset/analysis/ttWvsttbar_MC_newJEC_fixed/weights/";
    
    TString prefix = "TMVAClassification";
      
    TString methodName = TString("BDTG") + TString(" method");
    TString weightfile = dir + prefix + TString("_") + TString("BDTG") + TString(".weights.xml");
    //TString methodName = TString("BDT::BDT");
    //TString weightfile = TString("MVAtrainings/2017MC/ttVvsNPchargeMisID/trainingInSKlearn/bdt.weights.xml");
    readerTTWcsttbar->BookMVA( methodName, weightfile ); 
}

void fillBDTvariables(vector<Float_t> & varForBDT){

    usernJLoc = varForBDT.at(0);
    usernBLoc = varForBDT.at(1);
    userHTLoc = varForBDT.at(2);
    user_met = varForBDT.at(3);
    userminDeltaRlead = varForBDT.at(4);
    userminDeltaR = varForBDT.at(5);
    userll_deltaR = varForBDT.at(6);
    usermt = varForBDT.at(7);
    usermtlow = varForBDT.at(8);            
    userleadpt = varForBDT.at(9);
    usertrailpt = varForBDT.at(10);
    userleadingjetpt = varForBDT.at(11);
    usertrailjetpt = varForBDT.at(12);
    userleadeta = fabs(varForBDT.at(13));
    usertraileta = fabs(varForBDT.at(14));
    userchargeOfLeptons = varForBDT.at(15);
    usermll_ss = varForBDT.at(16);
    usermt2ll_ss = varForBDT.at(17);
}
*/

void setStackColors(Color_t & color, int sam){

    for(unsigned int i = 0; i < distribs.size(); i++){
        distribs[i].vectorHisto[sam].SetLineColor(color);
        distribs[i].vectorHisto[sam].SetFillColor(color);
        distribs[i].vectorHisto[sam].SetMarkerColor(color);
    }
}

double mt2ll(const TLorentzVector& l1, const TLorentzVector& l2, const TLorentzVector& metVec){
    return  asymm_mt2_lester_bisect::get_mT2(l1.M(), l1.Px(), l1.Py(), l2.M(), l2.Px(), l2.Py(), metVec.Px(), metVec.Py(), 0, 0);
}

// in readTreeSync.h the histos names are defined with input for th1 initialization

void initListToPrint(const std::string & selection){

  listToPrint["WZ"] = {"ptlead", "sublead", "trail", "njets", "nbjets", "mll", "ptZ", "ptNonZ", "mtW", "mll3e", "mll2e1mu", "mll1e2mu", "mll3mu", "met", "nPV", "mt_3m", "mt_2m1e",  "mt_1m2e", "mt_3e", "cosThetaStar", "flavour3L", "SRWZCR", "mlll", "etaLead", "etaSubl", "etaTrail"};
  listToPrint["Xgamma"] = {"ptlead", "sublead", "trail", "njets", "nbjets", "met", "nPV", "mlll", "flavour3L"};
  listToPrint["ttbar"] = {"ptlead", "sublead", "trail", "njets", "nbjets", "met", "nPV", "mlll", "flavour3L", "SRTTCR", "etaLead", "etaSubl", "etaTrail"};
  listToPrint["DY"] = {"ptlead", "sublead", "trail", "njets", "nbjets", "met", "nPV", "mll", "ptZ", "ptNonZ", "mtW", "mll3e", "mll2e1mu", "mll1e2mu", "mll3mu", "mt_3m", "mt_2m1e",  "mt_1m2e", "mt_3e", "mlll", "flavour3L"};
  listToPrint["ttW"] = {"ptlead", "sublead", "njets", "nbjets", "HT", "met", "nPV", "deltaR", "deltaRlead", "mtLeading", "mtTrailing", "leadJetPt", "trailJetPt", "etaLead", "etaSubl",     "mll_ss", "chargeOfLeptons", "ll_deltaR", "mt2ll_ss", "BDTpp", "BDTmm", "minMLeptonJet", "maxMLeptonJet", "minDeltaRLeptonJet", "maxDeltaRLeptonJet", "minDeltaPhiLeptonJet", "maxDeltaPhiLeptonJet", "minpTLeptonJet", "maxpTLeptonJet", "minMLeptonbJet", "maxMLeptonbJet", "minDeltaRLeptonbJet", "maxDeltaRLeptonbJet", "minDeltaPhiLeptonbJet", "maxDeltaPhiLeptonbJet", "minpTLeptonbJet", "maxpTLeptonbJet", "minMJetJet", "maxMJetJet", "minDeltaRJetJet", "maxDeltaRJetJet", "minDeltaPhiJetJet", "maxDeltaPhiJetJet", "minpTJetJet", "maxpTJetJet", "minDeltaPhiLeptonMET", "maxDeltaPhiLeptonMET", "minmTLeptonMET", "maxmTLeptonMET", "minpTLeptonMET", "maxpTLeptonMET", "minDeltaPhiJetMET", "maxDeltaPhiJetMET", "minmTJetMET", "maxmTJetMET", "minpTJetMET", "maxpTJetMET", "minDeltaPhiBJetMET", "maxDeltaPhiBJetMET", "minmTBJetMET", "maxmTBJetMET", "minpTBJetMET", "maxpTBJetMET"}; 
  listToPrint["ttWclean"] = {"ptlead", "sublead", "njets", "nbjets", "HT", "met", "nPV", "deltaR", "deltaRlead", "mtLeading", "mtTrailing", "leadJetPt", "trailJetPt", "etaLead", "etaSubl",     "mll_ss", "chargeOfLeptons", "ll_deltaR", "mt2ll_ss", "SR", "BDTpp", "BDTmm"}; 
  listToPrint["ZZ"] = {"ptlead", "sublead", "trail", "pt4th", "njets", "nbjets", "met", "nPV", "mll", "ptZ", "etaLead", "etaSubl", "etaTrail", "eta4th", "flavour4LZZ", "SRZZCR"};
  listToPrint["ttZ3L"] = {"ptlead", "sublead", "trail", "njets", "nbjets", "mll", "ptZ", "ptNonZ", "SR3L", "met", "cosThetaStar", "flavour3L"};
  listToPrint["ttZ3Lclean"] = {"ptlead", "sublead", "trail", "njets", "nbjets", "mll", "ptZ", "ptNonZ", "SR3L", "met", "cosThetaStar", "SRttZCleanPTZ", "SRttZCleanCosTheta", "flavour3L", "mll3e", "mll2e1mu", "mll1e2mu", "mll3mu"};
  listToPrint["ttZclean"] = {"ptlead", "sublead", "trail", "njets", "nbjets", "mll", "ptZ", "ptNonZ", "met", "cosThetaStar", "SRttZCleanPTZ", "SRttZCleanCosTheta", "flavour3L", "flavour4L", "flavour3L4L", "mllnoZcut"};
  listToPrint["ttZ4L"] = {"ptlead", "sublead", "trail", "pt4th", "njets", "nbjets", "met", "nPV", "mll", "ptZ", "etaLead", "etaSubl", "etaTrail", "eta4th", "SR4L", "cosThetaStar", "flavour4L"};
  listToPrint["tZq"] = {"ptlead", "sublead", "trail", "njets", "nbjets", "met", "nPV"};
  listToPrint["ttZ"] = {"SRallTTZ", "SR3L", "SR4L", "SRWZCR", "SRZZCR", "SRTTCR", "SR3L3m", "SR3L2m1e", "SR3L1m2e", "SR3L3e", "SRTTZ8SR3L", "flavour3L4L", "ptlead", "trail"};

  listToPrint["CTInMC"] = {"ptlead", "sublead", "trail", "njets", "nbjets", "mll", "met", "HT"};

  listToPrint["ttbar_emu"] = {"ptlead", "sublead", "etaLead", "etaSubl","njets", "nbjets", "met", "nPV","HT", "ptMuonPassedTight", "etaMuonPassedTight", "ptLepPassedLooseForEff", "etaLepPassedLooseForEff", "ptLepPassedTightForEff", "etaLepPassedTightForEff", "etaLepPassedLooseForEffLowPt", "etaLepPassedTightForEffLowPt", "etaLepPassedLooseForEffHighPt", "etaLepPassedTightForEffHighPt"};
  listToPrint["DYTo2L"] = {"ptlead", "sublead", "etaLead", "etaSubl","njets", "nbjets", "met", "nPV","HT"};

// for some processes the entries in the figNames are modified for process specific things.
  if(selection == "ttbar_emu"){
    figNames["ptlead"].nBins = 38;
    figNames["ptlead"].varMin = 10;
    figNames["ptlead"].varMax = 200;

    figNames["ptlead"].fancyName = "Muon p_{T} [GeV]";
    figNames["sublead"].fancyName = "Electron p_{T} [GeV]";

    figNames["etaLead"].fancyName = "Muon #eta";
    figNames["etaSubl"].fancyName = "Electron #eta";
  }

  if(selection == "DYTo2L"){
    figNames["etaLead"].nBins = 50;
    figNames["etaSubl"].nBins = 50;
  }

  if(selection == "ttZclean"){
    figNames["njets"].varMin = 1.5;
    figNames["njets"].nBins = 6;

    figNames["nbjets"].varMin = 0.5;
    figNames["nbjets"].nBins = 3;
  }

  if(selection == "ttZ3Lclean"){
    figNames["njets"].varMin = 2.5;
    figNames["njets"].nBins = 5;

    figNames["nbjets"].varMin = 0.5;
    figNames["nbjets"].nBins = 3;
  }

  if(selection == "ZZ"){
    figNames["njets"].varMax = 5.5;
    figNames["njets"].nBins = 6;

    figNames["nbjets"].varMax = 2.5;
    figNames["nbjets"].nBins = 3;

    figNames["mll"].varMin = 71.;
    figNames["mll"].varMax = 111.;
    figNames["mll"].nBins = 20;
  }

  if(listToPrint.find(selection) == listToPrint.end()){
      std::cerr << "control region selection is incorrect, please double check" << std::endl;
      exit(EXIT_FAILURE);
  }

}
