#ifndef Reweighter_H
#define Reweighter_H

//include c++ library classes

//include ROOT classes
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"

//include other parts of code 
#include "../interface/BTagCalibrationStandalone.h"
#include "../interface/Sample.h"

//Class storing scale-factor weights to be used in events
class Reweighter{
    public:
        Reweighter(const std::vector<Sample>&, const bool is2016, const int lepSel);
        ~Reweighter();

        // function for uncorrelated uncertainty calculation
        double uncorUncCalc(const double & centralOne, const double & centralTwo, const double & uncOne, const double & uncTwo) const;

        //pileup weight
        double puWeight(const double nTrueInt, const Sample&, const unsigned unc = 0) const;

        //b-tag weight
        double bTagWeight(const unsigned jetFlavor, const double jetPt, const double jetEta, const double jetCSV, const unsigned unc = 0) const;

        //b-tagging efficiency
        double bTagEff(const unsigned jetFlavor, const double jetPt, const double jetEta) const;

        //tracking + reconstruction weights
        double muonRecoWeight(const double eta, const unsigned unc) const;
        double electronRecoWeight(const double superClusterEta, const double pt, const unsigned unc) const;

        //lepton id + reconstruction weight
        /*
        // split into 2 functions, separately reco and tight, 14th of September, substitued with one SF, tight on top for reco
        double muonTightWeight(const double pt, const double eta, const unsigned unc) const{ 
            // no reco SF, email from Sergio and Sandro on 30th of July, 2018, they removed as well recommendations from twiki, varOne is always considered as 0
            //double varOne = TMath::Abs(muonRecoWeight(eta, 0) - muonRecoWeight(eta, 0));
            //std::cout << "muon sf with unc " << unc << " is " << muonTightIdWeight(pt,eta, unc) << std::endl;
            return muonTightIdWeight(pt,eta, unc);
        }
        */

        // split into 2 functions, separately reco and tight, 29th Aug 2018
        /*
        double electronTightWeight(const double pt, const double eta, const double superClusterEta, const unsigned unc) const{ 
            int var = unc == 2 ? -1 : int(unc);
            double varOne = (var != 0 ? fabs(electronRecoWeight(superClusterEta, pt, unc) - electronRecoWeight(superClusterEta, pt, 0)) : 0.);
            double varTwo = (var != 0 ? fabs(electronTightIdWeight(pt,superClusterEta, unc) - electronTightIdWeight(pt,superClusterEta, 0)) : 0.);
            // let's implement fully uncorrelated uncertainties
            double valueOne = electronRecoWeight(superClusterEta, pt, 0);
            double valueTwo = electronTightIdWeight(pt,superClusterEta, 0);
            //if(var == 2) std::cout << "electron tight id weight is " << valueOne*valueTwo + (var != 0 ? var * uncorUncCalc(valueOne, valueTwo, varOne, varTwo) : 0.) << std::endl;
            return valueOne*valueTwo + (var != 0 ? var * uncorUncCalc(valueOne, valueTwo, varOne, varTwo) : 0.);
        }
        */
        double muonTightWeightOnlyStat(const double pt, const double eta, const unsigned unc) const{ 
            return muonTightIdWeightOnlyStat(pt,eta, unc);
        }
        double muonTightWeightOnlySyst(const double pt, const double eta, const unsigned unc) const{ 
            return muonTightIdWeightOnlySyst(pt,eta, unc);
        }
        
        double electronTightWeightOnlyStat(const double pt, const double eta, const double superClusterEta, const unsigned unc) const{ 
            return electronTightIdWeightOnlyStat(pt,superClusterEta, unc);
        }
        double electronTightWeightOnlySyst(const double pt, const double eta, const double superClusterEta, const unsigned unc) const{ 
            return electronTightIdWeightOnlySyst(pt,superClusterEta, unc);
        }
        double muonLooseWeight(const double pt, const double eta, const unsigned unc) const{
            return  muonRecoWeight(eta, unc)*muonLooseIdWeight(pt,eta, unc);
        }

        double electronLooseWeight(const double pt, const double eta, const double superClusterEta, const unsigned unc) const{
            return electronRecoWeight(superClusterEta, pt, unc)*electronLooseIdWeight(pt,eta, unc);
        }
        
        // charge consistency weights (together with tight)
        double muonChargeConsWeight(const double pt, const double eta, const unsigned unc) const;
        double electronChargeConsWeight(const double pt, const double eta, const unsigned unc) const;

        //fakerates 
        double muonFakeRate(const double pt, const double eta, const unsigned unc = 0) const;
        double electronFakeRate(const double pt, const double eta, const unsigned unc = 0) const;

        // charge mis ID
        double electronCMIDRate(const double pt, const double eta, const unsigned unc = 0) const;

        // trigger SF
        double getTriggerSF(const double pt) const;

        //jet prefering probabilities
        double jetPrefiringProbability(const double pt, const double eta) const;

    private:
        //boolean flagging weights as 2016 or 2017
        bool is2016;
        int leptonSelection;

        //pu weights (one for every sample)
        std::map< std::string, std::vector< std::shared_ptr<TH1D> > > puWeights;

        //btag scale factors and efficiencies
        std::shared_ptr<BTagCalibration> bTagCalib;
        std::shared_ptr<BTagCalibrationReader> bTagCalibReader;
        std::shared_ptr<TH2D> bTagEffHist[3];

        //reconstruction scale factors
        std::shared_ptr<TGraph> muonRecoSF;
        std::shared_ptr<TH2D> electronRecoSF;
        std::shared_ptr<TH2D> electronRecoSF_pT0to20;
        std::shared_ptr<TH2D> electronRecoSF_pT20toInf;

        //muon id scale factors
        std::shared_ptr<TH2D> muonLooseToRecoSF;
        std::shared_ptr<TH2D> muonTightToLooseSF;
        std::shared_ptr<TH2D> muonTightToRecoSF;
        std::shared_ptr<TH2D> muonChargeConsToTightSF;
        std::shared_ptr<TH2D> muonTightToRecoSF_syst;
        std::shared_ptr<TH2D> muonTightToRecoSF_stat;

        //electron id scale factors
        std::shared_ptr<TH2D> electronLooseToRecoSF;
        std::shared_ptr<TH2D> electronTightToLooseSF;
        std::shared_ptr<TH2D> electronTightToRecoSF_syst;
        std::shared_ptr<TH2D> electronTightToRecoSF_stat;
        std::shared_ptr<TH2D> electronChargeConsToTightSF;


        //initialize all weight histograms
        void initialize2016Weights();
        void initialize2017Weights();

        //return jet flavor index 1 -> 0, 4 -> 1, 5 -> 2
        unsigned flavorInd(const unsigned jetFlavor) const{ 
            return 0 + (jetFlavor == 4) + 2*(jetFlavor == 5);
        }
        
        // keep order same as in other uncertainties
        const std::map<const unsigned, const int> var = {{0, 0}, {1, 1}, {2, -1}};

        //fake rate maps
        std::shared_ptr<TH2D> frMapEle[3];
        std::shared_ptr<TH2D> frMapMu[3];
        
        //charge mis ID rate maps
        std::shared_ptr<TH2D> CMIDMapEle[3];

        //jet prefiring probabilities
        std::shared_ptr<TH2D> prefiringMap;

        //loose id weights
        double muonLooseIdWeight(const double pt, const double eta, const unsigned unc) const;
        double electronLooseIdWeight(const double pt, const double eta, const unsigned unc) const;

        //tight id weights
        double muonTightIdWeightOnlySyst(const double pt, const double eta, const unsigned unc) const;
        double muonTightIdWeightOnlyStat(const double pt, const double eta, const unsigned unc) const;
        double electronTightIdWeightOnlySyst(const double pt, const double eta, const unsigned unc) const;
        double electronTightIdWeightOnlyStat(const double pt, const double eta, const unsigned unc) const;

        //read pu weights for a given list of samples
        void initializePuWeights(const std::vector< Sample >&); 

        //read b-tagging weights
        void initializeBTagWeights();

        //read electron id and reco weights
        void initializeElectronWeights();

        //read muon id and reco weights 
        void initializeMuonWeights();

        //initialize fake-rate
        void initializeFakeRate();

        // initialize charge mis ID 
        void initializeCMIDRate();

        //initialize all weights 
        void initializeAllWeights(const std::vector< Sample>&);


        //initialize jet prefiring probabilities
        void initializePrefiringProbabilities();
};
#endif
