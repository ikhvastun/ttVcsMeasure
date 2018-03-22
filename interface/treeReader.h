#ifndef treeReader_h
#define treeReader_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLorentzVector.h"

#include <TF1.h>
#include <TH1.h>


class treeReader {
    public :
        //Declare leaf types
        static const unsigned nL_max = 20;
        static const unsigned nJets_max = 20;
        static const unsigned gen_nL_max = 20;
        static const unsigned gen_nPh_max = 10;
        ULong_t       _runNb;
        ULong_t       _lumiBlock;
        ULong_t       _eventNb;
        UChar_t         _nVertex;
        Double_t        _weight;
        UChar_t         _nLheWeights;
        Double_t        _lheWeight[110];
        Float_t         _nTrueInt;
        Double_t        _gen_met;
        Double_t        _gen_metPhi;
        UChar_t         _gen_nL;
        Double_t        _gen_lPt[gen_nL_max];   
        Double_t        _gen_lEta[gen_nL_max];   
        Double_t        _gen_lPhi[gen_nL_max];   
        Double_t        _gen_lE[gen_nL_max];   
        UInt_t          _gen_lFlavor[gen_nL_max];   
        Int_t           _gen_lCharge[gen_nL_max];   
        Int_t           _gen_lMomPdg[gen_nL_max];   
        Bool_t          _gen_lIsPrompt[gen_nL_max];   
        Bool_t          _2017_e;
        Bool_t          _HLT_Ele35_WPTight_Gsf;
        Int_t           _HLT_Ele35_WPTight_Gsf_prescale;
        Bool_t          _HLT_Ele40_WPTight_Gsf;
        Int_t           _HLT_Ele40_WPTight_Gsf_prescale;
        Bool_t          _passTrigger_ee;
        Bool_t          _HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;
        Int_t           _HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_prescale;
        Bool_t          _HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
        Int_t           _HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_prescale;
        Bool_t          _HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
        Int_t           _HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale;
        Bool_t          _2017_eee;
        Bool_t          _HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
        Int_t           _HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale;
        Bool_t          _2017_em;
        Bool_t          _HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;
        Int_t           _HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_prescale;
        Bool_t          _HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
        Int_t           _HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale;
        Bool_t          _HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
        Int_t           _HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale;
        Bool_t          _2017_m;
        Bool_t          _HLT_IsoMu27;
        Int_t           _HLT_IsoMu27_prescale;
        Bool_t          _HLT_IsoMu30;
        Int_t           _HLT_IsoMu30_prescale;
        Bool_t          _2017_mee;
        Bool_t          _HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
        Int_t           _HLT_Mu8_DiEle12_CaloIdL_TrackIdL_prescale;
        Bool_t          _HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;
        Int_t           _HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_prescale;
        Bool_t          _passTrigger_mm;
        Bool_t          _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
        Int_t           _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_prescale;
        Bool_t          _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
        Int_t           _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_prescale;
        Bool_t          _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_Mass8;
        Int_t           _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_Mass8_prescale;
        Bool_t          _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
        Int_t           _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_prescale;
        Bool_t          _HLT_DoubleMu4_Mass8_DZ_PFHT350;
        Int_t           _HLT_DoubleMu4_Mass8_DZ_PFHT350_prescale;
        Bool_t          _HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;
        Int_t           _HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_prescale;
        Bool_t          _2017_mme;
        Bool_t          _HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
        Int_t           _HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale;
        Bool_t          _HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;
        Int_t           _HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_prescale;
        Bool_t          _2017_mmm;
        Bool_t          _HLT_TripleMu_10_5_5_DZ;
        Int_t           _HLT_TripleMu_10_5_5_DZ_prescale;
        Bool_t          _HLT_TripleMu_5_3_3_Mass3p8to60_DZ;
        Int_t           _HLT_TripleMu_5_3_3_Mass3p8to60_DZ_prescale;
        Bool_t          _TripleMu_12_10_5;
        Bool_t          _passMETFilters;
        Bool_t          _Flag_HBHENoiseFilter;
        Bool_t          _Flag_HBHENoiseIsoFilter;
        Bool_t          _Flag_EcalDeadCellTriggerPrimitiveFilter;
        Bool_t          _Flag_goodVertices;
        Bool_t          _Flag_eeBadScFilter;
        Bool_t          _Flag_globalTightHalo2016Filter;
        Bool_t          _Flag_BadPFMuonFilter;
        Bool_t          _Flag_BadChargedCandidateFilter;
        UChar_t         _nL;
        UChar_t         _nMu;
        UChar_t         _nEle;
        UChar_t         _nLight;
        UChar_t         _nTau;
        Double_t        _lPt[nL_max];   
        Double_t        _lEta[nL_max];   
        Double_t        _lEtaSC[nL_max];   
        Double_t        _lPhi[nL_max];   
        Double_t        _lE[nL_max];   
        UInt_t          _lFlavor[nL_max];   
        Int_t           _lCharge[nL_max];   
        Double_t        _dxy[nL_max];   
        Double_t        _dz[nL_max];   
        Double_t        _3dIP[nL_max];   
        Double_t        _3dIPSig[nL_max];   
        Float_t         _lElectronMva[nL_max];   
        Float_t         _lElectronMvaHZZ[nL_max];   
        Bool_t          _lElectronPassEmu[nL_max];   //[_nLight]
        Bool_t          _lElectronPassConvVeto[nL_max];   //[_nLight]
        Bool_t          _lElectronChargeConst[nL_max];   //[_nLight]
        UInt_t          _lElectronMissingHits[nL_max];   //[_nLight]  
        Double_t        _leptonMvaSUSY[nL_max];   //[_nLight]
        Double_t        _leptonMvaTTH[nL_max];   //[_nLight]
        Bool_t          _lHNLoose[nL_max];   
        Bool_t          _lHNFO[nL_max];   
        Bool_t          _lHNTight[nL_max];   
        Bool_t          _lEwkLoose[nL_max];   
        Bool_t          _lEwkFO[nL_max];   
        Bool_t          _lEwkTight[nL_max];   
        Bool_t          _lPOGVeto[nL_max];   
        Bool_t          _lPOGLoose[nL_max];   
        Bool_t          _lPOGMedium[nL_max];   
        Bool_t          _lPOGTight[nL_max];   
        Bool_t          _tauMuonVeto[nL_max];   
        Bool_t          _tauEleVeto[nL_max];   
        Bool_t          _decayModeFindingNew[nL_max];   
        Bool_t          _tauVLooseMvaNew[nL_max];   
        Bool_t          _tauLooseMvaNew[nL_max];   
        Bool_t          _tauMediumMvaNew[nL_max];   
        Bool_t          _tauTightMvaNew[nL_max];   
        Bool_t          _tauVTightMvaNew[nL_max];   
        Bool_t          _tauVTightMvaOld[nL_max];   
        Double_t        _relIso[nL_max];   
        Double_t        _relIso0p4Mu[nL_max];
        Double_t        _miniIso[nL_max];   
        Double_t        _miniIsoCharged[nL_max];   
        Double_t        _ptRel[nL_max];   
        Double_t        _ptRatio[nL_max];   
        Double_t        _closestJetCsvV2[nL_max];
        Double_t        _closestJetDeepCsv_b[nL_max];
        Double_t        _closestJetDeepCsv_bb[nL_max];
        UInt_t          _selectedTrackMult[nL_max];   
        Double_t        _lMuonSegComp[nL_max];   
        Double_t        _lMuonTrackPt[nL_max];   //[_nMu]
        Double_t        _lMuonTrackPtErr[nL_max];   //[_nMu]
        Bool_t          _lIsPrompt[nL_max];   
        Int_t           _lMatchPdgId[nL_max];   
        UChar_t         _nJets;
        Double_t        _jetPt[nJets_max];   
        Double_t        _jetPt_JECUp[nJets_max];   
        Double_t        _jetPt_JECDown[nJets_max];   
        Double_t        _jetPt_JERUp[nJets_max];   
        Double_t        _jetPt_JERDown[nJets_max];   
        Double_t        _jetEta[nJets_max];   
        Double_t        _jetPhi[nJets_max];   
        Double_t        _jetE[nJets_max];   
        Double_t        _jetCsvV2[nJets_max];   
        Double_t        _jetDeepCsv_udsg[nJets_max];   
        Double_t        _jetDeepCsv_b[nJets_max];   
        Double_t        _jetDeepCsv_c[nJets_max];   
        Double_t        _jetDeepCsv_bb[nJets_max];   
        UInt_t          _jetHadronFlavor[nJets_max];   
        UInt_t          _jetId[nJets_max];   
        
        Double_t        _met;
        Double_t        _metPhi;

        Double_t        _metJECDown;
        Double_t        _metJECUp;
        Double_t        _metJetResDown;
        Double_t        _metJetResUp;
        Double_t        _metUnclDown;
        Double_t        _metUnclUp;

        Double_t        _metPhiJECDown;
        Double_t        _metPhiJECUp;
        Double_t        _metPhiJetResDown;
        Double_t        _metPhiJetResUp;
        Double_t        _metPhiUnclDown;
        Double_t        _metPhiUnclUp; 

        Double_t        _rawmet;
        Double_t        _rawmetJECDown;
        Double_t        _rawmetJECUp;
        Double_t        _rawmetJetResDown;
        Double_t        _rawmetJetResUp;
        Double_t        _rawmetUnclDown;
        Double_t        _rawmetUnclUp;

        Double_t        _rawmetPhi;
        Double_t        _rawmetPhiJECDown;
        Double_t        _rawmetPhiJECUp;
        Double_t        _rawmetPhiJetResDown;
        Double_t        _rawmetPhiJetResUp;
        Double_t        _rawmetPhiUnclDown;
        Double_t        _rawmetPhiUnclUp;

        Double_t        _met_sm;                                                                              //met kinematics
        Double_t        _metJECDown_sm;
        Double_t        _metJECUp_sm;
        Double_t        _metJetResDown_sm;
        Double_t        _metJetResUp_sm;
        Double_t        _metUnclDown_sm;
        Double_t        _metUnclUp_sm;

        Double_t        _metPhi_sm;
        Double_t        _metPhiJECDown_sm;
        Double_t        _metPhiJECUp_sm;
        Double_t        _metPhiJetResDown_sm;
        Double_t        _metPhiJetResUp_sm;
        Double_t        _metPhiUnclDown_sm;
        Double_t        _metPhiUnclUp_sm; 

        Bool_t          _HLT_Ele27_WPTight_Gsf;
        Bool_t          _HLT_IsoMu24;
        Bool_t          _HLT_IsoTkMu24;
        

        //Constructor
        treeReader(TTree *tree = nullptr);

        //set up tree for reading and writing
        void initTree(TTree *tree, const bool isData = false);
        void setOutputTree(TTree*, const bool isData = false);

        //skim tree
        void skimTree(const std::string&, std::string outputDirectory = "", const bool isData = false);
        void combinePD(const std::vector<std::string>& datasets, std::string outputDirectory = "");

        //set up tree for analysis
        void readSamples(const std::string& list = ""); //read sample list from file
        void initSample();

        //functions to analyze tree
        void GetEntry(long unsigned entry);
        void Analyze();
        void Loop(const std::string& sample, const double xSection);

        //functions for event selection
        unsigned dilFlavorComb(const std::vector<unsigned>&);
        double coneCorr(const unsigned);
        void setConePt();
        bool lepIsGood(const unsigned);
        bool lepIsFOGood(const unsigned);
        bool lepIsTight(const unsigned);
        unsigned selectLep(std::vector<unsigned>&);
        unsigned selectFakeLep(std::vector<unsigned>&);
        unsigned tightLepCount(const std::vector<unsigned>&, const unsigned);
        bool passPtCuts2L(const std::vector<unsigned>&);
        bool passPtCuts3L(const std::vector<unsigned>&);
        bool jetIsClean(const unsigned, bool);
        bool jetIsGood(const unsigned, const unsigned ptCut = 25, const unsigned unc = 0, const bool clean = true, bool nonpromptSample = false);
        unsigned nJets(const unsigned unc, const bool clean, std::vector<unsigned>&, bool);
        double deltaRCalc(const std::vector<unsigned>& ind, unsigned & lept);
        bool bTaggedDeepCSV(const unsigned unc, const unsigned wp = 1);
        bool bTaggedCSVv2(const unsigned uncm, const unsigned wp = 1);
        unsigned nBJets(const unsigned unc = 0, const bool deepCSV = true, const bool clean = true, const unsigned wp = 1, bool nonpromptSample = false);
        double HTCalc(const std::vector<unsigned>& ind);
        double deltaMZ(const std::vector<unsigned>&, unsigned &, double & , double &, double &, double &);
        int getElectronNumber(const std::vector<unsigned>& ind);

        bool promptLeptons(const std::vector<unsigned>& ind);
        bool leptonIsPrompt(const unsigned& l);

        Color_t assignColor(std::string & name);

        double ptFake(double lpt, double ptratio, int flavour, double mvaTTHvalue, bool mediumIdPassed);

        std::vector<std::pair<double, unsigned>>  ptCorrV;

        std::vector<std::string> getNamesOfTheSample() {return namesOfTheSample;}
        std::vector<Color_t> getColsOfTheStack(){return colsOfTheStack;}
    private:
        TTree* fChain;                                                          //current Tree
        std::shared_ptr<TFile> sampleFile;                                      //current sample
        std::vector<std::tuple<std::string, std::string, double> > samples;     //list of samples
        std::vector<std::string> namesOfTheSample;
        std::vector<Color_t> colsOfTheStack;
        unsigned currentSample = 0;                                             //current index in list
        std::vector<unsigned> samplesOrder;
        std::vector<std::string> samplesOrderNames;
        int leptonSelection;
        bool isData = false;
        bool isDataPrompt = false;
        bool isDataNonprompt = false;
        double scale = 0;
        double weight = 1;                                                      //weight of given event
        unsigned long nEntries = 0;
        //const double dataLumi = 35.9;                                          //in units of 1/fb
        const double dataLumi = 41.9;           // 41.9                               //in units of 1/fb

        // List of branches
        TBranch        *b__runNb;   
        TBranch        *b__lumiBlock;   
        TBranch        *b__eventNb;   
        TBranch        *b__nVertex;   
        TBranch        *b__weight;   
        TBranch        *b__nLheWeights;   
        TBranch        *b__lheWeight;   
        TBranch        *b__nTrueInt;   
        TBranch        *b__gen_met;   
        TBranch        *b__gen_metPhi;   
        TBranch        *b__gen_nL;   
        TBranch        *b__gen_lPt;   
        TBranch        *b__gen_lEta;   
        TBranch        *b__gen_lPhi;   
        TBranch        *b__gen_lE;   
        TBranch        *b__gen_lFlavor;   
        TBranch        *b__gen_lCharge;   
        TBranch        *b__gen_lMomPdg;   
        TBranch        *b__gen_lIsPrompt;   
        TBranch        *b__2017_e;   
        TBranch        *b__HLT_Ele35_WPTight_Gsf;   
        TBranch        *b__HLT_Ele35_WPTight_Gsf_prescale;   
        TBranch        *b__HLT_Ele40_WPTight_Gsf;   
        TBranch        *b__HLT_Ele40_WPTight_Gsf_prescale;   
        TBranch        *b__passTrigger_ee;   
        TBranch        *b__HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;   
        TBranch        *b__HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_prescale;   
        TBranch        *b__HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;   
        TBranch        *b__HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_prescale;   
        TBranch        *b__HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   
        TBranch        *b__HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale;   
        TBranch        *b__2017_eee;   
        TBranch        *b__HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;   
        TBranch        *b__HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale;   
        TBranch        *b__2017_em;   
        TBranch        *b__HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;   
        TBranch        *b__HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_prescale;   
        TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   
        TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale;   
        TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   
        TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale;   
        TBranch        *b__2017_m;   
        TBranch        *b__HLT_IsoMu27;   
        TBranch        *b__HLT_IsoMu27_prescale;   
        TBranch        *b__HLT_IsoMu30;   
        TBranch        *b__HLT_IsoMu30_prescale;   
        TBranch        *b__2017_mee;   
        TBranch        *b__HLT_Mu8_DiEle12_CaloIdL_TrackIdL;   
        TBranch        *b__HLT_Mu8_DiEle12_CaloIdL_TrackIdL_prescale;   
        TBranch        *b__HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;   
        TBranch        *b__HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_prescale;   
        TBranch        *b__passTrigger_mm;   
        TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;   
        TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_prescale;   
        TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;   
        TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_prescale;   
        TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_Mass8;   
        TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_Mass8_prescale;   
        TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;   
        TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_prescale;   
        TBranch        *b__HLT_DoubleMu4_Mass8_DZ_PFHT350;   
        TBranch        *b__HLT_DoubleMu4_Mass8_DZ_PFHT350_prescale;   
        TBranch        *b__HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;   
        TBranch        *b__HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_prescale;   
        TBranch        *b__2017_mme;   
        TBranch        *b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL;   
        TBranch        *b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale;   
        TBranch        *b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;   
        TBranch        *b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_prescale;   
        TBranch        *b__2017_mmm;   
        TBranch        *b__HLT_TripleMu_10_5_5_DZ;   
        TBranch        *b__HLT_TripleMu_10_5_5_DZ_prescale;   
        TBranch        *b__HLT_TripleMu_5_3_3_Mass3p8to60_DZ;   
        TBranch        *b__HLT_TripleMu_5_3_3_Mass3p8to60_DZ_prescale;   
        TBranch        *b__TripleMu_12_10_5;   
        TBranch        *b__passMETFilters;   
        TBranch        *b__Flag_HBHENoiseFilter;   //!
        TBranch        *b__Flag_HBHENoiseIsoFilter;   //!
        TBranch        *b__Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
        TBranch        *b__Flag_goodVertices;   //!
        TBranch        *b__Flag_eeBadScFilter;   //!
        TBranch        *b__Flag_globalTightHalo2016Filter;   //!
        TBranch        *b__Flag_BadPFMuonFilter;   //!
        TBranch        *b__Flag_BadChargedCandidateFilter;   //!
        TBranch        *b__nL;   
        TBranch        *b__nMu;   
        TBranch        *b__nEle;   
        TBranch        *b__nLight;   
        TBranch        *b__nTau;   
        TBranch        *b__lPt;   
        TBranch        *b__lEta;   
        TBranch        *b__lEtaSC;   
        TBranch        *b__lPhi;   
        TBranch        *b__lE;   
        TBranch        *b__lFlavor;   
        TBranch        *b__lCharge;   
        TBranch        *b__dxy;   
        TBranch        *b__dz;   
        TBranch        *b__3dIP;   
        TBranch        *b__3dIPSig;   
        TBranch        *b__lElectronMva;   
        TBranch        *b__lElectronMvaHZZ;   
        TBranch        *b__lElectronPassEmu;   //!
        TBranch        *b__lElectronPassConvVeto;   //!
        TBranch        *b__lElectronChargeConst;   //!
        TBranch        *b__lElectronMissingHits;   //! 
        TBranch        *b__leptonMvaSUSY;   //!
        TBranch        *b__leptonMvaTTH;   //! 
        TBranch        *b__lHNLoose;   
        TBranch        *b__lHNFO;   
        TBranch        *b__lHNTight;   
        TBranch        *b__lEwkLoose;   
        TBranch        *b__lEwkFO;   
        TBranch        *b__lEwkTight;   
        TBranch        *b__lPOGVeto;   
        TBranch        *b__lPOGLoose;   
        TBranch        *b__lPOGMedium;   
        TBranch        *b__lPOGTight;   
        TBranch        *b__tauMuonVeto;   
        TBranch        *b__tauEleVeto;   
        TBranch        *b__decayModeFindingNew;   
        TBranch        *b__tauVLooseMvaNew;   
        TBranch        *b__tauLooseMvaNew;   
        TBranch        *b__tauMediumMvaNew;   
        TBranch        *b__tauTightMvaNew;   
        TBranch        *b__tauVTightMvaNew;   
        TBranch        *b__tauVTightMvaOld;   
        TBranch        *b__relIso;   
        TBranch        *b__relIso0p4Mu;
        TBranch        *b__miniIso;   
        TBranch        *b__miniIsoCharged;   
        TBranch        *b__ptRel;   
        TBranch        *b__ptRatio;   
        TBranch        *b__closestJetCsvV2;
        TBranch        *b__closestJetDeepCsv_b;
        TBranch        *b__closestJetDeepCsv_bb;
        TBranch        *b__selectedTrackMult;   
        TBranch        *b__lMuonSegComp;   
        TBranch        *b__lMuonTrackPt;   //!
        TBranch        *b__lMuonTrackPtErr;   //!
        TBranch        *b__lIsPrompt;   
        TBranch        *b__lMatchPdgId;   
        TBranch        *b__nJets;   
        TBranch        *b__jetPt;   
        TBranch        *b__jetPt_JECUp;   
        TBranch        *b__jetPt_JECDown;   
        TBranch        *b__jetPt_JERUp;   
        TBranch        *b__jetPt_JERDown;   
        TBranch        *b__jetEta;   
        TBranch        *b__jetPhi;   
        TBranch        *b__jetE;   
        TBranch        *b__jetCsvV2;   
        TBranch        *b__jetDeepCsv_udsg;   
        TBranch        *b__jetDeepCsv_b;   
        TBranch        *b__jetDeepCsv_c;   
        TBranch        *b__jetDeepCsv_bb;   
        TBranch        *b__jetHadronFlavor;   
        TBranch        *b__jetId;   
       
       TBranch        *b__met;   //!
       TBranch        *b__metJECDown;   //!
       TBranch        *b__metJECUp;   //!
       TBranch        *b__metJetResDown;   //!
       TBranch        *b__metJetResUp;   //!
       TBranch        *b__metUnclDown;   //!
       TBranch        *b__metUnclUp;   //!
       TBranch        *b__metPhi;   //!
       TBranch        *b__metPhiJECDown;   //!
       TBranch        *b__metPhiJECUp;   //!
       TBranch        *b__metPhiJetResDown;   //!
       TBranch        *b__metPhiJetResUp;   //!
       TBranch        *b__metPhiUnclDown;   //!
       TBranch        *b__metPhiUnclUp;   //!
    
       TBranch        *b__rawmetPhi;
       TBranch        *b__rawmetPhiJECDown;   //!
       TBranch        *b__rawmetPhiJECUp;   //!
       TBranch        *b__rawmetPhiJetResDown;   //!
       TBranch        *b__rawmetPhiJetResUp;   //!
       TBranch        *b__rawmetPhiUnclDown;   //!
       TBranch        *b__rawmetPhiUnclUp;   //!
       TBranch        *b__rawmet;   //!
       TBranch        *b__rawmetJECDown;   //!
       TBranch        *b__rawmetJECUp;   //!
       TBranch        *b__rawmetJetResDown;   //!
       TBranch        *b__rawmetJetResUp;   //!
       TBranch        *b__rawmetUnclDown;   //!
       TBranch        *b__rawmetUnclUp;   //!
       TBranch        *b__met_sm;   //!
       TBranch        *b__metJECDown_sm;   //!
       TBranch        *b__metJECUp_sm;   //!
       TBranch        *b__metJetResDown_sm;   //!
       TBranch        *b__metJetResUp_sm;   //!
       TBranch        *b__metUnclDown_sm;   //!
       TBranch        *b__metUnclUp_sm;   //!
       TBranch        *b__metPhi_sm;   //!
       TBranch        *b__metPhiJECDown_sm;   //!
       TBranch        *b__metPhiJECUp_sm;   //!
       TBranch        *b__metPhiJetResDown_sm;   //!
       TBranch        *b__metPhiJetResUp_sm;   //!
       TBranch        *b__metPhiUnclDown_sm;   //!
       TBranch        *b__metPhiUnclUp_sm;   //!

        TBranch        *b__HLT_Ele27_WPTight_Gsf;   //!
        TBranch        *b__HLT_Ele27_WPTight_Gsf_prescale;   //!
        TBranch        *b__HLT_IsoMu24;   //!
        TBranch        *b__HLT_IsoMu24_prescale;   //!
        TBranch        *b__HLT_IsoTkMu24;   //!
        TBranch        *b__HLT_IsoTkMu24_prescale;   //!
};
#endif
