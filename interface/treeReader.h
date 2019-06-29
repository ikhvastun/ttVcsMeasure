#ifndef treeReader_h
#define treeReader_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLorentzVector.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include <TF1.h>
#include <TH1.h>

#include "Reweighter.h"
#include "Sample.h"

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
        UChar_t         _nPsWeights;
        Double_t        _psWeight[14];
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
        Bool_t          _passTrigger_e;
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
        Bool_t          _passTrigger_eee;
        Bool_t          _HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
        Int_t           _HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale;
        Bool_t          _passTrigger_em;
        Bool_t          _HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;
        Int_t           _HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_prescale;
        Bool_t          _HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
        Int_t           _HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale;
        Bool_t          _HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
        Int_t           _HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale;
        Bool_t          _passTrigger_m;
        Bool_t          _HLT_IsoMu27;
        Int_t           _HLT_IsoMu27_prescale;
        Bool_t          _HLT_IsoMu30;
        Int_t           _HLT_IsoMu30_prescale;
        Bool_t          _passTrigger_mee;
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
        Bool_t          _passTrigger_mme;
        Bool_t          _HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
        Int_t           _HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale;
        Bool_t          _HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;
        Int_t           _HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_prescale;
        Bool_t          _passTrigger_emm;
        Bool_t          _passTrigger_eem;
        Bool_t          _passTrigger_mmm;
        Bool_t          _HLT_TripleMu_10_5_5_DZ;
        Int_t           _HLT_TripleMu_10_5_5_DZ_prescale;
        Bool_t          _HLT_TripleMu_5_3_3_Mass3p8to60_DZ;
        Int_t           _HLT_TripleMu_5_3_3_Mass3p8to60_DZ_prescale;
        Bool_t          _TripleMu_12_10_5;
        Bool_t          _passMETFilters;
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
        Float_t         _lElectronMvaFall17NoIso[nL_max];
        Bool_t          _lElectronPassEmu[nL_max];   //[_nLight]
        Bool_t          _lElectronPassConvVeto[nL_max];   //[_nLight]
        Bool_t          _lElectronChargeConst[nL_max];   //[_nLight]
        UInt_t          _lElectronMissingHits[nL_max];   //[_nLight]  
        Double_t        _leptonMvaSUSY[nL_max];   //[_nLight]
        Double_t        _leptonMvaTTH[nL_max];   //[_nLight]
        Double_t        _leptonMvatZqTTV[nL_max];   //[_nLight]
        Double_t        _leptonMvatZqTTV16[nL_max];   //[_nLight]
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
        UInt_t          _lProvenance[nL_max];   //[_nL]
        UInt_t          _lProvenanceCompressed[nL_max];
        UChar_t         _nJets;
        Double_t        _jetPt[nJets_max];   
        Double_t        _jetPt_JECUp[nJets_max];   
        Double_t        _jetPt_JECDown[nJets_max];   
        Double_t        _jetPt_JERUp[nJets_max];   
        Double_t        _jetPt_JERDown[nJets_max];   
        Double_t        _jetSmearedPt[nJets_max];
        Double_t        _jetSmearedPt_JECDown[nJets_max];
        Double_t        _jetSmearedPt_JECUp[nJets_max];
        Double_t        _jetSmearedPt_JERDown[nJets_max];
        Double_t        _jetSmearedPt_JERUp[nJets_max];
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
        Bool_t          _jetIsTight[nJets_max];   
        Double_t        _met;
        Double_t        _metJECDown;
        Double_t        _metJECUp;
        Double_t        _metUnclDown;
        Double_t        _metUnclUp;
        Double_t        _metPhi;
        Double_t        _metPhiJECDown;
        Double_t        _metPhiJECUp;
        Double_t        _metPhiUnclDown;
        Double_t        _metPhiUnclUp;  

        Bool_t          _HLT_Ele27_WPTight_Gsf;
        Bool_t          _HLT_IsoMu24;
        Bool_t          _HLT_IsoTkMu24;
        
        bool debug;
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
        //general function to read a list of samples
        void readSamples(const std::string&, std::vector<Sample>&);

        void initSample(std::string option);
        void initSample(const Sample&, std::string option);

        //functions to analyze tree
        void GetEntry(long unsigned entry);
        void Analyze();
        void AnalyzeFR(const std::vector<std::string>& fileToAnalyse, const std::string outputDir);
        void Analyze(const std::vector<std::string>& fileToAnalyse, const std::string option = "", const std::string selection = "", const std::string& sampleToDebug = "", long = -999);
        void GetEntry(const Sample&, long unsigned entry);
        void Loop(const std::string& sample, const double xSection);

        //functions for event selection
        unsigned dilFlavorComb(const std::vector<unsigned>&);
        double coneCorr(const unsigned);
        void setConePt();
        bool lepIsGoodFortZq(const unsigned);
        bool lepIsGood(const unsigned, const int);
        bool lepIsFOGoodFortZq(const unsigned);
        bool lepIsFOGood(const unsigned, const int);
        unsigned selectLep(std::vector<unsigned>&, const int);
        unsigned selectFakeLep(std::vector<unsigned>&, const int);
        unsigned selectLooseLep(std::vector<unsigned>&, const int);
        unsigned selectLepGoodForLeptonMVA(std::vector<unsigned>& ind);
        unsigned tightLepCount(const std::vector<unsigned>&, const unsigned);
        bool passPtCuts2LOF(const std::vector<unsigned>&);
        bool passPtCuts2LOSSF(const std::vector<unsigned>&);
        bool passPtCuts2L(const std::vector<unsigned>&);
        bool passPtCuts3L(const std::vector<unsigned>&);
        bool passPtCuts4L(const std::vector<unsigned>&);
        bool jetIsClean(const unsigned, const int);
        bool jetIsGood(const unsigned, const unsigned ptCut = 25, const unsigned unc = 0, const bool clean = true, bool is2017 = false, const double eta = 2.4);
        unsigned nJets(const unsigned unc, const bool clean, std::vector<unsigned>&, bool, const double eta = 2.4);
        bool isForwardJetPresent(bool is2017);
        unsigned nJetsNotB(const unsigned unc, const bool clean, std::vector<unsigned>&, const unsigned wp, bool);
        double deltaRCalc(const std::vector<unsigned>& ind, unsigned & lept, const bool deepCSV = true);
        bool bTaggedDeepCSV(const unsigned unc, const unsigned wp = 1);
        bool bTaggedCSVv2(const unsigned uncm, const unsigned wp = 1);
        //unsigned nBJets(const unsigned unc = 0, const bool deepCSV = true, const bool clean = true, std::vector<unsigned>& ind = vector<unsigned>(), const unsigned wp = 1, bool nonpromptSample = false);
        unsigned nBJets(const unsigned unc, const bool deepCSV, const bool clean, std::vector<unsigned>& ind, const unsigned wp, bool nonpromptSample);
        double HTCalc(const std::vector<unsigned>& ind);
        double deltaMZ(const std::vector<unsigned>&, unsigned &, double & , double &, double &, double &, std::vector<unsigned>&, TLorentzVector &, TLorentzVector &);
        bool invMassOfAny2Lbelow12GeV(const std::vector<unsigned>& ind);
        bool invMassOfAny2Lbelow20GeV(const std::vector<unsigned>& ind);
        int getElectronNumber(const std::vector<unsigned>& ind);
        bool elePassVLooseMvaIDSUSY(const unsigned ind);
        bool eleIsClean(const unsigned ind);
        bool lepIsLoose(const unsigned ind);
        bool lepIsGoodForLeptonMVA(const unsigned ind);
        bool passTTZSelection(const int, const double) const;
        bool passTTZCleanSelection(const int, const int, const double) const;
        bool passWZCRSelection(const int, const double) const;
        bool passttbarCRSelection(const int, const double, const double) const;
        bool passttbarCRintZqSelection(const int njets, const int nbjets, const double dMZ) const;
        bool passZGCRSelection(const double, const double) const;
        bool passDYCRSelection(const double, const double, const unsigned, const double, const double, const int, const int) const;
        double mtCalc(const TLorentzVector Vect, const double MET, const double MET_Phi) const;
        bool pass2Lpreselection(const int njets, const int nbjets, const std::vector<unsigned>& ind, const double met, const int nEle);
        bool pass2Lcleanpreselection(const int njets, const int nbjets, const std::vector<unsigned>& ind, const double met, const int nEle);
        bool passTTZ4LSelection(const std::vector<unsigned>& ind, std::vector<unsigned> indOf2LonZ, const int njets);
        bool passZZCRSelection(const std::vector<unsigned>& ind, std::vector<unsigned> indOf2LonZ, const int & njets);
        double SRIDTTZ(const std::vector<unsigned>& ind, std::vector<unsigned> indOf2LonZ, const int & njets, const int & nbjets, const double & dMZ, const double & mlll);
        double SRIDTTZ1L(const std::vector<unsigned>& ind, std::vector<unsigned> indOf2LonZ, const int & njets, const int & nbjets, const double & dMZ, const double & mlll);
        double SRIDTTCR(const int & njets, const int & nbjets, const double & dMZ, const double & mlll);
        double SRIDZZCR(const std::vector<unsigned>& ind, std::vector<unsigned> indOf2LonZ, const int & njets, const int & nbjets);
        double SRIDWZCR(const int & njets, const int & nbjets, const double & dMZ);
        double SRID3L(int & njets, int & nbjets, const double & dMZ);
        double SRID4L(int & njets, int & nbjets);
        double SRIDPTZ(const double & ptZ) const;
        double SRIDCosTheta(const double & cosTheta) const;
        double sumAllLeptonsCharge(const std::vector<unsigned>& ind);
        double SRID8SR3L(int & njets, int & nbjets, const double & dMZ);
        bool passTTZSRSelection(const std::vector<unsigned>& ind, std::vector<unsigned> indOf2LonZ, const int njets, const int nbjets, const double dMZ);
        TLorentzVector findBestNeutrinoAndTop(const TLorentzVector& wLep, const TLorentzVector& met, std::vector<unsigned>& taggedJetI, const std::vector<unsigned>& jetI, const std::vector<unsigned>&     bJetI, const TLorentzVector* jetV);
        std::pair<double, double> neutrinoPZ(const TLorentzVector& wLep, const TLorentzVector& met);

        bool promptLeptons(const std::vector<unsigned>& ind);
        bool leptonIsPrompt(const unsigned& l);
        bool noConversionInSelection(const std::vector<unsigned>& ind);
        bool leptonFromConversion(const unsigned & l);

        Color_t assignColor(const std::string & name);

        double ptFake(double lpt, double ptratio, int flavour, double mvaTTHvalue, bool mediumIdPassed);
        double ptFakeStIso(double lpt, int flavor, double isolation);
        bool twoLeptonsInEndcap(const std::vector<unsigned>& ind);
        double cosThetaStar(const TLorentzVector &, const TLorentzVector &);

        double puWeight(const unsigned unc = 0);
        double bTagWeight(const unsigned jetFlavor, const unsigned unc = 0);
        double bTagWeight_udsg(const unsigned unc = 0);
        double bTagWeight_c(const unsigned unc = 0);
        double bTagWeight_b(const unsigned unc = 0);
        double bTagWeight(const unsigned unc = 0);
        double leptonWeightOnlyStat(const unsigned unc = 0, const bool ignoreMuon = false, const bool ignoreEle = false);
        double leptonWeightOnlySyst(const unsigned unc = 0, const bool ignoreMuon = false, const bool ignoreEle = false);
        double leptonWeightOnlyReco(const unsigned unc = 0, const bool ignoreMuon = false, const bool ignoreEle = false);
        double triggerWeight();
        void initializeWeights();
        double sfWeight();
        double fakeRateWeight(const unsigned unc = 0);
        double CMIDRateWeight(const unsigned unc = 0);
        double jetPrefiringWeight();

        double largestAmongAll(const std::vector<double> & weights);
        double smallestAmongAll(const std::vector<double> & weights);
        void addVariablesToBDT();
        void fillBDTvariables(std::vector<Float_t> & varForBDT, int flavor);

        std::vector<std::pair<double, unsigned>>  ptCorrV;

        std::vector<std::string> getNamesOfTheFiles() {return namesOfTheFiles;}
        std::vector<std::string> getNamesOfTheProcesses() {return namesOfTheProcesses;}
        std::vector<Color_t> getColsOfTheStack(){return colsOfTheStack;}
    private:
        TTree* fChain;                                                          //current Tree
        std::shared_ptr<TFile> sampleFile;                                      //current sample
        std::vector<Sample> samples;
        //std::vector<std::tuple<std::string, std::string, double> > samples;     //list of samples
        std::vector<std::string> namesOfTheFiles;
        std::vector<std::string> namesOfTheProcesses;
        std::vector<Color_t> colsOfTheStack;
        //unsigned currentSampleIndex = 0;                                             //current index in list
        Sample currentSample; 
        int currentSampleIndex = -1;
        std::vector<unsigned> samplesOrder;
        std::vector<std::string> samplesOrderNames;
        int leptonSelection;
        double leptonMVAcut;
        double magicFactor;
        bool isData = false;
        bool isDataNonprompt = false;
        bool isChargeMisIDSample = false;
        bool is2017folder = false;
        double scale = 0;
        double sumSimulatedEventWeights = 0;
        double sumSimulatedEventWeightsScaleUp = 0;
        double sumSimulatedEventWeightsScaleDown = 0;
        double sumSimulatedEventWeightsISRScaleUp = 0;
        double sumSimulatedEventWeightsISRScaleDown = 0;
        double sumSimulatedEventWeightsFSRScaleUp = 0;
        double sumSimulatedEventWeightsFSRScaleDown = 0;
        double weight = 1;                                                      //weight of given event
        unsigned long nEntries = 0;
        double dataLumi = 41.5;                                          //in units of 1/fb
        int nonPromptSample = -999;
        int CMIDSample = -999;
        std::map<int, double> leptonMVAcutInAnalysis = {{2, 0.6}, {3, 0.4}, {4, -0.4}};
        std::map<int, double> magicFactorInAnalysis = {{2, 0.9}, {3, 0.85}};
        std::shared_ptr<Reweighter> reweighter;
        double crossSectionRatio[100][100];

        std::map<std::string, int> processToCounterMap;

        bool is2017() { return currentSample.is2017(); }
        bool is2016() { return currentSample.is2016(); }                  //if sample is not 2017 it is automatically 2016
        double lumi2016 = 35.9;
        double lumi2017 = 41.5;
        //bool isData() { return currentSample.isData(); }
        //bool isMC() { return currentSample.isMC(); } 
        Float_t user_pt, user_eta, user_trackMult,  user_miniIsoCharged, user_miniIsoNeutral, user_ptrel, user_ptratio, user_jetBtagCSV, user_sip3d, user_dxy, user_dz, user_segmComp, user_eleMVA, user_relIso;
        TMVA::Reader *readerLeptonMVAele = new TMVA::Reader( "!Color:!Silent" );
        TMVA::Reader *readerLeptonMVAmu = new TMVA::Reader( "!Color:!Silent" );

        bool isElectron(const unsigned leptonIndex) const { return (_lFlavor[leptonIndex] == 0); }
        bool isMuon(const unsigned leptonIndex) const { return (_lFlavor[leptonIndex] == 1); }

        // List of branches
        TBranch        *b__runNb;   
        TBranch        *b__lumiBlock;   
        TBranch        *b__eventNb;   
        TBranch        *b__nVertex;   
        TBranch        *b__weight;   
        TBranch        *b__nLheWeights;   
        TBranch        *b__lheWeight;   
        TBranch        *b__nPsWeights;
        TBranch        *b__psWeight;
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
        TBranch        *b__passTrigger_e;   
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
        TBranch        *b__passTrigger_eee;   
        TBranch        *b__passTrigger_eem;   
        TBranch        *b__passTrigger_emm;   
        TBranch        *b__HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;   
        TBranch        *b__HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale;   
        TBranch        *b__passTrigger_em;   
        TBranch        *b__HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;   
        TBranch        *b__HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_prescale;   
        TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   
        TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale;   
        TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   
        TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale;   
        TBranch        *b__passTrigger_m;   
        TBranch        *b__HLT_IsoMu27;   
        TBranch        *b__HLT_IsoMu27_prescale;   
        TBranch        *b__HLT_IsoMu30;   
        TBranch        *b__HLT_IsoMu30_prescale;   
        TBranch        *b__passTrigger_mee;   
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
        TBranch        *b__passTrigger_mme;   
        TBranch        *b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL;   
        TBranch        *b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale;   
        TBranch        *b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;   
        TBranch        *b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_prescale;   
        TBranch        *b__passTrigger_mmm;   
        TBranch        *b__HLT_TripleMu_10_5_5_DZ;   
        TBranch        *b__HLT_TripleMu_10_5_5_DZ_prescale;   
        TBranch        *b__HLT_TripleMu_5_3_3_Mass3p8to60_DZ;   
        TBranch        *b__HLT_TripleMu_5_3_3_Mass3p8to60_DZ_prescale;   
        TBranch        *b__TripleMu_12_10_5;   
        TBranch        *b__passMETFilters;   
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
        TBranch        *b__lElectronMvaFall17NoIso;
        TBranch        *b__lElectronPassEmu;   //!
        TBranch        *b__lElectronPassConvVeto;   //!
        TBranch        *b__lElectronChargeConst;   //!
        TBranch        *b__lElectronMissingHits;   //! 
        TBranch        *b__leptonMvaSUSY;   //!
        TBranch        *b__leptonMvaTTH;   //! 
        TBranch        *b__leptonMvatZqTTV;   //! 
        TBranch        *b__leptonMvatZqTTV16;   //! 
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
        TBranch        *b__lProvenance;   //!
        TBranch        *b__lProvenanceCompressed;   //!
        TBranch        *b__nJets;   
        TBranch        *b__jetPt;   
        TBranch        *b__jetPt_JECUp;   
        TBranch        *b__jetPt_JECDown;   
        TBranch        *b__jetPt_JERUp;   
        TBranch        *b__jetPt_JERDown;   
        TBranch        *b__jetSmearedPt;
        TBranch        *b__jetSmearedPt_JECDown;
        TBranch        *b__jetSmearedPt_JECUp;
        TBranch        *b__jetSmearedPt_JERDown;
        TBranch        *b__jetSmearedPt_JERUp;
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
        TBranch        *b__jetIsTight;   
        TBranch        *b__met;   
        TBranch        *b__metJECDown;   
        TBranch        *b__metJECUp;   
        TBranch        *b__metUnclDown;   
        TBranch        *b__metUnclUp;   
        TBranch        *b__metPhi;   
        TBranch        *b__metPhiJECDown;   
        TBranch        *b__metPhiJECUp;   
        TBranch        *b__metPhiUnclDown;   
        TBranch        *b__metPhiUnclUp;   

        TBranch        *b__HLT_Ele27_WPTight_Gsf;   //!
        TBranch        *b__HLT_Ele27_WPTight_Gsf_prescale;   //!
        TBranch        *b__HLT_IsoMu24;   //!
        TBranch        *b__HLT_IsoMu24_prescale;   //!
        TBranch        *b__HLT_IsoTkMu24;   //!
        TBranch        *b__HLT_IsoTkMu24_prescale;   //!
};
#endif
