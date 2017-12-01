#include <algorithm>
#include <vector>
#include <map>
#include <iomanip>

#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TApplication.h"
#include "TFrame.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "BTagCalibrationStandalone.h" 

#include "showHist.h"
#include "readTreeSync.h"
#include "Tools.h"
#include "Output.h"

#include "errors.h"

using namespace std;

Errors LastError::lasterror = Errors::UNKNOWN;
using Output::distribs;

void runttWvsTTbar()
{

  setTDRStyle(); 
  initdistribs();

    LastError::lasterror = Errors::UNKNOWN;

    vector<TH2D> fakeMaps;
    getFRmaps(fakeMaps);

    if(LastError::lasterror != Errors::OK){
     cout << "FR maps not found" << endl;
     return;
    }

    
    signalTree->Branch("nJLoc", &nJLoc, "nJLoc/I");
    signalTree->Branch("nBLoc", &nBLoc, "nBLoc/I");
    signalTree->Branch("HTLoc", &HTLoc, "HTLoc/D");
    signalTree->Branch("_met", &_met, "_met/D");

    signalTree->Branch("_weight", &_weightEventInTree, "_weight/D");

    signalTree->Branch("minDeltaR", &minDeltaR, "minDeltaR/D");
    signalTree->Branch("mt", &mtHighest, "mt/D");
    signalTree->Branch("mtlow", &mtLowest, "mtlow/D");

    signalTree->Branch("leadpt", &leadpt, "leadpt/D");
    signalTree->Branch("trailpt", &trailpt, "trailpt/D");
    signalTree->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/D");
    signalTree->Branch("trailJetPt", &trailJetPt, "trailJetPt/D");


    bkgTree->Branch("nJLoc", &nJLoc, "nJLoc/I");
    bkgTree->Branch("nBLoc", &nBLoc, "nBLoc/I");
    bkgTree->Branch("HTLoc", &HTLoc, "HTLoc/D");
    bkgTree->Branch("_met", &_met, "_met/D");

    bkgTree->Branch("_weight", &_weightEventInTree, "_weight/D");

    bkgTree->Branch("minDeltaR", &minDeltaR, "minDeltaR/D");
    bkgTree->Branch("mt", &mtHighest, "mt/D");
    bkgTree->Branch("mtlow", &mtLowest, "mtlow/D");

    bkgTree->Branch("leadpt", &leadpt, "leadpt/D");
    bkgTree->Branch("trailpt", &trailpt, "trailpt/D");
    bkgTree->Branch("leadingJetPt", &leadingJetPt, "leadingJetPt/D");  
    bkgTree->Branch("trailJetPt", &trailJetPt, "trailJetPt/D");


  readerEle->AddVariable( "LepGood_pt", &LepGood_pt );
  readerEle->AddVariable( "LepGood_eta", &LepGood_eta );
  readerEle->AddVariable( "LepGood_jetNDauChargedMVASel", &LepGood_jetNDauChargedMVASel );
  readerEle->AddVariable( "LepGood_miniRelIsoCharged", &LepGood_miniRelIsoCharged );
  readerEle->AddVariable( "LepGood_miniRelIsoNeutral", &LepGood_miniRelIsoNeutral );
  readerEle->AddVariable( "LepGood_jetPtRelv2", &LepGood_jetPtRelv2 );
  readerEle->AddVariable( "min(LepGood_jetPtRatiov2,1.5)", &LepGood_jetPtRatio );
  readerEle->AddVariable( "max(LepGood_jetBTagCSV,0)", &LepGood_jetBTagCSV );
  readerEle->AddVariable( "LepGood_sip3d", &LepGood_sip3d );
  readerEle->AddVariable( "log(abs(LepGood_dxy))", &LepGood_dxy );
  readerEle->AddVariable( "log(abs(LepGood_dz))", &LepGood_dz );
  readerEle->AddVariable( "LepGood_mvaIdSpring16HZZ", &LepGood_mvaIdSpring15 );

  readerEle->BookMVA( "BDTG method", "MVATTHxmlFiles/el_BDTG.weights.xml" );

  readerMu->AddVariable( "LepGood_pt", &LepGood_pt );
  readerMu->AddVariable( "LepGood_eta", &LepGood_eta );
  readerMu->AddVariable( "LepGood_jetNDauChargedMVASel", &LepGood_jetNDauChargedMVASel );
  readerMu->AddVariable( "LepGood_miniRelIsoCharged", &LepGood_miniRelIsoCharged );
  readerMu->AddVariable( "LepGood_miniRelIsoNeutral", &LepGood_miniRelIsoNeutral );
  readerMu->AddVariable( "LepGood_jetPtRelv2", &LepGood_jetPtRelv2 );
  readerMu->AddVariable( "min(LepGood_jetPtRatiov2,1.5)", &LepGood_jetPtRatio );
  readerMu->AddVariable( "max(LepGood_jetBTagCSV,0)", &LepGood_jetBTagCSV );
  readerMu->AddVariable( "LepGood_sip3d", &LepGood_sip3d );
  readerMu->AddVariable( "log(abs(LepGood_dxy))", &LepGood_dxy );
  readerMu->AddVariable( "log(abs(LepGood_dz))", &LepGood_dz );
  readerMu->AddVariable( "LepGood_segmentCompatibility", &LepGood_segmentCompatibility );

  readerMu->BookMVA( "BDTG method", "MVATTHxmlFiles/mu_BDTG.weights.xml" );
    

  for (int i = 0; i != distribsOrderMVATraining.size(); ++i) {

        int sam = distribsOrderMVATraining.at(i);

        if(sam == 0) continue;
        //if(sam == ttWsample) continue;

        //if(sam != 1) continue;
        //if(sam != 0 && sam != 1 && sam != 2 && sam != 7) continue;

        hfile[sam] = new TFile("/Users/illiakhvastunov/Desktop/CERN/MCsamples/80X/ttW_LeptonMVATTH/" + fileList[sam],"read"); // 
        
        hfile[sam]->cd("FakeElectrons");
        inputTreeFake[sam] = static_cast<TTree*>(hfile[sam]->Get("FakeElectrons/fakeTree"));
        
        initBranchAddresses(sam, inputTreeFake[sam]);

        _hCounter->Read("hCounter");
        Double_t scale = xSections[sam]*35900/(_hCounter->GetBinContent(1)); // 12.9, 16.9, 20.0, 27.5
        
        Long64_t nEntries = inputTreeFake[sam]->GetEntries();
        std::cout<<"Entries in "<<fileList[sam]<<" "<<nEntries<<std::endl;
        std::cout<<xSections[sam]<<" "<<_hCounter->GetBinContent(1)<<" "<<scale<<std::endl;

        double allEvents = 0;
        double negEvents = 0;
        
        for (Long64_t it=0; it!=nEntries; ++it) {

            inputTreeFake[sam]->GetEntry(it);

            //if(it > 10000) break;
            
            if (it%100000 == 0)
                cout<<'.'<<flush;

            //if(_eventNb != 455705) continue;
            bool printAddInfo = false;
            if(printAddInfo)
              cout << "NEW EVENT: " << _runNb << " " << _lumiBlock << " " << _eventNb << " " << _nEle + _nMu << endl;
            
            if (_nEle + _nMu < 2) continue;

            if(printAddInfo)
                cout << "Trigger decision: " << _triggers1l[0] << " " << _triggers1l[1] << " " << _triggers1l[3] << endl;

            if(!(_triggers1l[0] || _triggers1l[1] || _triggers1l[3] )) continue;

            nLoc = 0;
            nLocFake = 0;
            nLocLoose = 0;
            nLocEle = 0;
            nLocMu = 0;
            nJLoc = 0;
            nISR  = 0;
            nBLoc = 0;
            HTLoc = 0;

            nLocLoose = 0;

            for(int i = 0; i < _nMu + _nEle; i++){
                if(printAddInfo)
                  cout << "pt of the loose lepton: " << _lPt[i] << endl;
                if(_flavors[i] == 0 && _hitsNumber[i] > 1) continue;
                TLorentzVector ele;
                ele.SetPtEtaPhiE(_lPt[i], _lEta[i], _lPhi[i], _lE[i]);
                bool clean = true;
                for(int j = 0; j < _nMu; j++){
                    if(_flavors[i] == 1) continue;
                    TLorentzVector mu;
                    mu.SetPtEtaPhiE(_lPt[j], _lEta[j], _lPhi[j], _lE[j]);
                    if(printAddInfo)
                    cout << "delta R between electron and muon: " << ele.DeltaR(mu) << endl;
                    if (ele.DeltaR(mu) < 0.05)
                        clean = false;
                }
                if(!clean) continue;

                leptIndLoose[nLocLoose] = i;
                nLocLoose++;

            }

            if(printAddInfo)
            cout << "loose cleaned leptons: " << nLocLoose << endl;

            double deltaMZ = 999999.;
            double minll = 99999.;
            double mll = 99999.;
            int index3 = -1;
            double mtsmallest = 99999.;

            int nonPromptCounter = 0;

            int posCharge = 0;
            int negCharge = 0;

            bool chargeConsistencyEvent = true;

            for (int j=0; j!=nLocLoose; ++j) {

                int i = leptIndLoose[j];

                double leptMVAvalue = -999;

                if(_flavors[i] == 1 ){

                    LepGood_pt = _lPt[i];
                    LepGood_eta = _lEta[i];
                    LepGood_jetNDauChargedMVASel = _trackSelectionMultiplicity[i];
                    LepGood_miniRelIsoCharged = _miniisolationCharged[i][0];
                    LepGood_miniRelIsoNeutral = _miniisolation[i][0] - _miniisolationCharged[i][0];

                    //LepGood_miniRelIsoNeutral = 0.1328779906;
                    LepGood_jetPtRelv2 = _ptrel[i];
                    LepGood_jetPtRatio = TMath::Min(_ptratio[i],1.5);
                    LepGood_jetBTagCSV = TMath::Max(_closeJetCSVAll[i],0.);
                    LepGood_sip3d = _3dIPsig[i]; 
                    LepGood_dxy = TMath::Log(fabs(_ipPV[i]));
                    LepGood_dz = TMath::Log(fabs(_ipZPV[i]));
                    LepGood_segmentCompatibility = _muonSegmentComp[i]; 

                    if(printAddInfo)
                    cout << "Muon info (pt, eta, trackMult, isoCharge, isoNeutral, ptrel, ptratio, btag, sip3d, dxy, dz, muonSegm): " << LepGood_pt << " " << LepGood_eta << " " << LepGood_jetNDauChargedMVASel << " " << LepGood_miniRelIsoCharged << " " << LepGood_miniRelIsoNeutral << " " << LepGood_jetPtRelv2 << " " << LepGood_jetPtRatio << " " << LepGood_jetBTagCSV << " " << LepGood_sip3d << " " << LepGood_dxy << " " << LepGood_dz << " " << LepGood_segmentCompatibility << endl;

                    leptMVAvalue = readerMu->EvaluateMVA( "BDTG method" );

                }

                if(_flavors[i] == 0){

                    LepGood_pt = _lPt[i];
                    LepGood_eta = _lEta[i];
                    LepGood_jetNDauChargedMVASel = _trackSelectionMultiplicity[i];
                    LepGood_miniRelIsoCharged = _miniisolationCharged[i][0];
                    LepGood_miniRelIsoNeutral = _miniisolation[i][0] - _miniisolationCharged[i][0];
                    LepGood_jetPtRelv2 = _ptrel[i];
                    LepGood_jetPtRatio = TMath::Min(_ptratio[i],1.5);
                    LepGood_jetBTagCSV = TMath::Max(_closeJetCSVAll[i],0.);
                    LepGood_sip3d = _3dIPsig[i]; 
                    LepGood_dxy = TMath::Log(fabs(_ipPV[i]));
                    LepGood_dz = TMath::Log(fabs(_ipZPV[i]));
                    LepGood_mvaIdSpring15 = _mvaValue[i];

                    leptMVAvalue = readerEle->EvaluateMVA( "BDTG method" );
                }

                _leptMVAvalue[nLocFake] = leptMVAvalue;
                /*
                if (abs(lep.pdgId)!=13 or lep.mediumMuonId>0) and lep.mvaTTH > 0.90: return lep.pt
                else: return 0.90 * lep.pt / lep.jetPtRatiov2
                */

                double ptFakeTemp = ptFake(_lPt[i], _ptratio[i], _flavors[i], leptMVAvalue, _isloose[i]);

                if(printAddInfo){
                   cout << setprecision(3) << endl;
                   cout << "pt uncor/corr/ptratio/flavour: " << _lPt[i] << " " << ptFakeTemp << " " << _ptratio[i] << " " << _flavors[i] << endl;
                }



                if(ptFakeTemp < 10) continue;
                //if(heavyFlavourLeptons.find(i) != heavyFlavourLeptons.end()) continue;

                if(printAddInfo)
                cout << "PT\tETA\tPHI\tflavour\tiso\tisloose\t3dIP\tcharges\thitsNumb" << endl;

                if(printAddInfo)
                cout << _lPt[i] << "\t" << _lEta[i] << "\t" << _lPhi[i] << "\t" << _flavors[i] << "\t" << _miniisolation[i][0] << "\t" << _isloose[i] << "\t" << _3dIPsig[i] << "\t" << _charges[i] << "\t" << _hitsNumber[i] <<  endl;
          

                if(_flavors[i] == 1){
                  if(printAddInfo)
                  cout << "Muon delta pt/closest jet csv/muon segment compatibility/pt ratio: " << _muonDpt[i] << " " << _closeJetCSVAll[i] << " " << _muonSegmentComp[i] << " " << _ptratio[i] << endl;
                }
                

                if(_flavors[i] == 0){
                  if(printAddInfo)  
                    cout << "Electron info: " << _trigEmulator[i] << " " << _chargeConst[i] << " " << _vtxFitConversion[i] << " " << _hitsNumber[i] << endl;
                }

                if(_flavors[i] == 0 && !(_trigEmulator[i])) continue;
                if(_flavors[i] == 0 && _hitsNumber[i] != 0) continue;

                if(_closeJetCSVAll[i] > 0.8484) continue;

                if(_3dIPsig[i] > 8) continue;

                leptIndFake[nLocFake] = i;
                _isfakeLoc[nLocFake] = true;
                _istightLoc[nLocFake] = false;
                leptFakePtCorr[nLocFake] = ptFakeTemp;
                nLocFake++;
                
                if(_flavors[i] == 0 && _vtxFitConversion[i]) continue;
                if(_flavors[i] == 0 && !_chargeConst[i]) continue;

                if(_flavors[i] == 1 && !_isloose[i]) continue;
                if(_flavors[i] == 1 && _muonDpt[i] > 0.2) continue;
                
                if(printAddInfo)
                  cout << "Lepton passed preselection" << endl;

                if(_flavors[i] == 0) {

                      if(printAddInfo)  
                      cout << "Electron BDT value: " << leptMVAvalue << endl;

                      //myfile << _flavors[i] << endl;
                      //myfile << Form("%1d %9d %12d\t%5.1f %5.1f\t%1d %1.5f %1.5f\t%5.1f %1.3f %1.3f\t%2.3f %1.3f %1.3f\t%1.3f", _runNb, _lumiBlock, _eventNb, _lPt[i], _lEta[i], int(LepGood_jetNDauChargedMVASel), LepGood_miniRelIsoCharged, LepGood_miniRelIsoNeutral, LepGood_jetPtRelv2, LepGood_jetPtRatio, LepGood_jetBTagCSV, LepGood_sip3d, LepGood_dxy, LepGood_dz, LepGood_mvaIdSpring15) << endl;

                      if(leptMVAvalue > 0.5){

                          /*
                          if(!_vtxFitConversion[i] && _chargeConst[i]){

                            leptInd[nLoc] = i;
                            _istightLoc[nLocFake-1] = true;
                            nLoc++;

                          }
                          */
                          /*
                          if(!_vtxFitConversion[i]){

                           if(!_chargeConst[i]) {
                              chargeConsistencyEvent = false;
                              //continue;
                            }
                          */
                            leptInd[nLoc] = i;
                            _istightLoc[nLocFake-1] = true;
                            nLoc++;

                          //}

                      }
                      else if(!(_ptratio[i] > 0.5 && _closeJetCSVAll[i] < 0.3)){
                          
                          nLocFake--;
                          _isfakeLoc[nLocFake] = false;
                          _istightLoc[nLocFake] = false;
                            
                        }
                    }


                    if(_flavors[i] == 1){

                        if(printAddInfo)  
                        cout << "Muon BDT value: " << leptMVAvalue << endl;

                        //myfile << _flavors[i] << endl;
                        //myfile << Form("%1d %9d %12d\t%5.1f %5.1f\t%1d %1.5f %1.5f\t%5.1f %1.3f %1.3f\t%2.3f %1.3f %1.3f\t%1.3f", _runNb, _lumiBlock, _eventNb, LepGood_pt, LepGood_eta, int(LepGood_jetNDauChargedMVASel), LepGood_miniRelIsoCharged, LepGood_miniRelIsoNeutral, LepGood_jetPtRelv2, LepGood_jetPtRatio, LepGood_jetBTagCSV, LepGood_sip3d, LepGood_dxy, LepGood_dz, LepGood_segmentCompatibility) << endl;

                        if(leptMVAvalue > 0.5){

                            /*
                            if(!(_muonDpt[i] < 0.2)){
                              chargeConsistencyEvent = false;
                            }
                            */

                            //if(_isloose[i] && _muonDpt[i] < 0.2){

                              leptInd[nLoc] = i;
                              _istightLoc[nLocFake-1] = true;
                              nLoc++;
                            //}

                        }
                        else if(!(_ptratio[i] > 0.5 && _closeJetCSVAll[i] < 0.3 && _muonSegmentComp[i] > 0.3)){
                            nLocFake--;
                            _isfakeLoc[nLocFake] = false;
                            _istightLoc[nLocFake] = false;
                        }
                    }

            }

            if(printAddInfo){
            //cout << "tight leptons: " << nLocEle + nLocMu + nLocTau << " " << nLocEle << " " << nLocMu << " " << nLocTau << "; Non-prompt: " << nonPromptCounter << endl;
            cout << "tight leptons: " <<  nLoc <<  "; fake leptons: " << nLocFake << "; Pos and Neg: " << posCharge << " " << negCharge << endl;
            cout << "charge consistency conserved: " << chargeConsistencyEvent << endl;
            }

            //if (!chargeConsistencyEvent) continue;

            int samCategory = sam;

            if(nLoc != 2) continue;
            /*
            if (nLoc > 2) continue;
            if (nLoc < 1) continue;

            
            if(nLoc < 2 && nLocFake > 1){
                samCategory = nonpromptSample;
            }
            else if (nLoc == 2){
                samCategory = dataSample;
            }
            else
                continue;
            */
            minll = 99999.;
            
            for (int l0 = 0; l0<nLocFake; ++l0) {
                l0p4.SetPtEtaPhiE(leptFakePtCorr[leptIndFake[l0]] ,_lEta[leptIndFake[l0]],_lPhi[leptIndFake[l0]],_lE[leptIndFake[l0]] * leptFakePtCorr[leptIndFake[l0]] / _lPt[leptIndFake[l0]]);
                  
                for(int l1 = l0 + 1; l1 < nLocFake; ++l1){
                    if (_charges[leptIndFake[l0]] != _charges[leptIndFake[l1]]) {
                        
                        l1p4.SetPtEtaPhiE(leptFakePtCorr[leptIndFake[l1]],_lEta[leptIndFake[l1]],_lPhi[leptIndFake[l1]],_lE[leptIndFake[l1]] *  leptFakePtCorr[leptIndFake[l1]] / _lPt[leptIndFake[l1]]);
                        l1p4+=l0p4;
                        double mdiL = l1p4.M();
                        
                        if(printAddInfo)
                        cout << "mdiL: " << mdiL << endl;
                        if (_flavors[leptIndFake[l0]] == _flavors[leptIndFake[l1]] ) {
                            if (mdiL < minll) minll = mdiL;
                        }
                    }
                }
            }

            if(minll < 12) continue;

            double maxPt = 0.;
            double max2ndPt = 0.;

            int maxPtInd;
            int max2ndPtInd;

            vector<std::pair<double,int> > leptonPt;

            if(samCategory != nonpromptSample){
                for(int i = 0; i < nLoc; i++)
                    leptonPt.push_back(std::make_pair(_lPt[leptInd[i]], leptInd[i]));
            }
            else{
                for(int i = 0; i < nLocFake; i++)
                    leptonPt.push_back(std::make_pair(leptFakePtCorr[leptIndFake[i]], leptIndFake[i]));
            }


            if(printAddInfo)
            cout << "Vector size: " << leptonPt.size() << endl;
            
            sort (leptonPt.begin(), leptonPt.end(), comp);

            maxPt = leptonPt.at(0).first;
            max2ndPt = leptonPt.at(1).first;
            //minPt = leptonPt.at(2).first;

            maxPtInd = leptonPt.at(0).second;
            max2ndPtInd = leptonPt.at(1).second;
            //minPtInd = leptonPt.at(2).second;

            int maxInd[2] = {maxPtInd, max2ndPtInd};

            int chargeTwoLeptons = _charges[maxPtInd] + _charges[max2ndPtInd];
            
            if(chargeTwoLeptons != -2 && chargeTwoLeptons != 2) continue;


            nLocMu = (_flavors[maxPtInd] ? 1 : 0) + (_flavors[max2ndPtInd] ? 1 : 0); 
            nLocEle = 2 - nLocMu;
            
            

            bool isPromptMaxPt = _originReduced[maxPtInd] == 0;
            bool isPrompt2ndMaxPt = _originReduced[max2ndPtInd] == 0;

            bool isFakeMaxPt = _originReduced[maxPtInd] != 0;
            bool isFake2ndMaxPt = _originReduced[max2ndPtInd] != 0;

            bool prompt2leptons = isPromptMaxPt && isPrompt2ndMaxPt;
            bool fakeDecision = !prompt2leptons;
            
            //if(isFakeMaxPt && isFake2ndMaxPt) continue;
            
            if(sam == ttWsample && !prompt2leptons) continue;
            if(sam != ttWsample && prompt2leptons) continue;
            

            //if(sam == 1 && !prompt3leptons) continue;
            /*
            if (sam > 0 && sam < nSamples-1 && samCategory == 0 && !prompt3leptons)
                continue;
            */
            
            if(maxPt < 25 && _flavors[maxPtInd] == 1) continue;
            if(max2ndPt < 25 && _flavors[max2ndPtInd] == 1) continue;

            if(maxPt < 40 && _flavors[maxPtInd] == 0) continue;
            if(max2ndPt < 27 && _flavors[max2ndPtInd] == 0) continue;
            

            
            //double deltaRtrailLepLeadJet = -999.;
            double minDeltaRjetLep = 999.;
              
            if(printAddInfo)
            cout << "Leptons pt selection passed" << endl;

            int bQuarksCounter = 0;
            int cQuarksCounter = 0;
            int lightQuarksCounter = 0;

            double minDeltaRLoc = 9999.;

            for(int i = 0; i < _n_Jets; i++){
              if(printAddInfo)
                cout << "Jet info: " << _jetPt[i] << " " << _csv[i] << endl;
              if(_jetPt[i] > 30){

                vector<TLorentzVector> lep;
                for(unsigned int k = 0; k < nLocFake; k++){
                //for(unsigned int k = 0; k < nLocLoose; k++){
                  TLorentzVector lepTemp;
                  lepTemp.SetPtEtaPhiE(_lPt[leptIndFake[k]], _lEta[leptIndFake[k]], _lPhi[leptIndFake[k]], _lE[leptIndFake[k]]);
                  //lepTemp.SetPtEtaPhiE(_lPt[leptIndLoose[k]], _lEta[leptIndLoose[k]], _lPhi[leptIndLoose[k]], _lE[leptIndLoose[k]]);
                  lep.push_back(lepTemp);
                }

                TLorentzVector jet;
                jet.SetPtEtaPhiE(_jetPt[i], _jetEta[i], _jetPhi[i], _jetE[i]);

                bool jetClean = true;

                //for(int j = 0; j < nLocLoose; j++){
                for(int j = 0; j < nLocFake; j++){
                  if(printAddInfo)
                  cout << "delta R between jet and lepton: " << jet.DeltaR(lep.at(j)) << "; lepton pt is: " << lep.at(j).Pt() << endl;
                  if(jet.DeltaR(lep.at(j)) < minDeltaRjetLep) 
                    minDeltaRjetLep = jet.DeltaR(lep.at(j));
                  if(jet.DeltaR(lep.at(j)) < 0.4){
                    jetClean = false;
                    break;
                  }
                    
                }

                if(!jetClean) continue;

                double deltaR_jet_trailLep = jet.DeltaR(lep.at(1));
                if(deltaR_jet_trailLep < minDeltaRLoc)
                  minDeltaRLoc = deltaR_jet_trailLep;

                jetLooseInd[nJLoc] = i;
                nJLoc++;

                if(_csv[i] > 0.8484)
                  nBLoc++;

                HTLoc = HTLoc + _jetPt[i];

                if(fabs(_jetFlavour[i]) == 5) bQuarksCounter++;
                if(fabs(_jetFlavour[i]) == 4) cQuarksCounter++;
                if(fabs(_jetFlavour[i]) == 0) lightQuarksCounter++;

                

              }
            }

            minDeltaR = minDeltaRLoc;

            if(printAddInfo)
                cout << "Total number of jets: " << nJLoc << endl;


            double dataMCSF = 1.;
            double lepSF = 1.;


            double tempValue = (double) rand() / (RAND_MAX);
            int leptonFileDicision = -99;

            if(tempValue < 20./35.9)
                leptonFileDicision = 3;
            else
                leptonFileDicision = 5;  

            
            dataMCSF = 1;
            lepSF = 1;
            //dataMCSF *= lepSF;

            deltaMZ = 999999.;
            minll = 99999.;
            mll = 99999.;
            double mlll = 99999.;
            index3 = -1;
            mtsmallest = 99999.;

            double pt_Z = -99999.;
            

            for (int l0 = 0; l0<nLocFake; ++l0) {
                double ptfakeLep1 = leptFakePtCorr[leptIndFake[l0]]; // ptFake(_lPt[maxInd[l0]], _lEta[maxInd[l0]], _flavors[maxInd[l0]], _flavors[maxInd[l0]] ? _isolationDB[maxInd[l0]] : _isolation[maxInd[l0]]);
                //double ptfakeLep1 = _lPt[maxInd[l0]];
                l0p4.SetPtEtaPhiE(ptfakeLep1 ,_lEta[leptIndFake[l0]],_lPhi[leptIndFake[l0]],_lE[leptIndFake[l0]] * ptfakeLep1 / _lPt[leptIndFake[l0]]);
                  
                for(int l1 = l0 + 1; l1 < nLocFake; ++l1){
                    if (_charges[leptIndFake[l0]] != _charges[leptIndFake[l1]]) {
                        double ptfakeLep2 = leptFakePtCorr[leptIndFake[l1]]; //ptFake(_lPt[maxInd[l1]], _lEta[maxInd[l1]], _flavors[maxInd[l1]], _flavors[maxInd[l1]] ? _isolationDB[maxInd[l1]] : _isolation[maxInd[l1]]);
                        //double ptfakeLep2 = _lPt[maxInd[l1]];
                        l1p4.SetPtEtaPhiE(ptfakeLep2,_lEta[leptIndFake[l1]],_lPhi[leptIndFake[l1]],_lE[leptIndFake[l1]] *  ptfakeLep2 / _lPt[leptIndFake[l1]]);
                        l1p4+=l0p4;
                        double mdiL = l1p4.M();
                        
                        if(printAddInfo)
                        cout << "mdiL: " << mdiL << endl;
                        if (_flavors[leptIndFake[l0]] == _flavors[leptIndFake[l1]] ) {
                            if (mdiL < minll) minll = mdiL;
                            if (fabs(mdiL - 91.2) < deltaMZ) {
                                deltaMZ = fabs(mdiL - 91.2);
                                mll = mdiL;
                                pt_Z = l1p4.Pt();
                                /*
                                index3 = maxInd[2 - l0 - l1];
                                //double ptfakeLep3 = ptFake(_lPt[index3], _lEta[index3], _flavors[index3], _flavors[index3] ? _isolationDB[index3] : _isolation[index3]);
                                double ptfakeLep3 = _lPt[index3];
                                TLorentzVector l2p4;
                                l2p4.SetPtEtaPhiE(ptfakeLep3,_lEta[index3],_lPhi[index3],_lE[index3] *  ptfakeLep3 / _lPt[index3]);
                                l1p4+= l2p4;
                                mlll = l1p4.M();
                                
                                if(printAddInfo)
                                cout << "l0 and l1: " << l0 << " " << l1 << endl;
                                if(printAddInfo)
                                cout << "Index and deltaMZ: " << index3 << " " << deltaMZ << endl;
                                */
                            }
                        }
                    }
                }
            }


            if(printAddInfo)
            cout << "deltaMZ: " << deltaMZ << endl;
            
                          
            double btagSF_event = 1.;
            double btagSF_event_Up = 1.;
            double btagSF_event_Down = 1.;

            double btagSF_event_light = 1.;
            double btagSF_event_Up_light = 1.;
            double btagSF_event_Down_light = 1.;


            if(printAddInfo)
            cout << "MET Values/Number b quarks: " << _met << " " << nBLoc << endl;

            
            if(nJLoc < 2) continue;
            if(nLocEle == 2 && deltaMZ < 10) continue;
            if(_met < 30) continue;
            if(nBLoc < 1) continue;

            vector<std::pair<double,int> > jetPt;
            for(int i = 0; i < nJLoc; i++)
              jetPt.push_back(std::make_pair(_jetPt[jetLooseInd[i]], jetLooseInd[i]));

            sort (jetPt.begin(), jetPt.end(), comp); 

            if(nJLoc > 0)
                leadingJetPt = jetPt[0].first;
            if(nJLoc > 1)
                trailJetPt = jetPt[1].first;

            mtHighest = _mt[maxPtInd] > _mt[max2ndPtInd] ? _mt[maxPtInd] : _mt[max2ndPtInd];
            mtLowest = _mt[maxPtInd] < _mt[max2ndPtInd] ? _mt[maxPtInd] : _mt[max2ndPtInd];

            leadpt = maxPt;
            trailpt = max2ndPt;

            
            
            _weightEventInTree = scale * TMath::Abs(_weight);
            if(sam == ttWsample)
              signalTree->Fill();
            else
              bkgTree->Fill();
              

            double FRloc = 1.;
            int nFakeLepCounter = 0;

            if(samCategory == 0){ 

                FRloc = 0.;

                for(int i = 0; i < nLocFake; i++){
                  if(_istightLoc[leptIndFake[i]]) continue;
                  double ptf = leptFakePtCorr[leptIndFake[i]]; // ptFake(_lPt[maxInd[i]], _ptratio[maxInd[i]], _flavors[maxInd[i]], _mvaValue[maxInd[i]], _isloose[maxInd[i]]);      
                  double FRloc_loc = fakeMaps.at(_flavors[leptIndFake[i]]).GetBinContent(fakeMaps.at(_flavors[leptIndFake[i]]).FindBin(TMath::Min(ptf,ptBins[nPt-1]-1.), fabs(_lEta[leptIndFake[i]])));
                  if(printAddInfo)
                    cout << "FR from the map " << FRloc_loc << endl;
                  //FRloc += FRloc_loc/(1.-FRloc_loc);
                  FRloc += FRloc_loc;
                  nFakeLepCounter++;

                  if(printAddInfo)
                    cout << "FR loc is " << FRloc << endl;

                }
                
                FRloc *= TMath::Power(-1, nFakeLepCounter + 1);
            }

            if(printAddInfo)
            cout << "The b-tag SF: " << btagSF_event << endl;


            if(nLocMu == 2 && nLocEle == 0 && chargeTwoLeptons == -2)
              distribs[13].vectorHisto[samCategory].Fill(1.,scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            if(nLocMu == 1 && nLocEle == 1 && chargeTwoLeptons == -2)
              distribs[13].vectorHisto[samCategory].Fill(2.,scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            if(nLocMu == 0 && nLocEle == 2 && chargeTwoLeptons == -2)
              distribs[13].vectorHisto[samCategory].Fill(3.,scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            
           if(nLocMu == 2 && nLocEle == 0 && chargeTwoLeptons == 2)
              distribs[13].vectorHisto[samCategory].Fill(4.,scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            if(nLocMu == 1 && nLocEle == 1 && chargeTwoLeptons == 2)
              distribs[13].vectorHisto[samCategory].Fill(5.,scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            if(nLocMu == 0 && nLocEle == 2 && chargeTwoLeptons == 2)
              distribs[13].vectorHisto[samCategory].Fill(6.,scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            
            //distribs[17].vectorHisto[samCategory].Fill(SRID(nJLoc, nBLoc),scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
                
            distribs[7].vectorHisto[samCategory].Fill(TMath::Min(double(nJLoc),varMax[9]-0.1),scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            distribs[8].vectorHisto[samCategory].Fill(TMath::Min(double(nBLoc),varMax[8]-0.1),scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            
            distribs[6].vectorHisto[samCategory].Fill(TMath::Min(double(HTLoc),varMax[6]-0.1),scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            
            distribs[0].vectorHisto[samCategory].Fill(TMath::Min(maxPt,varMax[0]-0.1),scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            distribs[1].vectorHisto[samCategory].Fill(TMath::Min(max2ndPt,varMax[1]-0.1),scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            
            distribs[2].vectorHisto[samCategory].Fill(TMath::Min(mtHighest,varMax[2]-0.1),scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            distribs[3].vectorHisto[samCategory].Fill(TMath::Min(mtLowest,varMax[3]-0.1),scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            
            distribs[5].vectorHisto[samCategory].Fill(TMath::Min(_met,varMax[5]-0.1),scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            distribs[12].vectorHisto[samCategory].Fill(TMath::Min(double(_n_PV),varMax[12]-0.1),scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            //distribs[6].vectorHisto[samCategory].Fill(TMath::Min(deltaRmin,varMax[6]-0.1),scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);

            distribs[9].vectorHisto[samCategory].Fill(TMath::Min(_jetPt[jetLooseInd[0]],varMax[9]-0.001),scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            distribs[10].vectorHisto[samCategory].Fill(TMath::Min(_jetPt[jetLooseInd[1]],varMax[10]-0.001),scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            

            distribs[4].vectorHisto[samCategory].Fill(TMath::Min(minDeltaR,varMax[4]-0.001),scale*_weight*dataMCSF*FRloc*btagSF_event*btagSF_event_light);
            
            //if(nFakeLepCounter == 1){
            //if(sam == samCategory){
            allEvents += dataMCSF*FRloc*btagSF_event*btagSF_event_light;
             if(_weight < 0)
                negEvents += dataMCSF * FRloc * btagSF_event * btagSF_event_light;
            //}
  
        }
        
        //_weight = 1.;
        std::cout << "Scale and weight: " << scale << " " << _weight << std::endl; 
        std::cout << "Total nEv: " << (allEvents - 2 *negEvents) * scale * TMath::Abs(_weight) << " ; negEvents: " << negEvents * scale * TMath::Abs(_weight) << std::endl;

        weights[sam] = scale*TMath::Abs(_weight);

        cout << endl;
       
    }

    fileDummy->cd();
    signalTree->Write();
    bkgTree->Write();

    fileDummy->Close();


    //cout << "D weight is: " << h_nISRbefore->Integral() / h_nISRafter->Integral() << endl;


  // 0 - fakes, 1-2 ttZ, 3-WZ, 4-9 ttX, 10-21 rares, 22 - data

    TLegend* mtleg = new TLegend(0.25,0.89,0.95,0.77); 
    mtleg->SetNColumns(4);
    mtleg->SetFillColor(0);
    mtleg->SetFillStyle(0);
    mtleg->SetBorderSize(0);
    mtleg->SetTextFont(42);
    
    mtleg->AddEntry(&distribs[8].vectorHisto[nSamples-1],names[nSamples-1],"ep"); //data
    /*
    for (std::vector<int>::reverse_iterator it = distribsOrder.rbegin(); it != distribsOrder.rend(); it++) {

        int i = *it;  //  = {1, 2, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 0};
        //cout << i << " " << j << endl;
        if(i == 2) continue;
        if(i > 4 && i < 11) continue; // ttX
        //if(i > 4 && i < 7) continue; // ttX sym
        //if(i > 7 && i < 11) continue; // ttX as
        if(i > 11 && i < nSamples-1) continue;
        mtleg->AddEntry(&distribs[8].vectorHisto[i],names[i],"f");
    }
    */

    // 0 - fakes, 1-2 ttZ, 3,4,5 - WZ, 6-12 ttX, 13-23 rares, 24 - data
    
    for (std::vector<int>::reverse_iterator it = distribsOrder.rbegin(); it != distribsOrder.rend(); it++) {

        int i = *it;  //  = {1, 2, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 0};
        //cout << i << " " << j << endl;
        if(i == 5) continue;
        if(i > 9 && i < 20) continue;
        if(i > 20 && i < nSamples-1) continue;
        mtleg->AddEntry(&distribs[4].vectorHisto[i],names[i],"f");
    }
    
    //mtleg->AddEntry(&distribs[8].vectorHisto[nSamples], "uncertainties" ,"f");

    for (int i=0; i!=nVars; ++i)  {
        for(unsigned int sam = 0; sam != distribsOrder.size(); sam++)
          distribs[i].vectorHisto[uncertaintySample].Add(&distribs[i].vectorHisto[distribsOrder[sam]]);
        
        for (int ibin = 1; ibin!=nBins[i]+1; ++ibin) {
          //get syst. uncertainty band:
          double err = 0.;
          for (int sam=0; sam!=nSamples-1; ++sam) {
            if(distribs[i].vectorHisto[sam].GetBinContent(ibin) != 0){
                err += TMath::Power(distribs[i].vectorHisto[sam].GetBinContent(ibin) * systematics[sam], 2) + TMath::Power(distribs[i].vectorHisto[sam].GetBinError(ibin), 2);
                //err += TMath::Power(distribs[i].vectorHisto[sam].GetBinError(ibin), 2);
            }
            else{
                err += 0;
                //err += TMath::Power(weights[sam], 2);
            }
          }
            
          err = sqrt(err);
          distribs[i].vectorHisto[uncertaintySample].SetBinError(ibin, err); 
        }   

    }



    double scale_num = 1.6;
    
    TCanvas* plot[14];
      
    for(int i = 0; i < 14; i++){
      plot[i] = new TCanvas(Form("plot_%d", i),"",500,450);
    }

    plot[0]->cd();
    plot[0]->SetLogy();
    showHist(plot[0],distribs[6],"","H_{T} [GeV]","Events / " + std::to_string(int((varMax[6] - varMin[6])/nBins[6])) + " GeV",scale_num, mtleg);
    
    plot[1]->cd();
    plot[1]->SetLogy();
    showHist(plot[1],distribs[7],"","N_{jets}","Events",scale_num, mtleg);
        
    plot[2]->cd();
    showHist(plot[2],distribs[8],"","N_{b jets}","Events",scale_num, mtleg);
    
    plot[3]->cd();
    showHist(plot[3],distribs[0],"","Leading lepton p_{T} [GeV]","Events / "  + std::to_string(int((varMax[0] - varMin[0])/nBins[0])) + " GeV",scale_num, mtleg);
    
    plot[4]->cd();
    showHist(plot[4],distribs[1],"","Trailing lepton p_{T} [GeV]","Events / " + std::to_string(int((varMax[1] - varMin[1])/nBins[1])) + " GeV",scale_num, mtleg);
    
    plot[5]->cd();
    showHist(plot[5],distribs[5],"","E_{T}^{miss} [GeV]","Events / " + std::to_string(int((varMax[5] - varMin[5])/nBins[5])) + " GeV",scale_num, mtleg);

    plot[6]->cd();
    showHist(plot[6],distribs[2],"","M_{T}^{leading} [GeV]","Events / " + std::to_string(int((varMax[2] - varMin[2])/nBins[2])) + " GeV",scale_num, mtleg);

    plot[7]->cd();
    showHist(plot[7],distribs[3],"","M_{T}^{trailing} [GeV]","Events / " + std::to_string(int((varMax[3] - varMin[3])/nBins[3])) + " GeV",scale_num, mtleg);

    plot[8]->cd();
    showHist(plot[8],distribs[12],"","PU","Events",scale_num, mtleg);
    
    plot[9]->cd();
    showHist(plot[9],distribs[9],"","leading jet p_{T} [GeV]","Events / "  + std::to_string(int((varMax[9] - varMin[9])/nBins[9])) + " GeV",scale_num, mtleg);

    plot[10]->cd();
    showHist(plot[10],distribs[10],"","subleading jet p_{T} [GeV]","Events / "  + std::to_string(int((varMax[10] - varMin[10])/nBins[10])) + " GeV",scale_num, mtleg);

    plot[11]->cd();
    showHist(plot[11],distribs[11],"","BDT score","Events / "  + std::to_string(int((varMax[11] - varMin[11])/nBins[1])),scale_num, mtleg);

    plot[12]->cd();
    showHist(plot[12],distribs[13],"","flavour","Events" ,scale_num, mtleg);

    plot[13]->cd();
    showHist(plot[13],distribs[14],"","SR","Events" ,scale_num, mtleg);

    /*
    TFile * btagEffFile = new TFile("btagEff_medium.root", "RECREATE");
    TString letter[3] = {"L", "C", "B"};
    TCanvas *c1[3];
    for(unsigned int i = 0; i < 3; i++){
        c1[i] = new TCanvas(Form("c1_%d",i), Form("c1_%d",i));
        h2_effBtagSF[i][1]->Divide(h2_effBtagSF[i][0]);
        h2_effBtagSF[i][1]->Draw("text");
        h2_effBtagSF[i][1]->SetName("histo" + letter[i]);
        h2_effBtagSF[i][1]->Write();
    }
    btagEffFile->Close();
    */
    /*
    TFile *fileTTZ = new TFile("NonP_GH.root", "RECREATE");
    vector<int> writeTTZ = {8, 9, 3, 17, 28, 0, 1, 7, 4, 19, 31, 18, 32, 6, 29, 30, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46, 47, 49, 50, 51, 52, 53, 26, 27};
    vector<string> writeNamesTTZ = {"HT", "njets", "mll", "SR", "nbjets", "leading", "subleading", "met", "mt", "flavour", "pu", "pttrail", "invar3", "minDeltaR", "csv1st", "csv2nd", "csv3rd", "csv4th", "csv5th", "csv6th", "iso1Ele", "iso2Ele", "iso3Ele", "iso1Mu", "iso2Mu", "iso3Mu", "deltaRjetLepton", "ipPV", "ipZPV", "3Dsig", "isoEleH", "isoMuH", "ptjetlead", "ptjettrail"};
    for(unsigned int i = 0; i < writeTTZ.size(); i++){
        Double_t norm = distribs[writeTTZ.at(i)].vectorHisto[0].GetEntries();
        //distribs[writeTTZ.at(i)].vectorHisto[nSamples-1].Scale(1/norm);
        distribs[writeTTZ.at(i)].vectorHisto[0].SetName(writeNamesTTZ.at(i).c_str());
        distribs[writeTTZ.at(i)].vectorHisto[0].SetTitle(writeNamesTTZ.at(i).c_str());
        distribs[writeTTZ.at(i)].vectorHisto[0].Write();
    }
    fileTTZ->Close();
    return;
    */
/*
    TString folder = "plots/FullData/";

    plot[0]->SaveAs(folder + "HT.pdf");
    plot[1]->SaveAs(folder + "njets.pdf");
    plot[2]->SaveAs(folder + "mll.pdf");
    plot[3]->SaveAs(folder + "moneyPlot.pdf");
    plot[4]->SaveAs(folder + "nbjetsmedium.pdf");
    plot[5]->SaveAs(folder + "ptlead.pdf");
    plot[6]->SaveAs(folder + "ptsublead.pdf");
    plot[7]->SaveAs(folder + "met.pdf");
    plot[9]->SaveAs(folder + "flavour.pdf");
    plot[11]->SaveAs(folder + "pttrail.pdf");

    plot[22]->SaveAs(folder + "ptZ.pdf");
    plot[23]->SaveAs(folder + "jetptlead.pdf");
    plot[24]->SaveAs(folder + "jetptsublead.pdf");
    plot[30]->SaveAs(folder + "ptNonZ.pdf");


    plot[0]->SaveAs(folder + "HT.png");
    plot[1]->SaveAs(folder + "njets.png");
    plot[2]->SaveAs(folder + "mll.png");
    plot[3]->SaveAs(folder + "moneyPlot.png");
    plot[4]->SaveAs(folder + "nbjetsmedium.png");
    plot[5]->SaveAs(folder + "ptlead.png");
    plot[6]->SaveAs(folder + "ptsublead.png");
    plot[7]->SaveAs(folder + "met.png");
    plot[9]->SaveAs(folder + "flavour.png");
    plot[11]->SaveAs(folder + "pttrail.png");

    plot[22]->SaveAs(folder + "ptZ.png");
    plot[23]->SaveAs(folder + "jetptlead.png");
    plot[24]->SaveAs(folder + "jetptsublead.png");
    plot[30]->SaveAs(folder + "ptNonZ.png");
*/
    return;
    std::cout<<"Done 2"<<std::endl;

// ________________________________________

    // HAVE TO REDO IT COMPETELY

    //const int SRNumber = 15;

    int indSR[SRNumber];

    for(int i = 0; i < SRNumber; i++){
      indSR[i] = i+1;
    }

    double yields[10][SRNumber] = {0.}; // 0 - fakes, 1-2 ttZ, 3-WZ, 4-16 ttX (ttH-4, ttgamma - 5, tttt-6, ttW-7), 17-28 (11-ZZ) rares, 29 - data
    double yieldsTotal[10] = {0.};  // 0 - fakes, 1-2 ttZ, 3-WZ, 4-16 ttX (ttH-4, ttW - 5, tttt-6, ttgamma-7), 17-28 (11-ZZ) rares, 29 - data
    double yieldsErrs[10][SRNumber ] = {0.};
    double yieldsErrsSystematics[10][SRNumber ] = {0.};
    double yieldsTotalErrs[10] = {0.};
    double yieldsTotalErrsSystematics[10] = {0.};

    double yieldsTotalBackground = 0;
    double yieldsTotalBackgroundErrs = 0;

    double yieldsTotalTTW = 0;
    double yieldsTotalTTWErrs = 0;

    double yieldsTotalTTZ = 0;
    double yieldsTotalTTZErrs = 0;

    double yieldsTotalAll = 0;
    double yieldsTotalAllErrs = 0;

    double yieldsTotalObserved;

    int yInd[10] = {0,1,3,4,5,6,17,18,29,30};
    for (int i=0; i!=9; ++i) {
        for (int k=0; k!=SRNumber; ++k) {
            yields[i][k] = 0;
            yieldsErrs[i][k] = 0;
            yieldsErrsSystematics[i][k] = 0;
            yieldsTotalErrs[k] = 0; 
            yieldsTotalErrsSystematics[k] = 0; 
        }
    }

    
    for (int i=0; i!=9; ++i) {
        for (int j=yInd[i]; j!=yInd[i+1]; ++j) {

            for (int k=0; k!=SRNumber; ++k) {

                yields[i][k] += distribs[17].vectorHisto[j].GetBinContent(indSR[k]);
                yieldsErrs[i][k] += distribs[17].vectorHisto[j].GetBinError(indSR[k]) * distribs[17].vectorHisto[j].GetBinError(indSR[k]);
                yieldsErrsSystematics[i][k] += TMath::Power(distribs[17].vectorHisto[j].GetBinContent(indSR[k]) * systematics[j], 2);
                
                if(i == 8)
                  yieldsTotalObserved += distribs[17].vectorHisto[j].GetBinContent(indSR[k]);

                if(i < 8 && i != 1 && i != 3){
                  yields[9][k] += distribs[17].vectorHisto[j].GetBinContent(indSR[k]);
                  yieldsErrs[9][k] += distribs[17].vectorHisto[j].GetBinError(indSR[k]) * distribs[17].vectorHisto[j].GetBinError(indSR[k]); 
                  yieldsErrsSystematics[9][k] += TMath::Power(distribs[17].vectorHisto[j].GetBinContent(indSR[k]) * systematics[j], 2);
                  
                  yieldsTotalBackground += distribs[17].vectorHisto[j].GetBinContent(indSR[k]);
                  yieldsTotalBackgroundErrs += distribs[17].vectorHisto[j].GetBinError(indSR[k]) * distribs[17].vectorHisto[j].GetBinError(indSR[k]) + TMath::Power(distribs[17].vectorHisto[j].GetBinContent(indSR[k]) * systematics[j], 2); 
                } 

                if(i == 1){
                  yieldsTotalTTZ += distribs[17].vectorHisto[j].GetBinContent(indSR[k]);
                  yieldsTotalTTZErrs += distribs[17].vectorHisto[j].GetBinError(indSR[k]) * distribs[17].vectorHisto[j].GetBinError(indSR[k]) + TMath::Power(distribs[17].vectorHisto[j].GetBinContent(indSR[k]) * systematics[j], 2); 
                }

                if(i == 3){
                  yieldsTotalTTW += distribs[17].vectorHisto[j].GetBinContent(indSR[k]);
                  yieldsTotalTTWErrs += distribs[17].vectorHisto[j].GetBinError(indSR[k]) * distribs[17].vectorHisto[j].GetBinError(indSR[k]) + TMath::Power(distribs[17].vectorHisto[j].GetBinContent(indSR[k]) * systematics[j], 2); 
                }

                if(i < 8){
                    yieldsTotal[k] += distribs[17].vectorHisto[j].GetBinContent(indSR[k]);
                    yieldsTotalErrs[k] += distribs[17].vectorHisto[j].GetBinError(indSR[k]) * distribs[17].vectorHisto[j].GetBinError(indSR[k]);
                    yieldsTotalErrsSystematics[k] += TMath::Power(distribs[17].vectorHisto[j].GetBinContent(indSR[k]) * systematics[j], 2);
                
                    yieldsTotalAll += distribs[17].vectorHisto[j].GetBinContent(indSR[k]);
                    yieldsTotalAllErrs += distribs[17].vectorHisto[j].GetBinError(indSR[k]) * distribs[17].vectorHisto[j].GetBinError(indSR[k]) + TMath::Power(distribs[17].vectorHisto[j].GetBinContent(indSR[k]) * systematics[j], 2); 
                
                }
            }
        }
    }
    

    /*
    for (int i=0; i!=9; ++i) {
        for (int j=yInd[i]; j!=yInd[i+1]; ++j) {

            for (int k=0; k!=SRNumberZpt; ++k) {

                yields[i][k] += distribs[48].vectorHisto[j].GetBinContent(indSR[k]);
                yieldsErrs[i][k] += distribs[48].vectorHisto[j].GetBinError(indSR[k]) * distribs[48].vectorHisto[j].GetBinError(indSR[k]);
                yieldsErrsSystematics[i][k] += TMath::Power(distribs[48].vectorHisto[j].GetBinContent(indSR[k]) * systematics[j], 2);
                if(i < 8 && i != 1 && i != 3){
                  yields[9][k] += distribs[48].vectorHisto[j].GetBinContent(indSR[k]);
                  yieldsErrs[9][k] += distribs[48].vectorHisto[j].GetBinError(indSR[k]) * distribs[48].vectorHisto[j].GetBinError(indSR[k]); 
                  yieldsErrsSystematics[9][k] += TMath::Power(distribs[48].vectorHisto[j].GetBinContent(indSR[k]) * systematics[j], 2);
                } 

                if(i < 8){
                    yieldsTotal[k] += distribs[48].vectorHisto[j].GetBinContent(indSR[k]);
                    yieldsTotalErrs[k] += distribs[48].vectorHisto[j].GetBinError(indSR[k]) * distribs[48].vectorHisto[j].GetBinError(indSR[k]);
                    yieldsTotalErrsSystematics[k] += TMath::Power(distribs[48].vectorHisto[j].GetBinContent(indSR[k]) * systematics[j], 2);
                }
            }
        }
    }
    */


    
    string NameCategory[3] = {"0 b-tag jets", "1 b-tag jet", "$\\geq$ 2 b-tag jets"};
    string NameBKG[8] = {"nonprompt", "WZ", "ttX", "rare", "ttZ", "observed", "total", "total background"};
    string NameProcess[3] = {"2jets", "3jets", "4jets"};

    for(int cat = 0; cat < 3; cat++){

      ofstream tableBkg;
      tableBkg.open("tableFullAnalysis_CBshort_" + std::to_string(cat) + ".tex");
      tableBkg<<"\\begin{table}\n";
      tableBkg<<"\\begin{adjustbox}{width=1\\textwidth}\n";
      tableBkg<<"\\begin{tabular}{||c||c|c|c|c|c|c||}\\hline\n"; 
      tableBkg << std::fixed << setprecision(1) << "\n";

      tableBkg << "& \\multicolumn{3}{c|}{" + NameCategory[cat] << "} \\\\ \\hline \n";

      tableBkg << "Process";

      for(int num = 0; num < 3; num++){
        tableBkg << " & " << NameProcess[num] ;
      }

      tableBkg << "\\\\ \\hline\n";

      tableBkg << NameBKG[7] << " "; 
      for(int num = 3*cat; num < 3*(cat + 1); num++){
        tableBkg << " & $" << yields[9][num] << "\\pm"  << TMath::Sqrt(yieldsErrs[9][num] + yieldsErrsSystematics[9][num]) << "$  ";
      }

      tableBkg << " \\\\ \n";

      tableBkg << NameBKG[4] << " "; 
      for(int num = 3*cat; num < 3*(cat + 1); num++){
        tableBkg << " & $" << yields[1][num] << "\\pm"  << TMath::Sqrt(yieldsErrs[1][num] + yieldsErrsSystematics[1][num]) << "$  ";
      }

      tableBkg << " \\\\ \\hline \\hline \n";

      tableBkg << "ttW" << " "; 
      for(int num = 3*cat; num < 3*(cat + 1); num++){
        tableBkg << " & $" << yields[3][num] << "\\pm"  << TMath::Sqrt(yieldsErrs[3][num] + yieldsErrsSystematics[3][num]) << "$  ";
      }

      tableBkg << " \\\\ \\hline \\hline \n";

      tableBkg << NameBKG[6] << " "; 
      for(int num = 3*cat; num < 3*(cat + 1); num++){
        tableBkg << " & $" << yieldsTotal[num] << "\\pm" << TMath::Sqrt(yieldsTotalErrs[num] + yieldsTotalErrsSystematics[num]) << "$  ";
      }

      tableBkg << " \\\\ \\hline \\hline \n";

      tableBkg << NameBKG[5] << " "; 
      for(int num = 3*cat; num < 3*(cat + 1); num++){
        tableBkg << " & " << int(yields[8][num]); 
      }

      tableBkg << " \\\\ \\hline\n";

      tableBkg <<"\\end{tabular}\n";
      tableBkg <<"\\end{adjustbox}\n";
      tableBkg <<"\\end{table}\n";
      tableBkg.close();
    }

    /*
    // 0 - fakes, 1-2 ttZ, 3-WZ, 4-10 ttX (ttH-4, ttW-5), 11-22 (11-ZZ) rares, 23 - data
    string NameBKGextended[11] = {"nonprompt", "ttZ", "WZ", "ttH", "ttW", "ttX", "ZZ", "rare", "observed", "total", "total background"};

    for(int cat = 0; cat < 3; cat++){

      ofstream tableBkg;
      tableBkg.open("tableFullAnalysis_" + std::to_string(cat) + ".tex");
      tableBkg<<"\\begin{table}\n";
      tableBkg<<"\\begin{adjustbox}{width=1\\textwidth}\n";
      tableBkg<<"\\begin{tabular}{||c||c|c|c|c|c|c||}\\hline\n"; 
      tableBkg << std::fixed << setprecision(2) << "\n";

      tableBkg << "& \\multicolumn{3}{c|}{" + NameCategory[cat] << "} \\\\ \\hline \n";

      tableBkg << "Process";

      for(int num = 0; num < 3; num++){
        tableBkg << " & " << NameProcess[num] ;
      }

      tableBkg << "\\\\ \\hline\n";

      for(int bkg = 0; bkg < 8; bkg++){

        tableBkg << NameBKGextended[bkg] << " "; 
        for(int num = 3*cat; num < 3*(cat + 1); num++){
          tableBkg << " & $" << yields[bkg][num] << "\\pm"  << TMath::Sqrt(yieldsErrs[bkg][num]) << "\\pm"  << TMath::Sqrt(yieldsErrsSystematics[bkg][num]) << "$  ";
        }

         tableBkg << " \\\\ \n";
      }

      tableBkg << " \\hline \n";

      tableBkg << NameBKGextended[9] << " "; 
        for(int num = 3*cat; num < 3*(cat + 1); num++){
          tableBkg << " & $" << yieldsTotal[num] << "\\pm"  << TMath::Sqrt(yieldsTotalErrs[num]) << "\\pm"  << TMath::Sqrt(yieldsTotalErrsSystematics[num]) << "$  ";
        }

      tableBkg << " \\\\ \\hline \n";

      tableBkg << NameBKGextended[8] << " "; 
        for(int num = 3*cat; num < 3*(cat + 1); num++){
          tableBkg << " & $" << int(yields[8][num]) << "$  ";
        }

      tableBkg << " \\\\ \\hline\n";

      tableBkg <<"\\end{tabular}\n";
      tableBkg <<"\\end{adjustbox}\n";
      tableBkg <<"\\end{table}\n";
      tableBkg.close();
    }
*/



    ofstream tableBkgNew;
    tableBkgNew.open("tableFullAnalysisNew.tex");
    tableBkgNew<<"\\begin{table}\n";
    tableBkgNew<<"\\begin{adjustbox}{width=1\\textwidth}\n";
    tableBkgNew<<"\\begin{tabular}{ccccccc}\\hline\n"; 
    tableBkgNew<< std::fixed << setprecision(1) << "\n"; 
    tableBkgNew<< "\\Nbjets" << " & " << "\\Njets" << " & " << " Background" << " & " << "\\ttW" << " & " << "\\ttZ" << " & " << "Total" << " & " << "Observed \\\\ \\hline " << endl;

    int num = 0;

    TString nbjetsString[3] = {"\\multirow{3}{*}{= 0} ", "\\multirow{3}{*}{= 1} ", "\\multirow{3}{*}{ $>$ 1} "};
    TString njetsString[3] = {"= 2 ", "= 3 ", "$>$ 3"};


    for(int num = 0; num < SRNumber; num++){
          tableBkgNew << (num % 3 == 0 ? nbjetsString[int (num / 3)] : " ") << "&" << njetsString[num % 3] << " & " << yields[9][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[9][num] + yieldsErrsSystematics[9][num]) << " & " << yields[3][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[3][num] + yieldsErrsSystematics[3][num]) << "&" 
               << yields[1][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[1][num] + yieldsErrsSystematics[1][num]) << "&" << yieldsTotal[num] << "$\\pm$" << TMath::Sqrt(yieldsTotalErrs[num] + yieldsTotalErrsSystematics[num]) << "&" << int(yields[8][num]) <<  "\\\\ \n";
    }

    /*
    tableBkgNew << "\\multirow{3}{*}{= 0} & = 2 & " << yields[9][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[9][num] + yieldsErrsSystematics[3][num]) << " & " << yields[3][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[3][num] + yieldsErrsSystematics[3][num]) << "&" 
               << yields[1][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[1][num] + yieldsErrsSystematics[1][num]) << "&" << yieldsTotal[num] << "$\\pm$" << TMath::Sqrt(yieldsTotalErrs[num] + yieldsTotalErrsSystematics[num]) << "&" << int(yields[8][num]) <<  "\\\\ \n";
   
    num++;

    tableBkgNew << " & = 3 & " << yields[9][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[9][num] + yieldsErrsSystematics[9][num]) << " & " << yields[3][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[3][num] + yieldsErrsSystematics[3][num]) << "&" 
               << yields[1][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[1][num] + yieldsErrsSystematics[1][num]) << "&"<< yieldsTotal[num] << "$\\pm$" << TMath::Sqrt(yieldsTotalErrs[num] + yieldsTotalErrsSystematics[num]) << "&"  << int(yields[8][num]) <<  "\\\\ \n";
   

    num++;

    tableBkgNew << " & $>$ 3 & " << yields[9][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[9][num] + yieldsErrsSystematics[9][num]) << " & " << yields[3][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[3][num] + yieldsErrsSystematics[3][num]) << "&" 
               << yields[1][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[1][num] + yieldsErrsSystematics[1][num]) << "&" << yieldsTotal[num] << "$\\pm$" << TMath::Sqrt(yieldsTotalErrs[num] + yieldsTotalErrsSystematics[num]) << "&" << int(yields[8][num]) <<  "\\\\ \\hline \n";
   
    num++;

    tableBkgNew << "\\multirow{3}{*}{= 1} & = 2 & " << yields[9][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[9][num] + yieldsErrsSystematics[9][num]) << " & " << yields[3][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[3][num] + yieldsErrsSystematics[3][num]) << "&" 
               << yields[1][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[1][num] + yieldsErrsSystematics[1][num]) << "&" << yieldsTotal[num] << "$\\pm$" << TMath::Sqrt(yieldsTotalErrs[num] + yieldsTotalErrsSystematics[num]) << "&" << int(yields[8][num]) <<  "\\\\  \n";
   
    num++;

    tableBkgNew << " & = 3 & " << yields[9][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[9][num] + yieldsErrsSystematics[9][num]) << " & " << yields[3][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[3][num] + yieldsErrsSystematics[3][num]) << "&" 
               << yields[1][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[1][num] + yieldsErrsSystematics[1][num]) << "&" << yieldsTotal[num] << "$\\pm$" << TMath::Sqrt(yieldsTotalErrs[num] + yieldsTotalErrsSystematics[num]) << "&" << int(yields[8][num]) <<  "\\\\ \n";
   

    num++;

    tableBkgNew << " & $>$ 3 & " << yields[9][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[9][num] + yieldsErrsSystematics[9][num]) << " & " << yields[3][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[3][num] + yieldsErrsSystematics[3][num]) << "&" 
               << yields[1][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[1][num] + yieldsErrsSystematics[1][num]) << "&" << yieldsTotal[num] << "$\\pm$" << TMath::Sqrt(yieldsTotalErrs[num] + yieldsTotalErrsSystematics[num]) << "&" << int(yields[8][num]) <<  "\\\\ \\hline \n";
   
   
    num++;

    tableBkgNew << "\\multirow{3}{*}{$>$ 1} &=  2 & " << yields[9][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[9][num] + yieldsErrsSystematics[9][num]) << " & " << yields[3][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[3][num] + yieldsErrsSystematics[3][num]) << "&" 
               << yields[1][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[1][num] + yieldsErrsSystematics[1][num]) << "&" << yieldsTotal[num] << "$\\pm$" << TMath::Sqrt(yieldsTotalErrs[num] + yieldsTotalErrsSystematics[num]) << "&" << int(yields[8][num]) <<  "\\\\ \n";
   
    num++;

    tableBkgNew << " & 3 & " << yields[9][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[9][num] + yieldsErrsSystematics[9][num]) << " & " << yields[3][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[3][num] + yieldsErrsSystematics[3][num]) << "&" 
               << yields[1][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[1][num] + yieldsErrsSystematics[1][num]) << "&" << yieldsTotal[num] << "$\\pm$" << TMath::Sqrt(yieldsTotalErrs[num] + yieldsTotalErrsSystematics[num]) << "&" << int(yields[8][num]) <<  "\\\\ \n";
   

    num++;

    tableBkgNew << " & $>$ = 3 & " << yields[9][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[9][num] + yieldsErrsSystematics[9][num]) << " & " << yields[3][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[3][num] + yieldsErrsSystematics[3][num]) << "&" 
               << yields[1][num] << "$\\pm$"  << TMath::Sqrt(yieldsErrs[1][num] + yieldsErrsSystematics[1][num]) << "&" << yieldsTotal[num] << "$\\pm$" << TMath::Sqrt(yieldsTotalErrs[num] + yieldsTotalErrsSystematics[num]) << "&" << int(yields[8][num]) <<  "\\\\ \\hline \n";
   
    
    num++;
*/
    tableBkgNew << "\\multicolumn{2}{c}{Total} & " << yieldsTotalBackground << "$\\pm$"  << TMath::Sqrt(yieldsTotalBackgroundErrs) << " & " << yieldsTotalTTW << "$\\pm$"  << TMath::Sqrt(yieldsTotalTTWErrs) << "&" 
               << yieldsTotalTTZ << "$\\pm$"  << TMath::Sqrt(yieldsTotalTTZErrs) << "&" << yieldsTotalAll << "$\\pm$" << TMath::Sqrt(yieldsTotalAllErrs) << "&" << int(yieldsTotalObserved) <<  "\\\\ \\hline \n";
   
    

    tableBkgNew <<"\\end{tabular}\n";
    tableBkgNew <<"\\end{adjustbox}\n";
    tableBkgNew <<"\\end{table}\n";
    tableBkgNew.close();
    return;
    

    
    double tth_stat[SRNumberZpt];
    double ttz_stat[SRNumberZpt];
    double ttw_stat[SRNumberZpt];
    double ttxAs_stat[SRNumberZpt];
    double wz_stat[SRNumberZpt];
    double zz_stat[SRNumberZpt];
    double ttx_stat[SRNumberZpt];
    double rare_stat[SRNumberZpt];
    double fake_stat[SRNumberZpt];

    // 0 - fakes, 1-2 ttZ, 3-WZ, 4-16 ttX (ttH-4, ttgamma - 5, tttt-6, ttW-7), 17-28 (11-ZZ) rares, 29 - data
    TH1F* h_SR_yield[9]; 
    for (int i =0; i!=9; ++i) {
        TString name = Form("h_SR_yield_%d",i);
        h_SR_yield[i] = new TH1F(name,name,400,0,400);
        h_SR_yield[i]->Sumw2();
        for (int k=0; k!=SRNumberZpt; ++k) {
            h_SR_yield[i]->SetBinContent(k+1,yields[i][k]);
            h_SR_yield[i]->SetBinError(k+1,TMath::Sqrt(yieldsErrs[i][k]));
            double uncy;

            if(h_SR_yield[i]->GetBinContent(k+1) != 0)
              uncy = TMath::Abs(h_SR_yield[i]->GetBinError(k+1) / h_SR_yield[i]->GetBinContent(k+1));
            else
              uncy = 0.;

            if(i == 0)
              fake_stat[k] = 1 + uncy;
            if(i == 1)
              ttz_stat[k] = 1 + uncy;
            if(i == 2)
              wz_stat[k] = 1 + uncy;
            if(i == 3)
              tth_stat[k] = 1 + uncy;
            if(i == 4)
              ttw_stat[k] = 1 + uncy;
            if(i == 5)
              ttx_stat[k] = 1 + uncy;
            //if(i == 4)
            //  ttxAs_stat[k] = 1+uncy;
            if(i == 6)
              zz_stat[k] = 1 + uncy;
            if(i == 7)
              rare_stat[k] = 1 + uncy;
        }
    }
    
    
    TString datastr[SRNumberZpt];
    TString sttth[SRNumberZpt];
    TString stttw[SRNumberZpt];
    TString stfake[SRNumberZpt];
    TString stcharge[SRNumberZpt];
    TString stttz[SRNumberZpt];
    TString stttx[SRNumberZpt];
    TString stttxAs[SRNumberZpt];
    TString stwz[SRNumberZpt];
    TString stzz[SRNumberZpt];
    TString strare[SRNumberZpt];

    /*
    double ttw_jec[SRNumber]; // = {1.05, 1.02, 1.05, 1.05, 1.04, 1.02, 1.01, 1.02, 1.05, 1.04, 1.05, 1.02, 1.05, 1.05, 1.04, 1.02, 1.01, 1.02, 1.05, 1.04};

    double ttw_bl[SRNumber]; // = {1.01, 1.03, 1.01, 1.05, 1.01, 1.01, 1.03, 1.01, 1.05, 1.01, 1.01, 1.03, 1.01, 1.05, 1.01, 1.01, 1.03, 1.01, 1.05, 1.01};
    double ttw_bb[SRNumber]; // = {1.02, 1.03, 1.05, 1.07, 1.03, 1.04, 1.03, 1.06, 1.07, 1.05, 1.02, 1.03, 1.05, 1.07, 1.03, 1.04, 1.03, 1.06, 1.07, 1.05};

    double ttz_jec[SRNumber]; // = {1.07, 1.03, 1.08, 1.07, 1.02, 1.06, 1.01, 1.01, 1.05, 1.04, 1.07, 1.03, 1.08, 1.07, 1.02, 1.06, 1.01, 1.01, 1.05, 1.04};

    double ttz_bl[SRNumber]; // = {1.01, 1.03, 1.02, 1.05, 1.02, 1.01, 1.03, 1.01, 1.05, 1.01, 1.01, 1.03, 1.02, 1.05, 1.02, 1.01, 1.03, 1.01, 1.05, 1.01};
    double ttz_bb[SRNumber]; // = {1.02, 1.02, 1.05, 1.05, 1.03, 1.04, 1.02, 1.06, 1.07, 1.05, 1.02, 1.02, 1.05, 1.05, 1.03, 1.04, 1.02, 1.06, 1.07, 1.05};

    for(int i = 0; i < SRNumber; i++){
      ttw_jec[i] = 1.03;
      ttw_bl[i] = 1.03;
      ttw_bb[i] = 1.03;

      ttz_jec[i] = 1.03;
      ttz_bl[i] = 1.03;
      ttz_bb[i] = 1.03;
    }
    */

    
    /*
    process         ttZ     WZ      ttX     ttW     ttH    Fake rare    ZZ
    JES         lnN 1.0003  0.9916  1.0014  1.0114  0.9889  -   -   0.9853
    JES         lnN 0.9998  0.9888  1.0056  0.9863  1.0089  -   -   0.9784
    JES         lnN 1.0104  0.9760  1.0055  1.0339  1.0126  -   -   0.9693
    JES         lnN 0.9835  1.0045  1.0032  1.0014  0.9954  -   -   1.0004
    JES         lnN 1.0081  0.9852  1.0001  1.0040  0.9948  -   -   0.9995
    JES         lnN 1.0038  1.0125  1.0098  1.0069  1.0096  -   -   0.9922
    JES         lnN 0.9697  1.0175  1.0033  1.0030  0.9763  -   -   1.0031
    JES         lnN 0.9925  1.0025  0.9996  0.9989  1.0068  -   -   1.0150
    JES         lnN 1.0083  0.9888  1.0040  0.9933  1.0034  -   -   0.9617
    */
    /*
    double ttz_jec[SRNumber] = {1.0003, 0.9998, 1.0104, 0.9835, 1.0081, 1.0038, 0.9697, 0.9925, 1.0083};
    double ttw_jec[SRNumber] = {1.0014, 0.9863, 1.0339, 1.0014, 1.0040, 1.0069, 1.0030, 0.9989, 0.9933};
    double tth_jec[SRNumber] = {0.9889, 1.0089, 1.0126, 0.9954, 0.9948, 1.0096, 0.9763, 1.0068, 1.0034};
    double ttx_jec[SRNumber] = {1.0014, 1.0056, 1.0055, 1.0032, 1.0001, 1.0098, 1.0033, 0.9996, 1.0040};
    double wz_jec[SRNumber]  = {0.9916, 0.9888, 0.9760, 1.0045, 0.9852, 1.0125, 1.0175, 1.0025, 0.9888};

    double ttz_bl[SRNumber]  = {0.9987, 0.9978, 0.9962, 1.0003, 0.9996, 0.9979, 1.0044, 1.0050, 1.0038};
    double ttz_bb[SRNumber]  = {0.9804, 0.9726, 0.9635, 1.0102, 1.0040, 0.9979, 1.0357, 1.0348, 1.0344};

    double ttw_bl[SRNumber] = {0.9992, 0.9981, 0.9959, 1.0006, 0.9987, 1.0012, 0.9988, 1.0031, 1.0054};
    double ttw_bb[SRNumber] = {0.9669, 0.9696, 0.9641, 1.0031, 0.9975, 0.9970, 1.0343, 1.0306, 1.0290};

    double tth_bl[SRNumber] = {0.9988, 0.9979, 0.9962, 1.0011, 0.9991, 0.9982, 1.0000, 1.0063, 1.0070};
    double tth_bb[SRNumber] = {0.9766, 0.9716, 0.9642, 1.0062, 1.0028, 0.9973, 1.0353, 1.0341, 1.0330};

    double ttx_bl[SRNumber] = {0.9984, 0.9974, 0.9961, 0.9998, 0.9992, 0.9982, 1.0060, 1.0064, 1.0074};
    double ttx_bb[SRNumber] = {0.9804, 0.9760, 0.9700, 1.0124, 1.0082, 1.0027, 1.0359, 1.0353, 1.0336};

    double wz_bl[SRNumber] =  {0.9975, 0.9961, 0.9944, 1.0302, 1.0269, 1.0335, 1.0336, 1.0378, 1.0284};
    double wz_bb[SRNumber] =  {0.9985, 0.9979, 0.9971, 1.0250, 1.0232, 1.0183, 1.0462, 1.0419, 1.0358};
    */

    /*
    process         ttZ     WZ      ttX     ttW     ttH    Fake rare    ZZ
    btagl       lnN 0.9987  0.9975  0.9984  0.9992  0.9988  -   -   0.9975
    btagl       lnN 0.9978  0.9961  0.9974  0.9981  0.9979  -   -   0.9962
    btagl       lnN 0.9962  0.9944  0.9961  0.9959  0.9962  -   -   0.9946
    btagl       lnN 1.0003  1.0302  0.9998  1.0006  1.0011  -   -   1.0392
    btagl       lnN 0.9996  1.0269  0.9992  0.9987  0.9991  -   -   1.0353
    btagl       lnN 0.9979  1.0335  0.9982  1.0012  0.9982  -   -   1.0263
    btagl       lnN 1.0044  1.0336  1.0060  0.9988  1.0000  -   -   1.0276
    btagl       lnN 1.0050  1.0378  1.0064  1.0031  1.0063  -   -   1.0392
    btagl       lnN 1.0038  1.0284  1.0074  1.0054  1.0070  -   -   1.0142


    btagb       lnN 0.9804  0.9985  0.9804  0.9669  0.9766  -   -   0.9989
    btagb       lnN 0.9726  0.9979  0.9760  0.9696  0.9716  -   -   0.9982
    btagb       lnN 0.9635  0.9971  0.9700  0.9641  0.9642  -   -   0.9977
    btagb       lnN 1.0102  1.0250  1.0124  1.0031  1.0062  -   -   1.0152
    btagb       lnN 1.0040  1.0232  1.0082  0.9975  1.0028  -   -   1.0145
    btagb       lnN 0.9979  1.0183  1.0027  0.9970  0.9973  -   -   1.0140
    btagb       lnN 1.0357  1.0462  1.0359  1.0343  1.0353  -   -   1.0403
    btagb       lnN 1.0348  1.0419  1.0353  1.0306  1.0341  -   -   1.0300
    btagb       lnN 1.0344  1.0358  1.0336  1.0290  1.0330  -   -   1.0359
    */
    
    double wz_jec[SRNumber] = {1.03, 1.05, 1.06, 1.08, 1.02, 1.05, 1.07, 1.07};

    double wz_bl[SRNumber] = {1.02, 1.02, 1.03, 1.05, 1.04, 1.05, 1.05, 1.07};
    double wz_bb[SRNumber] = {1.00, 1.00, 1.00, 1.00, 1.01, 1.01, 1.01, 1.01};


    double ttx_jec[SRNumber] = {1.00, 1.03, 1.05, 1.06, 1.00, 1.02, 1.04, 1.07};

    double ttx_bl[SRNumber] = {1.01, 1.02, 1.03, 1.04, 1.01, 1.02, 1.02, 1.04};
    double ttx_bb[SRNumber] = {1.00, 1.00, 1.00, 1.00, 1.01, 1.02, 1.02, 1.02};


    double tth_jec[SRNumber] = {1.00, 0.97, 1.08, 1.01, 0.97, 1.00, 1.00, 1.06};

    double tth_bl[SRNumber] = {1.01, 1.01, 1.02, 1.04, 1.00, 1.01, 1.02, 1.04};
    double tth_bb[SRNumber] = {1.00, 1.00, 1.00, 1.00, 1.02, 1.02, 1.02, 1.02};

    
    double ttw_jec[SRNumber] = {1.00, 1.00, 1.02, 1.06, 0.99, 1.01, 1.04, 1.04};

    double ttw_bl[SRNumber] = {1.00, 1.01, 1.02, 1.04, 1.00, 1.01, 1.02, 1.03};
    double ttw_bb[SRNumber] = {1.00, 1.00, 1.00, 1.00, 1.02, 1.02, 1.02, 1.02};


    double ttz_jec[SRNumber] = {0.96, 0.98, 1.01, 1.04, 0.96, 0.98, 1.01, 1.04};

    double ttz_bl[SRNumber] = {1.01, 1.01, 1.02, 1.03, 1.01, 1.01, 1.02, 1.03};
    double ttz_bb[SRNumber] = {1.00, 1.00, 1.00, 1.00, 1.01, 1.02, 1.02, 1.02};

    
    for(int id = 0; id< SRNumber; id++){

        
        datastr[id]  = "B" + std::to_string(id); 

        stttw[id]  = "stttw" + std::to_string(id);
        sttth[id]  = "sttth" + std::to_string(id);
        stttz[id]  = "stttz" + std::to_string(id);
        stwz[id]  = "stwz" + std::to_string(id);
        stttx[id]  = "stttx" + std::to_string(id);
        stfake[id]  = "stfake" + std::to_string(id);
        strare[id]  = "strare" + std::to_string(id);
        stzz[id]  = "stzz" + std::to_string(id);
        stcharge[id] = "stcharge" + std::to_string(id);

        
        gSystem->Exec("rm datacards/" + datastr[id] + ".txt"); // delete previous tex file
        ofstream fileout;
        fileout .open ( "datacards/" + datastr[id] + ".txt", ios_base::app); // create a new tex file
        fileout << fixed << showpoint << setprecision(4);
        fileout << "#  " << datastr[id]  << endl;
        fileout << "imax 1  number of channels " <<  endl;
        fileout << "jmax 7  number of backgrounds " <<  endl;
        //fileout << "kmax 21  number of nuisance parameters (sources of systematical uncertainties) " <<  endl;
        //fileout << "kmax 19  number of nuisance parameters (sources of systematical uncertainties) " <<  endl;
        fileout << "kmax 23  number of nuisance parameters (sources of systematical uncertainties) " <<  endl;
        fileout << "----------- " <<  endl;
        fileout << "shapes * * FAKE" << endl;
        fileout << "----------- " <<  endl;
        fileout << "bin  " << '\t' << datastr[id] <<  endl;
        fileout << "observation  " << h_SR_yield[8]->GetBinContent(id+1) << endl;
        fileout << "----------- " <<  endl;
        fileout << "bin  " <<  '\t' << '\t' <<'\t' <<  datastr[id] << '\t' << datastr[id] <<  '\t' << datastr[id] <<  '\t' << datastr[id] <<  '\t' << datastr[id]  << '\t' << datastr[id] << '\t' << datastr[id] << '\t' << datastr[id]  << endl;
        fileout << "process  " <<  '\t' <<'\t' <<  "ttZ" << '\t' << "Fake" <<  '\t' << "WZ" <<  '\t' << "ttH" <<  '\t' << "ttW" << '\t' << "ttX" << '\t' << "ZZ" << '\t' << "rare" << endl;
        // 0 - fakes, 1 - charge mis-ID, 2-3 ttZ, 4-6-ttX, 7-WZ, 8-ZZ, 9-20 -rares, 21 - ttW
        fileout << "process  " <<  '\t' <<'\t' <<  "-1" << '\t' << "1" <<  '\t' << "2" <<  '\t' << "3" <<  '\t' << "4"  << '\t' << "5" << '\t' << "6" << '\t' << "7" << endl;
        //fileout << "process  " <<  '\t' <<'\t' <<  "0" << '\t' << "1" <<  '\t' << "2" <<  '\t' << "-1" <<  '\t' << "4"  << '\t' << "5" << '\t' << "6" << '\t' << "7" << endl;
        fileout << "rate  " <<  '\t' << '\t' << '\t' << h_SR_yield[1]->GetBinContent(id+1)  << '\t' << h_SR_yield[0]->GetBinContent(id+1) <<  '\t' << h_SR_yield[2]->GetBinContent(id+1) <<  '\t' << h_SR_yield[3]->GetBinContent(id+1) <<  '\t' << h_SR_yield[4]->GetBinContent(id+1) <<  '\t' << (h_SR_yield[5]->GetBinContent(id+1) > 0 ? h_SR_yield[5]->GetBinContent(id+1) : 0) <<  '\t' << (h_SR_yield[6]->GetBinContent(id+1) > 0 ? h_SR_yield[6]->GetBinContent(id+1) : 0) <<  '\t' << (h_SR_yield[7]->GetBinContent(id+1) > 0 ? h_SR_yield[7]->GetBinContent(id+1) : 0) << endl;
        fileout << "----------- " <<  endl;
        
        fileout << stttz[id] << '\t' << '\t'  << "    lnN"  << '\t' << ttz_stat[id] << '\t' << "-" <<  '\t' <<"-" <<  '\t' << "-" <<   '\t' << "-"  <<  '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" << endl;
        fileout << stfake[id] << '\t' << '\t'  << "    lnN"  << '\t'<< "-" << '\t' << fake_stat[id] <<  '\t' <<"-" <<  '\t' << "-" <<   '\t' << "-"  <<  '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" << endl;
        fileout << stwz[id] << '\t' << '\t' << "    lnN"   << '\t'<< "-" << '\t' << "-" <<  '\t' << wz_stat[id] <<  '\t' << "-" <<   '\t' << "-"  <<  '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" << endl;
        fileout << sttth[id] << '\t' << '\t' << "    lnN"  << '\t'<< "-" << '\t' << "-" <<  '\t' << "-" <<  '\t' << tth_stat[id] <<   '\t' << "-"  <<  '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" << endl;
        fileout << stttw[id] << '\t' << '\t' << "    lnN"  << '\t'<< "-" << '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" <<  '\t' << ttw_stat[id]  <<  '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" << endl;
        fileout << stttx[id] << '\t' << '\t' << "    lnN"  << '\t'<< "-" << '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" <<  '\t' <<  "-" <<  '\t' << ttx_stat[id] <<  '\t' << "-" <<  '\t' << "-" << endl;
        fileout << stzz[id] << '\t' << '\t' << "    lnN"  << '\t'<< "-" << '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" <<  '\t' <<  "-" <<  '\t' << "-" <<  '\t' << zz_stat[id] <<  '\t' << "-" << endl;
        fileout << strare[id] << '\t' << '\t' << "    lnN" << '\t' << "-" << '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" <<  '\t' <<  "-" <<  '\t' << "-" <<  '\t' << "-" <<  '\t' << rare_stat[id] << endl;
        fileout << endl;
        
        // prefit 
        /*    
        fileout << "lumi     " << '\t'  << '\t' << "lnN" << '\t' << 1.025 << '\t' <<  "-" <<  '\t' << 1.025 << '\t' << 1.025 << '\t' << 1.025  <<  '\t' << 1.026 << '\t' << 1.025 << '\t' << 1.025 << endl;

        fileout << "PU       " << '\t'  << '\t' << "lnN" << '\t' << 1.01        << '\t' <<  "-" << '\t' << 1.01        << '\t' << 1.01          << '\t' << 1.01       <<  '\t' << 1.01        << '\t' << 1.01 <<  '\t' << 1.01 << endl;
        fileout << "trigger  " << '\t'  << '\t' << "lnN" << '\t' << 1.04        << '\t' <<  "-" << '\t' << "-"         << '\t' << 1.04          << '\t' << 1.04       <<  '\t' << 1.04        << '\t' << 1.04 <<  '\t' << 1.04 << endl;
        fileout << "LeptonId " << '\t'  << '\t' << "lnN "<< '\t' << 1.05        << '\t' <<  "-" << '\t' << "-"         << '\t' << 1.05          << '\t' << 1.05       <<  '\t' << 1.05        << '\t' << 1.05 <<  '\t' << 1.05 << endl;
        fileout << "JES      " << '\t'  << '\t' << "lnN" << '\t' << ttz_jec[id] << '\t' <<  "-" << '\t' << wz_jec[id] << '\t' << tth_jec[id]   << '\t' << ttw_jec[id]<<  '\t' << ttx_jec[id] << '\t' << "-" <<  '\t' << "-" << endl;
        fileout << "JER      " << '\t'  << '\t' << "lnN" << '\t' << 1.01        << '\t' <<  "-" << '\t' << 1.01        << '\t' << 1.01          << '\t' << 1.01       <<  '\t' << "-"         << '\t' << "-" <<  '\t' << "-" << endl;
        fileout << "btagl    " << '\t'  << '\t' << "lnN" << '\t' << ttz_bl[id]  << '\t' <<  "-" << '\t' << wz_bl[id]  << "\t" << tth_bl[id]    << '\t' << ttw_bl[id] <<  '\t' << ttx_bl[id]  << '\t' << "-" <<  '\t' << "-"  << endl;
        fileout << "btagb    " << '\t'  << '\t' << "lnN" << '\t' << ttz_bb[id]  << '\t' <<  "-" << '\t' << wz_bb[id]  << "\t" << tth_bb[id]    << '\t' << ttw_bb[id] <<  '\t' << ttx_bb[id]  << '\t' << "-" <<  '\t' << "-"  << endl;
        fileout << endl;

        fileout << endl;
        fileout << "PDF      " << '\t'  << '\t' << "lnN" << '\t' << 1.01 << '\t' <<  "-" <<  '\t' << 1.01 <<  '\t' << 1.01  <<  '\t' << 1.01  <<  '\t' << 1.01 <<  '\t' << "-" <<  '\t' << "-" << endl;
        fileout << "Q2       " << '\t'  << '\t' << "lnN" << '\t' << 1.01 << '\t' <<  "-" <<  '\t' << 1.01 <<  '\t' << 1.01  <<  '\t' << 1.01  <<  '\t' << 1.01 <<  '\t' << "-" <<  '\t' << "-" << endl;
        fileout << endl;
        
        fileout << endl;
        fileout << "fake     " << '\t'  << '\t' << "lnN" << '\t' << "-" << '\t' << 1.30 << '\t' << "-"   << '\t' << "-"   <<  '\t' << "-"   << '\t' << "-"  << '\t' << "-"  <<  '\t' << "-"  << endl;
        fileout << "WZ       " << '\t'  << '\t' << "lnN" << '\t' << "-" << '\t' << "-"  << '\t' << WZbkg[id]   << '\t' << "-"   <<  '\t' << "-"   << '\t' << "-"  << '\t' << "-"  <<  '\t' << "-"  << endl;
        fileout << "ttX      " << '\t'  << '\t' << "lnN" << '\t' << "-" << '\t' << "-"  << '\t' << "-"   << '\t' << 1.11  <<  '\t' << 1.11  << '\t' << 1.11 << '\t' << "-"  <<  '\t' << "-"  << endl;
        fileout << "ZZ       " << '\t'  << '\t' << "lnN" << '\t' << "-" << '\t' << "-"  << '\t' << "-"   << '\t' << "-"   <<  '\t' << "-"   << '\t' << "-"  << '\t' << 1.20 <<  '\t' << "-"  <<  endl;
        fileout << "rare     " << '\t'  << '\t' << "lnN" << '\t' << "-" << '\t' << "-"  << '\t' << "-"   << '\t' << "-"   <<  '\t' << "-"   << '\t' << "-"  << '\t' << "-"  <<  '\t' << 1.50 << endl;
        fileout << endl;
        */
        
    }
    
    /*
    for(int id = 0; id< SRNumber; id++){

        
        datastr[id]  = "B" + std::to_string(id); 

        stttw[id]  = "stttw" + std::to_string(id);
        sttth[id]  = "sttth" + std::to_string(id);
        stttz[id]  = "stttz" + std::to_string(id);
        stwz[id]  = "stwz" + std::to_string(id);
        stttx[id]  = "stttx" + std::to_string(id);
        stttxAs[id]  = "stttxAs" + std::to_string(id);
        stfake[id]  = "stfake" + std::to_string(id);
        strare[id]  = "strare" + std::to_string(id);
        stzz[id]  = "stzz" + std::to_string(id);
        stcharge[id] = "stcharge" + std::to_string(id);

        
        gSystem->Exec("rm datacards/" + datastr[id] + ".txt"); // delete previous tex file
        ofstream fileout;
        fileout .open ( "datacards/" + datastr[id] + ".txt", ios_base::app); // create a new tex file
        fileout << fixed << showpoint << setprecision(2);
        fileout << "#  " << datastr[id]  << endl;
        fileout << "imax 1  number of channels " <<  endl;
        fileout << "jmax 6  number of backgrounds " <<  endl;
        //fileout << "kmax 21  number of nuisance parameters (sources of systematical uncertainties) " <<  endl;
        //fileout << "kmax 19  number of nuisance parameters (sources of systematical uncertainties) " <<  endl;
        fileout << "kmax 23  number of nuisance parameters (sources of systematical uncertainties) " <<  endl;
        fileout << "----------- " <<  endl;
        fileout << "shapes * * FAKE" << endl;
        fileout << "----------- " <<  endl;
        fileout << "bin  " << '\t' << datastr[id] <<  endl;
        fileout << "observation  " << h_SR_yield[7]->GetBinContent(id+1) << endl;
        fileout << "----------- " <<  endl;
        fileout << "bin  " <<  '\t' << '\t' <<'\t' <<  datastr[id] << '\t' << datastr[id] <<  '\t' << datastr[id] <<  '\t' << datastr[id]  << '\t' << datastr[id] << '\t' << datastr[id] << '\t' << datastr[id]  << endl;
        fileout << "process  " <<  '\t' <<'\t' <<  "ttZ" << '\t' << "Fake" <<  '\t' << "WZ" <<  '\t' << "ttX" << '\t' << "ttXas" << '\t' << "ZZ" << '\t' << "rare" << endl;
        // 0 - fakes, 1 - charge mis-ID, 2-3 ttZ, 4-6-ttX, 7-WZ, 8-ZZ, 9-20 -rares, 21 - ttW
        fileout << "process  " <<  '\t' <<'\t' <<  "-1" << '\t' << "1" <<  '\t' << "2" <<  '\t' << "3" <<  '\t' << "4"  << '\t' << "5" << '\t' << "6" << endl;
        //fileout << "process  " <<  '\t' <<'\t' <<  "0" << '\t' << "1" <<  '\t' << "2" <<  '\t' << "-1" <<  '\t' << "4"  << '\t' << "5" << '\t' << "6" << '\t' << "7" << endl;
        fileout << "rate  " <<  '\t' << '\t' << '\t' << h_SR_yield[1]->GetBinContent(id+1)  << '\t' << h_SR_yield[0]->GetBinContent(id+1) <<  '\t' << h_SR_yield[2]->GetBinContent(id+1) <<  '\t' << h_SR_yield[3]->GetBinContent(id+1) <<  '\t' << h_SR_yield[4]->GetBinContent(id+1) <<  '\t' << (h_SR_yield[5]->GetBinContent(id+1) > 0 ? h_SR_yield[5]->GetBinContent(id+1) : 0) <<  '\t' << (h_SR_yield[6]->GetBinContent(id+1) > 0 ? h_SR_yield[6]->GetBinContent(id+1) : 0) << endl;
        fileout << "----------- " <<  endl;
        
        fileout << stttz[id] << '\t' << '\t'  << "    lnN"  << '\t' << ttz_stat[id] << '\t' << "-" <<  '\t' <<"-" <<  '\t' << "-" <<   '\t' << "-"  <<  '\t' << "-" <<  '\t' << "-" << endl;
        fileout << stfake[id] << '\t' << '\t'  << "    lnN"  << '\t'<< "-" << '\t' << fake_stat[id] <<  '\t' <<"-" <<  '\t' << "-" <<   '\t' << "-"  <<  '\t' << "-" <<  '\t' << "-" << endl;
        fileout << stwz[id] << '\t' << '\t' << "    lnN"   << '\t'<< "-" << '\t' << "-" <<  '\t' << wz_stat[id] <<  '\t' << "-" <<   '\t' << "-"  <<  '\t' << "-" <<  '\t' << "-" << endl;
        fileout << stttx[id] << '\t' << '\t' << "    lnN"  << '\t'<< "-" << '\t' << "-" <<  '\t' << "-" <<  '\t' << ttx_stat[id] <<  '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" << endl;
        fileout << stttxAs[id] << '\t' << '\t' << "    lnN"  << '\t'<< "-" << '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" <<  '\t' << ttxAs_stat[id] <<  '\t' << "-" <<  '\t' << "-" << endl;
        fileout << stzz[id] << '\t' << '\t' << "    lnN"  << '\t'<< "-" << '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" <<  '\t' << zz_stat[id] <<  '\t' << "-" << endl;
        fileout << strare[id] << '\t' << '\t' << "    lnN" << '\t' << "-" << '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" <<  '\t' << "-" <<  '\t' << rare_stat[id] << endl;
        fileout << endl;
        
        // prefit 
    
        fileout << "lumi     " << '\t'  << '\t' << "lnN" << '\t' << 1.026 << '\t' <<  "-" <<  '\t' << 1.026 << '\t' << 1.026 << '\t' << 1.026 << '\t' << 1.026  <<  '\t' << 1.026 << endl;
        //fileout << "lumi.    " << '\t'  << '\t' << "lnN" << '\t' << 0.99 << '\t' <<  "-" <<  '\t' << "-" << '\t' << 0.99 << '\t' << 0.99  <<  '\t' << 0.99 << '\t' << 0.99 << '\t' << 0.99 << endl;
       
        fileout << "PU       " << '\t'  << '\t' << "lnN" << '\t' << 1.01        << '\t' <<  "-" << '\t' << 1.01        << '\t' << 1.01          << '\t' << 1.01          << '\t' << 1.01       <<  '\t' << 1.01        << endl;
        fileout << "trigger  " << '\t'  << '\t' << "lnN" << '\t' << 1.04        << '\t' <<  "-" << '\t' << "-"         << '\t' << 1.04          << '\t' << 1.04          << '\t' << 1.04       <<  '\t' << 1.04        << endl;
        fileout << "LeptonId " << '\t'  << '\t' << "lnN "<< '\t' << 1.05        << '\t' <<  "-" << '\t' << "-"         << '\t' << 1.05          << '\t' << 1.05          << '\t' << 1.05       <<  '\t' << 1.05        << endl;
        fileout << "JES      " << '\t'  << '\t' << "lnN" << '\t' << ttz_jec[id] << '\t' <<  "-" << '\t' << wz_jec[id]  << '\t' << ttx_jec[id] <<  '\t' << ttx_jec[id]    << '\t' << "-"        <<  '\t' << "-"         << endl;
        fileout << "JER      " << '\t'  << '\t' << "lnN" << '\t' << 1.01        << '\t' <<  "-" << '\t' << 1.01        << '\t' << 1.01          << '\t' << 1.01          <<  '\t' << "-"       <<  '\t' << "-"         << endl;
        fileout << "btagl    " << '\t'  << '\t' << "lnN" << '\t' << ttz_bl[id]  << '\t' <<  "-" << '\t' << wz_bl[id]   <<  '\t' << ttx_bl[id]  <<  '\t' << ttx_bl[id]    << '\t' << "-"        <<  '\t' << "-"         << endl;
        fileout << "btagb    " << '\t'  << '\t' << "lnN" << '\t' << ttz_bb[id]  << '\t' <<  "-" << '\t' << wz_bb[id]   <<  '\t' << ttx_bb[id]   <<  '\t' << ttx_bb[id]  << '\t' << "-"         <<  '\t' << "-"         << endl;
        fileout << endl;
        //fileout << "exp group = lumi. PU. trigger LeptonId JES JER btagl btagb"  << endl;
        //fileout << "bg group = lumi. PU. trigger LeptonId JES JER btagl btagb"  << endl;
        //fileout << "bgexp group = lumi. PU. trigger LeptonId JES JER btagl btagb"  << endl;
        fileout << endl;
        fileout << "PDF      " << '\t'  << '\t' << "lnN" << '\t' << 1.01 << '\t' <<  "-" <<  '\t' << 1.01 <<  '\t' << 1.01  <<  '\t' << 1.01  <<  '\t' << "-" <<  '\t' << "-" << endl;
        fileout << "Q2       " << '\t'  << '\t' << "lnN" << '\t' << 1.01 << '\t' <<  "-" <<  '\t' << 1.01 <<  '\t' << 1.01  <<  '\t' << 1.01  <<  '\t' << "-" <<  '\t' << "-" << endl;
        fileout << endl;
        //fileout << "sigtheo group = PDF Q2"  << endl;
        //fileout << "theo group = PDF Q2"  << endl;
        fileout << endl;
        fileout << "fake     " << '\t'  << '\t' << "lnN" << '\t' << "-" << '\t' << 1.30 << '\t' << "-"         << '\t' << "-"   <<  '\t' << "-"   << '\t' << "-"  << '\t'  << "-"  << endl;
        fileout << "WZ       " << '\t'  << '\t' << "lnN" << '\t' << "-" << '\t' << "-"  << '\t' << WZbkg[id]   << '\t' << "-"   <<  '\t' << "-"   << '\t' << "-"  << '\t'  << "-"  << endl;
        fileout << "ttX      " << '\t'  << '\t' << "lnN" << '\t' << "-" << '\t' << "-"  << '\t' << "-"         << '\t' << 1.11  <<  '\t' << "-"  << '\t' << "-" << '\t'  << "-"  << endl;
        fileout << "ttXas    " << '\t'  << '\t' << "lnN" << '\t' << "-" << '\t' << "-"  << '\t' << "-"         << '\t' << "-"   <<  '\t'  << 1.11 <<  '\t' << "-"  <<  '\t' << "-"  << endl;
        fileout << "ZZ       " << '\t'  << '\t' << "lnN" << '\t' << "-" << '\t' << "-"  << '\t' << "-"         << '\t' << "-"   <<  '\t' << "-"  <<  '\t' << 1.20 <<  '\t' << "-"  <<  endl;
        fileout << "rare     " << '\t'  << '\t' << "lnN" << '\t' << "-" << '\t' << "-"  << '\t' << "-"         << '\t'  << "-"  <<  '\t' << "-"  <<  '\t' << "-"  <<  '\t' << 1.50 << endl;
        fileout << endl;

        //fileout << "bgtheo group = fake charge ttX. WZ. rare."  << endl;
        //fileout << "theo group += fake charge ttX. WZ. rare."  << endl;
        //fileout << "bg group += fake charge ttX. WZ. rare."  << endl;
        
    }
    */
          std::cout << "datacard DONE" << std::endl;

    

}

int main(int argc, char *argv[]){

    TApplication *rootapp = new TApplication("example", &argc, argv);

    runttWvsTTbar();

    rootapp->Run();

    return 0;
}



