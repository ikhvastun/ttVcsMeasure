#include "THStack.h"

#include "../interface/treeReader.h"

using namespace std;

unsigned treeReader::dilFlavorComb(const std::vector<unsigned>& ind){
    unsigned flavCount[3] = {0,0,0};
    for(unsigned l = 0; l < 2; ++l) ++flavCount[_lFlavor[ind[l]]];
    if(flavCount[2] == 0){
        if(flavCount[1] == 0) return 0; //ee
        else if(flavCount[1] == 1) return 1; //em
        else return 2; //mm
    } else if(flavCount[2] == 1){
        if(flavCount[1] == 0) return 3; //et
        else return 4; //mt
    }
    return 5; //tt
}

double treeReader::coneCorr(const unsigned ind){
    return 1. + std::max(_relIso[ind] - 0.1, 0.);
}

void treeReader::setConePt(){
    for(unsigned l = 0; l < _nLight; ++l){
        double coneC = coneCorr(l);
        _lPt[l] *= coneC;
        _lE[l] *= coneC;
    }
}

/*
bool treeReader::lepIsGood(const unsigned l){
    if(!lepIsLoose(l)) return false;
    //if(_closestJetCsvV2[l] > 0.8484) return false;
    if(_lFlavor[l] == 0 && !_lElectronPassEmu[l]) return false;
    if(_lFlavor[l] == 0 && !eleIsClean(l)) return false;
    return true;
}
*/

bool treeReader::lepIsGood(const unsigned l){
    // what is used in leptonMVA analysis
    if(!lepIsFOGood(l)) return false;
    if(_lFlavor[l] == 1 && !_lPOGMedium[l]) return false;
    if(_leptonMvaTTH[l] < 0.9) return false;

    /*
    if(leptonSelection == 2){

        if(_lFlavor[l] == 1 && ((_lMuonTrackPtErr[l]/_lMuonTrackPt[l]) > 0.2)) return false;
    }
    */
    return true;
}

bool treeReader::lepIsFOGood(const unsigned l){
    // what is used in leptonMVA analysis
    if(!lepIsLoose(l)) return false;
    if(_closestJetCsvV2[l] > 0.8484) return false;

    if(_lFlavor[l] == 0 && !_lElectronPassEmu[l]) return false;
    if(_lFlavor[l] == 0 && !eleIsClean(l)) return false;

    /*
    if(leptonSelection == 2){
        if(_lFlavor[l] == 0 && _lElectronMissingHits[l] != 0) return false;

        if(_lFlavor[l] == 0 && !_lElectronPassConvVeto[l]) return false;
        if(_lFlavor[l] == 0 && !_lElectronChargeConst[l]) return false;

    }
    */
    
    if(_leptonMvaTTH[l] < 0.9){

        if(_lFlavor[l] == 0 && _lElectronMvaHZZ[l] <= 0.0 + (fabs(_lEta[l]) >= 1.479)*0.7) return false;
        
        if(_ptRatio[l] <= 0.5) return false;
        if(_closestJetCsvV2[l] >= 0.3) return false;
        if(_lFlavor[l] == 1 && _lMuonSegComp[l] <= 0.3) return false; 
    }

    return true;
}

bool treeReader::elePassVLooseMvaIDSUSY(const unsigned ind){
    if(_lFlavor[ind] != 0) return true;
    static const double gpCuts[3][2] = { {-0.48,-0.85}, {-0.67, -0.91}, {-0.49, -0.83} };
    static const double hzzCuts[3] = {0.46, -0.03, 0.06};
    unsigned eta = (fabs(_lEta[ind]) >= 0.8) + (fabs(_lEta[ind]) > 1.479);
    if(_lPt[ind] > 10){
        return _lElectronMva[ind] > std::min( gpCuts[eta][0], std::max(gpCuts[eta][1], gpCuts[eta][0] + (gpCuts[eta][1] - gpCuts[eta][0])*0.1*(_lPt[ind] - 15.) ) );
    } else{
        return _lElectronMvaHZZ[ind] > hzzCuts[eta];
    }
}

bool treeReader::eleIsClean(const unsigned ind){
    TLorentzVector ele;
    ele.SetPtEtaPhiE(_lPt[ind], _lEta[ind], _lPhi[ind], _lE[ind]);
    for(unsigned m = 0; m < _nMu; ++m){
        if(lepIsLoose(m)){
            TLorentzVector mu;
            mu.SetPtEtaPhiE(_lPt[m], _lEta[m], _lPhi[m], _lE[m]);
            if(ele.DeltaR(mu) < 0.05) return false;
        }
    }
    return true;
}

bool treeReader::lepIsLoose(const unsigned ind){
    // this is where potentially the discrepancy can come from
    //if(_closestJetCsvV2[ind] > 0.8484) return false;

    if(_lFlavor[ind] == 2) return false;  //don't consider taus here
    //if(_lPt[ind] <= 7 - 2*_lFlavor[ind]) return false;
    if(_lPt[ind] <= 10) return false;
    if(fabs(_lEta[ind]) >= (2.5 - 0.1*_lFlavor[ind])) return false;
    if(fabs(_dxy[ind]) >= 0.05) return false;
    if(fabs(_dz[ind]) >= 0.1) return false;
    if(_3dIPSig[ind] >= 8) return false;
    if(_miniIso[ind] >= 0.4) return false;
    if(_lFlavor[ind] == 1){
        //if(!_lPOGLoose[ind]) return false;
        if(!_lPOGMedium[ind]) return false;
    } else if(_lFlavor[ind] == 0){
        if(_lElectronMissingHits[ind] > 1) return false;
        if(!_lElectronPassEmu[ind]) return false;
        //if(!elePassVLooseMvaIDSUSY(ind)) return false;
    }
    return true;
}

bool treeReader::lepIsTight(const unsigned l){
    return _lEwkTight[l];
}

bool treeReader::lepIsFOGood_TTV(const unsigned l){

    //used for Cutbased WP

    if(_lPt[l] < 10) return false;

    if(leptonSelection == 2){
        if(_lFlavor[l] == 0 && _relIso[l] > 1.0) return false;
        if(_lFlavor[l] == 1 && _relIso0p4Mu[l] > 0.6) return false;
    }

    if(leptonSelection == 3){
        if(_lFlavor[l] == 0 && _relIso[l] > 0.7) return false;
        if(_lFlavor[l] == 1 && _relIso0p4Mu[l] > 0.7) return false;
    }

    if(_lFlavor[l] == 1 && !_lPOGMedium[l]) return false;

    if(fabs(_3dIPSig[l]) > 4) return false;

    if(_lFlavor[l] == 0){

       double looseMVAloose[3];

       if(leptonSelection == 2){
           looseMVAloose[0] = -0.56;
           looseMVAloose[1] = -0.72;
           looseMVAloose[2] = -0.1;
       }

       if(leptonSelection == 3){
           looseMVAloose[0] = 0.5;
           looseMVAloose[1] = 0.5;
           looseMVAloose[2] = 0.357;
       }

       int etaCategory = -1;

       if (TMath::Abs(_lEta[l]) < 0.8 ) {
           etaCategory = 0;
       } else if (TMath::Abs(_lEta[l]) < 1.479 ) {
           etaCategory = 1;
       } else {
           etaCategory = 2;
       }
       bool passedMVA = false;

       passedMVA = _lElectronMva[l] > looseMVAloose[etaCategory];

       if(!passedMVA) return false;

        if(leptonSelection == 2){

            if(!_lElectronPassConvVeto[l]) return false;
            if(!_lElectronChargeConst[l]) return false;
            if(_lElectronMissingHits[l] != 0) return false;

        }

        if(!_lElectronPassEmu[l]) return false;
    }

    return true;
}

bool treeReader::lepIsGood_TTV(const unsigned l){

    // used for cutBased analysis

    if(!lepIsFOGood_TTV(l)) return false;

    if(_lFlavor[l] == 0){

       double looseMVA[3];

       looseMVA[0] = 0.837;
       looseMVA[1] = 0.715;
       looseMVA[2] = 0.357;

       int etaCategory = -1;

       if (TMath::Abs(_lEta[l]) < 0.8 ) {
           etaCategory = 0;
       } else if (TMath::Abs(_lEta[l]) < 1.479 ) {
           etaCategory = 1;
       } else {
           etaCategory = 2;
       }
       bool passedMVA = false;

       passedMVA = _lElectronMva[l] > looseMVA[etaCategory];

       if (!passedMVA) return false;
       if (_relIso[l] > 0.1) return false;


    }

    if(leptonSelection == 2){
        // what's there originally
        if(_lFlavor[l] == 1 && _relIso0p4Mu[l] > 0.15) return false;
        //if(_lFlavor[l] == 1 && _relIso[l] > 0.1) return false;
    }

    if(leptonSelection == 3){
        if(_lFlavor[l] == 1 && _relIso0p4Mu[l] > 0.25) return false;
        //if(_lFlavor[l] == 1 && _relIso[l] > 0.1) return false;
    }

   return true;
}

unsigned treeReader::selectLep(std::vector<unsigned>& ind){
    //setConePt(); REMOVE CONE CORRECTION UNTIL MOVING TO FR
    unsigned lCount = 0;
    std::vector<std::pair<double, unsigned>> ptMap;
    for(unsigned l = 0; l < _nLight; ++l){
        //cout << "lepton info: " << _lPt[l] << " " << _lEwkLoose[l] << " " << _leptonMvaTTH[l] << endl;
        if(lepIsGood(l)){
            ++lCount;
            ptMap.push_back({_lPt[l], l});
        }
    }
    //if(lCount < 2) return 0;
    std::sort(ptMap.begin(), ptMap.end(), [](std::pair<double, unsigned>& p1, std::pair<double, unsigned>& p2){return p1.first > p2.first;} );
    for(unsigned l = 0; l < lCount; ++l){
        ind.push_back(ptMap[l].second);
    }
    return lCount;	
}

unsigned treeReader::selectLooseLep(std::vector<unsigned>& ind){
    //setConePt(); REMOVE CONE CORRECTION UNTIL MOVING TO FR
    unsigned lCount = 0;
    std::vector<std::pair<double, unsigned>> ptMap;
    for(unsigned l = 0; l < _nLight; ++l){
        //cout << "lepton info: " << _lPt[l] << " " << _lEwkLoose[l] << " " << _leptonMvaTTH[l] << endl;
        if(lepIsLoose(l)){
            ++lCount;
            ptMap.push_back({_lPt[l], l});
        }
    }
    //if(lCount < 2) return 0;
    std::sort(ptMap.begin(), ptMap.end(), [](std::pair<double, unsigned>& p1, std::pair<double, unsigned>& p2){return p1.first > p2.first;} );
    for(unsigned l = 0; l < lCount; ++l){
        ind.push_back(ptMap[l].second);
    }
    return lCount;  
}

double treeReader::ptFake(double lpt, double ptratio, int flavour, double mvaTTHvalue, bool mediumIdPassed) {

    /*
    if (abs(lep.pdgId)!=13 or lep.mediumMuonId>0) and lep.mvaTTH > 0.90: return lep.pt
    else: return 0.90 * lep.pt / lep.jetPtRatiov2
    */

    double ptf = lpt;
    if((flavour == 0 || (flavour == 1 && mediumIdPassed)) && mvaTTHvalue > 0.9)
      ptf = lpt;
    else
      ptf = 0.9 * lpt / ptratio;
    return ptf;
}

unsigned treeReader::tightLepCount(const std::vector<unsigned>& ind, const unsigned lCount){
    unsigned tightC = 0; 
    for(unsigned l = 0; l < lCount; ++l){
        if(lepIsTight(l)) ++tightC;
        else return tightC;
    }
    return tightC;
}

bool treeReader::passPtCuts3L(const std::vector<unsigned>& ind){
    
    std::vector<std::pair<double, unsigned>> ptMap;
    for(auto & i : ind){
        double ptcor = ptFake(_lPt[i], _ptRatio[i], _lFlavor[i], _leptonMvaTTH[i], _lPOGMedium[i]);
        ptMap.push_back({ptcor, i});
        
    }
    std::sort(ptMap.begin(), ptMap.end(), [](std::pair<double, unsigned>& p1, std::pair<double, unsigned>& p2){return p1.first > p2.first;} );

    //cout << "two leptons with pt: " << ptMap[0].first << " " << ptMap[1].first << endl;

    ptCorrV.clear();
    for(auto & i : ptMap)
    ptCorrV.push_back(i);
   
    if(ptMap[0].first < 40) return false;
    if(ptMap[1].first < 20) return false;
    if(ptMap[2].first < 10) return false;

    return true;
}

bool treeReader::passPtCuts2L(const std::vector<unsigned>& ind){
    
    std::vector<std::pair<double, unsigned>> ptMap;
    for(auto & i : ind){
        ptMap.push_back({_lPt[i], i});
        
    }
    std::sort(ptMap.begin(), ptMap.end(), [](std::pair<double, unsigned>& p1, std::pair<double, unsigned>& p2){return p1.first > p2.first;} );

    //cout << "two leptons with pt: " << ptMap[0].first << " " << ptMap[1].first << endl;

    ptCorrV.clear();
    for(auto & i : ptMap)
    ptCorrV.push_back(i);

    if(ptMap[0].first < 25) return false;
    if(ptMap[1].first < 15) return false;

    return true;
}

bool treeReader::jetIsClean(const unsigned ind, bool nonpromptSample){
    TLorentzVector jet;	
    jet.SetPtEtaPhiE(_jetPt[ind], _jetEta[ind], _jetPhi[ind], _jetE[ind]);
    for(unsigned l = 0; l < _nLight; ++l){
        if(!nonpromptSample && lepIsLoose(l)){ // cleaning with FO objects
            TLorentzVector lep;
            lep.SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
            //cout << "jet lepton cleaning is going on, delta R is: " << lep.DeltaR(jet) << endl;
            if(lep.DeltaR(jet) < 0.4) return false;
        }
    }
    return true;
}

bool treeReader::jetIsGood(const unsigned ind, const unsigned ptCut, const unsigned unc, const bool clean, bool nonpromptSample){
    //cout << "jet info (pt/eta/csv/unc/clean): " << _jetPt[ind] << " " << _jetEta[ind] << " " << _jetCsvV2[ind] << " " << unc << " " << clean << endl; 
    if(fabs(_jetEta[ind]) > 2.4) return false;
    //cout << "still here" << endl; 
    switch(unc){
        
        case 0: if(_jetPt[ind] < ptCut) return false; break;
        case 1: if(_jetPt_JECDown[ind] < ptCut) return false; break;
        case 2: if(_jetPt_JECUp[ind] < ptCut) return false; break;
        case 3: if(_jetPt_JERDown[ind] < ptCut) return false; break;
        case 4: if(_jetPt_JERUp[ind] < ptCut) return false; break;
        default: ;
    }
    //cout << "jet info (pt/eta/csv): " << _jetPt[ind] << " " << _jetEta[ind] << " " << _jetCsvV2[ind] << endl; 
    return !clean || jetIsClean(ind, nonpromptSample);
}

unsigned treeReader::nJets(const unsigned unc, const bool clean, std::vector<unsigned>& ind, bool nonpromptSample){
    unsigned nJets = 0;
    for(unsigned j = 0; j < _nJets; ++j){
        //cout << "jet info (pt/eta/csv): " << _jetPt[j] << " " << _jetEta[j] << " " << _jetCsvV2[j] << endl; 
        if(jetIsGood(j, 30, unc, clean, nonpromptSample)) {
            //cout << "this jet passes selection" << endl;
            ++nJets;
            ind.push_back(j);
        }
    }
    return nJets;
}

bool treeReader::bTaggedDeepCSV(const unsigned ind, const unsigned wp){
    static const double bTagWP[3] = {0.2219, 0.6324,  0.8958};
    return (_jetDeepCsv_b[ind] + _jetDeepCsv_bb[ind]) > bTagWP[wp];
}

bool treeReader::bTaggedCSVv2(const unsigned ind, const unsigned wp){
    static const double bTagWP[3] = {0.5426, 0.8484, 0.9535};
    return _jetCsvV2[ind] > bTagWP[wp];
}

unsigned treeReader::nBJets(const unsigned unc, const bool deepCSV, const bool clean, const unsigned wp, bool nonpromptSample){
    unsigned nbJets = 0;
    //cout << "nbjets calculation is starting here >>>>>>>>>>>>>>> " << endl;
    for(unsigned j = 0; j < _nJets; ++j){
        if(jetIsGood(j, 30, unc, clean, nonpromptSample)){
            
            if(deepCSV && bTaggedDeepCSV(j, wp)) ++nbJets;
            else if(bTaggedCSVv2(j, wp)) ++nbJets;
            
        }
    }
    return nbJets;
}

double treeReader::HTCalc(const std::vector<unsigned>& ind){
    double HT = 0;
    for(auto & i : ind)
        HT += _jetPt[i];
    return HT;
}

double treeReader::deltaRCalc(const std::vector<unsigned>& ind, unsigned & lept){
    TLorentzVector leptVec;
    leptVec.SetPtEtaPhiE(_lPt[lept], _lEta[lept], _lPhi[lept], _lE[lept]);

    double deltaRloc = 9999;
    for(auto & jetInd : ind){
        TLorentzVector jet;
        jet.SetPtEtaPhiE(_jetPt[jetInd], _jetEta[jetInd], _jetPhi[jetInd], _jetE[jetInd]);

        if(leptVec.DeltaR(jet) < deltaRloc)
            deltaRloc = leptVec.DeltaR(jet);
    }
    
    return deltaRloc;
}

double treeReader::deltaMZ(const std::vector<unsigned>& ind, unsigned & third, double & mll, double & ptZ, double & ptNonZ, double & phiZ){
                        
    TLorentzVector l0p4, l1p4;

    double deltaMZ = 999999.;

    for (unsigned l0 = 0; l0 < ind.size(); l0++) {
        double ptcor1 = ptFake(_lPt[ind.at(l0)], _ptRatio[ind.at(l0)], _lFlavor[ind.at(l0)], _leptonMvaTTH[ind.at(l0)], _lPOGMedium[ind.at(l0)]);
        //cout << "l0 is " << ind.at(l0) << endl;
        
        //l0p4.SetPtEtaPhiE(ptcor1 ,_lEta[ind.at(l0)],_lPhi[ind.at(l0)],_lE[ind.at(l0)] * ptcor1 / _lPt[ind.at(l0)] );          
        l0p4.SetPtEtaPhiE(_lPt[ind.at(l0)],_lEta[ind.at(l0)],_lPhi[ind.at(l0)],_lE[ind.at(l0)]);          
        for(unsigned l1 = l0; l1 < ind.size(); l1++){
            double ptcor2 = ptFake(_lPt[ind.at(l1)], _ptRatio[ind.at(l1)], _lFlavor[ind.at(l1)], _leptonMvaTTH[ind.at(l1)], _lPOGMedium[ind.at(l1)]);
            //cout << "l1 is " << ind.at(l1) << endl;
            
            if(ind.at(l0) == ind.at(l1)) continue;

            //cout << "info about 2 leptons(pt/flavor): " <<  _lPt[ind.at(l0)] << " " << _lFlavor[ind.at(l0)] << "; second lepton: " << _lPt[ind.at(l1)] << " " << _lFlavor[ind.at(l1)] << endl;
            //cout << "check the charge of 2 leptons: " << _lCharge[ind.at(l0)] << " " << _lCharge[ind.at(l1)] << " " << (_lCharge[ind.at(l0)] != _lCharge[ind.at(l1)]) << endl;

            if (_lCharge[ind.at(l0)] != _lCharge[ind.at(l1)]) {
                //l1p4.SetPtEtaPhiE(ptcor2 ,_lEta[ind.at(l1)],_lPhi[ind.at(l1)],_lE[ind.at(l1)] * ptcor2 / _lPt[ind.at(l1)]);
                l1p4.SetPtEtaPhiE(_lPt[ind.at(l1)],_lEta[ind.at(l1)],_lPhi[ind.at(l1)],_lE[ind.at(l1)]);
                l1p4+=l0p4;
                double mdiL = l1p4.M();
                
                //cout << "invariant mass of 2 leptons is: " << mdiL << endl;
                if (_lFlavor[ind.at(l0)] == _lFlavor[ind.at(l1)] ) {

                    //cout << "the leptons are of the same flavour with inv. mass: " << mdiL << endl;

                    if (fabs(mdiL - 91.2) < deltaMZ) {
                        deltaMZ = fabs(mdiL - 91.2);
                        mll = mdiL;
                        ptZ = l1p4.Pt();
                        phiZ = l1p4.Phi();
                        
                        if(leptonSelection == 3){
                            for(auto & lepThird : ind){
                                if(lepThird == ind.at(l0) || lepThird == ind.at(l1)) continue;
                                third = lepThird;
                                double ptcor3 = ptFake(_lPt[lepThird], _ptRatio[lepThird], _lFlavor[lepThird], _leptonMvaTTH[lepThird], _lPOGMedium[lepThird]);
                                ptNonZ = ptcor3;
                            }
                        }
                                
                    }
                }
            }
        }
    }

    return deltaMZ;
}

int treeReader::getElectronNumber(const std::vector<unsigned>& ind){
    int electronNumber = 0;

    for(auto & l : ind){
        if(_lFlavor[l] == 0) electronNumber++;
    }

    return electronNumber;
}

bool treeReader::promptLeptons(const std::vector<unsigned>& ind){
    bool allPrompt = true;

    for(auto & l : ind){
        //cout << "check this lepton: " << _lPt[l] << " " << _lEta[l] << " " << _lPhi[l] << endl;
        bool leptIsPr = leptonIsPrompt(l);
        allPrompt = allPrompt && leptIsPr; // _lIsPrompt[l];
    }
    
    return allPrompt;
}

bool treeReader::leptonIsPrompt(const unsigned & l){

    bool leptIsP = false;
    TLorentzVector l0reco;
    l0reco.SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
    //cout << "let's match, number of gen leptons: " << _gen_nL << endl;
    for(unsigned i = 0; i < _gen_nL; i++){
        TLorentzVector l0gen;
        //cout << "let's match: " << _gen_lPt[i] << " " << _gen_lEta[i] << " " << _gen_lPhi[i] << endl;
        l0gen.SetPtEtaPhiE(_gen_lPt[i], _gen_lEta[i], _gen_lPhi[i], _gen_lE[i]);
        //cout << "the deltaR is: " << l0reco.DeltaR(l0gen) << endl;
        if(l0reco.DeltaR(l0gen) < 1.0){
            //cout << "lepton is matched to (genFl/dR/genIsPr/momPDGid) " << _gen_lFlavor[i] << " " << l0reco.DeltaR(l0gen) << " " << _gen_lIsPrompt[i]  << " " << _gen_lMomPdg[i] << endl;
            if(_gen_lIsPrompt[i] && (_lFlavor[l] == _gen_lFlavor[i] || _gen_lFlavor[i] == 2)){ //  && _lFlavor[l] == _gen_lFlavor[i]
                leptIsP = true;
            }
        }
    }
    
    //cout << "what the matching function returns: " << leptIsP << endl;
    return leptIsP;
}

bool treeReader::leptonIsFromPromptTau(const unsigned & l){

    bool leptIsP = false;
    TLorentzVector l0reco;
    l0reco.SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
    //cout << "let's match, number of gen leptons: " << _gen_nL << endl;
    for(unsigned i = 0; i < _gen_nL; i++){
        TLorentzVector l0gen;
        //cout << "let's match: " << _gen_lPt[i] << " " << _gen_lEta[i] << " " << _gen_lPhi[i] << endl;
        l0gen.SetPtEtaPhiE(_gen_lPt[i], _gen_lEta[i], _gen_lPhi[i], _gen_lE[i]);
        //cout << "the deltaR is: " << l0reco.DeltaR(l0gen) << endl;
        if(l0reco.DeltaR(l0gen) < 0.1){
            //cout << "lepton is matched to (genFl/dR/genIsPr/momPDGid/provenance) " << _gen_lFlavor[i] << " " << l0reco.DeltaR(l0gen) << " " << _gen_lIsPrompt[i]  << " " << _gen_lMomPdg[i] << " " << _lProvenance[l]<< " " << _lProvenanceCompressed[l] << endl;
            if(_gen_lIsPrompt[i] && _gen_lFlavor[i] == 2){ //  && _lFlavor[l] == _gen_lFlavor[i]
                leptIsP = true;
            }
        }
    }
    
    //cout << "what the matching function returns: " << leptIsP << endl;
    return leptIsP;
}

Color_t treeReader::assignColor(std::string & name){

    if(name == "data") return kBlack;
    if(name == "nonprompt") return kBlue-9;
    if(name == "nonpromptData") return kBlue-9;
    if(name == "chargeMisID") return kMagenta-7;
    if(name == "ttZ") return 91;
    if(name == "ttW") return 98;
    if(name == "ttH") return kRed-10;
    if(name == "ttX") return kRed-10;
    if(name == "WZ") return 51;
    if(name == "ZZ") return 8;
    if(name == "rare") return 8;

    if(name == "DY") return kBlue-9;
    if(name == "Diboson") return 98; 
    if(name == "Triboson") return 8;
    if(name == "tight") return 91; 
    if(name == "loose") return 51; 
    if(name == "tW") return kRed-10;

    return kBlack;
}




