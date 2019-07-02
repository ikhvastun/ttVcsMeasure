#include "THStack.h"

#include "../interface/treeReader.h"

using namespace std;

bool treeReader::lepIsGoodFortZq(const unsigned l){
    // what is used in leptonMVA analysis
    //cout << "lepton info: " << _lPt[l] << " " << _lEta[l] << " " << _lFlavor[l] << " " << _leptonMvatZqTTV[l] << " " << _lPOGMedium[l] << " " << _lCharge[l] << endl;
    if(!lepIsFOGoodFortZq(l)) return false;
    if(_leptonMvatZqTTV[l] < 0.8) return false;

    return true;
}

bool treeReader::lepIsGood(const unsigned l, const int lepSel = 3){
    // what is used in leptonMVA analysis
    //cout << "lepton info: " << _lPt[l] << " " << _lEta[l] << " " << _lFlavor[l] << " " << _leptonMvatZqTTV[l] << " " << _lPOGMedium[l] << " " << _lCharge[l] << endl;
    if(!lepIsFOGood(l, lepSel)) return false;
    if(_leptonMvatZqTTV[l] < leptonMVAcutInAnalysis[lepSel]) return false;

    if(lepSel == 2){
        if(_lFlavor[l] == 1 && ((_lMuonTrackPtErr[l]/_lMuonTrackPt[l]) > 0.2)) return false;
        if(_lFlavor[l] == 0 && !_lElectronChargeConst[l]) return false;
        if(_lFlavor[l] == 0 && !_lElectronPassConvVeto[l]) return false;
        if(_lFlavor[l] == 0 && _lElectronMissingHits[l] != 0) return false;
    }

    //if(debug) cout << "lepton with flavour " << _lFlavor[l] << ", charge " << _lCharge[l] << ", pt " << _lPt[l] << ", eta " << _lEta[l] << ", phi " << _lPhi[l] << " and lepton MVA " << _leptonMvatZqTTV[l] << " is tight for " << lepSel << " selection" << endl;
    return true;
}

bool treeReader::lepIsFOGood(const unsigned l, const int lepSel = 3){

    if(_lPt[l] < 10) return false;
    if(!lepIsLoose(l)) return false;

    //if(debug) cout << "lepton with pt " << _lPt[l] << " is clean? " <<  eleIsClean(l) << endl;
    if(_lFlavor[l] == 0 && !eleIsClean(l)) return false;
    if(_lFlavor[l] == 1 && !_lPOGMedium[l]) return false;

    if(lepSel == 2){
        if(_lFlavor[l] == 1 && ((_lMuonTrackPtErr[l]/_lMuonTrackPt[l]) > 0.2)) return false;
        if(_lFlavor[l] == 0 && !_lElectronChargeConst[l]) return false;
        if(_lFlavor[l] == 0 && !_lElectronPassConvVeto[l]) return false;
        if(_lFlavor[l] == 0 && _lElectronMissingHits[l] != 0) return false;
    }

    //if(debug) cout << "info about FO leptons with pt " << _lPt[l] << " (closest jet deep csv): " << _closestJetDeepCsv_bb[l] + _closestJetDeepCsv_b[l]  << endl;
	 // here we tune the FO definition
	 // _ptRatio
	 // _closestJetDeepCsv_bb[l] + _closestJetDeepCsv_b[l]
	 // electronMVAvalue cut value
    if(_closestJetDeepCsv_bb[l] + _closestJetDeepCsv_b[l] > (currentSample.is2017() ? 0.8001 : 0.8958)) return false;

    if(lepSel != 4 && _leptonMvatZqTTV[l] < leptonMVAcutInAnalysis[leptonSelection]){
        if(_ptRatio[l] < (currentSample.is2017() ? (leptonSelection == 3 ? 0.4 : 0.5) : (leptonSelection == 3 ? 0.4 : 0.5))) return false;  // 0.4 original for 3L in 2016, 0.3 in 2017
        if(_closestJetDeepCsv_bb[l] + _closestJetDeepCsv_b[l] > (currentSample.is2017() ? (leptonSelection == 3 ? 0.5 : 0.2) : (leptonSelection == 3 ? 0.4 : 0.5))) return false; // 0.4 original one, medium WP :          (is2017 ? 0.4941 : 0.6324)
        double electronMVAvalue = currentSample.is2017() ? _lElectronMvaFall17NoIso[l] : _lElectronMva[l];
        if(_lFlavor[l] == 0 && electronMVAvalue < (currentSample.is2017() ? (leptonSelection == 3 ? -0.3 : -0.1) : (leptonSelection == 3 ? -0.1 : -0.6))    +    (fabs(_lEta[l]) >= 1.479)*(currentSample.is2017() ? (leptonSelection ==    3 ? 0.6 : 0.4) : (leptonSelection == 3 ? 0.8 : 0.4))) return false;
        // numbers for tZq
        /*
        if(_ptRatio[l] < (is2017 ? 0.6 : 0.6)) return false; // 0.3 for 2017
        if(_closestJetDeepCsv_bb[l] + _closestJetDeepCsv_b[l] > (is2017 ? 0.2 : 0.3)) return false;
        double electronMVAvalue = is2017 ? _lElectronMvaFall17NoIso[l] : _lElectronMva[l];
        if(_lFlavor[l] == 0 && electronMVAvalue < (is2017 ? 0.3 : 0.4) + (fabs(_lEta[l]) >= 1.479)*(is2017 ? 0.3 : 0.3)) return false;
        */
    }
    //cout << "info about FO leptons: " << _lPt[l] << " " << _lEta[l] << " " << _lFlavor[l] << " " << _lCharge[l] << " " << _leptonMvatZqTTV[l]  << endl;

    return true;
}

bool treeReader::lepIsFOGoodFortZq(const unsigned l){

    if(_lPt[l] < 10) return false;
    if(!lepIsLoose(l)) return false;

    //if(debug) cout << "lepton with pt " << _lPt[l] << " is clean? " <<  eleIsClean(l) << endl;
    if(_lFlavor[l] == 0 && !eleIsClean(l)) return false;
    if(_lFlavor[l] == 1 && !_lPOGMedium[l]) return false;

    //if(debug) cout << "info about FO leptons with pt " << _lPt[l] << " (closest jet deep csv): " << _closestJetDeepCsv_bb[l] + _closestJetDeepCsv_b[l]  << endl;
    if(_closestJetDeepCsv_bb[l] + _closestJetDeepCsv_b[l] > (currentSample.is2017() ? 0.8001 : 0.8958)) return false;

    if(_leptonMvatZqTTV[l] < 0.8){
        if(_ptRatio[l] < (currentSample.is2017() ? 0.6 : 0.6)) return false; // 0.3 for 2017
        if(_closestJetDeepCsv_bb[l] + _closestJetDeepCsv_b[l] > (currentSample.is2017() ? 0.3 : 0.4)) return false;
        double electronMVAvalue = currentSample.is2017() ? _lElectronMvaFall17NoIso[l] : _lElectronMva[l];
        if(_lFlavor[l] == 0 && electronMVAvalue < (currentSample.is2017() ? 0.3 : 0.4) + (fabs(_lEta[l]) >= 1.479)*(currentSample.is2017() ? 0.3 : 0.3)) return false;
    }
    return true;
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
    if(_lFlavor[ind] == 2) return false;  //don't consider taus here
    //if(debug) cout << "info about loose lepton with pt " << _lPt[ind] << ": " << (fabs(_dxy[ind]) >= 0.05) << " " << (fabs(_dz[ind]) >= 0.1) << " " << (_3dIPSig[ind] >= 8) << " " << (_miniIso[ind] >= 0.4) << " " << (_lElectronMissingHits[ind] > 1) << " " << (!_lElectronPassEmu[ind]) << " " << !_lPOGLoose[ind] << " " << !_lPOGMedium[ind] << endl;
    //cout << "lepton miniiso: " << _miniIso[ind] << endl;
    //if(_lPt[ind] <= 7 - 2*_lFlavor[ind]) return false;
    if(_lPt[ind] <= 10) return false;
    if(fabs(_lEta[ind]) >= (2.5 - 0.1*_lFlavor[ind])) return false;
    if(fabs(_dxy[ind]) >= 0.05) return false;
    if(fabs(_dz[ind]) >= 0.1) return false;
    if(_3dIPSig[ind] >= 8) return false;
    if(_miniIso[ind] >= 0.4) return false;
    if(_lFlavor[ind] == 1){
        if(!_lPOGLoose[ind]) return false;
    } else if(_lFlavor[ind] == 0){
        if(_lElectronMissingHits[ind] > 1) return false;
        if(!_lElectronPassEmu[ind]) return false;
    }
    return true;
}

bool treeReader::lepIsGoodForLeptonMVA(const unsigned ind){
    if(_lFlavor[ind] == 2) return false;  //don't consider taus here
    if(_lPt[ind] <= 10) return false;
    if(fabs(_lEta[ind]) >= (2.5 - 0.1*_lFlavor[ind])) return false;
    if(fabs(_dxy[ind]) >= 0.05) return false;
    if(fabs(_dz[ind]) >= 0.1) return false;
    if(_3dIPSig[ind] >= 8) return false;
    if(_miniIso[ind] >= 0.4) return false;
    if(_lFlavor[ind] == 1){
        if(!_lPOGLoose[ind]) return false;
        if(!_lPOGMedium[ind]) return false;
    } else if(_lFlavor[ind] == 0){
        if(_lElectronMissingHits[ind] > 1) return false;
        if(!_lElectronPassEmu[ind]) return false;
    }
    return true;
}

unsigned treeReader::selectLepGoodForLeptonMVA(std::vector<unsigned>& ind){
    unsigned lCount = 0;
    std::vector<std::pair<double, unsigned>> ptMap;
    for(unsigned l = 0; l < _nLight; ++l){
        if(lepIsGoodForLeptonMVA(l)){
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

unsigned treeReader::selectLep(std::vector<unsigned>& ind, const int lepSel = 3){
    unsigned lCount = 0;
    std::vector<std::pair<double, unsigned>> ptMap;
    for(unsigned l = 0; l < _nLight; ++l){
        if(lepIsGood(l, lepSel)){
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

unsigned treeReader::selectFakeLep(std::vector<unsigned>& ind, const int lepSel = 3){
    unsigned lCount = 0;
    std::vector<std::pair<double, unsigned>> ptMap;
    for(unsigned l = 0; l < _nLight; ++l){
        if(lepIsFOGood(l, lepSel)){
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

unsigned treeReader::selectLooseLep(std::vector<unsigned>& ind, const int lepSel = 3){
    unsigned lCount = 0;
    std::vector<std::pair<double, unsigned>> ptMap;
    for(unsigned l = 0; l < _nLight; ++l){
        if(_lPt[l] < 10) continue;
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

bool treeReader::passPtCuts4L(const std::vector<unsigned>& ind){
    
    std::vector<std::pair<double, unsigned>> ptMap;
    for(auto & i : ind){
        double ptcor = _lPt[i];
        ptMap.push_back({ptcor, i});
    }
    std::sort(ptMap.begin(), ptMap.end(), [](std::pair<double, unsigned>& p1, std::pair<double, unsigned>& p2){return p1.first > p2.first;} );

    ptCorrV.clear();
    for(auto & i : ptMap)
    ptCorrV.push_back(i);
    
    if(ptMap[0].first < 40) return false;
    if(ptMap[1].first < 10) return false;
    if(ptMap[2].first < 10) return false;
    if(ptMap[3].first < 10) return false;
    
    return true;
}

bool treeReader::passPtCuts3L(const std::vector<unsigned>& ind){
    
    std::vector<std::pair<double, unsigned>> ptMap;
    for(auto & i : ind){
        double ptcor = _leptonMvatZqTTV[i] > leptonMVAcutInAnalysis[leptonSelection] ? _lPt[i] : magicFactorInAnalysis[leptonSelection] * _lPt[i] / _ptRatio[i];
        ptMap.push_back({ptcor, i});
        
    }
    std::sort(ptMap.begin(), ptMap.end(), [](std::pair<double, unsigned>& p1, std::pair<double, unsigned>& p2){return p1.first > p2.first;} );

    ptCorrV.clear();
    for(auto & i : ptMap)
    ptCorrV.push_back(i);
    
    //if(ptMap[0].first > 100 && _lFlavor[ptMap[0].second] == 1) return false;

    if(ptMap[0].first < 40) return false;
    if(ptMap[1].first < 20) return false;
    if(ptMap[2].first < 10) return false;
    
    return true;
}

bool treeReader::passPtCuts2L(const std::vector<unsigned>& ind){
    
    std::vector<std::pair<double, unsigned>> ptMap;
    for(auto & i : ind){
        double ptcor = _leptonMvatZqTTV[i] > leptonMVAcutInAnalysis[leptonSelection] ? _lPt[i] : magicFactorInAnalysis[leptonSelection] * _lPt[i] / _ptRatio[i];
        //double ptcor = lepIsGood(i) ? _lPt[i] : magicFactor * _lPt[i] / _ptRatio[i];
        ptMap.push_back({ptcor, i});
    }
    std::sort(ptMap.begin(), ptMap.end(), [](std::pair<double, unsigned>& p1, std::pair<double, unsigned>& p2){return p1.first > p2.first;} );

    ptCorrV.clear();
    for(auto & i : ptMap)
        ptCorrV.push_back(i);

    if(ptMap[0].first < 25 && _lFlavor[ptMap[0].second] == 1) return false;
    if(ptMap[1].first < 25 && _lFlavor[ptMap[1].second] == 1) return false;

    if(ptMap[0].first < 40 && _lFlavor[ptMap[0].second] == 0) return false;
    if(ptMap[1].first < 27 && _lFlavor[ptMap[1].second] == 0) return false;
       
    return true;
}

bool treeReader::passPtCuts2LOF(const std::vector<unsigned>& ind){
    
    std::vector<std::pair<double, unsigned>> ptMap;
    for(auto & i : ind){
        double ptcor = _lPt[i];
        ptMap.push_back({ptcor, i});
    }
    std::sort(ptMap.begin(), ptMap.end(), [](std::pair<double, unsigned>& p1, std::pair<double, unsigned>& p2){return p1.first > p2.first;} );

    ptCorrV.clear();
    for(auto & i : ptMap)
        ptCorrV.push_back(i);

    if(ptMap[0].first < 10) return false;
    if(ptMap[1].first < 10) return false;
       
    return true;
}

bool treeReader::passPtCuts2LOSSF(const std::vector<unsigned>& ind){
    
    std::vector<std::pair<double, unsigned>> ptMap;
    for(auto & i : ind){
        double ptcor = _lPt[i];
        ptMap.push_back({ptcor, i});
    }
    std::sort(ptMap.begin(), ptMap.end(), [](std::pair<double, unsigned>& p1, std::pair<double, unsigned>& p2){return p1.first > p2.first;} );

    ptCorrV.clear();
    for(auto & i : ptMap)
        ptCorrV.push_back(i);

    if(ptMap[0].first < 40) return false;
    if(ptMap[1].first < 25) return false;
       
    return true;
}

bool treeReader::jetIsClean(const unsigned ind, const int lepSel){
    TLorentzVector jet;	
    jet.SetPtEtaPhiE(isData ? _jetPt[ind] : _jetSmearedPt[ind], _jetEta[ind], _jetPhi[ind], _jetE[ind] * _jetSmearedPt[ind] /  _jetPt[ind]); 
    //jet.SetPtEtaPhiE(_jetPt[ind], _jetEta[ind], _jetPhi[ind], _jetE[ind]); 
    for(unsigned l = 0; l < _nLight; ++l){
        //if(lepIsFOGoodFortZq(l)){ // cleaning with FO objects
        if(lepIsFOGood(l, lepSel)){ // cleaning with FO objects
            TLorentzVector lep;
            lep.SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
            //cout << "jet lepton cleaning is going on, delta R is: " << lep.DeltaR(jet) << endl;
            if(lep.DeltaR(jet) < 0.4) return false;
        }
    }
    return true;
}

bool treeReader::jetIsGood(const unsigned ind, const unsigned ptCut, const unsigned unc, const bool clean, bool is2017, const double eta){
    if(fabs(_jetEta[ind]) > eta) return false;
    if(is2017 && !_jetIsTight[ind]) return false;
    // temporary fix to sync with Daniel, 9 july 2018
    //if(!_jetIsTight[ind]) return false;
    switch(unc){
        case 0: if((isData ? _jetPt[ind] : _jetSmearedPt[ind]) < ptCut) return false; break;
        case 1: if(_jetSmearedPt_JECDown[ind] < ptCut) return false; break;
        case 2: if(_jetSmearedPt_JECUp[ind] < ptCut) return false; break;
        case 3: if(_jetSmearedPt_JERDown[ind] < ptCut) return false; break;
        case 4: if(_jetSmearedPt_JERUp[ind] < ptCut) return false; break;

        /*
        case 0: if((_jetPt[ind]) < ptCut) return false; break;
        case 1: if(_jetPt_JECDown[ind] < ptCut) return false; break;
        case 2: if(_jetPt_JECUp[ind] < ptCut) return false; break;
        case 3: if(_jetPt_JERDown[ind] < ptCut) return false; break;
        case 4: if(_jetPt_JERUp[ind] < ptCut) return false; break;
        */
        default: ;
    }
    return !clean || jetIsClean(ind, leptonSelection);
}

unsigned treeReader::nJets(const unsigned unc, const bool clean, std::vector<unsigned>& ind, bool is2017, const double eta){
    unsigned nJets = 0;
    for(unsigned j = 0; j < _nJets; ++j){
        //if(debug) cout << "jet with pt: " << _jetPt[j] << " " << _jetEta[j] << " " << _jetPhi[j] << ", and flavour: " << _jetHadronFlavor[j] << endl;
        if(jetIsGood(j, 30, unc, clean, is2017, eta)) {
            //if(debug) cout << "passed selection" << endl;
            ++nJets;
            ind.push_back(j);
        }
    }
    return nJets;
}

bool treeReader::isForwardJetPresent(bool is2017){
    bool presence = false;
    for(unsigned j = 0; j < _nJets; ++j){
        if(is2017 && !_jetIsTight[j]) continue;
        if(_jetPt[j] < 50) continue;
        if(TMath::Abs(_jetEta[j]) < 2.4) continue;
        if(!jetIsClean(j, leptonSelection)) continue;
        presence = true;
    }
    return presence;
}

unsigned treeReader::nJetsNotB(const unsigned unc, const bool clean, std::vector<unsigned>& ind, const unsigned btagWP, bool is2017){
    unsigned nJets = 0;
    for(unsigned j = 0; j < _nJets; ++j){
        if(jetIsGood(j, 30, unc, clean, is2017) && !bTaggedDeepCSV(j, btagWP)) {
            ++nJets;
            ind.push_back(j);
        }
    }
    return nJets;
}

bool treeReader::bTaggedDeepCSV(const unsigned ind, const unsigned wp){
    static const double bTagWP[3] = {(currentSample.is2017() ? 0.1522 : 0.2219), (currentSample.is2017() ? 0.4941 : 0.6324),  (currentSample.is2017() ? 0.8001 : 0.8958)};
    return (_jetDeepCsv_b[ind] + _jetDeepCsv_bb[ind]) > bTagWP[wp];
}

bool treeReader::bTaggedCSVv2(const unsigned ind, const unsigned wp){
    static const double bTagWP[3] = {0.5426, 0.8484, 0.9535};
    return _jetCsvV2[ind] > bTagWP[wp];
}

unsigned treeReader::nBJets(const unsigned unc, const bool deepCSV, const bool clean, std::vector<unsigned>& ind, const unsigned wp, bool is2017){
    unsigned nbJets = 0;
    //cout << "nbjets calculation is starting here >>>>>>>>>>>>>>> " << endl;
    for(unsigned j = 0; j < _nJets; ++j){
        if(jetIsGood(j, 30, unc, clean, is2017)){
            //if(fabs(_jetEta[j]) > 2.4) continue;
            if(deepCSV && bTaggedDeepCSV(j, wp)){ ++nbJets; ind.push_back(j);}
            else if(!deepCSV && bTaggedCSVv2(j, wp)){ ++nbJets; ind.push_back(j);}
            
        }
    }
    return nbJets;
}

double treeReader::HTCalc(const std::vector<unsigned>& ind){
    double HT = 0;
    for(auto & i : ind)
        HT += isData ? _jetPt[i] : _jetSmearedPt[i];
        //HT += _jetPt[i];
    return HT;
}

double treeReader::deltaRCalc(const std::vector<unsigned>& ind, unsigned & lept, const bool deepCSV){
    TLorentzVector leptVec;
    leptVec.SetPtEtaPhiE(_lPt[lept], _lEta[lept], _lPhi[lept], _lE[lept]);

    double deltaRloc = 9999;
    for(auto & jetInd : ind){
        //if(!((deepCSV && bTaggedDeepCSV(jetInd, 1)) || (!deepCSV && bTaggedCSVv2(jetInd, 1)))) continue;
        TLorentzVector jet;
        jet.SetPtEtaPhiE(_jetPt[jetInd], _jetEta[jetInd], _jetPhi[jetInd], _jetE[jetInd]);

        if(leptVec.DeltaR(jet) < deltaRloc)
            deltaRloc = leptVec.DeltaR(jet);
    }
    
    return deltaRloc;
}

bool treeReader::invMassOfAny2Lbelow12GeV(const std::vector<unsigned>& ind){
                        
    TLorentzVector l0p4, l1p4;
    bool belowMllheavyFlavourResonances = false;

    for (unsigned l0 = 0; l0 < ind.size(); l0++) {
        l0p4.SetPtEtaPhiE(_lPt[ind.at(l0)],_lEta[ind.at(l0)],_lPhi[ind.at(l0)],_lE[ind.at(l0)] );          
        for(unsigned l1 = l0+1; l1 < ind.size(); l1++){

            l1p4.SetPtEtaPhiE(_lPt[ind.at(l1)],_lEta[ind.at(l1)],_lPhi[ind.at(l1)],_lE[ind.at(l1)]);
            l1p4+=l0p4;
            double mdiL = l1p4.M();
                
            if (mdiL < 12) belowMllheavyFlavourResonances = true;       
        }
    }

    return belowMllheavyFlavourResonances;
}

bool treeReader::invMassOfAny2Lbelow20GeV(const std::vector<unsigned>& ind){
                        
    TLorentzVector l0p4, l1p4;
    bool belowMllheavyFlavourResonances = false;

    for (unsigned l0 = 0; l0 < ind.size(); l0++) {
        l0p4.SetPtEtaPhiE(_lPt[ind.at(l0)],_lEta[ind.at(l0)],_lPhi[ind.at(l0)],_lE[ind.at(l0)] );          
        for(unsigned l1 = l0+1; l1 < ind.size(); l1++){

            l1p4.SetPtEtaPhiE(_lPt[ind.at(l1)],_lEta[ind.at(l1)],_lPhi[ind.at(l1)],_lE[ind.at(l1)]);
            l1p4+=l0p4;
            double mdiL = l1p4.M();
                
            if (mdiL < 20) belowMllheavyFlavourResonances = true;       
        }
    }

    return belowMllheavyFlavourResonances;
}

double treeReader::deltaMZ(const std::vector<unsigned>& ind, unsigned & third, double & mll, double & ptZ, double & ptNonZ, double & mlll, std::vector<unsigned> & indLeptonOnZ, TLorentzVector & Zboson, TLorentzVector & lnegative){
                        
    TLorentzVector l0p4, l1p4, l2p4;

    double deltaMZ = 999999.;

    for (unsigned l0 = 0; l0 < ind.size(); l0++) {
        //double ptcor1 = _lPt[ind.at(l0)];
        double ptcor1 = leptonSelection != 4 ? (_leptonMvatZqTTV[ind.at(l0)] > leptonMVAcutInAnalysis[leptonSelection] ? _lPt[ind.at(l0)] : magicFactorInAnalysis[leptonSelection] * _lPt[ind.at(l0)] / _ptRatio[ind.at(l0)]) : _lPt[ind.at(l0)];
        
        l0p4.SetPtEtaPhiE(ptcor1 ,_lEta[ind.at(l0)],_lPhi[ind.at(l0)],_lE[ind.at(l0)] * ptcor1 / _lPt[ind.at(l0)] );          
        for(unsigned l1 = l0+1; l1 < ind.size(); l1++){

            //if(ind.at(l0) == ind.at(l1)) continue;
            //double ptcor2 = _lPt[ind.at(l1)];
            double ptcor2 = leptonSelection != 4 ? (_leptonMvatZqTTV[ind.at(l1)] > leptonMVAcutInAnalysis[leptonSelection] ? _lPt[ind.at(l1)] : magicFactorInAnalysis[leptonSelection] * _lPt[ind.at(l1)] / _ptRatio[ind.at(l1)]) : _lPt[ind.at(l1)];
            //cout << "l1 is " << ind.at(l1) << endl;

            if (_lCharge[ind.at(l0)] != _lCharge[ind.at(l1)]) {
                l1p4.SetPtEtaPhiE(ptcor2 ,_lEta[ind.at(l1)],_lPhi[ind.at(l1)],_lE[ind.at(l1)] * ptcor2 / _lPt[ind.at(l1)]);
                l1p4+=l0p4;
                double mdiL = l1p4.M();
                
                //cout << "invariant mass of 2 leptons is: " << mdiL << endl;
                if (_lFlavor[ind.at(l0)] == _lFlavor[ind.at(l1)] ) {

                    //cout << "the leptons are of the same flavour with inv. mass: " << mdiL << endl;
                    if (fabs(mdiL - 91.1876) < deltaMZ) {
                        deltaMZ = fabs(mdiL - 91.1876);
                        indLeptonOnZ.clear();
                        indLeptonOnZ.push_back(ind.at(l0));
                        indLeptonOnZ.push_back(ind.at(l1));
                        mll = mdiL;
                        ptZ = l1p4.Pt();
                        Zboson = l1p4;
                        lnegative = _lCharge[l0] == -1 ? l0p4 : (l1p4 - l0p4);
                        if(leptonSelection == 3){
                            for(auto & lepThird : ind){
                                if(lepThird == ind.at(l0) || lepThird == ind.at(l1)) continue;
                                third = lepThird;
                                //double ptcor3 = _lPt[lepThird];
                                double ptcor3 = leptonSelection != 4 ? (_leptonMvatZqTTV[lepThird] > leptonMVAcutInAnalysis[leptonSelection] ? _lPt[lepThird] : magicFactorInAnalysis[leptonSelection] * _lPt[lepThird] / _ptRatio[lepThird]) : _lPt[lepThird];
                                l2p4.SetPtEtaPhiE(ptcor3 ,_lEta[lepThird],_lPhi[lepThird],_lE[lepThird] * ptcor3 / _lPt[lepThird]);
                                ptNonZ = ptcor3;
                                mlll = (l1p4+l2p4).M();
                            }
                        }
                        if(leptonSelection == 4){
                            // here for 3rd lepton let's pick highest pt lepton
                            for(auto & lepThird : ind){
                                if(lepThird == ind.at(l0) || lepThird == ind.at(l1)) continue;
                                if(_lPt[lepThird] > ptNonZ)
                                    ptNonZ = _lPt[lepThird];
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
        //bool leptIsPr = leptonIsPrompt(l);
        bool leptIsPr = _lIsPrompt[l];
        allPrompt = allPrompt && leptIsPr; // _lIsPrompt[l];
    }
    
    return allPrompt;
}

bool treeReader::noConversionInSelection(const std::vector<unsigned>& ind){

    bool allGood = true;

    for(auto & l : ind){
        allGood = allGood && !leptonFromConversion(l);
    }
    
    return allGood;
}

bool treeReader::leptonFromConversion(const unsigned & l){
    return _lMatchPdgId[l] == 22;
}

// function defined by me, better to use CMSSW version
bool treeReader::leptonIsPrompt(const unsigned & l){

    bool leptIsP = false;
    TLorentzVector l0reco;
    l0reco.SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
    for(unsigned i = 0; i < _gen_nL; i++){
        TLorentzVector l0gen;
        //cout << "let's match: " << _gen_lPt[i] << " " << _gen_lEta[i] << " " << _gen_lPhi[i] << endl;
        l0gen.SetPtEtaPhiE(_gen_lPt[i], _gen_lEta[i], _gen_lPhi[i], _gen_lE[i]);
        //cout << "the deltaR is: " << l0reco.DeltaR(l0gen) << endl;
        if(l0reco.DeltaR(l0gen) < 0.1){
            if(_gen_lIsPrompt[i])
                leptIsP = true;
        }
    }
    
    //cout << "what the matching function returns: " << leptIsP << endl;
    return leptIsP;
}

Color_t treeReader::assignColor(const std::string & name){

    if(name == "data") return kBlack;
    if(name == "Nonprompt") return kBlue-9;
    if(name == "nonpromptData") return kBlue-9; // switch to kBlue-9
    if(name == "chargeMisID") return kMagenta-7;
    if(name == "chargeMisIDData") return kMagenta-7;
    if(name == "ttZ") return 91;
    if(name == "ttW") return 98;
    if(name == "ttH") return kRed-10;
    //if(name == "ttH") return kMagenta-7;
    if(name == "ttX") return kRed-10;
    if(name == "WZ") return 51;
    if(name == "ZZ") return kGreen+3;
    if(name == "rare") return 8;
    if(name == "Xgamma") return kGreen;
    if(name == "DY") return 8;
    if(name == "WJets") return kGreen;
    if(name == "VV") return 51;
    if(name == "tX") return 91;
    if(name == "TTto2L") return 98;
    if(name == "TTto1L") return kGreen+3;
    // ttX check
    if(name == "tWZ") return 51;
    if(name == "tZq") return kGreen+3;
    if(name == "tHq") return kGreen;
    if(name == "ttWZ") return 8;

    return kBlack;
}

bool treeReader::twoLeptonsInEndcap(const std::vector<unsigned>& ind){
   for(auto & i : ind){
    if(_lFlavor[i] == 0) continue;
    for(auto & l : ind){
     if(i == l) continue;
     if(_lFlavor[l] == 0) continue;
     if(fabs(_lEta[i]) > 1.6 && fabs(_lEta[l]) > 1.6) return true;
    }
   }
   return false;
}

double treeReader::cosThetaStar(const TLorentzVector & Z_tlv, const TLorentzVector & l_tlv){

    TVector3 Z, l;
    Z.SetPtEtaPhi( Z_tlv.Pt(), Z_tlv.Eta(), Z_tlv.Phi());
    l.SetPtEtaPhi( l_tlv.Pt(), l_tlv.Eta(), l_tlv.Phi());
    // get cos(theta) and the lorentz factor, calculate cos(theta*)
    //cout << "info for further calculations: " << Z_tlv.M() << " " << l_tlv.M() << endl;
    double cosTheta = Z*l / (sqrt(Z*Z) * sqrt(l*l));
    //cout << "cos theta in lab frame: " << cosTheta << endl;
    double gamma   = TMath::Sqrt( 1 + TMath::Power(Z_tlv.Pt()/Z_tlv.M(),2) * TMath::Power(TMath::CosH(Z_tlv.Eta()),2) );
    //cout << "gamma: " << gamma << endl;
    double beta    = sqrt( 1 - 1/TMath::Power(gamma,2) );
    //cout << "beta: " << beta << endl;
    return (-beta + cosTheta) / (1 - beta*cosTheta);

}

bool treeReader::passTTZSRSelection(const std::vector<unsigned>& ind, std::vector<unsigned> indOf2LonZ, const int njets, const int nbjets, const double dMZ){
    if(leptonSelection == 3){
        if(njets > 0 && passWZCRSelection(nbjets, dMZ)) return true;
        if(nbjets > 0 && passTTZSelection(njets, dMZ)) return true;
    }
    if(leptonSelection == 4){
        if(passTTZ4LSelection(ind, indOf2LonZ, njets)) return true;
    }
    return false;
}

bool treeReader::passTTZSelection(const int njets, const double dMZ) const{
    if(njets < 2) return false;
    if(dMZ > 10) return false; 
    return true;
}

bool treeReader::passTTZCleanSelection(const int njets, const int nbjets, const double dMZ) const{
    if(leptonSelection == 3){
        if(njets < 3) return false;
        if(nbjets < 1) return false;
        if(dMZ > 10) return false;
    }
    else if(leptonSelection == 4){
        if(njets < 2) return false;
        if(nbjets < 1) return false;
        if(dMZ > 20) return false;
    }
    return true;
}

bool treeReader::passWZCRSelection(const int nbjets, const double dMZ) const{
    if(nbjets != 0) return false;
    if(dMZ > 10) return false; 
    return true;
}


bool treeReader::passttbarCRSelection(const int nbjets, const double dMZ, const double mlll) const{
    //if(!(deltaMZ == 999999 || !((deltaMZ < 10) || (mlll < 105) || (nBLoc < 1)))) continue;
    if(dMZ < 10) return false;
    if(mlll < 101) return false;
    else if(dMZ > 10 && dMZ != 999999){
        if(nbjets == 0) return false;
    } 
    return true;
}

bool treeReader::passttbarCRintZqSelection(const int njets, const int nbjets, const double dMZ) const{
    if(dMZ != 999999) return false;
    if(nbjets != 1) return false;
    if(njets != 2 && njets != 3) return false;
    return true;
}

bool treeReader::passZGCRSelection(const double mlll, const double dMZ) const{
    if(dMZ < 10) return false;
    else if (mlll < 81 || mlll > 101) return false; 
    return true;
}

bool treeReader::passDYCRSelection(const double dMZ, const double ptNonZ, const unsigned third, const double met, const double metPhi, const int njets, const int nbjets) const{

    if(dMZ > 10) return false;
    if(met > 30) return false;
    if(njets > 1) return false;
    if(nbjets > 0) return false;

    TLorentzVector l0p4;
    l0p4.SetPtEtaPhiE(ptNonZ, _lEta[third], _lPhi[third], _lE[third] * ptNonZ / _lPt[third]);

    double mTthrid = mtCalc(l0p4, met, metPhi);

    if(mTthrid > 30) return false;
    return true;
}

bool treeReader::passZZCRSelection(const std::vector<unsigned>& ind, std::vector<unsigned> indOf2LonZ, const int & njets){

    //if(njets < 1) return false;
    double mll, ptZ, ptNonZ, mlll; 
    unsigned third = -9999;
    TLorentzVector Zboson, lnegative;
    std::vector<unsigned> vectorRemoved;
    vectorRemoved = ind;
    vectorRemoved.erase(std::remove(vectorRemoved.begin(), vectorRemoved.end(), indOf2LonZ.at(0)), vectorRemoved.end());
    vectorRemoved.erase(std::remove(vectorRemoved.begin(), vectorRemoved.end(), indOf2LonZ.at(1)), vectorRemoved.end());
    double seconddMZ = deltaMZ(vectorRemoved, third, mll, ptZ, ptNonZ, mlll, indOf2LonZ, Zboson, lnegative);
    if(seconddMZ > 20) return false;
    return true;
}

bool treeReader::passTTZ4LSelection(const std::vector<unsigned>& ind, std::vector<unsigned> indOf2LonZ, const int njets){

    if(njets < 2) return false;
    double mll, ptZ, ptNonZ, mlll; 
    unsigned third = -9999;
    TLorentzVector Zboson, lnegative;
    std::vector<unsigned> vectorRemoved;
    vectorRemoved = ind;
    vectorRemoved.erase(std::remove(vectorRemoved.begin(), vectorRemoved.end(), indOf2LonZ.at(0)), vectorRemoved.end());
    vectorRemoved.erase(std::remove(vectorRemoved.begin(), vectorRemoved.end(), indOf2LonZ.at(1)), vectorRemoved.end());
    double seconddMZ = deltaMZ(vectorRemoved, third, mll, ptZ, ptNonZ, mlll, indOf2LonZ, Zboson, lnegative);
    //if(debug) cout << "delta MZ of second pair is " << seconddMZ << endl;
    //if(debug && seconddMZ != 999999.) cout << "indexes of second pair leptons are " << indOf2LonZ.at(0) << " " << indOf2LonZ.at(1) << endl;
    if(seconddMZ < 20) return false; // remove second onZ pair
    return true;
}

double treeReader::SRIDTTZ(const std::vector<unsigned>& ind, std::vector<unsigned> indOf2LonZ, const int & njets, const int & nbjets, const double & dMZ, const double & mlll){
    
    if(leptonSelection == 3){
        if(passWZCRSelection(nbjets, dMZ)){
            if(njets == 0) return -999.;
            return njets < 4 ? njets-1 : 3; // SR 0, 1, 2, 3
        }
        else if(passTTZSelection(njets, dMZ)){
            if(nbjets < 1)  return -999.;
            int nbjetsInd = nbjets < 2 ? 0 : 1;
            int njetsInd  = njets < 5 ? njets-2 : 3;
            return 4 + nbjetsInd * 4 + njetsInd; // 4, 5, 6, 7, 8, 9, 10, 11
        }
    }
    else if(leptonSelection == 4){
        if(passTTZ4LSelection(ind, indOf2LonZ, njets)){
            int nbjetsInd = nbjets < 1 ? 0 : 1;
            return 12 + nbjetsInd; // 12, 13, 14 // nj1nb1, nj2nb0, nj2nb1
        }
    }
    return -999;

}

double treeReader::SRIDTTZ1L(const std::vector<unsigned>& ind, std::vector<unsigned> indOf2LonZ, const int & njets, const int & nbjets, const double & dMZ, const double & mlll){
    
    if(leptonSelection == 3){
        if(nbjets == 0){
            if(njets == 0) return -999.;
            return njets < 4 ? njets-1 : 3; // SR 0, 1, 2, 3
        }
        else if(nbjets != 0){
            if(nbjets < 1)  return -999.;
            int nbjetsInd = nbjets < 2 ? 0 : 1;
            int njetsInd  = njets < 5 ? njets-2 : 3;
            return 4 + nbjetsInd * 4 + njetsInd; // 4, 5, 6, 7, 8, 9, 10, 11
        }
    }
    return -999;

}

double treeReader::SRID8SR3L(int & njets, int & nbjets, const double & dMZ) {

    if(leptonSelection == 3){
        if(passTTZSelection(njets, dMZ)){
            if(nbjets < 1)  return -999.;
            int nbjetsInd = nbjets < 2 ? 0 : 1;
            int njetsInd = njets < 5 ? njets-2 : 3;
            return nbjetsInd * 4 + njetsInd; // 0, 1, 2, 3, 4, 5, 6, 7
        }
    }
    return -999;
}

double treeReader::SRID3L(int & njets, int & nbjets, const double & dMZ) {
    if(leptonSelection == 3){
        if(passWZCRSelection(nbjets, dMZ)){
            if(njets == 0) return -999.;
            return njets < 4 ? njets-1 : 3; // SR 0, 1, 2, 3
        }
        else if(passTTZSelection(njets, dMZ)){
            if(nbjets < 1)  return -999.;
            int nbjetsInd = nbjets < 2 ? 0 : 1;
            int njetsInd = njets < 5 ? njets-2 : 3;
            return 4 + nbjetsInd * 4 + njetsInd; // 4, 5, 6, 7, 8, 9, 10, 11
        }
    }
    return -999;
}

double treeReader::SRIDPTZ(const double & ptZ) const{
    if(ptZ < 75) return 0.;
    else if (ptZ < 150) return 1.;
    else if (ptZ < 250) return 2.;
    else return 3.;
}

double treeReader::SRIDCosTheta(const double & cosTheta) const{
    if(cosTheta < -0.5) return 0.;
    else if (cosTheta < 0.) return 1.;
    else if (cosTheta < 0.5) return 2.;
    else return 3.;
}

double treeReader::SRID4L(int & njets, int & nbjets){
    if(nbjets == 0) return 0.;
    else if(nbjets > 0) return 1.;
}

double treeReader::SRIDWZCR(const int & njets, const int & nbjets, const double & dMZ){
    
    if(leptonSelection == 3){
        if(passWZCRSelection(nbjets, dMZ)){
            if(njets == 0) return -999.;
            return njets < 4 ? njets-1 : 3; // SR 0, 1, 2, 3
        }
    }
    return -999;
}

double treeReader::SRIDZZCR(const std::vector<unsigned>& ind, std::vector<unsigned> indOf2LonZ, const int & njets, const int & nbjets){
    
    if(leptonSelection == 4){
        if(passZZCRSelection(ind, indOf2LonZ, njets)){
            if(njets < 1) return -999.;
            int nbjetsInd = nbjets < 1 ? 0 : 1;
            int njetsInd = njets < 2 ? 0 : 1;
            return nbjetsInd * 2 + njetsInd; // 4, 5, 6, 7
        }
    }
    return -999;

}

double treeReader::SRIDTTCR(const int & njets, const int & nbjets, const double & dMZ, const double & mlll){
    if(leptonSelection == 3){
        if(passttbarCRSelection(nbjets, dMZ, mlll)){
            //if(njets < 2) return -999.;
            int nbjetsInd = nbjets < 2 ? nbjets : 2;
            int njetsInd = njets < 3 ? 0 : (njets == 3 ? 1 : 2);
            return nbjetsInd * 3 + njetsInd; // 8, 9, 10, 11, 12, 13, 14, 15, 16
        }
    }
    return -999;

}


bool treeReader::pass2Lpreselection(const int njets, const int nbjets, const std::vector<unsigned>& ind, const double met, const int nEle){
    if(njets < 2) return false;
    if(nbjets < 1) return false;

    TLorentzVector l0p4, l1p4;
    double lepCorrLead = _leptonMvatZqTTV[ind.at(0)] > leptonMVAcutInAnalysis[2] ? _lPt[ind.at(0)] : magicFactorInAnalysis[2] * _lPt[ind.at(0)] / _ptRatio[ind.at(0)];
    double lepCorrSubLead = _leptonMvatZqTTV[ind.at(1)] > leptonMVAcutInAnalysis[2] ? _lPt[ind.at(1)] : magicFactorInAnalysis[2] * _lPt[ind.at(1)] / _ptRatio[ind.at(1)];
    l0p4.SetPtEtaPhiE(lepCorrLead, _lEta[ind.at(0)], _lPhi[ind.at(0)], _lE[ind.at(0)] * lepCorrLead / _lPt[ind.at(0)]);
    l1p4.SetPtEtaPhiE(lepCorrSubLead, _lEta[ind.at(1)], _lPhi[ind.at(1)], _lE[ind.at(1)] * lepCorrSubLead / _lPt[ind.at(1)]);

    double ele_mll = (l0p4+l1p4).M();

    if(ele_mll < 12) return false;
    if(ele_mll > 76 && ele_mll < 106 && nEle == 2) return false;
    if(met < 30) return false;
    return true;
}

bool treeReader::pass2Lcleanpreselection(const int njets, const int nbjets, const std::vector<unsigned>& ind, const double met, const int nEle){
    if(njets < 3) return false;
    if(nbjets < 2) return false;

    TLorentzVector l0p4, l1p4;
    double lepCorrLead = _leptonMvatZqTTV[ind.at(0)] > leptonMVAcutInAnalysis[2] ? _lPt[ind.at(0)] : magicFactorInAnalysis[2] * _lPt[ind.at(0)] / _ptRatio[ind.at(0)];
    double lepCorrSubLead = _leptonMvatZqTTV[ind.at(1)] > leptonMVAcutInAnalysis[2] ? _lPt[ind.at(1)] : magicFactorInAnalysis[2] * _lPt[ind.at(1)] / _ptRatio[ind.at(1)];
    l0p4.SetPtEtaPhiE(lepCorrLead, _lEta[ind.at(0)], _lPhi[ind.at(0)], _lE[ind.at(0)] * lepCorrLead / _lPt[ind.at(0)]);
    l1p4.SetPtEtaPhiE(lepCorrSubLead, _lEta[ind.at(1)], _lPhi[ind.at(1)], _lE[ind.at(1)] * lepCorrSubLead / _lPt[ind.at(1)]);

    double ele_mll = (l0p4+l1p4).M();

    if(ele_mll < 12) return false;
    if(ele_mll > 76 && ele_mll < 106 && nEle == 2) return false;
    if(met < 30) return false;
    return true;
}

double treeReader::mtCalc(const TLorentzVector Vect, const double MET, const double MET_Phi) const{

    double MT=sqrt(2* Vect.Pt() * MET * ( 1 - (TMath::Cos(Vect.Phi() - MET_Phi )) ) );
    return MT;
}

double treeReader::sumAllLeptonsCharge(const std::vector<unsigned>& ind){
    double sum = 0.;
    for(auto & i : ind)
        sum += _lCharge[i];
    return sum;
}

double treeReader::largestAmongAll(const std::vector<double> & weights){

    double max_value = weights.front() ;
    for( std::size_t i = 1 ; i < weights.size() ; ++i ) 
        if( weights[i] > max_value ) max_value = weights[i] ;
    return max_value;

}

double treeReader::smallestAmongAll(const std::vector<double> & weights){

    double min_value = weights.front() ;
    for( std::size_t i = 1 ; i < weights.size() ; ++i ) 
        if( weights[i] < min_value ) min_value = weights[i] ;
    return min_value;

}

void treeReader::addVariablesToBDT(){

    readerLeptonMVAele->AddVariable( "pt", &user_pt);
    readerLeptonMVAele->AddVariable( "eta", &user_eta);
    readerLeptonMVAele->AddVariable( "trackMultClosestJet", &user_trackMult);
    readerLeptonMVAele->AddVariable( "miniIsoCharged", &user_miniIsoCharged);
    readerLeptonMVAele->AddVariable( "miniIsoNeutral", &user_miniIsoNeutral);
    readerLeptonMVAele->AddVariable( "pTRel", &user_ptrel);
    readerLeptonMVAele->AddVariable( "ptRatio", &user_ptratio);
    readerLeptonMVAele->AddVariable( "relIso", &user_relIso);
    readerLeptonMVAele->AddVariable( "deepCsvClosestJet", &user_jetBtagCSV);
    readerLeptonMVAele->AddVariable( "sip3d", &user_sip3d);
    readerLeptonMVAele->AddVariable( "dxy", &user_dxy);
    readerLeptonMVAele->AddVariable( "dz", &user_dz);
    if(currentSample.is2017())
        readerLeptonMVAele->AddVariable( "electronMvaFall17NoIso", &user_eleMVA);
    else
        readerLeptonMVAele->AddVariable( "electronMvaSpring16GP", &user_eleMVA);

    readerLeptonMVAmu->AddVariable( "pt", &user_pt);
    readerLeptonMVAmu->AddVariable( "eta", &user_eta);
    readerLeptonMVAmu->AddVariable( "trackMultClosestJet", &user_trackMult);
    readerLeptonMVAmu->AddVariable( "miniIsoCharged", &user_miniIsoCharged);
    readerLeptonMVAmu->AddVariable( "miniIsoNeutral", &user_miniIsoNeutral);
    readerLeptonMVAmu->AddVariable( "pTRel", &user_ptrel);
    readerLeptonMVAmu->AddVariable( "ptRatio", &user_ptratio);
    readerLeptonMVAmu->AddVariable( "relIso", &user_relIso);
    readerLeptonMVAmu->AddVariable( "deepCsvClosestJet", &user_jetBtagCSV);
    readerLeptonMVAmu->AddVariable( "sip3d", &user_sip3d);
    readerLeptonMVAmu->AddVariable( "dxy", &user_dxy);
    readerLeptonMVAmu->AddVariable( "dz", &user_dz);
    readerLeptonMVAmu->AddVariable( "segmentCompatibility", &user_segmComp);

    TString methodName = TString("BDTG method");
    TString weightfile;
    if(currentSample.is2017())
        weightfile = "/user/ikhvastu/Work/ttVcsMeasure/checkLeptonMVAvar/MVAtrainings/2017MC/WillemFiles/BDTG/TMVAClassification_BDTG_200Cuts_Depth4_baggedGrad_1000trees_shrinkage0p1_electron_2017.weights.xml";
    else
        weightfile = "/user/ikhvastu/Work/ttVcsMeasure/checkLeptonMVAvar/MVAtrainings/2016MC/WillemFiles/BDTG/toBeTested/TMVAClassification_BDTG_200Cuts_Depth4_baggedGrad_1000trees_shrinkage0p1_electron.weights.xml";
    readerLeptonMVAele->BookMVA( methodName, weightfile );

    methodName = TString("BDTG method");
    if(currentSample.is2017())
        weightfile = "/user/ikhvastu/Work/ttVcsMeasure/checkLeptonMVAvar/MVAtrainings/2017MC/WillemFiles/BDTG/TMVAClassification_BDTG_200Cuts_Depth4_baggedGrad_1000trees_shrinkage0p1_muon_2017.weights.xml";
    else
        weightfile = "/user/ikhvastu/Work/ttVcsMeasure/checkLeptonMVAvar/MVAtrainings/2016MC/WillemFiles/BDTG/toBeTested/TMVAClassification_BDTG_200Cuts_Depth4_baggedGrad_1000trees_shrinkage0p1_muon.weights.xml";
    readerLeptonMVAmu->BookMVA( methodName, weightfile );

}

void treeReader::fillBDTvariables(std::vector<Float_t> & varForBDT, int flavor){

    user_pt =  varForBDT.at(0);
    user_eta = varForBDT.at(1);
    user_trackMult = varForBDT.at(2);
    user_miniIsoCharged = varForBDT.at(3);
    user_miniIsoNeutral = varForBDT.at(4);
    user_ptrel = varForBDT.at(5);
    user_ptratio = varForBDT.at(6);
    user_relIso = varForBDT.at(7);
    user_jetBtagCSV = varForBDT.at(8);
    user_sip3d = varForBDT.at(9);
    user_dxy = varForBDT.at(10);
    user_dz = varForBDT.at(11);
    if(flavor == 0)
       user_eleMVA = varForBDT.at(12);
    if(flavor == 1)
       user_segmComp = varForBDT.at(12);
}


TLorentzVector treeReader::findBestNeutrinoAndTop(const TLorentzVector& wLep, const TLorentzVector& met, std::vector<unsigned>& taggedJetI, const std::vector<unsigned>& jetI, const std::vector<unsigned>& bJetI, const TLorentzVector* jetV){
    static const double mTop = 173.1;
    //set up tagged jets vector
    taggedJetI = {99, 99};
    //Check if there are any jets
    if(jetI.size() == 0) return met;    //dont attempt to find neutrino Z since we cant fix the ambiguity here
    std::pair<double, double> pzSol = neutrinoPZ(wLep, met);
    TLorentzVector neutrinoPlus(met.Px(), met.Py(), pzSol.first, sqrt(met.Px()*met.Px() + met.Py()*met.Py() + pzSol.first*pzSol.first) );
    TLorentzVector neutrinoMin(met.Px(), met.Py(), pzSol.second, sqrt(met.Px()*met.Px() + met.Py()*met.Py() + pzSol.second*pzSol.second) );
    bool plus = true;  //true if + solution is chosen, otherwise false
    if(bJetI.size() == 0){
        taggedJetI[0] = jetI[0];
    } else if(bJetI.size() == 1){
        taggedJetI[0] = bJetI[0];
    } else{
        double minDiff = 999999.;
        for(unsigned b = 0; b < bJetI.size(); ++b){
            double diffPlus = fabs( (neutrinoPlus + wLep + jetV[bJetI[b]] ).M() - mTop);
            double diffMin = fabs( (neutrinoMin + wLep + jetV[bJetI[b]] ).M() - mTop);
            if(diffPlus < minDiff && diffPlus < diffMin){
                plus = true;
                minDiff = diffPlus;
                taggedJetI[0] = bJetI[b];
            } else if(diffMin < minDiff){
                plus = false;
                minDiff = diffMin;
                taggedJetI[0] = bJetI[b];
            }
        }
    }
    if(bJetI.size() < 2){
        if(fabs( (neutrinoPlus + wLep + jetV[taggedJetI[0]] ).M() - mTop)  < fabs( (neutrinoMin + wLep + jetV[taggedJetI[0]] ).M() - mTop) ){
            plus = true;
        } else{
            plus = false;
        }
    }
    //find recoiling jet ( := leading jet not considered the b-jet from top)
    for(unsigned j = 0; j < jetI.size(); ++j){
        if(jetI[j] != taggedJetI[0]){
            taggedJetI[1] = jetI[j];
            break;                              //jet collection is assumed oredered
        }
    }
    if(plus) return neutrinoPlus;//( neutrinoPlus + wLep + jetV[taggedJetI[0]] ).M();
    else     return neutrinoMin;//( neutrinoMin  + wLep + jetV[taggedJetI[0]] ).M();
}

std::pair<double, double> treeReader::neutrinoPZ(const TLorentzVector& wLep, const TLorentzVector& met){
    static const double mW = 80.385;
    double mSquared = 0.5*mW*mW + wLep.Px()*met.Px() + wLep.Py()*met.Py();
    double preFac = mSquared/(wLep.Pt() * wLep.Pt());
    double term2 = wLep.P()*sqrt( std::max(0., 1 - met.Pt()*met.Pt()*wLep.Pt()*wLep.Pt()/(mSquared*mSquared) ) );
    return {preFac*(wLep.Pz() + term2), preFac*(wLep.Pz() - term2)};
}
