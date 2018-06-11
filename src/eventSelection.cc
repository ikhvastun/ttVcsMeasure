#include "THStack.h"

#include "../interface/treeReader.h"

using namespace std;

bool treeReader::lepIsGood(const unsigned l){
    // what is used in leptonMVA analysis
    //cout << "lepton info: " << _lPt[l] << " " << _lEta[l] << " " << _lFlavor[l] << " " << _leptonMvatZqTTV[l] << " " << _lPOGMedium[l] << " " << _lCharge[l] << endl;
    if(!lepIsFOGood(l)) return false;
    if(_leptonMvatZqTTV[l] < leptonMVAcut) return false;

    if(leptonSelection == 2){
        if(_lFlavor[l] == 1 && ((_lMuonTrackPtErr[l]/_lMuonTrackPt[l]) > 0.2)) return false;
        if(_lFlavor[l] == 0 && !_lElectronChargeConst[l]) return false;
        if(_lFlavor[l] == 0 && !_lElectronPassConvVeto[l]) return false;
        if(_lFlavor[l] == 0 && _lElectronMissingHits[l] != 0) return false;
    }

    return true;
}

bool treeReader::lepIsFOGood(const unsigned l){

    if(_lPt[l] < 10) return false;

    if(!lepIsLoose(l)) return false;
    if(_lFlavor[l] == 0 && !eleIsClean(l)) return false;
    if(_lFlavor[l] == 1 && !_lPOGMedium[l]) return false;

    if(_closestJetDeepCsv_bb[l] + _closestJetDeepCsv_b[l] > (is2017 ? 0.8001 : 0.8958)) return false;

    if(_leptonMvatZqTTV[l] < leptonMVAcut){
        if(_ptRatio[l] < (is2017 ? (leptonSelection == 3 ? 0.3 : 0.4) : (leptonSelection == 3 ? 0.4 : 0.5))) return false;
        if(_closestJetDeepCsv_bb[l] + _closestJetDeepCsv_b[l] > (is2017 ? (leptonSelection == 3 ? 0.2 : 0.4) : (leptonSelection == 3 ? 0.4 : 0.3))) return false;
        double electronMVAvalue = is2017 ? _lElectronMvaFall17NoIso[l] : _lElectronMva[l];
        if(_lFlavor[l] == 0 && electronMVAvalue < (is2017 ? (leptonSelection == 3 ? 0.0 : 0.4) : (leptonSelection == 3 ? -0.1 : 0.3))    +    (fabs(_lEta[l]) >= 1.479)*(is2017 ? (leptonSelection == 3 ? 0.3 : 0.4) : (leptonSelection == 3 ? 0.8 : 0.5))) return false;
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
    //cout << "info about loose lepton with pt " << _lPt[ind] << ": " << (fabs(_dxy[ind]) >= 0.05) << " " << (fabs(_dz[ind]) >= 0.1) << " " << (_3dIPSig[ind] >= 8) << " " << (_miniIso[ind] >= 0.4) << " " << (_lElectronMissingHits[ind] > 1) << " " << (!_lElectronPassEmu[ind]) << endl;
    //cout << "lepton miniiso: " << _miniIso[ind] << endl;
    if(_lPt[ind] <= 7 - 2*_lFlavor[ind]) return false;
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

unsigned treeReader::selectLep(std::vector<unsigned>& ind){
    unsigned lCount = 0;
    std::vector<std::pair<double, unsigned>> ptMap;
    for(unsigned l = 0; l < _nLight; ++l){
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

unsigned treeReader::selectFakeLep(std::vector<unsigned>& ind){
    unsigned lCount = 0;
    std::vector<std::pair<double, unsigned>> ptMap;
    for(unsigned l = 0; l < _nLight; ++l){
        if(lepIsFOGood(l)){
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
        double ptcor = _leptonMvatZqTTV[i] > leptonMVAcut ? _lPt[i] : magicFactor * _lPt[i] / _ptRatio[i];
        ptMap.push_back({ptcor, i});
        
    }
    std::sort(ptMap.begin(), ptMap.end(), [](std::pair<double, unsigned>& p1, std::pair<double, unsigned>& p2){return p1.first > p2.first;} );

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
        double ptcor = _leptonMvatZqTTV[i] > leptonMVAcut ? _lPt[i] : magicFactor * _lPt[i] / _ptRatio[i];
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

bool treeReader::jetIsClean(const unsigned ind){
    TLorentzVector jet;	
    jet.SetPtEtaPhiE(_jetPt[ind], _jetEta[ind], _jetPhi[ind], _jetE[ind]);
    for(unsigned l = 0; l < _nLight; ++l){
        if(lepIsFOGood(l)){ // cleaning with FO objects
            TLorentzVector lep;
            lep.SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
            //cout << "jet lepton cleaning is going on, delta R is: " << lep.DeltaR(jet) << endl;
            if(lep.DeltaR(jet) < 0.4) return false;
        }
    }
    return true;
}

bool treeReader::jetIsGood(const unsigned ind, const unsigned ptCut, const unsigned unc, const bool clean, bool is2017){
    if(fabs(_jetEta[ind]) > 2.4) return false;
    if(is2017 && !_jetIsTight[ind]) return false;
    switch(unc){
        
        case 0: if(_jetPt[ind] < ptCut) return false; break;
        case 1: if(_jetPt_JECDown[ind] < ptCut) return false; break;
        case 2: if(_jetPt_JECUp[ind] < ptCut) return false; break;
        case 3: if(_jetPt_JERDown[ind] < ptCut) return false; break;
        case 4: if(_jetPt_JERUp[ind] < ptCut) return false; break;
        default: ;
    }
    return !clean || jetIsClean(ind);
}

unsigned treeReader::nJets(const unsigned unc, const bool clean, std::vector<unsigned>& ind, bool is2017){
    unsigned nJets = 0;
    for(unsigned j = 0; j < _nJets; ++j){
        if(jetIsGood(j, 30, unc, clean, is2017)) {
            ++nJets;
            ind.push_back(j);
        }
    }
    return nJets;
}

unsigned treeReader::nJetsNotB(const unsigned unc, const bool clean, std::vector<unsigned>& ind, const unsigned btagWP, bool nonpromptSample){
    unsigned nJets = 0;
    for(unsigned j = 0; j < _nJets; ++j){
        if(jetIsGood(j, 30, unc, clean, nonpromptSample) && !bTaggedDeepCSV(j, btagWP)) {
            ++nJets;
            ind.push_back(j);
        }
    }
    return nJets;
}

bool treeReader::bTaggedDeepCSV(const unsigned ind, const unsigned wp){
    static const double bTagWP[3] = {(is2017 ? 0.1522 : 0.2219), (is2017 ? 0.4941 : 0.6324),  (is2017 ? 0.8001 : 0.8958)};
    return (_jetDeepCsv_b[ind] + _jetDeepCsv_bb[ind]) > bTagWP[wp];
}

bool treeReader::bTaggedCSVv2(const unsigned ind, const unsigned wp){
    static const double bTagWP[3] = {0.5426, 0.8484, 0.9535};
    return _jetCsvV2[ind] > bTagWP[wp];
}

unsigned treeReader::nBJets(const unsigned unc, const bool deepCSV, const bool clean, const unsigned wp, bool is2017){
    unsigned nbJets = 0;
    //cout << "nbjets calculation is starting here >>>>>>>>>>>>>>> " << endl;
    for(unsigned j = 0; j < _nJets; ++j){
        if(jetIsGood(j, 30, unc, clean, is2017)){
            
            if(deepCSV && bTaggedDeepCSV(j, wp)) ++nbJets;
            else if(!deepCSV && bTaggedCSVv2(j, wp)) ++nbJets;
            
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

double treeReader::deltaMZ(const std::vector<unsigned>& ind, unsigned & third, double & mll, double & ptZ, double & ptNonZ, double & mlll, std::vector<unsigned> & indLeptonOnZ){
                        
    TLorentzVector l0p4, l1p4, l2p4;

    double deltaMZ = 999999.;

    for (unsigned l0 = 0; l0 < ind.size(); l0++) {
        //double ptcor1 = _lPt[ind.at(l0)];
        double ptcor1 = leptonSelection != 4 ? (_leptonMvatZqTTV[ind.at(l0)] > leptonMVAcut ? _lPt[ind.at(l0)] : magicFactor * _lPt[ind.at(l0)] / _ptRatio[ind.at(l0)]) : _lPt[ind.at(l0)];
        
        l0p4.SetPtEtaPhiE(ptcor1 ,_lEta[ind.at(l0)],_lPhi[ind.at(l0)],_lE[ind.at(l0)] * ptcor1 / _lPt[ind.at(l0)] );          
        for(unsigned l1 = l0+1; l1 < ind.size(); l1++){

            //if(ind.at(l0) == ind.at(l1)) continue;
            //double ptcor2 = _lPt[ind.at(l1)];
            double ptcor2 = leptonSelection != 4 ? (_leptonMvatZqTTV[ind.at(l1)] > leptonMVAcut ? _lPt[ind.at(l1)] : magicFactor * _lPt[ind.at(l1)] / _ptRatio[ind.at(l1)]) : _lPt[ind.at(l1)];
            //cout << "l1 is " << ind.at(l1) << endl;

            if (_lCharge[ind.at(l0)] != _lCharge[ind.at(l1)]) {
                l1p4.SetPtEtaPhiE(ptcor2 ,_lEta[ind.at(l1)],_lPhi[ind.at(l1)],_lE[ind.at(l1)] * ptcor2 / _lPt[ind.at(l1)]);
                l1p4+=l0p4;
                double mdiL = l1p4.M();
                
                //cout << "invariant mass of 2 leptons is: " << mdiL << endl;
                if (_lFlavor[ind.at(l0)] == _lFlavor[ind.at(l1)] ) {

                    //cout << "the leptons are of the same flavour with inv. mass: " << mdiL << endl;
                    if (fabs(mdiL - 91.2) < deltaMZ) {
                        deltaMZ = fabs(mdiL - 91.2);
                        indLeptonOnZ.clear();
                        indLeptonOnZ.push_back(ind.at(l0));
                        indLeptonOnZ.push_back(ind.at(l1));
                        mll = mdiL;
                        ptZ = l1p4.Pt();
                        if(leptonSelection == 3 || leptonSelection == 4){
                            for(auto & lepThird : ind){
                                if(lepThird == ind.at(l0) || lepThird == ind.at(l1)) continue;
                                third = lepThird;
                                //double ptcor3 = _lPt[lepThird];
                                double ptcor3 = leptonSelection != 4 ? (_leptonMvatZqTTV[lepThird] > leptonMVAcut ? _lPt[lepThird] : magicFactor * _lPt[lepThird] / _ptRatio[lepThird]) : _lPt[lepThird];
                                l2p4.SetPtEtaPhiE(ptcor3 ,_lEta[lepThird],_lPhi[lepThird],_lE[lepThird] * ptcor3 / _lPt[lepThird]);
                                ptNonZ = ptcor3;
                                mlll = (l1p4+l2p4).M();
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

Color_t treeReader::assignColor(std::string & name){

    if(name == "data") return kBlack;
    if(name == "Nonprompt") return kBlue-9;
    if(name == "nonpromptData") return kBlue-9;
    if(name == "chargeMisID") return kMagenta-7;
    if(name == "ttZ") return 91;
    if(name == "ttW") return 98;
    if(name == "ttH") return kRed-10;
    if(name == "ttX") return kRed-10;
    if(name == "WZ") return 51;
    if(name == "ZZ") return kGreen+3;
    if(name == "rare") return 8;
    if(name == "Z#gamma") return kGreen;
    if(name == "DY") return kBlue-9;

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
