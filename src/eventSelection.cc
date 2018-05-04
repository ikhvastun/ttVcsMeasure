#include "THStack.h"

#include "../interface/treeReader.h"

using namespace std;

bool treeReader::lepIsGood(const unsigned l){
    // what is used in leptonMVA analysis
    if(!lepIsFOGood(l)) return false;
    if(_leptonMvatZqTTV[l] < leptonMVAcut) return false;

    if(leptonSelection == 2){
        if(_lFlavor[l] == 1 && ((_lMuonTrackPtErr[l]/_lMuonTrackPt[l]) > 0.2)) return false;
    }

    return true;
}

bool treeReader::lepIsFOGood(const unsigned l){

    if(_lPt[l] < 10) return false;

    if(!lepIsLoose(l)) return false;
    if(_lFlavor[l] == 0 && !eleIsClean(l)) return false;
    if(_lFlavor[l] == 1 && !_lPOGMedium[l]) return false;

    if(leptonSelection == 2){
        if(_lFlavor[l] == 0 && !_lElectronChargeConst[l]) return false;
    }

    if(_closestJetDeepCsv_bb[l] + _closestJetDeepCsv_b[l] > (is2017 ? 0.8001 : 0.8958)) return false;

    if(_leptonMvatZqTTV[l] < leptonMVAcut){
        if(_ptRatio[l] < (is2017 ? 0.3 : 0.4)) return false;
        if(_closestJetDeepCsv_bb[l] + _closestJetDeepCsv_b[l] > (is2017 ? 0.2 : 0.4)) return false;
        double electronMVAvalue = is2017 ? _lElectronMvaFall17NoIso[l] : _lElectronMva[l];
        if(_lFlavor[l] == 0 && electronMVAvalue < (is2017 ? 0.0 : -0.1) + (fabs(_lEta[l]) >= 1.479)*(is2017 ? 0.3 : 0.8)) return false;
        //if(_lFlavor[l] == 0 && (leptonSelection == 2 ? (_lElectronMvaFall17NoIso[l] < 0.2 + (fabs(_lEta[l]) >= 1.479)*0.5) : (_lElectronMvaFall17NoIso[l] < 0.0 + (fabs(_lEta[l]) >= 1.479)*0.7))) return false;
    }

    return true;
}

bool treeReader::eventChargeConsistent(const std::vector<unsigned>& ind){
    for(auto & l : ind){
        if(_lFlavor[l] == 0 && !_lElectronChargeConst[l]) return false;
        if(_lFlavor[l] == 1 && ((_lMuonTrackPtErr[l]/_lMuonTrackPt[l]) > 0.2)) return false;
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
    std::sort(ptMap.begin(), ptMap.end(), [](std::pair<double, unsigned>& p1, std::pair<double, unsigned>& p2){return p1.first > p2.first;} );
    for(unsigned l = 0; l < lCount; ++l){
        ind.push_back(ptMap[l].second);
    }
    return lCount;	
}

unsigned treeReader::selectFOLep(std::vector<unsigned>& ind){
    unsigned lCount = 0;
    std::vector<std::pair<double, unsigned>> ptMap;
    for(unsigned l = 0; l < _nLight; ++l){
        if(lepIsFOGood(l)){
            ++lCount;
            ptMap.push_back({_lPt[l], l});
        }
    }
    if(lCount < 2) return 0;
    std::sort(ptMap.begin(), ptMap.end(), [](std::pair<double, unsigned>& p1, std::pair<double, unsigned>& p2){return p1.first > p2.first;} );
    for(unsigned l = 0; l < lCount; ++l){
        ind.push_back(ptMap[l].second);
    }
    return lCount;  
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
   
    if(ptMap[0].first < 10) return false;
    if(ptMap[1].first < 10) return false;
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

    if(ptMap[0].first < 10) return false;
    if(ptMap[1].first < 10) return false;

    return true;
}

bool treeReader::jetIsClean(const unsigned ind, bool nonpromptSample){
    TLorentzVector jet;	
    jet.SetPtEtaPhiE(_jetPt[ind], _jetEta[ind], _jetPhi[ind], _jetE[ind]);
    for(unsigned l = 0; l < _nLight; ++l){
        //if(lepIsFOGood(l)){
        if(lepIsLoose(l)){
            TLorentzVector lep;
            lep.SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
            if(lep.DeltaR(jet) < 0.4) return false;
        }
    }
    return true;
}

bool treeReader::jetIsGood(const unsigned ind, const unsigned ptCut, const unsigned unc, const bool clean, bool nonpromptSample){
    if(fabs(_jetEta[ind]) > 2.4) return false;
    switch(unc){
        
        case 0: if(_jetPt[ind] < ptCut) return false; break;
        case 1: if(_jetPt_JECDown[ind] < ptCut) return false; break;
        case 2: if(_jetPt_JECUp[ind] < ptCut) return false; break;
        case 3: if(_jetPt_JERDown[ind] < ptCut) return false; break;
        case 4: if(_jetPt_JERUp[ind] < ptCut) return false; break;
        default: ;
    }
    return !clean || jetIsClean(ind, nonpromptSample);
}

unsigned treeReader::nJets(const unsigned unc, const bool clean, std::vector<unsigned>& ind, bool nonpromptSample){
    unsigned nJets = 0;
    for(unsigned j = 0; j < _nJets; ++j){
        if(jetIsGood(j, 30, unc, clean, nonpromptSample)) {
            ++nJets;
            ind.push_back(j);
        }
    }
    return nJets;
}

bool treeReader::bTaggedDeepCSV(const unsigned ind, const unsigned wp){
    //static const double bTagWP[3] = {0.2219, 0.6324,  0.8958};
    static const double bTagWP[3] = {0.1522, 0.4941,  0.8001};
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
        double ptcor1 = _leptonMvatZqTTV[ind.at(l0)] > leptonMVAcut ? _lPt[ind.at(l0)] : magicFactor * _lPt[ind.at(l0)] / _ptRatio[ind.at(l0)];
        
        l0p4.SetPtEtaPhiE(ptcor1 ,_lEta[ind.at(l0)],_lPhi[ind.at(l0)],_lE[ind.at(l0)] * ptcor1 / _lPt[ind.at(l0)] );          
        for(unsigned l1 = l0; l1 < ind.size(); l1++){
            double ptcor2 = _leptonMvatZqTTV[ind.at(l1)] > leptonMVAcut ? _lPt[ind.at(l1)] : magicFactor * _lPt[ind.at(l1)] / _ptRatio[ind.at(l1)];
            
            if(ind.at(l0) == ind.at(l1)) continue;

            if (_lCharge[ind.at(l0)] != _lCharge[ind.at(l1)]) {
                l1p4.SetPtEtaPhiE(ptcor2 ,_lEta[ind.at(l1)],_lPhi[ind.at(l1)],_lE[ind.at(l1)] * ptcor2 / _lPt[ind.at(l1)]);
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
                                //double ptcor3 = ptFake(_lPt[lepThird], _ptRatio[lepThird], _lFlavor[lepThird], _leptonMvaTTH[lepThird], _lPOGMedium[lepThird]);
                                double ptcor3 = _leptonMvatZqTTV[lepThird] > leptonMVAcut ? _lPt[lepThird] : magicFactor * _lPt[lepThird] / _ptRatio[lepThird];
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
        bool leptIsPr = leptonIsPrompt(l);
        allPrompt = allPrompt && leptIsPr; 
    }
    
    return allPrompt;
}

bool treeReader::noConversionInSelection(const std::vector<unsigned>& ind){

    bool allGood = true;

    for(auto & l : ind){
        allGood = allGood && !(_lMatchPdgId[l] == 22);
    }

    return allGood;
}

// test function that checks the gen lepton in cone 0.1 around reco lepton, not used anymore
bool treeReader::leptonIsPrompt(const unsigned & l){

    bool leptIsP = false;
    TLorentzVector l0reco;
    l0reco.SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
    for(unsigned i = 0; i < _gen_nL; i++){
        TLorentzVector l0gen;
        l0gen.SetPtEtaPhiE(_gen_lPt[i], _gen_lEta[i], _gen_lPhi[i], _gen_lE[i]);
        if(l0reco.DeltaR(l0gen) < 0.1){
            if(_gen_lIsPrompt[i] && (_lFlavor[l] == _gen_lFlavor[i] || _gen_lFlavor[i] == 2)){
                leptIsP = true;
            }
        }
    }
    return leptIsP;
}

bool treeReader::leptonIsFromPromptTau(const unsigned & l){

    bool leptIsP = false;
    TLorentzVector l0reco;
    l0reco.SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
    for(unsigned i = 0; i < _gen_nL; i++){
        TLorentzVector l0gen;
        l0gen.SetPtEtaPhiE(_gen_lPt[i], _gen_lEta[i], _gen_lPhi[i], _gen_lE[i]);
        if(l0reco.DeltaR(l0gen) < 0.1){
            if(_gen_lIsPrompt[i] && _gen_lFlavor[i] == 2){
                leptIsP = true;
            }
        }
    }
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
    if(name == "tight") return kBlack; 
    if(name == "loose") return 51; 
    if(name == "tW") return kRed-10;

    return kBlack;
}
