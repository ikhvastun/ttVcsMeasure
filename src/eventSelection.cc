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

bool treeReader::lepIsGood(const unsigned l){
    // what is used in leptonMVA analysis

    if(!_lPOGTight[l]) return false;
    if(_lFlavor[l] == 1 && _relIso0p4Mu[l] > 0.25) return false;

    return true;
}


bool treeReader::lepIsFOGood(const unsigned l){
    // what is used in leptonMVA analysis

    if(!_lEwkLoose[l]) return false;
    if(_lFlavor[l] == 0 && !_lElectronPassEmu[l]) return false;
    if(_closestJetCsvV2[l] > 0.8484) return false;

    // this should be done with electronMVAHZZ (0, 0, 0.7) cuts
    if(_lFlavor[l] == 0){
        
        if(fabs(_lEta[l]) < 0.8 && _lElectronMvaHZZ[l] <  0)         return false;
        else if (fabs(_lEta[l]) < 1.479 && _lElectronMvaHZZ[l] <  0)   return false;
        else if (fabs(_lEta[l]) < 2.5 && _lElectronMvaHZZ[l] <  0.7) return false;

    }
    
        
    if(leptonSelection == 2){
        if(_lFlavor[l] == 0 && _lElectronMissingHits[l] != 0) return false;
    }

    
    if(_leptonMvaTTH[l] < 0.9){

        if(_lFlavor[l] == 1 && (_ptRatio[l] > 0.5 && _closestJetCsvV2[l] < 0.3 && _lMuonSegComp[l] > 0.3)) return true;
        else if (_lFlavor[l] == 0 && (_ptRatio[l] > 0.5 && _closestJetCsvV2[l] < 0.3)) return true;
        else return false;
    }
    
    
    return true;
}

bool treeReader::lepIsTight(const unsigned l){
    return _lEwkTight[l];
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

unsigned treeReader::selectFakeLep(std::vector<unsigned>& ind){
    //setConePt(); REMOVE CONE CORRECTION UNTIL MOVING TO FR
    unsigned lCount = 0;
    std::vector<std::pair<double, unsigned>> ptMap;
    for(unsigned l = 0; l < _nLight; ++l){
        //cout << "lepton info: " << _lPt[l] << " " << _lEwkLoose[l] << " " << _leptonMvaTTH[l] << endl;
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
        double ptcor = ptFake(_lPt[i], _ptRatio[i], _lFlavor[i], _leptonMvaTTH[i], _lPOGMedium[i]);
        ptMap.push_back({ptcor, i});
        
    }
    std::sort(ptMap.begin(), ptMap.end(), [](std::pair<double, unsigned>& p1, std::pair<double, unsigned>& p2){return p1.first > p2.first;} );

    //cout << "two leptons with pt: " << ptMap[0].first << " " << ptMap[1].first << endl;

    ptCorrV.clear();
    for(auto & i : ptMap)
    ptCorrV.push_back(i);

    if(ptMap[0].first < 25 && _lFlavor[ptMap[0].second] == 1) return false;
    if(ptMap[1].first < 25 && _lFlavor[ptMap[1].second] == 1) return false;

    if(ptMap[0].first < 40 && _lFlavor[ptMap[0].second] == 0) return false;
    if(ptMap[1].first < 27 && _lFlavor[ptMap[1].second] == 0) return false;
       
    return true;
}

bool treeReader::jetIsClean(const unsigned ind, bool nonpromptSample){
    TLorentzVector jet;	
    jet.SetPtEtaPhiE(_jetPt[ind], _jetEta[ind], _jetPhi[ind], _jetE[ind]);
    for(unsigned l = 0; l < _nLight; ++l){
        if((nonpromptSample && lepIsFOGood(l)) || (!nonpromptSample && lepIsGood(l))){ // cleaning with FO objects
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

double treeReader::deltaMZ(const std::vector<unsigned>& ind, unsigned & third, double & mll, double & ptZ, double & ptNonZ){
                        
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
    if(name == "Top") return 91; 

    return kBlack;
}




