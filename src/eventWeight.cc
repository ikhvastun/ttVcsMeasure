//include c++ library classes
#include <iostream>

//include other parts of code
#include "../interface/treeReader.h"

//pu SF 
double treeReader::puWeight(const unsigned unc){
    if(debug) std::cout << "pile up true int " << _nTrueInt << "; pu weight: " << reweighter->puWeight(_nTrueInt, currentSample, unc) << std::endl;
    return reweighter->puWeight(_nTrueInt, currentSample, unc);
}

//b-tagging SF for given flavor
double treeReader::bTagWeight(const unsigned jetFlavor, const unsigned unc){
    //WARNING: reactivate this code once the b-tag efficiencies have been computed 
    double pMC = 1.;
    double pData = 1.;
    for(unsigned j = 0; j < _nJets; ++j){
        if(_jetHadronFlavor[j] == jetFlavor){
            //QUESTION: should JEC and b-tag weights also be varied up and down at the same time when computing systematics?
            if(jetIsGood(j, 30., 0, true, currentSample.is2017())){
                double sf = reweighter->bTagWeight(_jetHadronFlavor[j], _jetSmearedPt[j], _jetEta[j], _closestJetDeepCsv_b[j] + _closestJetDeepCsv_bb[j], unc);
                double eff = reweighter->bTagEff(_jetHadronFlavor[j], _jetSmearedPt[j], _jetEta[j]);
                //if(debug) std::cout << "jet with pt: " << _jetPt[j] << "; has btagEffcalc: " << eff << "; and centralValue: " << sf << "; csv value: " << _closestJetDeepCsv_b[j] + _closestJetDeepCsv_bb[j] << std::endl;
                if(bTaggedDeepCSV(j, 1)){
                    pMC *= eff;
                    pData *= eff*sf;
                } else {
                    pMC *= (1 - eff);
                    pData *= (1 - eff*sf);
                }
            }
        }
    }
    return pData/pMC;
}

//light flavor b-tagging SF
double treeReader::bTagWeight_udsg(const unsigned unc){
    return bTagWeight(0, unc);
}

//heavy flavor b-tagging SF
double treeReader::bTagWeight_c(const unsigned unc){
    return bTagWeight(4, unc);
}

//beauty flavor b-tagging SF
double treeReader::bTagWeight_b(const unsigned unc){
    return bTagWeight(5, unc);
}

//total b-tagging SF
double treeReader::bTagWeight(const unsigned unc){
    if(debug) std::cout << "btag light: " << bTagWeight_udsg(unc) << "; btag c: " << bTagWeight_c(unc) << "; btag b: " << bTagWeight_b(unc) << std::endl; 
    return bTagWeight_udsg(unc)*bTagWeight_c(unc)*bTagWeight_b(unc); 
}

//total lepton SF
// first only Syst uncertainty
double treeReader::leptonWeightOnlySyst(const unsigned unc, const bool ignoreMuon, const bool ignoreEle){
    double sf = 1.;
    for(unsigned l = 0; l < _nLight; ++l){
        if( lepIsGood(l, leptonSelection) ){
        //if( lepIsGoodFortZq(l) ){
            if( isMuon(l) && !ignoreMuon){
                sf *= reweighter->muonTightWeightOnlySyst(_lPt[l], _lEta[l], unc);
            } else if( isElectron(l) && !ignoreEle){
                sf *= reweighter->electronTightWeightOnlySyst(_lPt[l], _lEta[l], _lEtaSC[l], unc);
            }
            if(debug) std::cout << "sf after the lepton with pt " << _lPt[l] << " and flavour " << _lFlavor[l] << " is " << sf << std::endl;
        } else if( lepIsLoose(l) ){
            if( isMuon(l) && !ignoreMuon ){
                sf *= reweighter->muonLooseWeight(_lPt[l], _lEta[l], unc);
            } else if( isElectron(l) && !ignoreEle ){
                sf *= reweighter->electronLooseWeight(_lPt[l], _lEta[l], _lEtaSC[l], unc);
            }
        }  
    }
    return sf;
}

double treeReader::leptonWeightOnlyStat(const unsigned unc, const bool ignoreMuon, const bool ignoreEle){
    double sf = 1.;
    for(unsigned l = 0; l < _nLight; ++l){
        if( lepIsGood(l, leptonSelection) ){
        //if( lepIsGoodFortZq(l) ){
            if( isMuon(l) && !ignoreMuon){
                sf *= reweighter->muonTightWeightOnlyStat(_lPt[l], _lEta[l], unc);
            } else 
            if( isElectron(l) && !ignoreEle){
                sf *= reweighter->electronTightWeightOnlyStat(_lPt[l], _lEta[l], _lEtaSC[l], unc);
            }
            if(debug) std::cout << "sf after the lepton with pt " << _lPt[l] << " and flavour " << _lFlavor[l] << " is " << sf << std::endl;
        }else if( lepIsLoose(l) ){
            if( isMuon(l) && !ignoreMuon){
                sf *= reweighter->muonLooseWeight(_lPt[l], _lEta[l], unc);
            } else if( isElectron(l) && !ignoreEle){
                sf *= reweighter->electronLooseWeight(_lPt[l], _lEta[l], _lEtaSC[l], unc);
            }
        } 
    }
    return sf;
}


double treeReader::leptonWeightOnlyReco(const unsigned unc, const bool ignoreMuon, const bool ignoreEle){
    double sf = 1.;
    for(unsigned l = 0; l < _nLight; ++l){
        if( lepIsGood(l, leptonSelection) ){
        //if( lepIsGoodFortZq(l) ){
            /*
            if( isMuon(l) ){
                sf *= reweighter->muonTightWeight(_lPt[l], _lEta[l], unc);
            } else 
            */
            if( isElectron(l) && !ignoreEle){
                sf *= reweighter->electronRecoWeight(_lEtaSC[l], _lPt[l], unc);
            }
            if(debug) std::cout << "reco sf after the lepton with pt " << _lPt[l] << " and flavour " << _lFlavor[l] << " is " << sf << std::endl;
        } 
    }
    return sf;
}

// trigger SF
double treeReader::triggerWeight(){
    return reweighter->getTriggerSF(ptCorrV[0].first);
}

//check if scale-factors have to be initialized, and do so if needed
void treeReader::initializeWeights(){
    static bool weightsAre2016 = is2016();
    bool firstTime = ( reweighter.use_count() == 0 );
    bool changedEra = ( weightsAre2016 != is2016() );
    if( firstTime || changedEra){
        weightsAre2016 = is2016();
        reweighter.reset(new Reweighter(samples, is2016(), leptonSelection ) );
    } 
}
    
double treeReader::sfWeight(){
    initializeWeights();
    double sf = puWeight();
    sf *= bTagWeight();
    //sf *= leptonWeightOnlyStat() * leptonWeightOnlySyst() * leptonWeightOnlyReco();
    sf *= leptonWeightOnlySyst() * leptonWeightOnlyReco(); // consider only one because central value is the same for stat and syst
    sf *= fakeRateWeight();
    sf *= triggerWeight();
    if( sf == 0){
        std::cerr << "Error: event sf is zero! This has to be debugged! evNumber: " << _eventNb << "; sf(pu, btag, lep(Reco, Stat and Syst), fr, trigger): " << puWeight() << " " << bTagWeight() << " " << leptonWeightOnlyReco()  << " " << leptonWeightOnlyStat()  << " " << leptonWeightOnlySyst() << " " << fakeRateWeight() << " " << triggerWeight() << std::endl;
    } else if( std::isnan(sf) ){
        std::cerr << "Error: event sf is nan! This has to be debugged" << std::endl;
    }
    return sf;
}


//fake rate
double treeReader::fakeRateWeight(const unsigned unc){
    initializeWeights();
    double sf = -1.;
    for(unsigned l = 0; l < _nLight; ++l){
        if(lepIsFOGood(l, leptonSelection) && !lepIsGood(l, leptonSelection) ){
            double fr = 1.;
            if( isMuon(l) ){
                fr = reweighter->muonFakeRate(magicFactorInAnalysis[leptonSelection] * _lPt[l] / _ptRatio[l], _lEta[l], unc);
                //std::cout << "fr for muon is with pt and eta: " << magicFactorInAnalysis[leptonSelection] * _lPt[l] / _ptRatio[l] << " " << _lEta[l] << " " << fr << std::endl;
            } else if( isElectron(l) ){
                fr = reweighter->electronFakeRate(magicFactorInAnalysis[leptonSelection] * _lPt[l] / _ptRatio[l], _lEta[l], unc);
                //std::cout << "fr for electron is with pt and eta: " << magicFactorInAnalysis[leptonSelection] * _lPt[l] / _ptRatio[l] << " " << _lEta[l] << " " << fr << std::endl;
            }
            sf *= -fr/(1 - fr);
        }
    }
    if(!isData) sf *= -1;
    if(debug) std::cout << "fr from class is " << sf << std::endl;
    return sf;
}

//charge mis ID rate
double treeReader::CMIDRateWeight(const unsigned unc){
    initializeWeights();
    double fr = 0.;
    for(unsigned l = 0; l < _nLight; ++l){
        if(lepIsGood(l, leptonSelection) && _lFlavor[l] == 0 ){
            fr += reweighter->electronCMIDRate(_lPt[l], _lEta[l], unc);
        }
    }
    return fr;
}

//jet prefiring probability 
double treeReader::jetPrefiringWeight(){
    double sf = 1.;
    for(unsigned j = 0; j < _nJets; ++j){
        if( jetIsGood(j, 30., 0, true, currentSample.is2017()) ){
            double prefiringProbability = reweighter->jetPrefiringProbability( _jetPt[j], _jetEta[j] );
            sf *= ( 1. - prefiringProbability );
        }
    }
    return sf;
}
