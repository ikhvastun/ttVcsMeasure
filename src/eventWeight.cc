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
            if(jetIsGood(j, 30., 0, true, is2017) && fabs(_jetEta[j]) < 2.4){
                double sf = reweighter->bTagWeight(_jetHadronFlavor[j], _jetPt[j], _jetEta[j], _closestJetDeepCsv_b[j] + _closestJetDeepCsv_bb[j], unc);
                double eff = reweighter->bTagEff(_jetHadronFlavor[j], _jetPt[j], _jetEta[j]);
                //std::cout << "jet with pt: " << _jetPt[j] << "; has btagEffcalc: " << eff << "; and centralValue: " << sf << std::endl;
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
double treeReader::leptonWeight(const unsigned unc){
    double sf = 1.;
    for(unsigned l = 0; l < _nLight; ++l){
        if( lepIsGood(l, leptonSelection) ){
            if( isMuon(l) ){
                sf *= reweighter->muonTightWeight(_lPt[l], _lEta[l], unc);
            } else if( isElectron(l) ){
                sf *= reweighter->electronTightWeight(_lPt[l], _lEta[l], _lEtaSC[l], unc);
            }
        } 
        /*
        else if( lepIsLoose(l) ){
            if( isMuon(l) ){
                sf *= reweighter->muonLooseWeight(_lPt[l], _lEta[l], unc);
            } else if( isElectron(l) ){
                sf *= reweighter->electronLooseWeight(_lPt[l], _lEta[l], _lEtaSC[l], unc);
            }
        }
        */
        //std::cout << "sf after lepton with pt: " << _lPt[l] << " and flavour " << _lFlavor[l] << "is " << sf << std::endl;
    }
    if(debug) std::cout << "lepton SF is " << sf << std::endl;
    return sf;
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
    sf *= leptonWeight();
    sf *= fakeRateWeight();
    if( sf == 0){
        std::cerr << "Error: event sf is zero! This has to be debugged! evNumber: " << _eventNb << "; sf(pu, btag, lep): " << puWeight() << " " << bTagWeight() << " " << leptonWeight() << std::endl;
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
                fr = reweighter->muonFakeRate(_lPt[l], _lEta[l], unc);
            } else if( isElectron(l) ){
                fr = reweighter->electronFakeRate(_lPt[l], _lEta[l], unc);
            }
            sf *= -fr/(1 - fr);
        }
    }
    if(!isData) sf *= -1;
    if(debug) std::cout << "fr from class is " << sf << std::endl;
    return sf;
}
