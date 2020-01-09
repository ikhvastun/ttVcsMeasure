#!/bin/bash

#include bash functiosn to set up CMSSW
source setCMSSW.sh

cwd=$(pwd)                                          #current working directory needed to locate code 

skimSample(){                                           #function to skim one sample
    name="${1%/*}"                                      #remove everything before the last "/" in the path to the sample

    if [[ $1 = *"_realistic_v10"* ]]; then
        name="${name}_realistic_v10"
    elif [[ $1 = *"_realistic_v14"* ]]; then
        name="${name}_realistic_v14"
    elif [[ $1 = *"_realistic_v11"* ]]; then
        name="${name}_realistic_v11"
    fi

    if [[ $1 = *"_RECOPF"* ]]; then
        name="${name}_RECOPF"
    fi

    if [[ $1 = *"MiniAOD2016"* ]] || [[ $1 = *"Run2016"* ]]; then # added now
        name="${name}_MiniAOD2016"                                # added now
    elif [[ $1 = *"MiniAOD2017"* ]] || [[ $1 = *"Run2017"* ]]; then # added now
        name="${name}_MiniAOD2017"                                # added now
    elif [[ $1 = *"MiniAOD2018"* ]] || [[ $1 = *"Run2018"* ]]; then # added now
        name="${name}_MiniAOD2018"                                # added now
    elif [[ $1 = *"Fall17"* ]] || [[ $1 = *"Run2017"* ]]; then
        name="${name}_Fall17"
    else 
        name="${name}_Summer16"
    fi
    echo "$name"
    outputDir=~/Work/ntuples_temp_2L_ttZ_${name}
    if [ ! -d "$outputDir" ]; then                      #make output directory if it doesn't exist 
        mkdir -p $outputDir
    fi
    submit=~/skimJob.sh
    makeSubmit $submit $2                               #make temporary submission script

    count=0                                             #file counter
    echo "${1}"
    files=${1}/*/*/*root
    for f in $files
#        do if(( ($((count)) == 0) )); then
#            /user/mniedzie/ttZ/ttVcsMeasure/skimTree $f
#        fi
#        count=$((count+1))
        do if (( $count % 50 == 0)); then
            submitJob $submit "12:00:00"
            makeSubmit $submit $2
        fi
        #filename=${f##*/}                               
        filename=${f///}
        filename=${filename%.*}
        echo "${cwd}/../skimTree $f $outputDir/ > ${outputDir}/${filename}_log.txt 2> ${outputDir}/${filename}_err.txt" >> $submit
        count=$((count+1))
    done
    submitJob $submit "12:00:00"
    rm $submit                                          #remove temporary submit file
}

baseFolder=/pnfs/iihe/cms/store/user/wverbeke/heavyNeutrino/
#baseFolder=/pnfs/iihe/cms/store/user/tutran/heavyNeutrino/
#baseFolder=/pnfs/iihe/cms/store/user/ikhvastu/heavyNeutrino/
cd $baseFolder
#working copies
#foldersMC=*/*ewkino2016MCList-v16p2                              #add suffix for newer versions
#ufoldersMC17=*/*ewkino2017MCList-v16p2                              #add suffix for newer versions
#foldersData=*/*2016LeptonicDataList_v16p2
#foldersData17=*/*2017LeptonicDataList_v16p1
# old one for 2016, new one is v20p1 both for MC and data

#foldersMC=TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/*leptonMvaTrainingList-v5
#foldersMC=TTJets_SingleLeptFromT*TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*leptonMvaTrainingList-v5
#foldersMC=TTJets_SingleLeptFromT*TuneCUETP8M1*/*ewkino2016MCList-v16p2
#foldersMC=TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*ewkino2016MCList-v16p2
#foldersMC=ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/*ewkino2016MCList-v16p2
#foldersMC=DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/*ewkino2016MCList-v20p1
#foldersMC=TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_Moriond2017-v1_ewkino2016MCList-v20p1
#foldersMC=TTJets_SingleLeptFromT*TuneCUETP8M1_13TeV-madgraphMLM-pythia8/*ewkino2016MCList-v21
#foldersMC=*/*ewkino2016MCList-v27
#foldersMC=TTV_MC18_v2/QCD_Pt-15to20*/*
#foldersMC=TTV_MC18_v2/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/*
#foldersMC=WJets*/*ttVMCList-v15
#foldersMC=TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/*colorReconSamples-v1
#foldersMC=TT_TuneCUETP8M2T4_erdON_13TeV-powheg-pythia8/crab_Moriond2017-v1_colorReconSamples-v1
#foldersMC=TT_TuneCUETP8M2T4_GluonMoveCRTune_13TeV-powheg-pythia8/crab_Moriond2017-v1_colorReconSamples-v1
#foldersMC=TT_TuneCUETP8M2T4_QCDbasedCRTune_erdON_13TeV-powheg-pythia8/crab_Moriond2017-v1_colorReconSamples-v1
#foldersMC=TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1_ttV_94X_MCList-v6
#foldersMC=DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/*v10*_ewkino2017MCList-v25
#foldersMC=*/*QCD_2016List_v3
#foldersMC=WJets*/*92X*ttV_94X_MCList-v7
#foldersMC=*/*ewkino2016MCList-v22p2

#foldersData=DoubleMuon/*2017LeptonicDataList_ReReco_v9
#foldersData=SingleElectron/*2017LeptonicDataList_ReReco_v9
#foldersData=DoubleEG/*2017LeptonicDataList_v22p3
#foldersData=SingleElectron/*2016LeptonicDataList_v16p2
#foldersData=DoubleMuon/*2016LeptonicDataList_v9
#foldersData=*/*2017LeptonicDataList_ReReco_v11
#foldersData=DoubleMuon/*2016LeptonicDataList_v20p1
#foldersData=*/*07Aug17*2016LeptonicDataList_v27
#foldersData=*/*2017LeptonicDataList_ReReco_ZMET_v4

foldersMC=*/*dilepton_MC_2018_v5
#foldersMC=*/*dilepton_MC_2018_v7
foldersData=*/*dilepton_data_2018_v5

declare -a array=("TTZToLLNuNu_M-10_T"
                  "ttHToNonbb_M125_Tu"
                  "tZq_ll_4f_ckm_NLO_"
                  "ST_tWll_5f_LO_Tune"
                  "TTWJetsToLNu_TuneC"
                  "THQ_4f_Hincl_13TeV"
                  "THW_5f_Hincl_13TeV"
                  "TTZToLL_M-1to10_Tu"
                  "TTTT_TuneCP5_13TeV"
                  "TTWW_TuneCP5_13TeV"
                  "TTWZ_TuneCP5_13TeV"
                  "TTZZ_TuneCP5_13TeV"
                  "WZTo3LNu_TuneCP5_13TeV-amcat"
                  "DYJetsToLL_M-50_TuneCP5_13TeV-amc"
                  "TTGJets_TuneCP5_13"
                  "ZZTo4L_13TeV_powhe"
                  "GluGluToContinToZZ"
                  "GluGluHToZZTo4L_M1"
                  "VBF_HToZZTo4L_M125"
                  "WplusH_HToZZTo4L_M"
                  "WminusH_HToZZTo4L_"
                  "ZH_HToZZ_4LFilter_"
                  "WWW_4F_TuneCP5_13T"
                  "WWZ_4F_TuneCP5_13T"
                  "WZZ_TuneCP5_13TeV-"
                  "ZZZ_TuneCP5_13TeV-"
                  "WWTo2L2Nu_NNPDF31_"
)
declare -a array2L=("TTZToLLNuNu_M-10_TuneCP5"
																		"DYJetsToLL_M-50_TuneCP5"
																		"TTTo2L2Nu_TuneCP5"
)

declare -a data_array=("MuonEG"
                       "SingleMuon"
)
for d in $foldersMC
do                #skim all samples
    for n in "${array2L[@]}"
    do
        if [[ $d == *$n* ]] 
        then 
            skimSample $d $baseFolder
        fi
    done
done

#for d in $foldersData                #skim all samples
#    do skimSample $d $baseFolder
#done

#for d in $foldersData
#do                #skim all samples
#    for n in "${data_array[@]}"
#    do
#        if [[ $d == *$n* ]] 
#        then 
#            skimSample $d $baseFolder
#        fi
#    done
#done
