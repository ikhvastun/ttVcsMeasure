#!/bin/bash

#include functions to export CMSSW in submit
source setCMSSW.sh
outputFolder=ntuples_ttz_2L_ttZ_2018_v2

fillJob(){
    name=${1%/}
    name=${name##*/}
    name=${name#ntuples_temp_2L_ttZ_v2*}
    name=${name}_v2.root
    if [ ! -d ~/Work/${outputFolder} ]
        then mkdir ~/Work/${outputFolder}
    fi
#    echo "if [ -f ~/Work/$outputFolder/$name ]" >> $2
#    echo "    then rm ~/Work/$outputFolder/$name" >> $2
#    echo "fi" >> $2
    echo "hadd ~/Work/$outputFolder/$name ${1}/*root" >> $2
#    echo "rm -r $1" >> $2
}

mergeTuple(){
    > mergeJob.sh
    fillJob $1 mergeJob.sh 
    bash mergeJob.sh
    rm mergeJob.sh
}

submitMergeTuple(){
    > mergeJob.sh
    setCMSSW mergeJob.sh
    fillJob $1 mergeJob.sh
    submitJob mergeJob.sh "12:00:00"
    rm mergeJob.sh
}

skimmedTuples=~/Work/ntuples_temp_2L_ttZ_v2*
mergeTuples(){
    for d in $skimmedTuples
        do mergeTuple $d
        #rm -r $d
    done
}

submitMergeTuples(){
    for d in $skimmedTuples
        do submitMergeTuple $d
    done
}
