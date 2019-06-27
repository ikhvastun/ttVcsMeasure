#!/bin/bash

# 
#outputDir=/user/ikhvastu/Work/ttVcsMeasure/ttVselection/output/
outputDir=$PWD/../output/
ttVdir=$PWD/../
mkdir -p ${outputDir}
CMSSW="CMSSW_9_4_4"

if [ -z "$1" ] && [ -z "$2" ]
  then
    echo "please specify dataset and control region name"
    echo "dataset from data/samples/ folder"
    echo "possible control regions: ttZ3L, WZ, Zgamma, ttbar and DY"
    exit 1
fi

sampleFile=$1
selection=$2

year="2017"
if [[ $sampleFile = *","* ]]; then
    year="comb"
elif [[ $sampleFile = *"2017"* ]]; then
    year="2017"
else 
    year="2016"
fi

echo "#!/bin/bash              

cd $CMSSW_BASE/src
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval \`scram runtime -sh\`

echo "Job started..."

cd ${ttVdir}

./main $sampleFile runFullSelection selection:$selection > ${outputDir}/out_log_${selection}_${year}.txt 2> ${outputDir}/out_err_${selection}_${year}.txt

echo "...Job Ended!"

" > submit.sh ;

qsub submit.sh -l walltime=12:00:00;
#echo 'output files are', out_err_${selection}_${year}
#echo 'everything is fine'
rm submit.sh

#cd /user/${USER}/${CMSSW}/src
