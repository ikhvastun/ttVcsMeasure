#!/bin/bash

workFolder=$(pwd)
if [ ! -d "$workFolder/output" ]; then                      #make output directory if it doesn't exist
    mkdir -p $workFolder/output
fi
echo "currently in $workFolder"

outputDir=$workFolder/output
CMSSW="CMSSW_9_2_4"

echo "#!/bin/bash              

cd /user/${USER}/${CMSSW}/src
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval \`scram runtime -sh\`

echo "Job started..."

cd $workFolder
cd .. 

${workFolder}/../main > ${outputDir}/out_log.txt 2> ${outputDir}/out_err.txt

echo "...Job Ended!"

" > submit.sh ;

qsub submit.sh -l walltime=01:00:00;
rm submit.sh
