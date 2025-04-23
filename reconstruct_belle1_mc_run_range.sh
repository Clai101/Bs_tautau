#!/bin/bash
EXPERIMENT=$1
RUN_START=$2
RUN_END=$3
TY=$4
DATA_TYPE=$5
STREAM=$6

WORKDIR=~/basf2/LC_eff/mc_${EXPERIMENT}_${RUN_START}_${RUN_END}_${TY}_${DATA_TYPE}_${STREAM}

URL="http://bweb3/montecarlo.php?ex=${EXPERIMENT}&rs=${RUN_START}&re=${RUN_END}&ty=${TY}&dt=${DATA_TYPE}&bl=caseB&dv=zfserv&st=${STREAM}"

rm -rf ${WORKDIR}
mkdir -p ${WORKDIR}

cp run_efficenty_mc.py ${WORKDIR}

cd ${WORKDIR}
OUTPUT_FILE=${EXPERIMENT}_${RUN_START}_${RUN_END}.root

echo "basf2 -l error run_efficenty_mc.py \"${URL}\" ${OUTPUT_FILE} &> my_output_hack.log" > job_script
chmod 755 job_script
bsub -q s -e error.log -o output.log ./job_script
#./job_script -e error.log -o output.log
