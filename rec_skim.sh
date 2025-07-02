#!/bin/bash
EXPERIMENT=$1
RUN_START=$2
RUN_END=$3
TY=$4
DATA_TYPE=$5
STREAM=$6
TIME=$7

WORKDIR=/group/belle2/users2022/matrk/B_Dlnu/mc_${EXPERIMENT}_${RUN_START}_${RUN_END}_${TY}_${DATA_TYPE}_${STREAM}
WORkCODE="Bs_my_skim.py" 
absolute_path="$(pwd)/${WORkCODE}"
OUTPUT_FILE="mc"

URL="http://bweb3/montecarlo.php?ex=${EXPERIMENT}&rs=${RUN_START}&re=${RUN_END}&ty=${TY}&dt=${DATA_TYPE}&bl=caseB&dv=zfserv&st=${STREAM}"

rm -rf ${WORKDIR}
mkdir -p ${WORKDIR}
cd ${WORKDIR}

echo "basf2 -l error $absolute_path \"${URL}\" ${OUTPUT_FILE}.root ${TY}_${STREAM} &> ${OUTPUT_FILE}.log" > job_script
chmod 755 job_script
bsub -q ${TIME} ./job_script