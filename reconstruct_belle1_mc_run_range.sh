#!/bin/bash
EXPERIMENT=$1
RUN_START=$2
RUN_END=$3
TY=$4
DATA_TYPE=$5
STREAM=$6

WORKDIR=~/B_tautau/run_Bs/Gen_MC
WORkCODE="Bs_2tau_mc_FEI.py" 
absolute_path="$(cd "$(dirname "$0")" && pwd)/${WORkCODE}"
OUTPUT_FILE="mc_${EXPERIMENT}_${RUN_START}_${RUN_END}_${TY}_${DATA_TYPE}_${STREAM}"

URL="http://bweb3/montecarlo.php?ex=${EXPERIMENT}&rs=${RUN_START}&re=${RUN_END}&ty=${TY}&dt=${DATA_TYPE}&bl=caseB&dv=zfserv&st=${STREAM}"

cd ${WORKDIR}

echo "basf2 -l error $absolute_path \"${URL}\" ${OUTPUT_FILE}.root ${TY}_${STREAM} &> ${OUTPUT_FILE}.log"
chmod 755 job_script
bsub -q l ./job_script
