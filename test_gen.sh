#!/bin/bash


EXPERIMENT=69
RUN_START=1097
RUN_END=1097
TY=evtgen-bsbs
DATA_TYPE=5S_onresonance
STREAM=1


WORkCODE="Bs_2tau_mc_FEI.py" 
absolute_path="$(pwd)/${WORkCODE}"
OUTPUT_FILE="mc_gen"

URL="http://bweb3/montecarlo.php?ex=${EXPERIMENT}&rs=${RUN_START}&re=${RUN_END}&ty=${TY}&dt=${DATA_TYPE}&bl=caseB&dv=zfserv&st=${STREAM}"



echo "basf2 -l error $absolute_path \"${URL}\" ${OUTPUT_FILE}.root ${TY}_${STREAM} &> ${OUTPUT_FILE}.log" > job_script
chmod 755 job_script
./job_script
rm job_script
