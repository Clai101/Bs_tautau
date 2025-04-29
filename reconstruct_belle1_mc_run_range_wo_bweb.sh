#!/bin/bash
EXPERIMENT=$1
RUN_START=$2
RUN_END=$3
TY=$4
DATA_TYPE=$5
STREAM=$6

WORKDIR=/home/belle2/matrk/B_tautau/run_Bs/gen_mc_${EXPERIMENT}_${RUN_START}_${RUN_END}_${TY}_${DATA_TYPE}_${STREAM}

rm -rf ${WORKDIR}
mkdir -p ${WORKDIR}

OUTPUT_FILE_txt=${EXPERIMENT}_${RUN_START}_${RUN_END}_script

cp /home/belle2/matrk/B_tautau/run_Bs/run.py ${WORKDIR}
cp ~/hse_against_kekcc/running_scripts/${TY}_${STREAM}/${OUTPUT_FILE_txt} ${WORKDIR}
cd ${WORKDIR}

chmod 755 ${OUTPUT_FILE_txt}
bsub -q l -e error.log -o output.log ./${OUTPUT_FILE_txt}
