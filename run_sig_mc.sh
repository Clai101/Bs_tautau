
#!/bin/bash

EXPERIMENT=$1

WORKDIR=~/belle2_MC/mc_${EXPERIMENT}

URL="/gpfs/group/belle2/dataprod/MC/MC15ri/ccbar/sub00/mdst_000${EXPERIMENT}_prod00024786_task10020000${EXPERIMENT}.root"
rm -rf ${WORKDIR}
mkdir -p ${WORKDIR}

#cp run.py channels.py backward_compatibility_layer.py ${WORKDIR}
cp run.py ${WORKDIR}

cd ${WORKDIR}
OUTPUT_FILE=${EXPERIMENT}.root

for file in /gpfs/home/belle2/matrk/Extend/gsim_out/mdst/*.mdst; 
do
  echo "basf2 -l error run.py \"${URL}\" ${OUTPUT_FILE} &> my_output_hack.log" > job_script
  chmod 755 job_script
  bsub -q s -e error.log -o output.log ./job_script

  echo "$file"
done