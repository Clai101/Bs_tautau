
#!/bin/bash


for file in /gpfs/home/belle2/matrk/Extend/gsim_out/mdst/*.mdst; 
do


  filename=$(basename "$file")
  name_no_ext="${filename%.mdst}"

  WORKDIR=/gpfs/home/belle2/matrk/Extend/sig_mc/mc_${name_no_ext}

  rm -rf ${WORKDIR}
  mkdir -p ${WORKDIR}

  cp Bs_2tau_mc.py ${WORKDIR}

  cd ${WORKDIR}
  OUTPUT_FILE=${name_no_ext}.root
  echo "basf2 -l error Bs_2tau_mc.py \"${file}\" ${OUTPUT_FILE} &> mc.log" > job_script
  chmod 755 job_script
  bsub -q s ./job_script

  echo "$file"
done