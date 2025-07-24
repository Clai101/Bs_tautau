
#!/bin/bash


for file in /gpfs/home/belle2/matrk/Extend/gsim_out/mdst/bszeroeff/*.mdst; 
do


  filename=$(basename "$file")
  name_no_ext="${filename%.mdst}"

  WORKDIR=/gpfs/home/belle2/matrk/Extend/eff/mc_${name_no_ext}

  rm -rf ${WORKDIR}
  mkdir -p ${WORKDIR}

  cp eff0.py ${WORKDIR}
  cp eff1.py ${WORKDIR}

  cd ${WORKDIR}
  echo "basf2 -l error eff0.py \"${file}\" mc0.root Sig_mc &> mc0.log" > job_script
  echo "basf2 -l error eff1.py \"${file}\" mc1.root Sig_mc &> mc1.log" > job_script
  chmod 755 job_script
  bsub -q s ./job_script

  echo "$file"
done