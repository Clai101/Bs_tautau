
#!/bin/bash


for file in /gpfs/home/belle2/matrk/Extend/gsim_out/mdst/bszeroeff/*.mdst; 
do


  filename=$(basename "$file")
  name_no_ext="${filename%.mdst}"

  WORKDIR=/gpfs/home/belle2/matrk/Extend/eff/mc_${name_no_ext}

  rm -rf ${WORKDIR}
  mkdir -p ${WORKDIR}

  cp eff0.py ${WORKDIR}
  cp Bs_2tau_mc_FEI.py ${WORKDIR}

  cd ${WORKDIR}
  echo "basf2 -l error eff0.py \"${file}\" mc0.root Sig_mc &> mc0.log" > job_script0
  echo "basf2 -l error Bs_2tau_mc_FEI.py \"${file}\" mc1.root Sig_mc &> mc1.log" > job_script1
  chmod 755 job_script0
  chmod 755 job_script1
  bsub -q s ./job_script0
  bsub -q s ./job_script1

  echo "$file"
done