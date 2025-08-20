
#!/bin/bash


for file in /gpfs/home/belle2/matrk/Extend/gsim_out/mdst/bszero_lepton/*.mdst; 
do


  filename=$(basename "$file")
  name_no_ext="${filename%.mdst}"

  WORKDIR=/gpfs/home/belle2/matrk/Extend/lep_mc/mc_${name_no_ext}

  rm -rf ${WORKDIR}
  mkdir -p ${WORKDIR}

  cp Bs_2tau_mc_FEI.py ${WORKDIR}

  cd ${WORKDIR}
  echo "basf2 -l error Bs_2tau_mc_FEI.py \"${file}\" mc.root Sig_mc &> mc.log" > job_script
  chmod 755 job_script
  bsub -q s ./job_script

  echo "$file"
done