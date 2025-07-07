#!/bin/bash

file="/home/belle2/matrk/Extend/gsim_out/mdst/bszerotautau/evtgen_exp_69_bszerotautau-148.mdst"


OUTPUT_FILE="mc_sig.root"
echo "basf2 -l error Bs_2tau_mc_FEI.py \"${file}\" ${OUTPUT_FILE} Sig_mc &> mc.log" > job_script
chmod 755 job_script
./job_script
rm job_script

