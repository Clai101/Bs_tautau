#!/bin/bash

cd reconstruction/
export USE_GRAND_REPROCESS_DATA=1
export BELLE_LEVEL=b20090127_0910
export BELLE_DEBUG=debug
export BELLE_MSG_MAX_SHLVL=1
export ROOTSYS=/sw/belle/cern/root_v5.34.36
export ANALYSIS_TOP_DIR=./
source /cvmfs/belle.cern.ch/tools/b2setup
export B2BII=true

basf2 run.py "http://bweb3/mdst.php?ex=39&rs=0&re=12&skm=HadronBorJ&dt=on_resonance&bl=caseB" 39_0_12.root