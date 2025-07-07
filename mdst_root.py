#!/usr/bin/env python3

##########################################################################
# basf2 (Belle II Analysis Software Framework)                           #
# Author: The Belle II Collaboration                                     #
#                                                                        #
# See git log for contributors and copyright holders.                    #
# This file is licensed under LGPL-3.0, see LICENSE.md.                  #
##########################################################################

import sys
import os

import basf2 as b2
from modularAnalysis import *
import b2biiConversion
import mdst
import ROOT
import vertex
import fei
from ROOT import Belle2



if (len(sys.argv) < 3):
    print('Usage: B_converted_apply.py input output')
    exit(1)
input_file  = sys.argv[1]
output_file = sys.argv[2]
skim_id = sys.argv[3]

print(input_file)
print(output_file)

# Create path
path = b2.create_path()
b2.register_module("EnableMyVariable")
b2.register_module("EnableMyMetaVariable")
b2.register_module("SkimFiles")

if skim_id == "Sig_mc":
    if ('evtgen_exp_53' in input_file):
        b2.conditions.prepend_testing_payloads('/home/belle/yasaveev/bb/bin/53_localdb/database.txt')

os.environ['PGUSER'] = 'g0db'
b2biiConversion.convertBelleMdstToBelleIIMdst(input_file, path=path)
setAnalysisConfigParams({'mcMatchingVersion': 'Belle'}, path)

mdst.add_mdst_output(path, mc=True, filename="mdst.root")
mdst.add_udst_output(path, output_file="udst.root", mc=True)




