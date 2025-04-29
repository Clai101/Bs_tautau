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
import ROOT
import vertex
from ROOT import Belle2

Lamc_m = 2.28646
Lamc_25_m = 2.5925
Lamc_26_m = 2.628
Lam_m = 1.115683
D_0_m = 1.86483
D_p_m = 1.86966
D_st_p_m = 2.01026
D_st_0_m = 2.00685
Pi_p_m = 0.13957
Pi_0_m = 0.13498
D_st_D_dif = 0.142014
K_s_m = 0.497611
K_p_m = 0.493677
Bs_m = 5.36693
tau_m = 1.77693
mu_m = 0.1056583755
D_s_m = 1.96835

# arguments
if (len(sys.argv) < 2):
    print('Usage: B_converted_apply.py input output')
    exit(1)

input_file  = sys.argv[1]
output_file = sys.argv[2]

# Create path
path = b2.create_path()

b2.register_module("EnableMyVariable")
b2.register_module("EnableMyMetaVariable")

# Load input ROOT files
os.environ['PGUSER'] = 'g0db'
b2biiConversion.convertBelleMdstToBelleIIMdst(input_file, path=path)
setAnalysisConfigParams({'mcMatchingVersion': 'Belle'}, path)

#Part 1

#FSP
fillParticleList('pi+:alle','abs(dr) < 0.5 and abs(dz) < 2', path = path)
fillParticleList('K+:alle','abs(dr) < 0.5 and abs(dz) < 2 and atcPIDBelle(3,2) > 0.6', path=path)

path.add_module('MCMatcherParticles', listName='pi+:alle', looseMCMatching=True)
path.add_module('MCMatcherParticles', listName='K+:alle', looseMCMatching=True)

copyParticles('pi0:alle','pi0:mdst', path=path)
copyParticles('K_S0:alle','K_S0:mdst', path=path)

copyList('K_L0:alle','K_L0:mdst',path=path)
applyCuts('K_L0:alle','klmClusterBelleTrackFlag == 0',path=path)
applyEventCuts('[formula(nParticlesInList(K_L0:alle))<=0]', path=path)

applyCuts('K_S0:alle',f'abs(InvM - {K_s_m}) < 0.015 and goodBelleKshort == 1',path=path)
applyCuts('pi0:alle',f'abs(InvM - {Pi_0_m}) < 0.015 and daughter(0, E) > 0.05 and daughter(1, E) > 0.05',path=path)

vertex.kFit('pi0:alle', conf_level = 0.0, fit_type = 'mass',path=path)
vertex.kFit('K_S0:alle', conf_level = 0.0, fit_type = 'massvertex',path=path)

#Intermediate particles
reconstructDecay('D_s+:1 -> K+:alle K-:alle pi+:alle', f'abs(InvM - {D_s_m}) < 0.015', 1, path=path)
copyLists('D_s+:alle',['D_s+:1', ], path=path)

path.add_module('MCMatcherParticles', listName='D_s+:alle', looseMCMatching=True)
vertex.kFit('D_s+:alle', conf_level = 0.0, fit_type = 'massvertex',path=path)

reconstructDecay('B_s0:1 -> D_s-:alle pi+:alle', f'5.25 <= InvM  and  InvM <= 5.51', 1, path=path)
copyLists('B_s0:alle',['B_s0:1', ], path=path)

from variables import variables as vm

vm.addAlias('pcm','useCMSFrame(p)')

path.add_module('MCMatcherParticles', listName='B_s0:alle', looseMCMatching=True)
variablesToNtuple('B_s0:alle', ['M', 'pcm', 'isSignal'],
                     treename='BS', filename='Bs_2tau_sig_MC.root', path=path)
#vertex.kFit('B_s0:alle', conf_level = 0.0, fit_type = 'massvertex', path=path)

#Process 1000 events
print(path)
#b2.process(path, max_event=100000)
b2.process(path)
print(b2.statistics)