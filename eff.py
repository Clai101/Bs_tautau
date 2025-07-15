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
Ds_D_dif = 0.142014
K_s_m = 0.497611
K_p_m = 0.493677

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
#fillParticleList('pi+:skim','pseudo_skim_data_hadron_y4s == 1',path=path)
#applyEventCuts('[nParticlesInList(pi+:skim)!=0]', path=path)
setAnalysisConfigParams({'mcMatchingVersion': 'Belle'}, path)

# To properly read the Belle database the user name is set to g0db
#os.environ['PGUSER'] = 'g0db'

#b2.register_module("EnableMyVariable")
#b2.register_module("EnableMyMetaVariable")

# Load input ROOT files
#b2biiConversion.convertBelleMdstToBelleIIMdst(input_file, enableNisKsFinder=False, path=path)
#inputMdstList(environmentType='default', filelist=[input_file], path=path)


#FSP
fillParticleList('e+:alle','abs(dz) < 2 and dr < 0.5 and eIDBelle > 0.9', path = path)
fillParticleList('mu+:alle','abs(dz) < 2 and dr < 0.5 and muIDBelle > 0.9', path=path)

applyCuts('mu+:alle','isSignal == 1', path=path)
applyCuts('e+:alle','isSignal == 1', path=path)

#reconstructDecay('rho+:1 -> pi+:alle pi0:alle', f'0.52 < InvM and InvM < 1.06', 1, path=path)
#copyLists('rho+:alle',['rho+:1', ], path=path)

#tau_dec = ["e+:alle", "mu+:alle", "pi+:alle", "rho+:alle", "pi+:alle pi+:alle pi-:alle", "pi+:alle gamma:alle"]
tau_dec = ["e+:alle", "mu+:alle"]

for i, dec in enumerate(tau_dec):
    reconstructDecay(f'tau+:{int(i)} -> {dec}', '', int(i), path=path)


copyLists('tau+:alle',[f"tau+:{int(i)}" for i in range(len(tau_dec))], path=path)
applyCuts('tau+:alle','isSignalAcceptMissing == 1', path=path)

variablesToNtuple('tau+:alle', ['pcm', ], treename='tau', filename=output_file, path=path)



fillParticleListFromMC('nu_e:mc_alle', '', path=path)
fillParticleListFromMC('nu_mu:mc_alle', '', path=path)
fillParticleListFromMC('nu_tau:mc_alle', '', path=path)
fillParticleListFromMC('mu+:mc_alle', '', path=path)
fillParticleListFromMC('e+:mc_alle', '', path=path)
reconstructMCDecay('tau+:mc_1 -> e+:mc_alle nu_e:mc_alle anti-nu_tau:mc_alle', '', 1, path=path)
reconstructMCDecay('tau+:mc_2 -> e+:mc_alle nu_e:mc_alle anti-nu_tau:mc_alle', '', 2, path=path)

copyLists('tau+:mc_alle',["tau+:mc_1", "tau+:mc_2"], path=path)

variablesToNtuple('tau+:mc_alle', ['pcm', ], treename='tau_mc', filename=output_file, path=path)

reconstructDecay('B_s0:tautau -> tau+:alle tau-:alle', '', 1, path=path)
applyCuts('B_s0:tautau','isSignalAcceptMissing == 1', path=path)

reconstructMCDecay('B_s0:mc -> tau+:mc_alle tau-:mc_alle', '', 1, path=path)

fillParticleListFromMC('B_s0:mc_only', '', path=path)

from variables import variables as vm

vm.addAlias('pcm','useCMSFrame(p)')
vm.addAlias('cos_Theta_Lc','useCMSFrame(cosTheta)')



# Ntuples
variablesToNtuple('B_s0:tautau', ['pcm', ], treename='BS', filename=output_file, path=path)
variablesToNtuple('B_s0:mc', ['pcm', ], treename='BS_MC', filename=output_file, path=path)
variablesToNtuple('B_s0:mc_only', ['pcm', ], treename='Lc_MC_only', filename=output_file, path=path)



#Process 1000 events
print(path)
#b2.process(path, max_event=100000)
b2.process(path)
print(b2.statistics)



