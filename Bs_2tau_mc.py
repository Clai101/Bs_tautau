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
from pathlib import Path
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

"""
# arguments
if (len(sys.argv) < 2):
    print('Usage: B_converted_apply.py input output')
    exit(1)
input_file  = sys.argv[1]
output_file = sys.argv[2]
"""

# Create path
path = b2.create_path()
b2.register_module("EnableMyVariable")
b2.register_module("EnableMyMetaVariable")

# Load input ROOT files


directory = Path("/gpfs/home/belle2/matrk/Extend/gsim_out/mdst")

mdst_files = list(directory.glob("*.mdst"))

full_paths = [str(file.resolve()) for file in mdst_files]

os.environ['PGUSER'] = 'g0db'
b2biiConversion.convertBelleMdstToBelleIIMdst(full_paths[:5], path=path)
setAnalysisConfigParams({'mcMatchingVersion': 'Belle'}, path)

#Part 1

#FSP
fillParticleList('pi+:alle','abs(dr) < 0.5 and abs(dz) < 2', path = path)
fillParticleList('K+:alle','abs(dr) < 0.5 and abs(dz) < 2 and atcPIDBelle(3,2) > 0.6', path=path)

path.add_module('MCMatcherParticles', listName='pi+:alle', looseMCMatching=True)
path.add_module('MCMatcherParticles', listName='K+:alle', looseMCMatching=True)

copyParticles('pi0:alle','pi0:mdst', path=path)
copyParticles('K_S0:alle','K_S0:mdst', path=path)

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

path.add_module('MCMatcherParticles', listName='B_s0:alle', looseMCMatching=True)
variablesToNtuple('B_s0:alle', ['M', 'pcm', 'isSignal'],
                     treename='BS', filename='Bs_2tau_sig_MC.root', path=path)
#vertex.kFit('B_s0:alle', conf_level = 0.0, fit_type = 'massvertex', path=path)

#applyEventCuts('[formula(nParticlesInList(B_s0:alle) + nParticlesInList(anti-B_s0:alle)) == 1]', path=path)

# Part 2

#FSP +

copyList('K_L0:alle','K_L0:mdst',path=path)

copyParticles('gamma:alle','gamma:mdst',path=path)
applyCuts('gamma:alle','goodBelleGamma == 1 and clusterBelleQuality == 0',path=path)


fillParticleList('e+:alle','abs(dz) < 2 and dr < 0.5 and eIDBelle > 0.9', path = path)
fillParticleList('mu+:alle','abs(dz) < 2 and dr < 0.5 and muIDBelle > 0.9', path=path)


reconstructDecay('rho+:1 -> pi+:alle pi0:alle', f'0.52 < InvM and InvM < 1.06', 1, path=path)
copyLists('rho+:alle',['rho+:1', ], path=path)

tau_dec = ["e+:alle", "mu+:alle", "pi+:alle", "rho+:alle", "pi+:alle pi+:alle pi-:alle", "pi+:alle gamma:alle"]
for i, dec in enumerate(tau_dec):
    reconstructDecay(f'tau+:{int(i)} -> {dec}', '', int(i), path=path)

copyLists('tau+:alle',[f"tau+:{int(i)}" for i in range(len(tau_dec))], path=path)
path.add_module('MCMatcherParticles', listName='tau+:alle', looseMCMatching=True)


reconstructDecay('B_s0:tautau -> tau+:alle tau-:alle', '', 1, ignoreIfTooManyCandidates = False, path=path)
reconstructDecay('B_s0:tau -> tau+:alle', '', 2, ignoreIfTooManyCandidates = False, allowChargeViolation = True, path=path)
copyLists('B_s0:tauonic',['B_s0:tautau', 'B_s0:tau'], path=path)

reconstructDecay('Upsilon(5S):alle0 -> B_s0:alle B_s0:tauonic', '', 1, path=path)
reconstructDecay('Upsilon(5S):alle1 -> B_s0:alle anti-B_s0:tauonic', '', 2, path=path)
copyLists('Upsilon(5S):alle',['Upsilon(5S):alle0', 'Upsilon(5S):alle1'], path=path)

applyCuts('Upsilon(5S):alle', 'N_tracks_in_ROE == 0', path=path)

#applyEventCuts('[formula(nParticlesInList(Upsilon(5S):alle)) == 1]', path=path)

from variables import variables as vm

vm.addAlias('pcm','useCMSFrame(p)')
vm.addAlias('p0','daughter(0,useCMSFrame(p))')
vm.addAlias('M0','daughter(0,M)')
vm.addAlias('missedE', 'formula(Ecms - useCMSFrame(E))')

vm.addAlias('ecm','useCMSFrame(E)')

vm.addAlias('recM2', 'formula((beamE - E)**2 - (beamPx - px)**2 - (beamPy - py)**2 - (beamPz - pz)**2)')
vm.addAlias('idec0', 'daughter(1, daughter(0, extraInfo(decayModeID)))')
vm.addAlias('idec1', 'daughter(1, daughter(1, extraInfo(decayModeID)))')
vm.addAlias('is0', 'daughter(0, isSignal)')
vm.addAlias('N_KL', 'nParticlesInList(K_L0:alle)')

vm.addAlias('lost_nu_0', 'formula(daughter(1, daughter(0, genNMissingDaughter(12))) + daughter(1, daughter(0, genNMissingDaughter(14))) + daughter(1, daughter(0, genNMissingDaughter(16))))')
vm.addAlias('Miss_id_0', 'daughter(1, daughter(0, isSignalAcceptMissing))')
vm.addAlias('lost_gamma_0', 'daughter(1, daughter(0, genNMissingDaughter(22)))')
vm.addAlias('lost_pi_0', 'daughter(1, daughter(0, genNMissingDaughter(211)))')
vm.addAlias('lost_K_0', 'daughter(1, daughter(0, genNMissingDaughter(321)))')

vm.addAlias('lost_nu_1', 'formula(daughter(1, daughter(1, genNMissingDaughter(12))) + daughter(1, daughter(1, genNMissingDaughter(14))) + daughter(1, daughter(1, genNMissingDaughter(16))))')
vm.addAlias('Miss_id_1', 'daughter(1, daughter(1, isSignalAcceptMissing))')
vm.addAlias('lost_gamma_1', 'daughter(1, daughter(1, genNMissingDaughter(22)))')
vm.addAlias('lost_pi_1', 'daughter(1, daughter(1, genNMissingDaughter(211)))')
vm.addAlias('lost_K_1', 'daughter(1, daughter(1, genNMissingDaughter(321)))')

vm.addAlias('p_tau_d_1_0', 'daughter(1, daughter(1,  daughter(0, p)))')
vm.addAlias('p_tau_d_1_1', 'daughter(1, daughter(1,  daughter(1, p)))')
vm.addAlias('p_tau_d_1_2', 'daughter(1, daughter(1,  daughter(2, p)))')

vm.addAlias('p_tau_d_0_0', 'daughter(1, daughter(0,  daughter(0, p)))')
vm.addAlias('p_tau_d_0_1', 'daughter(1, daughter(0,  daughter(1, p)))')
vm.addAlias('p_tau_d_0_2', 'daughter(1, daughter(0,  daughter(2, p)))')

vm.addAlias('p_tau_d_1_0', 'daughter(1, daughter(0,  daughter(0, p)))')
vm.addAlias('p_tau_d_1_1', 'daughter(1, daughter(0,  daughter(1, p)))')
vm.addAlias('p_tau_d_1_2', 'daughter(1, daughter(0,  daughter(2, p)))')

vm.addAlias('p_tau_dd_0_0', 'daughter(1, daughter(0,  daughter(0, daughter(0, p))))')
vm.addAlias('p_tau_dd_0_1', 'daughter(1, daughter(0,  daughter(1, daughter(1, p))))')

vm.addAlias('p_tau_dd_1_0', 'daughter(1, daughter(1,  daughter(0, daughter(0, p))))')
vm.addAlias('p_tau_dd_1_1', 'daughter(1, daughter(1,  daughter(0, daughter(1, p))))')

# Ntuples
variablesToNtuple('Upsilon(5S):alle', ['missedE','M0', 'p0', 'recM2', 'idec0', 'idec1', 'totalEnergyMC', 
                                       'E_gamma_in_ROE', 'N_tracks_in_ROE', 'is0', 'lost_nu_0', 'lost_gamma_0', 
                                       'lost_pi_0', 'lost_K_0', 'Miss_id_0', 'lost_nu_1', 'lost_gamma_1', 'lost_pi_1', 
                                       'lost_K_1', 'Miss_id_1', 
                                       'p_tau_d_1_0', 'p_tau_d_1_1', 'p_tau_d_1_2', 
                                       'p_tau_d_0_0', 'p_tau_d_0_1', 'p_tau_d_0_2',
                                       'p_tau_dd_0_0', 'p_tau_dd_0_1', 'p_tau_dd_1_0', 'p_tau_dd_1_1'],
                     treename='Y5S', filename='Bs_2tau_sig_MC.root', path=path)

#Process 1000 events
print(path)
#b2.process(path, max_event=100000)
b2.process(path)
print(b2.statistics)
