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
import fei
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

print('\n\n\n_________________________________________0_________________________________________\n\n\n')
# arguments


if (len(sys.argv) < 2):
    print('Usage: B_converted_apply.py input output')
    exit(1)
input_file  = sys.argv[1]
output_file = sys.argv[2]

print(input_file)
print(output_file)
print('\n\n\n_________________________________________1_________________________________________\n\n\n')

# Create path
path = b2.create_path()
b2.register_module("EnableMyVariable")
b2.register_module("EnableMyMetaVariable")
b2.register_module("SkimFiles")


print('\n\n\n_________________________________________2_________________________________________\n\n\n')

# Load input ROOT files
# Load input ROOT files
os.environ['PGUSER'] = 'g0db'
b2biiConversion.convertBelleMdstToBelleIIMdst(input_file, path=path)
setAnalysisConfigParams({'mcMatchingVersion': 'Belle'}, path)

print('\n\n\n_________________________________________3_________________________________________\n\n\n')


#Fei
skim_id = sys.argv[3]
print(skim_id)
if (skim_id[-1] == '0' or skim_id[-1] == '1' or skim_id[-1] == '2'):
    b2.conditions.prepend_testing_payloads('/home/belle/yasaveev/bb/fei/training_results/250623_all/localdb/database.txt')
else:
    b2.conditions.prepend_testing_payloads('/home/belle/yasaveev/bb/fei/training_results/060123_all/localdb/database.txt')

print('\n\n\n_________________________________________4_________________________________________\n\n\n')
#SKIM
fillParticleList('pi+:skim','pseudo_skim_y5s_' + skim_id[7:] + ' == 1',path=path)
applyEventCuts('[nParticlesInList(pi+:skim)!=0]', path=path)

print('\n\n\n_________________________________________5_________________________________________\n\n\n')

particles = fei.get_channels()
configuration = fei.config.FeiConfiguration(prefix='FEI_TEST', training=False, monitor=False, cache=0)
feistate = fei.get_path(particles, configuration)
path.add_path(feistate.path)

print('\n\n\n_________________________________________6_________________________________________\n\n\n')


rankByHighest('B_s0:generic', 'extraInfo(SignalProbability)', numBest=1, outputVariable='iCand', path=path)
path.add_module('MCMatcherParticles', listName='B_s0:generic', looseMCMatching=True)
applyEventCuts('[nParticlesInList(B_s0:generic)!=0]', path=path)

print('\n\n\n_________________________________________7_________________________________________\n\n\n')


#Part 1

#FSP
copyLists('B_s0:alle',['B_s0:generic', ], path=path)
path.add_module('MCMatcherParticles', listName='B_s0:alle', looseMCMatching=True)

# Part 2

#FSP +
fillParticleList('pi+:alle','abs(dr) < 0.5 and abs(dz) < 2', path = path)
path.add_module('MCMatcherParticles', listName='pi+:alle', looseMCMatching=True)

copyParticles('pi0:alle','pi0:mdst', path=path)
copyParticles('gamma:alle','gamma:mdst',path=path)
applyCuts('pi0:alle',f'abs(InvM - {Pi_0_m}) < 0.015 and daughter(0, E) > 0.05 and daughter(1, E) > 0.05',path=path)
applyCuts('gamma:alle','goodBelleGamma == 1 and clusterBelleQuality == 0',path=path)

vertex.kFit('pi0:alle', conf_level = 0.0, fit_type = 'mass',path=path)


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

copyList('K_L0:alle','K_L0:mdst',path=path)
applyCuts('K_L0:alle','klmClusterBelleTrackFlag == 0',path=path)

print('\n\n\n_________________________________________8_________________________________________\n\n\n')


from variables import variables as vm

vm.addAlias('pcm','useCMSFrame(p)')
vm.addAlias('p0','daughter(0,useCMSFrame(p))')
vm.addAlias('theta_Bs','daughter(0,useCMSFrame(cosTheta))')
vm.addAlias('M0','daughter(0,M)')
vm.addAlias('missedE', 'formula(Ecms - useCMSFrame(E))')

vm.addAlias('ecm','useCMSFrame(E)')

vm.addAlias('recM2', 'formula((beamE - E)**2 - (beamPx - px)**2 - (beamPy - py)**2 - (beamPz - pz)**2)')
vm.addAlias('idec0', 'daughter(1, daughter(0, extraInfo(decayModeID)))')
vm.addAlias('idec1', 'daughter(1, daughter(1, extraInfo(decayModeID)))')
vm.addAlias('is0', 'daughter(0, isSignal)')
vm.addAlias('N_KL', 'nParticlesInList(K_L0:alle)')

vm.addAlias('lost_nu_0', 'daughter(1, daughter(0, genNMissingDaughter(18)))')
vm.addAlias('Miss_id_0', 'daughter(1, daughter(0, isSignalAcceptMissing))')
vm.addAlias('lost_gamma_0', 'daughter(1, daughter(0, genNMissingDaughter(22)))')
vm.addAlias('lost_pi_0', 'daughter(1, daughter(0, genNMissingDaughter(211)))')
vm.addAlias('lost_K_0', 'daughter(1, daughter(0, genNMissingDaughter(321)))')

vm.addAlias('lost_nu_1', 'daughter(1, daughter(1, genNMissingDaughter(18)))')
vm.addAlias('Miss_id_1', 'daughter(1, daughter(1, isSignalAcceptMissing))')
vm.addAlias('lost_gamma_1', 'daughter(1, daughter(1, genNMissingDaughter(22)))')
vm.addAlias('lost_pi_1', 'daughter(1, daughter(1, genNMissingDaughter(211)))')
vm.addAlias('lost_K_1', 'daughter(1, daughter(1, genNMissingDaughter(321)))')
vm.addAlias('Bs_lik', 'daughter(0, extraInfo(SignalProbability))')

print('\n\n\n_________________________________________9_________________________________________\n\n\n')


for tau_index in [0, 1]: 
    for d_index in range(3):  
        alias_name = f'theta_tau_d_{tau_index}_{d_index}'
        daughter_path = f'daughter(1, daughter({tau_index}, daughter({d_index}, useCMSFrame(cosTheta))))'
        vm.addAlias(alias_name, daughter_path)

# p_tau_dd_{tau_index}_{d1}_{d2} — если есть глубже: дочка от дочки
for tau_index in [0, 1]:
        for d2 in range(2):
            alias_name = f'theta_tau_dd_{tau_index}_0_{d2}'
            daughter_path = f'daughter(1, daughter({tau_index}, daughter(0, daughter({d2}, useCMSFrame(cosTheta)))))'
            vm.addAlias(alias_name, daughter_path)

print('\n\n\n_________________________________________10_________________________________________\n\n\n')


# Ntuplestau_dec = ["e+:alle", "mu+:alle", "pi+:alle", "rho+:alle", "pi+:alle pi+:alle pi-:alle", "pi+:alle gamma:alle"]

for tau_index in [0, 1]:  # для двух τ
    # last point z
    vm.addAlias(f'tau_last_z_{tau_index}', f'daughter(1, daughter({tau_index}, z))')
    vm.addAlias(f'tau_last_r_{tau_index}', f'daughter(1, daughter({tau_index}, r))')

print('\n\n\n_________________________________________11_________________________________________\n\n\n')

variablesToNtuple('Upsilon(5S):alle', ['missedE','M0', 'p0', 'recM2', 'idec0', 'idec1', 'totalEnergyMC', 
                                       'E_gamma_in_ROE', 'N_tracks_in_ROE', 'is0', 'lost_nu_0', 'lost_gamma_0', 
                                       'lost_pi_0', 'lost_K_0', 'Miss_id_0', 'lost_nu_1', 'lost_gamma_1', 'lost_pi_1', 
                                       'lost_K_1', 'Miss_id_1',],
                     treename='Y5S', filename=output_file, path=path)

"""variablesToNtuple('Upsilon(5S):alle', ['N_KL', 'idec0', 'idec1', 'totalEnergyMC', 'E_gamma_in_ROE', 
                                       'N_tracks_in_ROE', 'Bs_lik', 'is0', 'lost_nu_0', 'lost_gamma_0', 'lost_pi_0', 'lost_K_0', 
                                       'Miss_id_0', 'lost_nu_1', 'lost_gamma_1', 'lost_pi_1', 'lost_K_1', 'Miss_id_1',
                                        #Full-event парметры
                                        'missedE','M0', 'p0', 'recM2',


                                        # Последние точки трека
                                        'tau_last_z_0', 'tau_last_z_1',
                                        'tau_last_r_0', 'tau_last_r_1',

                                        # Углы дочерей tau
                                        'theta_tau_d_0_0', 'theta_tau_d_0_1', 'theta_tau_d_0_2',
                                        'theta_tau_d_1_0', 'theta_tau_d_1_1', 'theta_tau_d_1_2',

                                        # Углы дочерей rho
                                        'theta_tau_dd_0_0_0', 'theta_tau_dd_0_0_1',
                                        'theta_tau_dd_1_0_0', 'theta_tau_dd_1_0_1'
                                       ],
                     treename='Y5S', filename=output_file, path=path)"""

#Process 1000 events
print(path)
#b2.process(path, max_event=100000)
b2.process(path)
print(b2.statistics)