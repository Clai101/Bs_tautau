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

#Fei

print(skim_id)
if (skim_id[-1] == '0' or skim_id[-1] == '1' or skim_id[-1] == '2'):
    b2.conditions.prepend_testing_payloads('/home/belle/yasaveev/bb/fei/training_results/250623_all/localdb/database.txt')
else:
    b2.conditions.prepend_testing_payloads('/home/belle/yasaveev/bb/fei/training_results/060123_all/localdb/database.txt')

if skim_id != "Sig_mc":
    fillParticleList('pi+:skim','pseudo_skim_y5s_' + skim_id[7:] + ' == 1',path=path)
    applyEventCuts('[nParticlesInList(pi+:skim)!=0]', path=path)


particles = fei.get_channels()
configuration = fei.config.FeiConfiguration(prefix='FEI_TEST', training=False, monitor=False, cache=0)
feistate = fei.get_path(particles, configuration)
path.add_path(feistate.path)

rankByHighest('B_s0:generic', 'extraInfo(SignalProbability)', numBest=1, outputVariable='iCand', path=path)
path.add_module('MCMatcherParticles', listName='B_s0:generic', looseMCMatching=True)
applyEventCuts('[nParticlesInList(B_s0:generic)!=0]', path=path)

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
buildEventShape(path=path)

#applyEventCuts('[formula(nParticlesInList(Upsilon(5S):alle)) == 1]', path=path)

copyList('K_L0:alle','K_L0:mdst',path=path)
applyCuts('K_L0:alle','klmClusterBelleTrackFlag == 0',path=path)


from variables import variables as vm

__Alias_names = list()

def add_aliases(*args, **kwargs):
    global __Alias_names
    if args:
        __Alias_names.append(args[0])
    else:
        raise ValueError("No alias name provided")
    vm.addAlias(*args, **kwargs)

#Beam
add_aliases('pcm','useCMSFrame(p)')
add_aliases('ecm','useCMSFrame(E)')
add_aliases('missedE', 'formula(Ecms - useCMSFrame(E))')
add_aliases('recM2_Ups', 'formula((beamE - E)**2 - (beamPx - px)**2 - (beamPy - py)**2 - (beamPz - pz)**2)')

#Ups
add_aliases('pmiss','formula(((beamPx - px)*v*2 + (beamPy - py)**2 + (beamPz - pz)**2)**0.5)')
add_aliases('thetamiss','formula((beamPz - pz) / ((beamPx - px)**2 + (beamPy - py)**2 + (beamPz - pz)**2)**0.5)')
add_aliases('fox','foxWolframR2')
add_aliases('asymmetry', '''formula( 
            (
                daughter(1, daughter(0, daughter(0, pz))) - daughter(1, daughter(1, daughter(1, pz)))
            ) / (
                daughter(1, daughter(0, daughter(0, pz))) + daughter(1, daughter(1, daughter(1, pz)))
            ) 
            ) 
            )''')

#Bs_tag
add_aliases('p0','daughter(0,useCMSFrame(p))')
add_aliases('theta_Bs','daughter(0,useCMSFrame(cosTheta))')
add_aliases('M0','daughter(0,M)')
add_aliases('recM2_Bs', 'formula((beamE - daughter(0,E))**2 - (beamPx - daughter(0,px))**2 - (beamPy - daughter(0,py))**2 - (beamPz - daughter(0,pz))**2)')


add_aliases('idec0', 'daughter(1, daughter(0, extraInfo(decayModeID)))')
add_aliases('idec1', 'daughter(1, daughter(1, extraInfo(decayModeID)))')
add_aliases('is0', 'daughter(0, isSignal)')
add_aliases('N_KL', 'nParticlesInList(K_L0:alle)')

add_aliases('lost_nu_0', '''formula(
            daughter(1, daughter(0, genNMissingDaughter(12))) + 
            daughter(1, daughter(0, genNMissingDaughter(14))) + 
            daughter(1, daughter(0, genNMissingDaughter(16)))
            )''')
add_aliases('Miss_id_0', 'daughter(1, daughter(0, isSignalAcceptMissing))')
add_aliases('lost_gamma_0', 'daughter(1, daughter(0, genNMissingDaughter(22)))')
add_aliases('lost_pi_0', 'daughter(1, daughter(0, genNMissingDaughter(211)))')
add_aliases('lost_K_0', 'daughter(1, daughter(0, genNMissingDaughter(321)))')

add_aliases('lost_nu_1', '''formula(
            daughter(1, daughter(1, genNMissingDaughter(12))) + 
            daughter(1, daughter(1, genNMissingDaughter(14))) + 
            daughter(1, daughter(1, genNMissingDaughter(16)))
            )''')
add_aliases('Miss_id_1', 'daughter(1, daughter(1, isSignalAcceptMissing))')
add_aliases('lost_gamma_1', 'daughter(1, daughter(1, genNMissingDaughter(22)))')
add_aliases('lost_pi_1', 'daughter(1, daughter(1, genNMissingDaughter(211)))')
add_aliases('lost_K_1', 'daughter(1, daughter(1, genNMissingDaughter(321)))')
add_aliases('Bs_lik', 'daughter(0, extraInfo(SignalProbability))')

copyParticles('K_S0:good','K_S0:mdst',path=path)
applyCuts('K_S0:good','extraInfo(ksnbStandard) == 1',path=path)
reconstructDecay('Upsilon(5S):Ks -> B_s0:generic K_S0:good',' ',path=path)

add_aliases('N_KS','nParticlesInList(Upsilon(5S):Ks)')

for tau_index in [0, 1]: 
    alias_name = f'theta_tau_{tau_index}'
    daughter_path = f'daughter(1, daughter({tau_index}, useCMSFrame(cosTheta)))'
    add_aliases(alias_name, daughter_path)
    __Alias_names.append(alias_name)
    
    alias_name = f'p_tau_{tau_index}'
    daughter_path = f'daughter(1, daughter({tau_index}, p))'
    add_aliases(alias_name, daughter_path)
    __Alias_names.append(alias_name)

for tau_index in [0, 1]: 
    alias_name = f'tau_d_{tau_index}_0'
    daughter_path = f'daughter(1, daughter({tau_index}, daughter(0, M)))'
    add_aliases(alias_name, daughter_path)
    __Alias_names.append(alias_name)

for tau_index in [0, 1]:
    add_aliases(f'tau_last_z_{tau_index}', f'daughter(1, daughter({tau_index}, daughter(0, dz)))')
    add_aliases(f'tau_last_r_{tau_index}', f'daughter(1, daughter({tau_index}, daughter(0, dr)))')
    __Alias_names.append(f'tau_last_z_{tau_index}')
    __Alias_names.append(f'tau_last_r_{tau_index}')

for tau_ind in [0, 1]:
    for hypo1 in [0, 1, 2, 4]:  
        for hypo2 in [0, 1, 2, 4]:  
            expr = f'daughter(1, daughter({tau_ind}, daughter(0, atcPIDBelle({hypo1}, {hypo2}))))'
            alias_name = f'PID_{hypo1}_vs_{hypo2}_tau{tau_ind}'
            __Alias_names.append(alias_name)
            add_aliases(alias_name, expr)

applyCuts('Upsilon(5S):alle', '[[[idec0 == 1] or [idec0 == 0]] and [[idec1 == 1] or [idec1 == 0]] and [is0 == 1] and [E_gamma_in_ROE < 0.02]]', path=path)

applyEventCuts('nParticlesInList(Upsilon(5S):alle) > 0', path=path)

mdst.add_mdst_output(filename = 'MCPrint_'+output_file, path=path)


from variables.MCGenTopo import mc_gen_topo
add_module('Upsilon(5S):alle', mc_gen_topo(20000), 'MCGenTopo', output_file, path=path)


#Process 1000 events
print(path)
#b2.process(path, max_event=100000)
b2.process(path)

# Print out the summary
print(b2.statistics)

# Invoke the TopoAna program
os.system('/gpfs/home/belle2/matrk/TopoAna/topoana-master/bin/topoana.exe /gpfs/home/belle2/matrk/TopoAna/topoana-master/test/topoana.card')

print(b2.statistics)