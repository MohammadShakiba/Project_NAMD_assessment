import sys
import cmath
import math
import os
import h5py
import matplotlib.pyplot as plt   # plots
import numpy as np
import time
import warnings
import argparse

from liblibra_core import *
import util.libutil as comn
from libra_py import units
import models
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.dynamics.tsh.plot as tsh_dynamics_plot
import libra_py.data_savers as data_savers
import libra_py.models.GLVC as GLVC

parser = argparse.ArgumentParser()
parser.add_argument('--iter', default=0, type=int)
parser.add_argument('--ntraj', default=10, type=int)
parser.add_argument('--istate', default=1, type=int)
parser.add_argument('--rep', default=1, type=int)
parser.add_argument('--elec_int', default=2, type=int)
parser.add_argument('--mqcmethod', default=1, type=int)
parser.add_argument('--tshmethod', default='fssh', type=str)
parser.add_argument('--dectime', default='', type=str)
parser.add_argument('--nsteps', default=100000, type=int)
parser.add_argument('--dt', default=1.0, type=float)
parser.add_argument('--temperature', default=300.0, type=float)
parser.add_argument('--reorg_energy', default=0.0125, type=float)
parser.add_argument('--A', default=10.0, type=float)
parser.add_argument('--eps', default=0.1, type=float)
parser.add_argument('--ssy', default=False, type=bool)
parser.add_argument('--wpwidth', default=1.0, type=float)
args = parser.parse_args()

s = units.inv_cm2Ha
freqs = np.array([129.743484,263.643268,406.405874,564.000564,744.784527,
                  961.545250,1235.657107,1606.622430,2157.558683,3100.000000])*s
freqs = freqs.tolist()
g_vals = [0.000563,0.001143,0.001762,0.002446,0.003230,0.004170,0.005358,0.006966,0.009356,0.013444]
eps = 0.02
v0 = 0.001
m = 1836.0 # 1 amu
temperature = args.temperature # K
model_params = {"model": 0, "model0": 0,"nstates": 2, "num_osc": 10, "coupling_scaling": [1.0, -1.0], 
                "omega": [freqs]*2, "coupl": [g_vals]*2, "mass": [m]*10, "Ham": [[eps,-v0],[-v0,-eps]],
                "beta": 1/(units.kB*temperature)}


############################################################
##### 2. Choose the Nonadiabatic Dynamics Methodology 
############################################################
NSTATES = 2
A = MATRIX(10,1); 
for i in range(10):
    A.set(i, 0, 1/args.A) #; A.set(1,1, args.A)
# Initialize the parameters in dyn_general, these parameters will
# be updated for each mqc method below
dyn_general = { "nsteps":args.nsteps, "ntraj":args.ntraj, "nstates":NSTATES,
                "dt": args.dt, "num_electronic_substeps":1, "isNBRA":0, "is_nbra":0,
                "progress_frequency":0.1, "which_adi_states":range(NSTATES), "which_dia_states":range(NSTATES),      
                "mem_output_level":4,
                #"properties_to_save":[ "timestep", "time","se_pop_adi", "se_pop_dia", "sh_pop_adi", "sh_pop_dia",
                #                       "ave_decoherence_rates" ],  
                "properties_to_save":[ "time", "sh_pop_adi", "sh_pop_dia", "se_pop_adi", "se_pop_dia"],
                "prefix":"", "prefix2":"",
                "ham_update_method":1, "ham_transform_method":1, "time_overlap_method":1, 
                "nac_update_method":1, "hvib_update_method":1, "force_method":1, "rep_force":1,
                "hop_acceptance_algo":20, "momenta_rescaling_algo":201, 
                "use_Jasper_Truhlar_criterion":1, "tsh_method":0, "decoherence_algo":-1,
                "decoherence_times_type":-1, "dephasing_informed":0, 
                "decoherence_rates":MATRIX(2,2), 
                "ave_gaps":MATRIX(2,2), 
                "do_ssy":0, "rep_tdse":args.rep, "electronic_integrator":args.elec_int, 
                "state_tracking_algo":-1, "do_phase_correction":1, 
                "A_val_schw": args.A, "eps_val_edc": args.eps,
                "schwartz_decoherence_inv_alpha": A,
                "decoherence_C_param": 1.0, "decoherence_eps_param":0.1
              }
# The assessing method
# The MQC method:
#   1: no decoherence
#   2: BCSH
#   3: ID-A
#   4: SHXF -> with external parameters of: wp_width, coherence_threshold, project_out_aux
#   5: SDM  -> with external parameters of: decoherence_times_type, decoherence_C_param, decoherence_eps_param, schwartz_decoherence_inv_alpha
#   6: DISH  -> with external parameters of: decoherence_times_type, decoherence_C_param, decoherence_eps_param, schwartz_decoherence_inv_alpha
mqcmethod = args.mqcmethod
# The TSH mehod: fssh, fssh2, gfsh
tshmethod = args.tshmethod
# Decohrence time method: schw1, schw2, and edc
dectime = args.dectime
# =========== No decoherence
if mqcmethod==1:
    folder_prefix = ''
# ============ BCSH, ID-A, and SHXF
elif mqcmethod==2:
    dyn_general.update({"decoherence_algo": 3}) # BCSH
    folder_prefix = 'BCSH'
elif mqcmethod==3:
    dyn_general.update({"decoherence_algo": 1}) # ID-A
    folder_prefix = 'ID-A'
elif mqcmethod==4:
    dyn_general.update({ "decoherence_algo":5}) # SHXF
    #=========== XF setting ==============================
    dyn_general.update({"wp_width":0.3, "coherence_threshold": 0.01 })
    dyn_general.update({"project_out_aux": 1})
    folder_prefix = 'SHXF'
# ============ SDM
elif mqcmethod==5:
    dyn_general.update({ "decoherence_algo":0}) # msdm
    if dectime=='edc':
        folder_prefix = 'SDM_EDC'
        dyn_general.update( { "decoherence_times_type":1, "decoherence_C_param": 1.0, "decoherence_eps_param": args.eps } )  # EDC + default params
    elif dectime=='schw2':
        folder_prefix = 'SDM_SCHW2'
        dyn_general.update( { "decoherence_times_type":3, "schwartz_decoherence_inv_alpha":A } ) # Schwartz version 2
    elif dectime=='gu_franco':
        folder_prefix = 'SDM_Gu_Franco' 
        dyn_general.update({"Temperature": args.temperature, "reorg_energy": args.reorg_energy, "decoherence_times_type": 5 })
# ============ DISH_REV23
elif mqcmethod==6:
    dyn_general.update({ "decoherence_algo":7}) # new DISH
    if dectime=='edc':
        folder_prefix = 'DISH_REV23_EDC'
        dyn_general.update( { "decoherence_times_type":1, "decoherence_C_param": 1.0, "decoherence_eps_param": args.eps } )  # EDC + default params
    elif dectime=='schw1':
        folder_prefix = 'DISH_REV23_SCHW1'
        dyn_general.update( { "decoherence_times_type":2, "schwartz_decoherence_inv_alpha":A } ) # Schwartz version 1
    elif dectime=='schw2':
        folder_prefix = 'DISH_REV23_SCHW2'
        dyn_general.update( { "decoherence_times_type":3, "schwartz_decoherence_inv_alpha":A } ) # Schwartz version 2
    elif dectime=='gu_franco':
        folder_prefix = 'DISH_REV23_Gu_Franco'
        dyn_general.update({"Temperature": args.temperature, "reorg_energy": args.reorg_energy, "decoherence_times_type": 5 })
#=========== Phase correction of SSY =================
if args.ssy:
    dyn_general.update({"do_ssy": 1}) 
else:
    dyn_general.update({"do_ssy": 0}) 
# ============ The TSH method used for the MQC method 
if tshmethod=='fssh':    
    dyn_general.update({"tsh_method": 0}) # FSSH
elif tshmethod=='fssh2':    
    dyn_general.update({"tsh_method": 7}) # FSSH2
elif tshmethod=='gfsh':
    dyn_general.update({"tsh_method": 9}) # GFSH original

############################################################
### 3. Choose the initial conditions: Nuclear and Electronic
############################################################
# Nuclear initial conditions 
beta = model_params["beta"]
_nosc = ndof = model_params["num_osc"]
_nst = NSTATES 

# The minimum of the position and the Wigner sampling momenta
q_init = []
p_init = []
for i in range(len(freqs)):
    r_i = g_vals[i]/(model_params["mass"][i]*freqs[i]*freqs[i])
    p_i = np.sqrt(model_params["mass"][i]*freqs[i])
    q_init.append(r_i)
    p_init.append(p_i)
    print(f'for frequency {i}, the initial conditions are R={r_i}, P={p_i}')

# 1/2 * mass * beta * w**2 = 1/2 * (1/sigma_q**2) => sigma_q = 1/(w * sqrt(mass * beta))
# 1/2 * (beta/mass) = 1/2 * (1/sigma_p**2) => sigma_p = 1/(sqrt(beta/mass))
# force_constant = mass * w**2
nucl_params = { "ndof":ndof,
                #"q":[0.0 for i in range(ndof)],
                # Average of the position for each dof
                "q": q_init,
                # Average of the momentum for each dof
                "p": [0.0 for i in range(ndof)],
                "mass":list(model_params["mass"]),
                "force_constant":[ model_params["mass"][i]*model_params["omega"][0][i]**2 for i in range(_nosc) ],
                "q_width":[ 1.0/( np.sqrt(model_params["mass"][i] * beta) * model_params["omega"][0][i]) for i in range(_nosc)],
                #"q_width":[ 1.0/np.sqrt(2*model_params["mass"][i] * model_params["omega"][0][i]) for i in range(_nosc)],
                "p_width":[ 1.0/( np.sqrt(beta/model_params["mass"][i])) for i in range(_nosc)],
                #"p_width":[ np.sqrt(model_params["mass"][i] * model_params["omega"][0][i] / 2 ) for i in range(_nosc)],
                "init_type": 4
              }
# Electronic initial conditions - start on the second diabatic state (starting from 0 = ground)
##########
istate = args.istate
##########
istates = [0.0 for i in range(NSTATES)]
istates[istate] = 1.0     
print(istates)
elec_params = {"verbosity":2, "init_dm_type":0,"ndia":NSTATES, "nadi":NSTATES, 
               "rep": args.rep, "init_type":0, "istates":istates, "istate":istate     }


def compute_model(q, params, full_id):
    model = params["model"]
    res = None
    if model==0:
        res = GLVC.GLVC(q, params, full_id)
    else:
        pass
    return res


def potential(q, params):
    full_id = Py2Cpp_int([0,0]) 
    
    return compute_model(q, params, full_id)

if mqcmethod==4:
    # Set a Gaussian for the quantum momentum calculation
    WP_W = MATRIX(ndof,1); 
    for idof in range(ndof):
        # You can test a bunch of width for this Gaussian, and this is not the case that this width must follow the initial wavepacket width.
        # Sometimes this value is set to the width from nuclear sampling, but sometimes other options are taken such as \sigma / 10.0
        # The smaller width, the larger decoherencce, that is, decoherence correction \propto \sigma^{-2}
        WP_W.set(idof,0, nucl_params["q_width"][idof] )
    
    dyn_general.update({"wp_width":WP_W*args.wpwidth, "coherence_threshold": 0.01 })
    
    #dyn_general["properties_to_save"] += ["p_quant", "VP", "q_aux", "p_aux", "nab_phase", "wp_width"]   # This is usually for debugging. Add this if you are interested in checking the XF variables.


#============ Types of surface hopping acceptance and momenta rescaling opntions =================

############################################################
### 4. Run the calculations and save the results into the HDF5 files
############################################################
x1 = time.time()
dyn_params = dict(dyn_general)

# Setting up the folders names
pref = tshmethod + F"_{folder_prefix}_ntraj_{args.ntraj}_iter_{args.iter}_dt_{args.dt}"
if dectime=='edc':
    pref = pref + F'_eps_param_{args.eps}'
elif 'schw' in dectime:
    pref = pref + F'_A_{args.A}'
elif mqcmethod==4:
    pref = pref + F'_wpwidth_scale_{args.wpwidth}'
elif mqcmethod==8 or 'gu_franco' in dectime:
    pref = pref + F'_temperature_{args.temperature}_reorg_energy_{args.reorg_energy}'

if args.ssy:
    pref = pref + '_SSY'

pref = pref + f'_istate_{args.istate}_rep_{args.rep}_elec_int_{args.elec_int}'

dyn_params.update({ "prefix":pref, "prefix2":pref })
dyn_general.update({ "prefix":pref, "prefix2":pref })
print(F"Computing {pref}")    

rnd = Random()
res = tsh_dynamics.generic_recipe(dyn_general, compute_model, model_params, elec_params, nucl_params, rnd)
print(F'ELAPSED TIME for {pref}', time.time()-x1,' seconds') 


