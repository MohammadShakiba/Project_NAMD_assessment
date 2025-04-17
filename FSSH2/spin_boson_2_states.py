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

from recipes import fssh, fssh2, fssh3, gfsh, ida, sdm, mfsd, shxf, bcsh, sdm_schw1, sdm_schw2, shxf_edc, shxf_schw1, shxf_schw2, mfsd_edc, mfsd_schw2, dish_rev2023_edc, dish_rev2023_schw1, dish_rev2023_schw2

import libra_py.models.GLVC as GLVC

parser = argparse.ArgumentParser()
parser.add_argument('--iter', default=0, type=int)
parser.add_argument('--ntraj', default=10, type=int)
parser.add_argument('--tshmethod', default=1, type=int)
parser.add_argument('--nsteps', default=100000, type=int)
parser.add_argument('--dt', default=1.0, type=float)
parser.add_argument('--A', default=10.0, type=float)
parser.add_argument('--eps', default=0.1, type=float)
args = parser.parse_args()

s = units.inv_cm2Ha
freqs = np.array([129.743484,263.643268,406.405874,564.000564,744.784527,
                  961.545250,1235.657107,1606.622430,2157.558683,3100.000000])*s
freqs = freqs.tolist()
g_vals = [0.000563,0.001143,0.001762,0.002446,0.003230,0.004170,0.005358,0.006966,0.009356,0.013444]
eps = 0.02
v0 = 0.001
m = 1836.0 # 1 amu
temperature = 300.0 # K
model_params = {"model": 0, "model0": 0,"nstates": 2, "num_osc": 10, "coupling_scaling": [1.0, -1.0], 
                "omega": [freqs]*2, "coupl": [g_vals]*2, "mass": [m]*10, "Ham": [[eps,-v0],[-v0,-eps]],
                "beta": 1/(units.kB*temperature)}


############################################################
##### 2. Choose the Nonadiabatic Dynamics Methodology 
############################################################
NSTATES = 2
dyn_general = { "nsteps":args.nsteps, "ntraj":args.ntraj, "nstates":NSTATES,
                "dt": 1.0, "num_electronic_substeps":1, "isNBRA":0, "is_nbra":0,
                "progress_frequency":0.1, "which_adi_states":range(NSTATES), "which_dia_states":range(NSTATES),      
                "mem_output_level":3,
                #"properties_to_save":[ "timestep", "time","se_pop_adi", "se_pop_dia", "sh_pop_adi", "sh_pop_dia"],
                "properties_to_save":[ "sh_pop_adi"],
                "prefix":"adiabatic_md", "prefix2":"adiabatic_md"
              }

tshmethod = args.tshmethod
if tshmethod==1:
    fssh.load(dyn_general) # FSSH
    tsh_method = 'FSSH'
elif tshmethod==2:
    fssh2.load(dyn_general) # FSSH2
    tsh_method = 'FSSH2'
elif tshmethod==3:
    bcsh.load(dyn_general) # BCSH
    tsh_method = 'BCSH'
elif tshmethod==4:
    gfsh.load(dyn_general) # GFSH
    tsh_method = 'GFSH'
elif tshmethod==5:
    ida.load(dyn_general) # ID-A
    tsh_method = 'ID-A'
elif tshmethod==6:
    mfsd.load(dyn_general) # MFSD SCHW1
    tsh_method = 'MFSD'
    A = MATRIX(2,2); A.set(0, 0, args.A); A.set(1,1, args.A)
    dyn_general.update( { "decoherence_times_type":2, "schwartz_decoherence_inv_alpha":A } ) # Schwartz version 1 - only this option is possible for MFSD!

elif tshmethod==7:
    sdm.load(dyn_general) # sdm EDC
    tsh_method = 'SDM'
    dyn_general.update( { "decoherence_times_type":1, "decoherence_C_param": 1.0, "decoherence_eps_param": args.eps } )  # EDC + default params
elif tshmethod==8:
    shxf.load(dyn_general) # shxf, infinite decoherence times
    tsh_method = 'SHXF'
elif tshmethod==9:
    gfsh.load(dyn_general) # GFSH
    tsh_method = 'GFSH'
elif tshmethod==10:
    fssh3.load(dyn_general) # FSSH3
    tsh_method = 'FSSH3'
elif tshmethod==11:
    shxf_edc.load(dyn_general) # SHXF EDC
    tsh_method = 'SHXF_EDC'
    dyn_general.update( { "decoherence_times_type":1, "decoherence_C_param": 1.0, "decoherence_eps_param": args.eps } )  # EDC + default params
elif tshmethod==12:
    shxf_schw1.load(dyn_general) # SHXF SCHW1
    tsh_method = 'SHXF_SCHW1'
elif tshmethod==13:
    shxf_schw2.load(dyn_general) # SHXF SCHW2
    tsh_method = 'SHXF_SCHW2'
elif tshmethod==14:
    sdm_schw1.load(dyn_general) # SDM SCHW1
    tsh_method = 'SDM_SCHW1'
    A = MATRIX(2,2); A.set(0, 0, args.A); A.set(1,1, args.A)
    dyn_general.update( { "decoherence_times_type":2, "schwartz_decoherence_inv_alpha":A } ) # Schwartz version 1 - only this option is possible for MFSD!

elif tshmethod==15:
    sdm_schw2.load(dyn_general) # SDM SCHW2
    tsh_method = 'SDM_SCHW2'
    A = MATRIX(2,2); A.set(0, 0, args.A); A.set(1,1, args.A)
    dyn_general.update( { "decoherence_times_type":3, "schwartz_decoherence_inv_alpha":A } ) # Schwartz version 2

elif tshmethod==16:
    mfsd_edc.load(dyn_general) # MFSD EDC
    tsh_method = 'MFSD_EDC'
    dyn_general.update( { "decoherence_times_type":1, "decoherence_C_param": 1.0, "decoherence_eps_param": args.eps } )  # EDC + default params
elif tshmethod==17:
    mfsd_schw2.load(dyn_general) # MFSD SCHW2
    tsh_method = 'MFSD_SCHW2'
    A = MATRIX(2,2); A.set(0, 0, args.A); A.set(1,1, args.A)
    dyn_general.update( { "decoherence_times_type":3, "schwartz_decoherence_inv_alpha":A } ) # Schwartz version 2

elif tshmethod==18:
    dish_rev2023_edc.load(dyn_general) # DISH rev23 EDC
    tsh_method = 'DISH_REV23_EDC'
    dyn_general.update( { "decoherence_times_type":1, "decoherence_C_param": 1.0, "decoherence_eps_param": args.eps } )  # EDC + default params
elif tshmethod==19:
    dish_rev2023_schw1.load(dyn_general) # DISH rev23 SCHW1
    tsh_method = 'DISH_REV23_SCHW1'
    A = MATRIX(2,2); A.set(0, 0, args.A); A.set(1,1, args.A)
    dyn_general.update( { "decoherence_times_type":2, "schwartz_decoherence_inv_alpha":A } ) # Schwartz version 1 - only this option is possible for MFSD!

elif tshmethod==20:
    dish_rev2023_schw2.load(dyn_general) #DISH rev23 SCHW2
    tsh_method = 'DISH_REV23_SCHW2'
    A = MATRIX(2,2); A.set(0, 0, args.A); A.set(1,1, args.A)
    dyn_general.update( { "decoherence_times_type":3, "schwartz_decoherence_inv_alpha":A } ) # Schwartz version 2


#=========== Phase correction of SSY =================
#dyn_general.update({"do_ssy": 0}) 
dyn_general.update({"do_ssy": 1 }) 
dyn_general.update({"tsh_method": 7}) # FSSH2
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
                "p":[0.0 for i in range(ndof)],
                "mass":list(model_params["mass"]),
                "force_constant":[ model_params["mass"][i]*model_params["omega"][0][i]**2 for i in range(_nosc) ],
                "q_width":[ 1.0/( np.sqrt(model_params["mass"][i] * beta) * model_params["omega"][0][i]) for i in range(_nosc)  ],
                "p_width":[ 1.0/( np.sqrt(beta/model_params["mass"][i])) for i in range(_nosc)  ],
                "init_type": 4
              }
dt = 1.0 # atomic unit
print(F"dt = {dt}")

# Electronic initial conditions - start on the second diabatic state (starting from 0 = ground)
##########
istate = 1
##########
istates = [0.0 for i in range(NSTATES)]
istates[istate] = 1.0     
print(istates)
elec_params = {"verbosity":2, "init_dm_type":0,"ndia":NSTATES, "nadi":NSTATES, 
               "rep":1, "init_type":0, "istates":istates, "istate":istate     }


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

if tshmethod==8:
    # Set a Gaussian for the quantum momentum calculation
    WP_W = MATRIX(ndof,1); 
    for idof in range(ndof):
        # You can test a bunch of width for this Gaussian, and this is not the case that this width must follow the initial wavepacket width.
        # Sometimes this value is set to the width from nuclear sampling, but sometimes other options are taken such as \sigma / 10.0
        # The smaller width, the larger decoherencce, that is, decoherence correction \propto \sigma^{-2}
        WP_W.set(idof,0, nucl_params["q_width"][idof] )
    
    dyn_general.update({"wp_width":WP_W, "coherence_threshold": 0.01 })
    
    #dyn_general["properties_to_save"] += ["p_quant", "VP", "q_aux", "p_aux", "nab_phase", "wp_width"]   # This is usually for debugging. Add this if you are interested in checking the XF variables.


#============ Types of surface hopping acceptance and momenta rescaling opntions =================

############################################################
### 4. Run the calculations and save the results into the HDF5 files
############################################################
x1 = time.time()
dyn_params = dict(dyn_general)
#pref = F"{tsh_method}_ntraj_{args.ntraj}_iter_{args.iter}_A_1.0_eps_param_0.05" 
#pref = F"{tsh_method}_ntraj_{args.ntraj}_iter_{args.iter}_A_5.0_eps_param_0.2" 
if tshmethod in [6,14,15,17,19,20]:
    #pref = F"{tsh_method}_ntraj_{args.ntraj}_iter_{args.iter}_A_{args.A}" 
    pref = F"{tsh_method}_ntraj_{args.ntraj}_iter_{args.iter}_A_{args.A}_SSY" 
elif tshmethod in [18,16,11,7]:
    #pref = F"{tsh_method}_ntraj_{args.ntraj}_iter_{args.iter}_eps_param_{args.eps}" 
    pref = F"{tsh_method}_ntraj_{args.ntraj}_iter_{args.iter}_eps_param_{args.eps}_SSY" 
else:
    #pref = F"{tsh_method}_ntraj_{args.ntraj}_iter_{args.iter}"
    pref = F"{tsh_method}_ntraj_{args.ntraj}_iter_{args.iter}_SSY"
#pref = F"{tsh_method}_ntraj_{args.ntraj}_iter_{args.iter}_A_1.0_eps_param_0.1_SSY" 
dyn_params.update({ "prefix":pref, "prefix2":pref })
print(F"Computing {pref}")    

rnd = Random()
res = tsh_dynamics.generic_recipe(dyn_params, compute_model, model_params, elec_params, nucl_params, rnd)
print('ELAPSED TIME', time.time()-x1,' seconds')


