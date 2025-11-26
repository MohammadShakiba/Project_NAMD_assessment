# Project NAMD assessment


This repository contains the files for assessing different mixed-quantum classical methods. 
For trajectory surface hopping (TSH) methods, we adopt FSSH, FSSH-2, and GFSH methods. The decoherece correction methods are
SDM, DISH, ID-A, BCSH, and SHXF. For SDM and DISH methods, difference decoherence time methods are adopted: Schwartz, EDC, and Gu-Franco.
The pythonic file that contain the code for running TSH calculations is `spin_boson_2_states.py`. 
Each have different properties that is passed by the user and will result in an overall **66** different MQC methods.

```python
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
```

Since we aim to run multiple runs, we need to submit these files on different compute nodes. `submit_template.slm` is a sample file used to submit a job on UB CCR. The script to submit all jobs are `run_all*.sh` files. The description of different methods to be used are shown in the pythonic file. 


The results can be plotted using the `plot_results_no_ssy.ipynb` and `plot_results_with_ssy.ipynb` Jupyter files.


For more details on how to run these calculations see [this file](1_NAMD/README.md).
