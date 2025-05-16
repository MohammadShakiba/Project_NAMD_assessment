# Project_NAMD_assessment


This repository contains the files for assessing different mixed-quantum classical methods. For trajectory surface hopping (TSH) methods, we adopt FSSH, FSSH-2, and GFSH methods. The pythonic file that contain the code for running TSH calculations is `spin_boson_2_states.py`. 
Each have different properties that is passed by user:

```python
parser.add_argument('--iter', default=0, type=int)
parser.add_argument('--ntraj', default=10, type=int)
parser.add_argument('--mqcmethod', default=1, type=int)
parser.add_argument('--tshmethod', default='fssh', type=str)
parser.add_argument('--dectime', default='', type=str)
parser.add_argument('--nsteps', default=100000, type=int)
parser.add_argument('--dt', default=1.0, type=float)
parser.add_argument('--temperature', default=300.0, type=float)
parser.add_argument('--A', default=10.0, type=float)
parser.add_argument('--eps', default=0.1, type=float)
parser.add_argument('--ssy', default=False, type=bool)
parser.add_argument('--wpwidth', default=1.0, type=float)
```

Since we aim to run multiple runs, we need to submit these files on different compute nodes. `submit_template.slm` is a sample file used to submit a job on UB CCR. The script to submit all jobs are `run_all*.sh` files. Different methods to be used are shown in the pythonic files. 

