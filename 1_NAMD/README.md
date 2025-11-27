# Spin–Boson (2-State) Computational Workflow

This small workflow automates large sets of spin–boson surface-hopping simulations on an HPC cluster using a single Python driver script and two shell wrappers.

## Files

### `spin_boson_2_states.py`
Main simulation driver.  
It:

- Parses command-line arguments (e.g. `--ntraj`, `--tshmethod`, `--mqcmethod`, `--dectime`, `--A`, `--eps`, `--reorg_energy`, `--temperature`, `--ssy`, `--wpwidth`, etc.).
- Sets up a 2-state GLVC spin–boson model with 10 harmonic bath modes and their couplings.
- Builds the general dynamics dictionary (`dyn_general`) controlling:
  - Number of trajectories (`ntraj`), time step (`dt`), number of steps (`nsteps`),
  - TSH algorithm (`tshmethod`: `fssh`, `fssh2`, `gfsh`),
  - Decoherence / MQC correction (`mqcmethod`: 1–6) and decoherence-time scheme (`dectime`: `edc`, `schw1`, `schw2`, `gu_franco`),
  - Phase correction (`do_ssy`), electronic integrator (`elec_int`), representation (`rep`), etc.
- Constructs initial nuclear and electronic conditions and calls the Libra surface-hopping driver (`tsh_dynamics.generic_recipe`).
- Automatically builds an output folder prefix from the chosen methods and parameters (e.g. `fssh_BCSH_ntraj_1000_iter_0_dt_5.0_...`) and writes results (HDF5 + text) there.

You can also run this script directly for one-off calculations (see **Running a single calculation** below).

---

### `run_all.sh`
Batch launcher **without** SSY phase correction.

It is a pure Bash script that:

- Defines common simulation parameters at the top:
  - `ntraj`, `dt`, `nstep`, `elec_int`, `rep`, `istate`.
- Loops over:
  - TSH methods: `fssh`, `fssh2`, `gfsh`,
  - MQC methods: `mqcmethod = 1–6` (with subsets depending on the block),
  - Decoherence-time options (`dectime = edc, schw1, schw2, gu_franco`),
  - Additional parameters such as `eps`, `A`, and `reorg_energy` for the different decoherence models,
  - Wave packet width scaling `wpwidth` for `mqcmethod = 4` (SHXF).
- For each parameter combination, it:
  1. Replaces the `python ...` line in `submit_template.slm` **in place** using `sed -i`, inserting the appropriate `python spin_boson_2_states.py ...` command with all arguments.
  2. Submits the job via `sbatch submit_template.slm`.

This script **does not** pass `--ssy`, so `spin_boson_2_states.py` runs with its default `ssy=False` (no SSY correction).

---

### `run_all_ssy.sh`
Batch launcher **with** SSY phase correction turned on.

It has the same structure and parameter sweeps as `run_all.sh`, but in every `sed` substitution the command includes `--ssy True`. This sets `args.ssy = True` inside `spin_boson_2_states.py`, which flips `do_ssy` on and activates the SSY phase correction in the dynamics.

Use this script when you explicitly want to generate SSY-corrected trajectories for the same parameter sets as `run_all.sh`.

---

## Requirements

- **Python environment** with:
  - `liblibra_core` and `libra_py` (Libra library),
  - `numpy`, `h5py`, `matplotlib`, `argparse`.
- **HPC scheduler**: the scripts assume **SLURM**, using:
  - `sbatch submit_template.slm`
- A Slurm submit template file:
  - `submit_template.slm` must exist and contain a line with a `python ...` command to be replaced by `sed`.  
    The `sed -i "s/python.*/python spin_boson_2_states.py .../g"` pattern assumes there is exactly one line starting with `python`.

---

## Running a single calculation

You can run the Python script directly (e.g. on a login or compute node, respecting your cluster policies). A minimal example:

```bash
python spin_boson_2_states.py \
  --iter 0 \
  --ntraj 1000 \
  --nsteps 100000 \
  --dt 1.0 \
  --tshmethod fssh \
  --mqcmethod 1 \
  --istate 1 \
  --rep 1 \
  --elec_int 2

```

---


