#!/bin/bash
#for j in 6 14 15 17; do
#for j in 6 14 15 17 19 20; do
for j in 19 20; do
    echo "Submitting jobs for j $j ##############"
    for i in $(seq 0 49); do
        for A in 0.01 0.1 1.0 10.0 20.0; do
            sed -i "s/python.*/python spin_boson_2_states.py --iter $i --tshmethod $j --ntraj 100 --nsteps 2100000 --A $A/g" submit_template.slm
            echo "tshmethod $j iter $i ntraj 100 nsteps 2100000"
            sbatch submit_template.slm
        done
    done
done
