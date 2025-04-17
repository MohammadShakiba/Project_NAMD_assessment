#!/bin/bash
for j in 7 11 16 18; do
    echo "Submitting jobs for j $j ##############"
    for i in $(seq 0 49); do
        for eps in 0.01 0.05 0.1 0.2 0.4; do
            sed -i "s/python.*/python spin_boson_2_states.py --iter $i --tshmethod $j --ntraj 100 --nsteps 2100000 --eps $eps/g" submit_template.slm
            echo "tshmethod $j iter $i ntraj 100 nsteps 2100000"
            sbatch submit_template.slm
        done
    done
done
