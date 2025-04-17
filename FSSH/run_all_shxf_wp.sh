#!/bin/bash
for j in 0.2 0.3 0.4 ; do
    echo "Submitting jobs for j $j ##############"
    for i in $(seq 0 49); do
        sed -i "s/python.*/python spin_boson_2_states_shxf.py --iter $i --tshmethod 8 --ntraj 100 --nsteps 2100000 --wpwidth $j /g" submit_template.slm
        echo "tshmethod 8 iter $i ntraj 100 nsteps 2100000 wpwidth $j"
        sbatch submit_template.slm
    done
done
