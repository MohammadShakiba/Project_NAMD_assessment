#!/bin/bash
ntraj=1000
dt=5.0
elec_int=2
nstep=420000
rep=1
istate=1
# For no decoherence, BCSH, and ID-A
for tshmethod in fssh fssh2 gfsh; do 
  for mqcmethod in 1 2 3 ; do
    sed -i "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --ntraj $ntraj --nsteps $nstep --dt $dt --istate $istate --rep $rep --elec_int $elec_int  /g" submit_template.slm
    echo "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --ntraj $ntraj --nsteps $nstep --dt $dt --istate $istate --rep $rep --elec_int $elec_int  /g" 
    sbatch submit_template.slm 
  done
done
# For SHXF with different wave packet width scale, mqcmethod = 4
for tshmethod in fssh fssh2 gfsh; do
  for wpwidth in 0.05 0.1 0.25 0.5 1.0 2.0 4.0 ; do
    sed -i "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod 4 --ntraj $ntraj --nsteps $nstep --dt $dt --wpwidth $wpwidth --istate $istate --rep $rep --elec_int $elec_int  /g" submit_template.slm
    echo "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod 4 --ntraj $ntraj --nsteps $nstep --dt $dt --wpwidth $wpwidth --istate $istate --rep $rep --elec_int $elec_int  /g"
    sbatch submit_template.slm
  done
done
# For SDM and DISH_REV23 with EDC decoherence time, mqcmethod = 5 or 6
for tshmethod in fssh fssh2 gfsh; do
  for mqcmethod in 5 6 ; do
    for eps in 0.01 0.05 0.1 0.2 0.4 1.0 5.0 10.0 50.0 100.0 500.0 1000.0; do
      sed -i "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime edc --ntraj $ntraj --nsteps $nstep --dt $dt --eps $eps --istate $istate --rep $rep --elec_int $elec_int   /g" submit_template.slm
      echo "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime edc --ntraj $ntraj --nsteps $nstep --dt $dt --eps $eps --istate $istate --rep $rep --elec_int $elec_int   /g"
      sbatch submit_template.slm
    done 
  done
done
# For DISH_REV23 with Schwartz 1 decoherence time, mqcmethod = 6
for tshmethod in fssh fssh2 gfsh; do
  for mqcmethod in 6 ; do
    for A in 0.0001 0.001 0.01 0.1 1.0 10.0 100.0 1000.0 10000.0 100000.0 1000000.0 10000000.0 100000000.0 1000000000.0; do
      sed -i "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime schw1 --ntraj $ntraj --nsteps $nstep --dt $dt --A $A --istate $istate --rep $rep --elec_int $elec_int  /g" submit_template.slm
      echo "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime schw1 --ntraj $ntraj --nsteps $nstep --dt $dt --A $A --istate $istate --rep $rep --elec_int $elec_int  /g"
      sbatch submit_template.slm
    done
  done
done

# For SDM and DISH_REV23 with Schwartz 2 decoherence time, mqcmethod = 5 or 6
for tshmethod in fssh fssh2 gfsh; do
  for mqcmethod in 5 6 ; do
    for A in 0.0001 0.001 0.01 0.1 1.0 10.0 100.0 1000.0 10000.0 100000.0 1000000.0 10000000.0 100000000.0 1000000000.0; do
      sed -i "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime schw2 --ntraj $ntraj --nsteps $nstep --dt $dt --A $A --istate $istate --rep $rep --elec_int $elec_int  /g" submit_template.slm
      echo "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime schw2 --ntraj $ntraj --nsteps $nstep --dt $dt --A $A --istate $istate --rep $rep --elec_int $elec_int  /g"
      sbatch submit_template.slm
    done
  done
done

# For SDM and DISH_REV23 with Gu-Franco decoherence time
for tshmethod in fssh fssh2 gfsh; do
  for mqcmethod in 5 6 ; do
    for reorg_energy in 0.000000000125 0.0000125 0.000125 0.0125 1.25 12.5 125.0 1250.0 ; do
      sed -i "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime gu_franco --ntraj $ntraj --nsteps $nstep --dt $dt --istate $istate --rep $rep --elec_int $elec_int --reorg_energy $reorg_energy /g" submit_template.slm
      echo "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime gu_franco --ntraj $ntraj --nsteps $nstep --dt $dt --istate $istate --rep $rep --elec_int $elec_int --reorg_energy $reorg_energy /g" 
      sbatch submit_template.slm
    done
  done
done



