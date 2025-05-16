#!/bin/bash
ntraj=2000
nstep=105000
# For no decoherence, BCSH, and ID-A
for tshmethod in fssh fssh2 gfsh; do 
  for mqcmethod in 1 2 3 ; do
    sed -i "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --ntraj $ntraj --nsteps $nstep --dt 20.0 --ssy True /g" submit_template.slm
    echo "python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --ntraj $ntraj --nsteps $nstep --dt 20.0 /g" submit_template.slm
    sbatch submit_template.slm 
  done
done

# For SHXF with different wave packet width scale
for tshmethod in fssh fssh2 gfsh; do
  #for wpwidth in 0.05 0.1 0.25 0.5 1.0 2.0 4.0 ; do
  for wpwidth in 0.0001 0.001 0.005 ; do
    sed -i "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod 4 --ntraj $ntraj --nsteps $nstep --dt 20.0 --wpwidth $wpwidth --ssy True /g" submit_template.slm
    echo "python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod 4 --ntraj $ntraj --nsteps $nstep --dt 20.0 --wpwidth $wpwidth /g" submit_template.slm
    sbatch submit_template.slm
  done
done
 
 # For SDM, MFSD, and DISH_REV23 with EDC decoherence time
 for tshmethod in fssh fssh2 gfsh; do
   for mqcmethod in 5 6 7 ; do
     for eps in 0.01 0.05 0.1 0.2 0.4 1.0 5.0 10.0 20.0 40.0 80.0; do
       sed -i "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime edc --ntraj $ntraj --nsteps $nstep --dt 20.0 --eps $eps --ssy True /g" submit_template.slm
       echo "python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime edc --ntraj $ntraj --nsteps $nstep --dt 20.0 --eps $eps /g" submit_template.slm
       sbatch submit_template.slm
     done 
   done
 done
 
 # For SDM, MFSD, and DISH_REV23 with Schwartz 1 decoherence time
 for tshmethod in fssh fssh2 gfsh; do
   for mqcmethod in 5 6 7 ; do
     for A in 0.0001 0.001 0.01 0.1 1.0 10.0 ; do
       sed -i "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime schw1 --ntraj $ntraj --nsteps $nstep --dt 20.0 --A $A --ssy True/g" submit_template.slm
       echo "python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime schw1 --ntraj $ntraj --nsteps $nstep --dt 20.0 --A $A /g" submit_template.slm
       sbatch submit_template.slm
     done
   done
 done


# For SDM, MFSD, and DISH_REV23 with Schwartz 2 decoherence time
  for tshmethod in fssh fssh2 gfsh; do
    for mqcmethod in 5 6 7 ; do
      for A in 0.0001 0.001 ; do
      #for A in 0.01 0.1 1.0 10.0 20.0; do
        sed -i "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime schw2 --ntraj $ntraj --nsteps $nstep --dt 20.0 --A $A --ssy True/g" submit_template.slm
        echo "python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime schw2 --ntraj $ntraj --nsteps $nstep --dt 20.0 --A $A /g" submit_template.slm
        sbatch submit_template.slm
      done
    done
  done

# For SDM, MFSD, and DISH_REV23 with Gu-Franco decoherence time
for tshmethod in fssh fssh2 gfsh; do
  for mqcmethod in 5 6 7 ; do
    sed -i "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime gu_franco --ntraj $ntraj --nsteps $nstep --dt 20.0 --ssy True/g" submit_template.slm
    echo "python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod $mqcmethod --dectime gu_franco --ntraj $ntraj --nsteps $nstep --dt 20.0 /g" submit_template.slm
      sbatch submit_template.slm
  done
done

# For the new decoherence algorithm proposed by Alexey with Gu-Franco decoherence time 
for tshmethod in fssh fssh2 gfsh; do
  sed -i "s/python.*/python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod 8 --dectime gu_franco --ntraj $ntraj --nsteps $nstep --dt 20.0 --ssy True/g" submit_template.slm
  echo "python spin_boson_2_states.py --iter 0 --tshmethod $tshmethod --mqcmethod 8 --dectime gu_franco --ntraj $ntraj --nsteps $nstep --dt 20.0 /g" submit_template.slm
  sbatch submit_template.slm
done


