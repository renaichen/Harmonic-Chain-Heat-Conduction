#! bin/bash

valuable=()

# for j in `seq 0.01 0.01 0.1`
# do
#     valuable+=($j)
# done

valuable+=(0.01 0.1 1. 10.)

len=${#valuable[@]}

for i in `seq 1 $len`
do
    # qsub ~/Harmonic-Chain-Heat-Conduction/diatomic_omega_pin_toLndCompare/src/qsub_dia.sh ${xleft[$i-1]} 
    qsub qsub_submit.sh ${valuable[$i-1]} 
   # sleep 60 ## seconds
done

