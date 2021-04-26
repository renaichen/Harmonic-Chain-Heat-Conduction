#! bin/bash

valuable=()

# for j in `seq 0.001 0.1 20.0`
# do
#     valuable+=($j)
# done

for j in `seq 20.1 0.1 50.0`
do
    valuable+=($j)
done

# valuable+=(0.1 0.5 1.0 2)

len=${#valuable[@]}

for i in `seq 1 $len`
do
    # qsub ~/Harmonic-Chain-Heat-Conduction/diatomic_omega_pin_toLndCompare/src/qsub_dia.sh ${xleft[$i-1]} 
    qsub qsub_submit.sh ${valuable[$i-1]} 
   sleep 60 ## seconds
done

