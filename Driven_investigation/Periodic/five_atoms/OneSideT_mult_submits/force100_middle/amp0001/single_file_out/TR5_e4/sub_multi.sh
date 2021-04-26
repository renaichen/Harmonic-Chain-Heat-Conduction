#! bin/bash

valuable=()

for j in `seq 0.00001 0.1 40.`
do
    valuable+=($j)
done

# for j in `seq 0.01 0.01 0.1`
# do
#     valuable+=($j)
# done

# for j in `seq 0.1 0.1 1.0`
# do
#     valuable+=($j)
# done

# valuable+=(0.5 1.0)

len=${#valuable[@]}

for i in `seq 1 $len`
do
    # qsub ~/Harmonic-Chain-Heat-Conduction/diatomic_omega_pin_toLndCompare/src/qsub_dia.sh ${xleft[$i-1]} 
    qsub qsub_submit.sh ${valuable[$i-1]} 
   sleep 10 ## seconds
done

