#! bin/bash

valuable=()
# for i in `seq 1 1 5`
# do
#     xleft+=($i)
# done

# for j in `seq 0.1 0.1 3.0`
# do
#     valuable+=($j)
# done

# valuable+=(6 9)
valuable+=(9)

len=${#valuable[@]}

for i in `seq 1 $len`
do
    # qsub ~/Harmonic-Chain-Heat-Conduction/diatomic_omega_pin_toLndCompare/src/qsub_dia.sh ${xleft[$i-1]} 
    qsub qsub_submit.sh ${valuable[$i-1]} 
   # sleep 18m
done

