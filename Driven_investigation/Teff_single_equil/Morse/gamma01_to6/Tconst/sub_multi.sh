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

valuable+=(0.1 0.2 0.3 0.4 0.5 0.7 1 1.4 1.8 2.3 3 4 5 6)
# valuable+=(1)

len=${#valuable[@]}

for i in `seq 1 $len`
do
    # qsub ~/Harmonic-Chain-Heat-Conduction/diatomic_omega_pin_toLndCompare/src/qsub_dia.sh ${xleft[$i-1]} 
    qsub qsub_submit.sh ${valuable[$i-1]} 
   # sleep 18m
done

