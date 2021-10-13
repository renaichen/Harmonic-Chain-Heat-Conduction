#! bin/bash

valuable=()
# for i in `seq 1 1 5`
# do
#     valuable+=($i)
# done

for j in `seq 0.03756 0.07512 1.1268`
do
    valuable+=($j)
done

# valuable+=(0.05 1)
# valuable+=(0.5)

len=${#valuable[@]}

for i in `seq 1 $len`
do
    # qsub ~/Harmonic-Chain-Heat-Conduction/diatomic_omega_pin_toLndCompare/src/qsub_dia.sh ${xleft[$i-1]} 
    qsub qsub_submit.sh ${valuable[$i-1]} 
    echo $i
   # sleep 18m
done

