#! bin/bash

xleft=()
# for i in `seq 1 1 5`
# do
#     xleft+=($i)
# done

# for j in `seq 0.1 0.1 2.0`
# do
#     xleft+=($j)
# done

# xleft+=(0 0.001 0.01 0.1 0.5 0.8 1 1.5 2 5 10 30 50 100 200 300 400 500)
# xleft+=(0.001 0.01 0.1 0.5 0.8 1.5 5 10 100 200)
xleft+=(0.01 0.1 0.5 0.8 1.5 5 10 100)
# xleft+=(50 100)
# xleft=10

# xright=()
# for j in `seq 54.5 -0.1 53.3`
# do
    # xright+=($j)
# done

len=${#xleft[@]}

for i in `seq 1 $len`
do
    qsub ~/Harmonic-Chain-Heat-Conduction/diatomic_omega_pin_toLndCompare/src/qsub_dia.sh ${xleft[$i-1]} 
    # qsub qsub_single.sh ${xleft[$i-1]} 
   # sleep 18m
#    bash justecho.sh ${xleft[$i-1]} ${xright[$i-1]}
done

