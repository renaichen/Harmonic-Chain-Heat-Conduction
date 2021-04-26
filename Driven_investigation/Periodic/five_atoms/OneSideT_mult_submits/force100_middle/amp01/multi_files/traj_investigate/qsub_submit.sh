#! /bin/bash
#$ -N Renai-muti
#$ -V
#$ -cwd
#$ -pe openmp 12

python3 harmonic_heat_multi_traj.py $1
# python3 harmonic_heat_multi_traj.py
