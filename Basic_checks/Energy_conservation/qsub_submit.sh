#! /bin/bash
#$ -N Renai-muti
#$ -V
#$ -cwd
#$ -pe openmp 2

python3 harmonic_heat_multi_traj.py
