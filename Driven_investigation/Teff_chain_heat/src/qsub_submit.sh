#! /bin/bash
#$ -N Renai-muti
#$ -V
#$ -cwd
#$ -pe openmp 24

# python3 ../src/harmonic_heat_multi_traj.py
python3 ../src/harmonic_multi_Teff.py
