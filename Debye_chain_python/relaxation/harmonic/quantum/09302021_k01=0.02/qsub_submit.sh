#! /bin/bash
#$ -N Debye-Quantum
#$ -V
#$ -cwd
#$ -pe openmp 24

python3 Relaxation_harmonic_Debye_Helium.py $1
# python3 Relaxation_harmonic_Debye_Helium.py
