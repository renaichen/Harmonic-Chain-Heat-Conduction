#! /bin/bash
#$ -N Debye-Classical
#$ -V
#$ -cwd
#$ -pe openmp 8

python3 Relaxation_harmonic_Debye_Helium.py $1
# python3 Relaxation_harmonic_Debye_Helium.py