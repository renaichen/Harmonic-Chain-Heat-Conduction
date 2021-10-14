#! /bin/bash
#$ -N Debye-Classical
#$ -V
#$ -cwd
#$ -pe openmp 48

python3 HC_Debye_Helium.py $1
