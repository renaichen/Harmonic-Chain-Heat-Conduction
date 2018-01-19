#! /bin/bash
#$ -N Renai-muti-tests
#$ -V
#$ -cwd
#$ -pe openmp 2

python multi_traj_test.py
