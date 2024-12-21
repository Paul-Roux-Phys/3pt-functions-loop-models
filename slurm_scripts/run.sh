#!/bin/bash
#
#SBATCH -p C-Short   # partition 
#SBATCH -N 2 # number of nodes 
#SBATCH -n 16 # number of cores 
#SBATCH --mem 128 # memory pool for all cores 
#SBATCH -t 1-00:00 # time (D-HH:MM) 
#SBATCH -o output/check_2pt_1_0 # STDOUT 
#SBATCH -e error/check_2pt_1_0 # STDERR 

cd ../julia_scripts
julia -t auto check_2pt.jl