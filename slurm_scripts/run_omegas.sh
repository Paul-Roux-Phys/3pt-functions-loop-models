#!/bin/bash
#
#SBATCH -p C-Short   # partition 
#SBATCH -N 1 # number of nodes 
#SBATCH -n 16 # number of cores 
#SBATCH --mem 128 # memory pool for all cores 
#SBATCH -t 1-00:00 # time (D-HH:MM) 
#SBATCH -o output/omega_11_10_11 # STDOUT 
#SBATCH -e error/omega_11_10_11 # STDERR 

cd ../On_loops/cpp_program
julia -t auto omega_3pt.jl