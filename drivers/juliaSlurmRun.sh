#!/bin/bash
#SBATCH -A uci128
#SBATCH --job-name="jOpt"
#SBATCH --output="output/jOpt.%j.%N.out"
#SBATCH --partition=gpu-shared
#SBATCH --gres=gpu:2
#SBATCH --nodes=1
#SBATCH --export=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user=crackauc@uci.edu
#SBATCH --ntasks-per-node=12
#SBATCH -t 0:10:00
export SLURM_NODEFILE=`generate_pbs_nodefile`
/home/crackauc/julia-a2f713dea5/bin/julia --machinefile $SLURM_NODEFILE /home/crackauc/projectCode/ImprovedSRK/Optimization/cometDriver.jl
