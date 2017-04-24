#!/bin/bash

#$ -N srkIm
#$ -q math
#$ -pe openmp 64
#$ -cwd            		# run the job out of the current directory
#$ -m beas
#$ -ckpt blcr
#$ -o output/
#$ -e output/
#$ -l mem_free=490G
module load Mathematica/10.0
math -noprompt -run "<<FullImplicitSolve.m"
