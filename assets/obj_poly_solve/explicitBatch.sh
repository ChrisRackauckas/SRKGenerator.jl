#!/bin/bash

#$ -N srkEx
#$ -q math
#$ -pe openmp 64
#$ -cwd            		# run the job out of the current directory
#$ -m beas
#$ -ckpt blcr
#$ -o output/
#$ -e output/
module load Mathematica/10.0
math -noprompt -run "<<ExplicitSolve.m"
