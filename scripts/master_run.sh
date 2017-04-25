#!/bin/bash
for i in {1..50}
  sbatch scripts/slurm_run.sh -6 3
done
