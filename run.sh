#!/bin/bash

#SBATCH --job-name=ecmc_sc
#SBATCH --output=%x.o
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --partition=cpu_med

#Source necessary modules
source modules_load.sh

# Run MPI script
srun --unbuffered build/gauge_single_core inputs/testsc.txt
