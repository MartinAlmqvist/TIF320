#!/usr/bin/env bash
#SBATCH -A C3SE2023-2-18 # Project
#SBATCH -J GPAWTask1MartinLovisa # Name of the job
#SBATCH -N 1 # Use 1 node
#SBATCH -n 16 # 
#SBATCH -t 01:00:00 # Maximum time
#SBATCH -o std.out # stdout goes to this file
#SBATCH -e err.out # stderr goes to this file

# Unload all modules, to make sure no incompatibilities arise
module purge

# Load desired modules
module load GPAW/22.8.0-foss-2022a
module load ASE/3.22.1-foss-2022a

# Run program
srun gpaw python test.py
