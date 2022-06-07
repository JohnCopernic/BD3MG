#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --job-name=JobCVN
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --cpus-per-task=1



# Load necessary modules
module purge
module load anaconda3/2020.02/gcc-9.2.0

# Run python script
echo "Running python script"

python3 BD3MG/v3_3MG.py


date