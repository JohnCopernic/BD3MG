#!/bin/bash
#SBATCH -p cpu_long
#SBATCH --time=24:00:00
#SBATCH --job-name=BP3MG_long
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --cpus-per-task=1



# Load necessary modules
module purge
module load anaconda3/2020.02/gcc-9.2.0

# Run python script
echo "Running python script"

python3 BD3MG/restore_BP3MG.py


date