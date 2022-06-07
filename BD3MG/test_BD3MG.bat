#!/bin/bash
#SBATCH -p cpu_med
#SBATCH --time=04:00:00
#SBATCH --job-name=JobCVN
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1



::Load necessary modules
module purge
module load anaconda3/2020.02/gcc-9.2.0

::Run python script
echo "Running python script"

python3 BD3MG/test_BD3MG.py


date