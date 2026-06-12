#!/bin/bash -l
# Standard output and error:
#SBATCH -D ./
# Job name
#SBATCH -J job_name
#
#SBATCH --ntasks=1
#SBATCH --constraint="apu"
#
# --- default case: use a single GPU on a shared node ---
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --mem=110000
#
#SBATCH --time=24:00:00

# Load compiler and MPI modules
module purge
module load gcc/15 openmpi/5.0 rocm/7.2 openmpi_gpu/5.0 gromacs/2026.1
srun gmx_mpi mdrun  -deffnm prod -cpi prod.cpt -pin on  -ntomp 24

