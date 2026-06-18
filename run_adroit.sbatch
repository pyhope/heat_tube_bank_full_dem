#!/bin/bash
#SBATCH --job-name=he3-dem-base
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
#SBATCH --constraint=ice
#SBATCH --output=slurm-%j.out

LAMMPS_BIN=/path/to/your/lammps/executable

module purge
module load modules_for_lammps

srun "$LAMMPS_BIN" -in in.heat_tube_bank
