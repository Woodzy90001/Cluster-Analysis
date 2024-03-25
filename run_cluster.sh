#!/bin/bash --login
#SBATCH --account=group-code
#SBATCH --partition=work
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=PIMD-runs

lammpsexe=lammps-directory

srun -N 1 -n 1 -c 1 -m block:block:block ${lammpsexe} -in PIMD_$1.in > PIMD_$1.out
