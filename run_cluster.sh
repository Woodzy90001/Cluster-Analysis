#!/bin/bash --login
#SBATCH --account=p87
#SBATCH --partition=work
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=PIMD-test

lammpsexe=/software/projects/p87/cwood1/LAMMPS/mylammps/PIMD-build/lmp

srun -N 1 -n 1 -c 1 -m block:block:block ${lammpsexe} -in PIMD_$1.in > PIMD_$1.out