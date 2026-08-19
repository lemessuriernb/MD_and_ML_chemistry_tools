#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=ucb-general
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=2
#SBATCH --job-name=cp2k-l1
#SBATCH --constraint=ib
#SBATCH --mail-type=END
#SBATCH --mail-user=nale8203@colorado.edu
#SBATCH --output=cp2k.%j.out
cd $SLURM_SUBMIT_DIR

module purge

module load gcc/11.2.0
module load openmpi/4.1.1
module load cp2k/2023.1

export SLURM_EXPORT_ENV=ALL

NAME=md_pyrene_1w

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}


mpirun -np $SLURM_NTASKS cp2k.psmp -o $NAME.out -i $NAME.inp
