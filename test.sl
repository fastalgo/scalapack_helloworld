#!/bin/bash -l

#SBATCH --partition debug
#SBATCH --nodes 3
#SBATCH --time=00:03:00

cd $SLURM_SUBMIT_DIR
srun -n 3 ./execute;
