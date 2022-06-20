#!/bin/sh
# slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=test
#SBATCH --time 0:10:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=12
#SBATCH --cpus-per-task=1

# Place [budget code] below with your budget code
#SBATCH --acount=c01-eng
#SBATCH --partition=standard
#SBATCH --qos=standard

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1

# Launch the parallel job
#   Using 512 MPI processes and 128 MPI processes per node
#   srun picks up the distribution from the sbatch options

srun --distribution=block:block --hint=nomultithread ./../../bin/CHAPSim < input_chapsim.ini



