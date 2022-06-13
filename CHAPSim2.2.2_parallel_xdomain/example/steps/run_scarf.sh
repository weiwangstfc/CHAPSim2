#!/bin/sh -l
#SBATCH -p scarf ## which queue to run n
#SBATCH -J test
#SBATCH -o std%J.out
#SBATCH -e std%J.err
## #SBATCH --exclusive  ## exclusive node use, to avoid sharing with others. Only for parallel jobs.
#SBATCH -t 0:10:00
###SBATCH -C "[scarf15|scarf16|scarf17|scarf18]"
#SBATCH -n 1 ## number of processors to run on

# load module
#module load intel_mpi

# cd $LS_SUBCWD

# count how many processors are allocated
NP=0
for TOKEN in $LSB_HOSTS
do
   ((NP++))
done

mpirun -np $NP ../bin/CHAPSim <input_chapsim.ini  ## > OUTPUT_$(date +%Y-%m-%d_%H.%M).log
