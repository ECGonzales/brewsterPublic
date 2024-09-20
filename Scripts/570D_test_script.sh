#!/bin/bash
#SBATCH -J G570D_test
#SBATCH --output=G570_test.o%j
#SBATCH --error=G570_test.e%j
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --partition=cpucluster
#SBATCH --time=00:30:00
#SBATCH --mail-user=egonzales@sfsu.edu
#SBATCH --mail-type=ALL

source ~/.bashrc
module load mpi/openmpi-x86_64


declare -xr WDIR="/Users/912032971/Brewster_Ydwarfs"

declare PATH=${PATH}:${WDIR}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${WDIR}:~/

# Active python 3 environment
source activate retrievals

time_start=`date '+%T%t%d_%h_06'`
  
echo ------------------------------------------------------
echo -n 'Job is running on node '$SLURM_JOB_NODELIST
echo ------------------------------------------------------
echo SBATCH: sbatch is running on $SLURM_SUBMIT_HOST
echo SBATCH: originating queue is $SLURM_SUBMIT_PARTITION
echo SBATCH: executing queue is $SLURM_JOB_PARTITION
echo SBATCH: working directory is $SLURM_SUBMIT_DIR
echo SBATCH: job identifier is $SLURM_JOBID
echo SBATCH: job name is $SLURM_JOB_NAME
echo SBATCH: node file is $SLURM_JOB_NODELIST
echo SBATCH: current home directory is $SLURM_SUBMIT_HOME
echo SBATCH: PATH = $SLURM_SUBMIT_PATH
echo ------------------------------------------------------


cd ${WDIR}

mpirun python G570D_test.py > /Users/912032971/Brewster_Ydwarfs/Results/G570D/G570D_test.log

time_end=`date '+%T%t%d_%h_06'`
echo Started at: $time_start
echo Ended at: $time_end
echo ------------------------------------------------------
echo Job ends

