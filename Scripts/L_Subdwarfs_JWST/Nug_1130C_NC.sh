#!/bin/csh
#PBS -S /bin/tcsh
#PBS -N Wolf1130C_McNug
#PBS -m abe
#PBS -l select=4:ncpus=24:mpiprocs=24:model=has
#PBS -l walltime=10:00:00
#PBS -k oe
#PBS -r n
#PBS -q long
#PBS -W group_list=s2429
#PBS -M egonzales@sfsu.edu

source /usr/share/Modules/init/csh
module purge
module load mpi-hpe/mpt.2.28_25Apr23_rhel87
module load comp-intel/2020.4.304

setenv MPI_REQUEST_MAX 512
setenv MPI_SHEPHERD true
setenv MPI_BUFS_PER_PROC 512

setenv WDIR /home5/egonza65/retrievals/Arcana_Subdwarfs
setenv PYTHONPATH ${WDIR}
#setenv PATH ${PATH}:${WDIR}
setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:${WDIR}"

source /swbuild/analytix/tools/miniconda3_220407/etc/profile.d/conda.csh
conda activate brewster376


# Debug output to check paths
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
echo "PYTHONPATH: $PYTHONPATH"

#setenv OMP_NUM_THREADS 20

unlimit stacksize

limit coredumpsize 0

set time_start=`date '+%T%t%d_%h_06'`
  
echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------




cd ${WDIR}


mpiexec -np 96 python 1130C_NC_mcnuggets.py > /nobackupp27/egonza65/Results/LSubdwarfs_JWST/1130C_NC_mcnuggets.log

set time_end=`date '+%T%t%d_%h_06'`
echo Started at: $time_start
echo Ended at: $time_end
echo ------------------------------------------------------
echo Job ends

