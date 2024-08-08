#!/bin/bash
#PBS -N ADCIRC_sub_job
#PBS -A UCDV0023
#PBS -l walltime=00:30:00
#PBS -q main
#PBS -l job_priority=regular
#PBS -j oe
#PBS -k eod
#PBS -l select=1:ncpus=128:mpiprocs=128

source ../../scripts/set_env.sh
mpiexec -n 128 ../../work/padcirc