#!/bin/bash
#PBS -N MPAS_240k_uniform_run_job
#PBS -A UCDV0023
#PBS -l walltime=00:10:00
#PBS -q main
#PBS -l job_priority=regular
#PBS -j oe
#PBS -k eod
#PBS -l select=1:ncpus=64:mpiprocs=64:mem=235GB

source ../../scripts/set_MPAS_env_intel.sh
time ./run_script.sh