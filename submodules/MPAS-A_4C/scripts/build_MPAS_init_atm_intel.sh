#!/bin/bash
#PBS -N mpas_build_job
#PBS -A UCDV0023
#PBS -l walltime=0:10:00
#PBS -q preempt
#PBS -l job_priority=regular
#PBS -j oe
#PBS -k eod
#PBS -l select=1:ncpus=16:mpiprocs=16

export PRECISION=double

if [ -z ${TUNING_REPO_ROOT} ]; then
  TUNING_REPO_ROOT=$(realpath $(dirname $(realpath $0))"/..")
fi

cd $TUNING_REPO_ROOT
source ./scripts/set_MPAS_env_intel.sh
make clean CORE=init_atmosphere
make -j16 ifort USE_PIO2=true CORE=init_atmosphere PRECISION="${PRECISION}"