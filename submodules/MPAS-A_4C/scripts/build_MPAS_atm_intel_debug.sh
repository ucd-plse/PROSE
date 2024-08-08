#!/bin/bash
#PBS -N mpas_build_job
#PBS -A UCDV0023
#PBS -l walltime=0:10:00
#PBS -q develop
#PBS -l job_priority=economy
#PBS -j oe
#PBS -k eod
#PBS -l select=1:ncpus=16:mpiprocs=16

export PRECISION=double

if [ -z ${TUNING_REPO_ROOT} ]; then
  TUNING_REPO_ROOT=$(realpath $(dirname $(realpath $0))"/..")
fi

cd $TUNING_REPO_ROOT
source ./scripts/set_MPAS_env_intel.sh
make clean CORE=atmosphere
make -j16 ifort DEBUG=true USE_PIO2=true CORE=atmosphere PRECISION="${PRECISION}"