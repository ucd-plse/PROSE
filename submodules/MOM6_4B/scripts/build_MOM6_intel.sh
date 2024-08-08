#!/bin/bash
#PBS -N mom6_build_job
#PBS -A UCDV0023
#PBS -l walltime=0:15:00
#PBS -q preempt
#PBS -l job_priority=economy
#PBS -j oe
#PBS -k eod
#PBS -l select=1:ncpus=8:mpiprocs=8

if [ -z ${TUNING_REPO_ROOT} ]; then
  TUNING_REPO_ROOT=$(realpath $(dirname $(realpath $0))"/..")
fi

source ${TUNING_REPO_ROOT}/scripts/set_env.sh

FMS_BLD_DIR="${TUNING_REPO_ROOT}"/build/intel/shared/repro
MOM6_BLD_DIR="${TUNING_REPO_ROOT}"/build/intel/ocean_only/repro
INTEL_MKMF_TEMPLATE="${TUNING_REPO_ROOT}"/src/mkmf/templates/derecho-intel.mk

TARGET_MOM6_OBJ=$1

set -v 
set -e # abort if any command fails

# build FMS:
if test -f "$FMS_BLD_DIR/libfms.a"; then
    echo "FMS already built."
else
  mkdir -p $FMS_BLD_DIR
  cd $FMS_BLD_DIR
  rm -f path_names *.?90 *.mod *.rmod *.o
  ../../../../src/mkmf/bin/list_paths ../../../../src/FMS
  ../../../../src/mkmf/bin/mkmf -t $INTEL_MKMF_TEMPLATE -p libfms.a -c "-Duse_libMPI -Duse_netCDF -DSPMD" path_names
  make clean
  make NETCDF=3 REPRO=1 libfms.a
fi

# build MOM6:
mkdir -p $MOM6_BLD_DIR
cd $MOM6_BLD_DIR
rm -f path_names *.?90 *.mod *.rmod *.o
../../../../src/mkmf/bin/list_paths ./ ../../../../src/MOM6/{config_src/dynamic,config_src/solo_driver,src/{*,*/*}}/
../../../../src/mkmf/bin/mkmf -t $INTEL_MKMF_TEMPLATE -o '-I../../shared/repro' -p MOM6 -l '-L../../shared/repro -lfms' -c '-Duse_libMPI -Duse_netCDF -DSPMD' path_names
make NETCDF=3 REPRO=1 $TARGET_MOM6_OBJ

echo "Successfully built MOM6 with INTEL."