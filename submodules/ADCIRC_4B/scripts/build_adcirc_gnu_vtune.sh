#!/bin/bash

TUNING_REPO_ROOT=$(realpath $(dirname $(realpath $0))/.. )
ADCIRC_SRC_DIR="${TUNING_REPO_ROOT}"/src
ADCIRC_BUILD_DIR="${TUNING_REPO_ROOT}"/work

module load gcc

cd $ADCIRC_BUILD_DIR
make clean
make clobber
export compiler=gnu
export NETCDFHOME=$NETCDF
export NETCDF=enable
export VTUNE=enable

set -e # abort if any command fails
make padcirc
echo "Successfully built padcirc with gnu and VTune flags"
