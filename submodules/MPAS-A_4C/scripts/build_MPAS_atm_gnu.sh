#!/bin/bash

export PRECISION=double

cd $(dirname "$0")/../
source ./scripts/set_MPAS_env_gnu.sh
make clean CORE=atmosphere
make -j16 gfortran CORE=atmosphere PRECISION="${PRECISION}"