#!/bin/bash
source ${PROSE_REPO_PATH}/scripts/set_env.sh
module load ncarenv/23.06
module load intel-classic/2023.0.0
module load cray-mpich/8.1.25
module load parallel-netcdf/1.12.3
module load parallelio/2.6.1
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/glade/u/home/jdvanover/gptl-8.1.1/install-intel/lib