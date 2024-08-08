#!/bin/bash

# the MPAS tutorial at https://www2.mmm.ucar.edu/projects/mpas/tutorial/Boulder2019/index.html recommends the following
# module load ncarenv/1.0
# module load ncarbinlibs/1.1
# module load gnu/8.3.0
# module load ncarcompilers/1.0
# module load netcdf/4.7.0
# module load mpich/3.3.1
# module load pnetcdf/1.11.2
# module load pio/1.9.23
# module load ncl/6.6.2
# module load ncview/2.1.7
# module load metis/5.1.0
# module load hdf5/1.10.5

module reset
module purge

echo "## LOADING MPAS GNU DEPENDENCIES"

# below is the closest match I could find on the current cheyenne system
module load ncarenv/1.3
module load gnu/8.3.0
module load ncarcompilers/0.5.0
module load netcdf/4.7.1
    # module load netcdf/4.6.3
module load mpt/2.22
module load pnetcdf/1.12.1
    # module load pnetcdf/1.11.1
module load pio/1.10.1
module load ncl/6.6.2
module load ncview/2.1.7

# LMOD reports that hdf5 cannot be loaded concurrently with netcdf
# module load hdf5/1.10.5

module list