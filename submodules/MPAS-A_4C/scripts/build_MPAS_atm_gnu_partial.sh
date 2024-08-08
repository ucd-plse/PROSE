#!/bin/bash

export REPO_ROOT=/glade/work/jdvanover/latest_and_greatest/MPAS/MPAS-Model
export PIO_INCLUDE=/glade/u/apps/ch/opt/pio/1.10.1/mpt/2.22/gnu/8.3.0//include
export NETCDF_INCLUDE=/glade/u/apps/ch/opt/netcdf/4.7.1/gnu/8.3.0//include
export PNETCDF_INCLUDE=/glade/u/apps/ch/opt/pnetcdf/1.12.1/mpt/2.22/gnu/8.3.0//include

cd $(dirname "$0")/../
source ./scripts/set_MPAS_env_gnu.sh

set -v
set -e

cd src/framework
mpif90 -D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=unknown -O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 -c mpas_log.F90 -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I../external/esmf_time_f90
mpif90 -D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=unknown -O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 -c mpas_pool_routines.F90 -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I../external/esmf_time_f90
mpif90 -D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=v7.3-39-g9496c86b-dirty -O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 -c mpas_dmpar.F90 -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I../external/esmf_time_f90
mpif90 -D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=unknown -O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 -c mpas_timekeeping.F90 -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I../external/esmf_time_f90
ar -ru libframework.a mpas_kind_types.o mpas_framework.o mpas_timer.o mpas_timekeeping.o mpas_constants.o mpas_attlist.o mpas_hash.o mpas_sort.o mpas_block_decomp.o mpas_block_creator.o mpas_dmpar.o mpas_abort.o mpas_decomp.o mpas_threading.o mpas_io.o mpas_io_streams.o mpas_bootstrapping.o mpas_io_units.o mpas_stream_manager.o mpas_stream_list.o mpas_forcing.o mpas_c_interfacing.o random_id.o pool_hash.o mpas_derived_types.o mpas_domain_routines.o mpas_field_routines.o mpas_pool_routines.o xml_stream_parser.o regex_matching.o mpas_field_accessor.o mpas_log.o ../external/ezxml/ezxml.o
cd ..
ln -sf framework/libframework.a libframework.a


cd operators
mpif90 -D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=unknown -O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 -c mpas_vector_reconstruction.F90 -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I../framework -I../external/esmf_time_f90
ar -ru libops.a mpas_vector_operations.o mpas_matrix_operations.o mpas_tensor_operations.o mpas_rbf_interpolation.o mpas_vector_reconstruction.o mpas_spline_interpolation.o mpas_tracer_advection_helpers.o mpas_tracer_advection_mono.o mpas_tracer_advection_std.o mpas_geometry_utils.o
cd ..
ln -sf operators/libops.a libops.a

cd core_atmosphere/physics
mpif90 -D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=v7.3-39-g9496c86b-dirty -Dmpas  -O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 -c mpas_atmphys_todynamics.F90 -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I./physics_wrf -I.. -I../../framework -I../../external/esmf_time_f90
ar -ru libphys.a mpas_atmphys_camrad_init.o mpas_atmphys_control.o mpas_atmphys_driver.o mpas_atmphys_driver_cloudiness.o mpas_atmphys_driver_convection.o mpas_atmphys_driver_gwdo.o mpas_atmphys_driver_lsm.o mpas_atmphys_driver_microphysics.o mpas_atmphys_driver_oml.o mpas_atmphys_driver_pbl.o mpas_atmphys_driver_radiation_lw.o mpas_atmphys_driver_radiation_sw.o mpas_atmphys_driver_sfclayer.o mpas_atmphys_finalize.o mpas_atmphys_init.o mpas_atmphys_init_microphysics.o mpas_atmphys_interface.o mpas_atmphys_landuse.o mpas_atmphys_lsm_noahinit.o mpas_atmphys_manager.o mpas_atmphys_o3climatology.o mpas_atmphys_packages.o mpas_atmphys_rrtmg_lwinit.o mpas_atmphys_rrtmg_swinit.o mpas_atmphys_todynamics.o mpas_atmphys_update_surface.o mpas_atmphys_update.o mpas_atmphys_vars.o
cd ..

cd libphys
ar -x ../physics/libphys.a
cd ..

cd dynamics
mpif90 -D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=v7.3-39-g9496c86b-dirty -DDO_PHYSICS -O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 -c mpas_atm_boundaries.F90 -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I.. -I../../framework -I../../operators -I../physics -I../physics/physics_wrf -I../../external/esmf_time_f90
mpif90 -D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=v7.3-39-g9496c86b-dirty -DDO_PHYSICS -O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 -c mpas_atm_iau.F90 -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I.. -I../../framework -I../../operators -I../physics -I../physics/physics_wrf -I../../external/esmf_time_f90
mpif90 -D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=v7.3-39-g9496c86b-dirty -DDO_PHYSICS -O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 -c mpas_atm_time_integration.F90 -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I.. -I../../framework -I../../operators -I../physics -I../physics/physics_wrf -I../../external/esmf_time_f90
cd ..
mpif90 -D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=v7.3-39-g9496c86b-dirty -DDO_PHYSICS -O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 -c mpas_atm_core.F90 -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I${PIO_INCLUDE} -I${NETCDF_INCLUDE} -I${PNETCDF_INCLUDE} -I./inc -I../framework -I../operators -I./physics -I./dynamics -I./diagnostics -I./physics/physics_wrf -I../external/esmf_time_f90
ar -ru libdycore.a mpas_atm_core.o mpas_atm_core_interface.o mpas_atm_dimensions.o mpas_atm_threading.o dynamics/*.o libphys/*.o diagnostics/*.o
cd ..
ln -sf core_atmosphere/libdycore.a libdycore.a

cd driver
rm -f *.o *.mod *.f90
rm -f *.i
mpif90 -D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=v7.3-39-g9496c86b-dirty -O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 -c mpas_subdriver.F -I/glade/u/apps/ch/opt/pio/1.10.1/mpt/2.22/gnu/8.3.0//include -I/glade/u/apps/ch/opt/netcdf/4.7.1/gnu/8.3.0//include -I/glade/u/apps/ch/opt/pnetcdf/1.12.1/mpt/2.22/gnu/8.3.0//include -I/glade/u/apps/ch/opt/pio/1.10.1/mpt/2.22/gnu/8.3.0//include -I/glade/u/apps/ch/opt/netcdf/4.7.1/gnu/8.3.0//include -I/glade/u/apps/ch/opt/pnetcdf/1.12.1/mpt/2.22/gnu/8.3.0//include -I../framework -I../operators -I../core_atmosphere -I../external/esmf_time_f90
mpif90 -D_MPI -DCORE_ATMOSPHERE -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_EXE_NAME=atmosphere_model -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=v7.3-39-g9496c86b-dirty -O3 -ffree-line-length-none -fconvert=big-endian -ffree-form -fdefault-real-8 -fdefault-double-8 -c mpas.F -I/glade/u/apps/ch/opt/pio/1.10.1/mpt/2.22/gnu/8.3.0//include -I/glade/u/apps/ch/opt/netcdf/4.7.1/gnu/8.3.0//include -I/glade/u/apps/ch/opt/pnetcdf/1.12.1/mpt/2.22/gnu/8.3.0//include -I/glade/u/apps/ch/opt/pio/1.10.1/mpt/2.22/gnu/8.3.0//include -I/glade/u/apps/ch/opt/netcdf/4.7.1/gnu/8.3.0//include -I/glade/u/apps/ch/opt/pnetcdf/1.12.1/mpt/2.22/gnu/8.3.0//include -I../framework -I../operators -I../core_atmosphere -I../external/esmf_time_f90
cd ..
mpif90 -O3 -o atmosphere_model driver/*.o -L. -ldycore -lops -lframework -L/glade/u/apps/ch/opt/pio/1.10.1/mpt/2.22/gnu/8.3.0//lib -lpio -L/glade/u/apps/ch/opt/netcdf/4.7.1/gnu/8.3.0//lib -lnetcdff -lnetcdf -L/glade/u/apps/ch/opt/pnetcdf/1.12.1/mpt/2.22/gnu/8.3.0//lib -lpnetcdf -I./external/esmf_time_f90 -L./external/esmf_time_f90 -lesmf_time
cd ..

if [ -e src/atmosphere_model ]; then mv src/atmosphere_model .; fi

cd src/core_atmosphere
if [ ! -e ${REPO_ROOT}/default_inputs ]; then mkdir ${REPO_ROOT}/default_inputs; fi
cp default_inputs/* ${REPO_ROOT}/default_inputs/.
cd ${REPO_ROOT}/default_inputs; for FILE in `ls -1`; do if [ ! -e ../$FILE ]; then cp $FILE ../.; fi; done
