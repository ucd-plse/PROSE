module write_ocean_obs_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains a set of dummy interfaces for compiling the MOM6 DA
! driver code. These interfaces are not finalized and will be replaced by supported
! interfaces at some later date.
!
! 3/22/18
! matthew.harrison@noaa.gov
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 use mpp_io_mod, only : fieldtype, axistype, mpp_open,&
      MPP_OVERWR, MPP_NETCDF, MPP_MULTI, MPP_SINGLE,&
      mpp_write_meta, mpp_write, mpp_close
 use mpp_mod, only : mpp_error, FATAL
 use ocean_da_types_mod, only : ocean_profile_type
 use time_manager_mod!, only : time_type, get_time, set_date, operator ( - )

 implicit none

 private

 public :: open_profile_file, write_profile, close_profile_file, &
      write_ocean_obs_init

#include <netcdf.inc>

contains

function open_profile_file(name, nvar, grid_lon, grid_lat,thread,fset)

  character(len=*), intent(in) :: name
  integer, intent(in), optional :: nvar
  real, dimension(:), optional, intent(in) :: grid_lon, grid_lat
  integer, intent(in), optional :: thread, fset

  integer :: open_profile_file

  open_profile_file=-1
end function open_profile_file


subroutine write_profile(unit,profile)

integer, intent(in) :: unit
type(ocean_profile_type), intent(in) :: profile


return
end subroutine write_profile

subroutine close_profile_file(unit)

  integer, intent(in) :: unit

  return
end subroutine close_profile_file

subroutine write_ocean_obs_init()

  return

end subroutine write_ocean_obs_init

end module write_ocean_obs_mod


