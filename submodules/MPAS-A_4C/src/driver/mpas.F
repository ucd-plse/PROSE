! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
program mpas

   use mpas_subdriver

   implicit none

#ifdef GPTL
  include 'gptl.inc'
  integer :: gptl_ret, gptl_i
  gptl_ret = gptlsetoption(gptlabort_on_error, 1)                  ! Abort on GPTL error
  gptl_ret = gptlsetoption(gptloverhead, 0)                        ! Turn off overhead estimate
  gptl_ret = gptlsetutr(gptlnanotime)                              ! Set underlying timer
  gptl_ret = gptlsetoption(gptlprint_method, gptlfull_tree)        ! Print full tree
  gptl_ret = gptlsetoption (gptlsync_mpi, 1)                       ! separates true MPI time from time spent waiting
  gptl_ret = gptlinitialize()                                      ! Initialize GPTL
#endif

   call mpas_init()

   call mpas_run() 

   call mpas_finalize()

#ifdef GPTL
  do gptl_i=0,63
        gptl_ret = gptlpr(gptl_i)
  end do
  gptl_ret = gptlfinalize()
#endif


   stop

end program mpas
