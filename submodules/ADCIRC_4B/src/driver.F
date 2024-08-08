!******************************************************************************
!******************************************************************************
!*    Non-ESMF ADCIRC DRIVER
!******************************************************************************
!******************************************************************************
      PROGRAM ADCIRC
!
      USE ADCIRC_Mod, ONLY : ADCIRC_Init, ADCIRC_Run, ADCIRC_Final
!
#ifdef CSWAN
!asey 090302: Include the following routines for coupling to unstructured SWAN.
      USE Couple2swan, ONLY: PADCSWAN_INIT, PADCSWAN_FINAL
      USE SIZES, ONLY: MNPROC,MYPROC
#endif
!
      IMPLICIT NONE
!
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
      CALL ADCIRC_Init
!
#ifdef CSWAN
!asey 090302: Allow SWAN to initialize stuff before the start
!             of the time step loop.  This subroutine is inside
!             the 'couple2swan.F' src file.
!asey 110518: Added this IF statement.
      IF(MYPROC.LT.MNPROC)THEN
         CALL PADCSWAN_INIT
      ENDIF
#endif
!
      CALL ADCIRC_Run
!
#ifdef CSWAN
!asey 090302: Let SWAN clean up stuff.
!asey 110518: Added this IF statement.
      IF(MYPROC.LT.MNPROC)THEN
         CALL PADCSWAN_FINAL
      ENDIF
#endif
!
      CALL ADCIRC_Final
!
#ifdef GPTL
      do gptl_i=0,127
            gptl_ret = gptlpr(gptl_i)
      end do
      gptl_ret = gptlfinalize()
#endif
      CALL EXIT(0)
      END PROGRAM ADCIRC
