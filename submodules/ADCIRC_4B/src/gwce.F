! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
! 
!                       M O D U L E   G W C E
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!
! Generalized Wave Continuity Equation (GWCE) sets up and solves the 
! continuity (conservation of mass) equation. 
!
!******************************************************************************
! Logical Variable List (default value .FALSE., set in global.f)              *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     C2DDI            - 2D Depth Integrated model run                        *
!     C3D              - 3D model run                                         *
!     C3DDSS           - Stress form of the 3D momentum equations             *
!     C3DVS            - Velocity form of the 3D momentum equations           *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     C2D_BTrans       - Include a 2D baroclinic transport calculation        *
!     C2D_PTrans       - Include a 2D passive transport calculation           *
!     C3D_BTrans       - Include a 3D baroclinic transport calculation        *
!                        (used in 3D subroutines only)                        *
!     C3D_PTrans       - Include a 3D passive transport calculation           *
!                        (used in 3D subroutines only)                        *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     CBaroclinic      - Include baroclinic terms                             *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     CGWCE_New        - New ADCIRC GWCE formulation (old algorithm, new code)*
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     CGWCE_Lump       - Lump the GWCE matrix                                 *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     GWCE_New() and GWCE_New_pc() only:                                      *
!     CGWCE_LS_KGQ     - Kolar-Gray, flux-based, lateral stress formulation   *
!                        in the GWCE (same as original formulation)           *
!     CGWCE_LS_2PartQ  - 2 Part, flux-based, lateral stress formulation       *
!                        in the GWCE                                          *
!     CGWCE_LS_2PartV  - 2 Part, velocity-based, lateral stress formulation   *
!                        in the GWCE                                          *
!     CGWCE_LS_2PartSQ - 2 Part, flux-based, symmetric lateral stress         *
!                        formulation in the GWCE                              *
!     CGWCE_LS_2PartSV - 2 Part, velocity-based, symmetric lateral stress     *
!                        formulation in the GWCE                              *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     GWCE_New() and GWCE_New_pc() only:                                      *
!     CGWCE_Advec_NC   - Non-conservative advection formulation in the GWCE   *
!                        (same as original formulation)                       *
!     CGWCE_Advec_C1   - Use conservative advection formulation 1 in the GWCE *
!     CGWCE_Advec_C2   - Use conservative advection formulation 2 in the GWCE *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     CME_Orig         - Original Momentum Eq. formulation                    *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     CPRECOR          - Use the predictor-corrector algorithm for GWCE       *
!                        and momentum equations (package deal)                *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!     CSmag_Eh         - Use Smagorinski, spatially varying, vertically       *
!                        integrated lateral viscosity coefficient             *
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
!******************************************************************************

      module gwce
      use sizes, only : sz
      use global, only : DEBUG, INFO, WARNING, ERROR, logMessage, 
     &   setMessageSource, unsetMessageSource, allMessage
      implicit none

      real(kind=8),allocatable :: lsxx(:),lsxy(:),lsyx(:),lsyy(:)
      real(kind=8),allocatable :: etas0(:)      
      real(kind=8),allocatable :: eta0(:)
      real(kind=8) :: convcr
      integer :: isldia
      integer :: itmax
      real(kind=8) :: a00,b00,c00
      real(kind=8) :: dt2
      real(kind=8) :: ga00 

!     Arrays used by JCG iterative solver
      integer :: nw
      integer, allocatable :: iwksp(:)
      real(kind=8), allocatable :: wksp(:)
      integer, allocatable :: iparm(:)
      real(kind=8), allocatable :: rparm(:)
      integer :: numitr ! number of iterations per timestep 

!jjw...added 2 lines
!.....for matrix conditioning,global to allow for reading in from hotstart file
      REAL(kind=8) :: EP !jgf45.08 EP global for predictor-corrector
      !
      ! for managing datasets for aperiodic elevation specified boundaries (fort.19)
      real(kind=8) :: etime1, etime2, etiminc 

      ! jgf52.30.04: Allow the inverted barometer condition to be 
      ! activated on elevation-specified boundaries so that low pressure
      ! systems crossing the boundary won't trigger numerical instabilities.
      ! This is normally .false. but can be set to .true. in the fort.15
      ! in the metControl namelist. 
      logical :: invertedBarometerOnElevationBoundary = .false.

      contains

!----------------------------------------------------------------------
!           S U B R O U T I N E   I N I T  G  W  C  E
!----------------------------------------------------------------------
!     Allocates arrays and initializes constants.
!----------------------------------------------------------------------
      subroutine initGWCE()
      use sizes, only : mnp
      use global, only : dt, g, ilump, alloc_main11, alloc_main11_lumped
      use mesh, only : np
      use itpackv, only : dfault
      implicit none

      allocate(lsxx(mnp),lsxy(mnp),lsyx(mnp),lsyy(mnp))
      allocate(etas0(mnp))
      allocate(eta0(mnp))
      allocate( iparm(12),rparm(12) )
      allocate( iwksp(3*mnp),wksp(4*mnp+400) )
      dt2=dt*2.e0
      ga00=g*a00

!     jgf48.4619: Accommodate Seizos changes for explicit solve
      IF (ILump.eq.0) THEN ! default, fully consistent LHS
!        allocate arrays needed by GWCE matrix and iterative solver
         call alloc_main11()
!        initialize parameter arrays needed by iterative solver
         CALL DFAULT(IPARM,RPARM)
         IPARM(1)=ITMAX
         IPARM(2)=ISLDIA
         !cms51.06: moved opening of fort.33 to openLogFile sub in global.F
         IPARM(4)=33
         RPARM(1)=CONVCR
         NW = 4*NP + 4*ITMAX
      ELSE ! lumped LHS
         call alloc_main11_lumped()
      ENDIF

!----------------------------------------------------------------------
      end subroutine initGWCE
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!     S U B R O U T I N E   S O L V E   G   W   C   E
!----------------------------------------------------------------------
!     Two mutually exclusive possibilities are (1)
!     the conventional GWCE or (2) the predictor-corrector GWCE.
!----------------------------------------------------------------------
      subroutine solveGWCE(it, timeloc, timeh)
      use sizes, only : mnp
      use global, only : CGWCE_New, CPRECOR, ETA2, QX2, QY2, UU2, VV2
#ifdef CMPI
     & , dumy1, dumy2
      use messenger
#endif
      use momentum, only : mom_eqs_new_nc
      use nodalattributes, only : loadEleSlopeLim
      implicit none
      integer, intent(in) :: it ! time step 
      real(8), intent(in) :: timeloc ! simulation time in seconds
      real(8), intent(in) :: timeh

!       -------------PREDICTOR-CORRECTOR-GWCE---------------------------
!kmd    Changed to include the predictor-corrector algorithm
      IF(CPRECOR) THEN
         CALL GWCE_New(IT,TimeLoc,TimeH)
!...tcm added call for slope limiting         
         IF (LoadEleSlopeLim) THEN
            call check_slopes(it,TimeLoc)
            call apply_slope_limits(ETA2,MNP)
#ifdef CMPI
            CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif
         ENDIF         
         CALL Mom_Eqs_New_NC()
!...  If running in parallel, update velocities & fluxes on all processors
#ifdef CMPI
         CALL UPDATER(UU2,VV2,DUMY1,2)
         CALL UPDATER(QX2,QY2,DUMY1,2)
#endif
         CALL GWCE_New_pc(IT,TimeLoc,TimeH)
!...tcm added call to slope limiter         
         IF (LoadEleSlopeLim) THEN
            call check_slopes(it,TimeLoc)
            call apply_slope_limits(ETA2,MNP)
#ifdef CMPI
            CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif
         ENDIF         
      ENDIF
!       -------------PREDICTOR-CORRECTOR-GWCE---------------------------

!                         -OR-

!      -------------------GWCE------------------------------------------
      IF(CGWCE_New) THEN
         CALL GWCE_New(IT,TimeLoc,TimeH)
!.... tcm added call to slope limiter         
         IF (LoadEleSlopeLim) THEN
            call check_slopes(it,TimeLoc)
            call apply_slope_limits(ETA2,MNP)
#ifdef CMPI
            CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif
         ENDIF         
      ENDIF
!     -------------------GWCE------------------------------------------
!----------------------------------------------------------------------
      end subroutine solveGWCE
!----------------------------------------------------------------------


!******************************************************************************
!                                                                             *
!    Subroutine to compute the elevation using the GWCE formluation           *
!    Re-written to conform to the ADCIRC Theory Report                        *
!                                                                             *
!                            r.l.  06/22/2005                                 *
!******************************************************************************
      SUBROUTINE GWCE_New(IT,TimeLoc,TimeH)
!
#ifdef IEEE_DEBUG
      USE, INTRINSIC :: IEEE_ARITHMETIC
#endif
      USE SIZES, ONLY : sz, myproc
      USE GLOBAL, ONLY : GWCE_LV, COEF, COEFD, OBCCOEF, ETAS, ETA1, ETA2,
     &   ESBIN1, ESBIN2, C2DDI, C3D, CBaroclinic, CGWCE_Advec_C1,
     &   CGWCE_Advec_C2, CGWCE_Advec_NC, CGWCE_LS_2PartQ, CGWCE_LS_2PartSQ,
     &   CGWCE_LS_2PartV, CGWCE_LS_2PartSV, CGWCE_LS_KGQ, CPRECOR,
     &   CSmag_Eh, Smag_Comp_Flag, CTIP, DT, FluxSettlingIT, G, IFNLCAT,
     &   IFNLCT, IFNLFA, NBFR, NODECODE, NRS, NSCREEN, NCCHANGE, IDEN,
     &   NWS, RAMPELEV, RAMPMETE, RES_BC_FLAG, RHOWAT0, screenUnit,
     &   SMAG_LOWER_LIM, SMAG_UPPER_LIM, TK, NOFF, UU1, VV1, QX1, QY1,
     &   PR1, TiP1, BSX1, BSY1, DUU1, DUV1, DVV1, WSX1, WSY1, VIDBCPDXOH,
     &   VIDBCPDYOH, CORIF, PER, FACE, FF, EFA, EMO, AMIG, PR2, LNM_BC,
     &   QN1, EN1, QN0, ElevDisc, QN2, EN0, iLump, BCFLAG_LNM, EN2,
     &   TKM, NPERSEG, NTIP, NNPERBC, IPERCONN, VIDispDXOH, Bd, Ad, Cs2,
     &   VIDispDYOH, IFSFM, IFNL_HDP, CGWCE_HDP, CAliDisp, 
     &   usingDynamicWaterLevelCorrection, dynamicWaterLevelCorrection1
#ifdef CMPI
     &    , dumy1, dumy2, rnp_global
#endif
      USE MESH!, ONLY : NE, NP, DP, NNeigh, NeiTab, TotalArea, FDXE,
      USE MESH, ONLY : nm_mesh => NM
      USE BOUNDARIES, ONLY : NETA, NFLUXF, NFLUXB, NFLUXGBC, NFLUXIB, 
     &   NFLUXRBC, NOPE, NVEL, NBD, NBV, LBCODEI, BndLen2O3, NPEBC
      USE ITPACKV
      USE NodalAttributes, ONLY :
     &     LoadGeoidOffset, GeoidOffset, EVM, NOLIBF,
     &     TAU0VAR, HighResTimeVaryingTau0, FullDomainTimeVaryingTau0,
     &     CalculateTimeVaryingTau0, LoadAdvectionState, advectlocal,
!... DW
     &     LoadAbsLayerSigma, absorblayer_sigma_eta,
     &     absorblayer_sigma_mnx, absorblayer_sigma_mny
      USE SPONGELAYER
!... DW

#ifdef CMPI
      USE MESSENGER
#endif
      USE SUBDOMAIN, ONLY : subdomainOn, enforceBN, enforceECB, 
     &  checkChange, enforceEOB, enforceEIB, enforceGWCELVOB
      USE WIND, ONLY: PRBCKGRND_MH2O 

      IMPLICIT NONE

      INTEGER IE, JN, IJ, I, J                           !local loop counters
      INTEGER IT
      INTEGER NM1, NM2, NM3, NMI1, NMI2, NMI3, NMJ1, NMJ2, NMJ3
      INTEGER NC1, NC2, NC3, NCEle, NCI, NCJ
      INTEGER NCyc
      INTEGER NBDI
      INTEGER OnDiag, OffDiag

      LOGICAL  DIE
      REAL(kind=8) A00pB00
      REAL(kind=8) BCXAvg, BCYAvg
      REAL(kind=8) BndLenO6NC    !BNDLEN2O3NC, NCBND need to be removed from global.f and put in original GWCE subroutine
      REAL(kind=8) BSXN1, BSXN2, BSXN3, BSYN1, BSYN2, BSYN3, BSXAvg, BSYAvg
      REAL(kind=8) CorifAvg
      REAL(kind=8) DPAvg, GDPAvgOAreaIE4
      REAL(kind=8) DPAvgS0, GDPAvgOAreaIE4_S0
      REAL(kind=8) DispX, DispY, DispXAvg, DispYAvg
      REAL(kind=8) E0N1, E0N2, E0N3, E0XGrad2A, E0YGrad2A
      REAL(kind=8) E1N1, E1N2, E1N3, E1XGrad2A, E1YGrad2A
      REAL(kind=8) E1N1SQ, E1N2SQ, E1N3SQ
      REAL(kind=8) ESN1, ESN2, ESN3, ESAvg
      REAL(kind=8) EVMN1, EVMN2, EVMN3, EVMXGrad, EVMYGrad, EVMAvgODT
      REAL(kind=8) EVMEle, EVMSmag
      REAL(kind=8) GA00DPAvgOAreaIE4
      REAL(kind=8) GHAvg, GHAvgOAreaIE2, GOAreaIE4
      REAL(kind=8) H1N1, H1N2, H1N3, HAvg, H1, H2
      REAL(kind=8) H2OTotalArea
      REAL(kind=8) LSXXGradA, LSXYGradA, LSYXGradA, LSYYGradA
      REAL(kind=8) LSXXEle, LSXYEle, LSYXEle, LSYYEle
      REAL(kind=8) MsFacR, MsFacLOnDiag, MsFacLOffDiag
      REAL(kind=8) MX, MY, MXAvg, MYAvg
      REAL(kind=8) JXAvg, JYAvg
      REAL(kind=8) Pr1N1, Pr1N2, Pr1N3
      REAL(kind=8) Off1N1, Off1N2, Off1N3 ! water level offsets
      REAL(kind=8) QX1N1, QX1N2, QX1N3, QY1N1, QY1N2, QY1N3, QX1Avg, QY1Avg
      REAL(kind=8) SFacAvg, SFmxAvg, SFmyAvg, sfdxfac, sfdyfac 
      REAL(kind=8) SFCXAvg, SFCYAvg, YCFACAvg, SFCXADJ, SFCYADJ, TANPHIAvg
      REAL(kind=8) T0N1,T0N2, T0N3
      REAL(kind=8) Tau0Avg, Tau0QXAvg, Tau0QYAvg
      REAL(kind=8) Tau0XGrad2A, Tau0YGrad2A, Tau0SpaVar
      REAL(kind=8) TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      REAL(kind=8) TiPN1, TiPN2, TiPN3
      REAL(kind=8) U1N1,U1N2,U1N3, U1Avg
      REAL(kind=8) V1N1,V1N2,V1N3, V1Avg
      REAL(kind=8) WSXAvg, WSYAvg
      REAL(8) AreaIE, AreaIE2, AreaIE4
      REAL(kind=8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      REAL(kind=8) FDXX1, FDXX2, FDXX3, FDYY1, FDYY2, FDYY3
      REAL(8) TimeLoc, TimeH
      REAL(kind=8) HH1 !jgf46.02 Added for Katrina.
      REAL(kind=8) RDIAG ! jgf48.4619 Seizo parameter for fully explicit mode
      INTEGER :: j12, j13, j21, j23, j31, j32
      REAL(kind=8) :: arg, argj
      REAL(kind=8) :: celerity
      REAL(kind=8) :: etratio
      INTEGER :: ier
      INTEGER :: nbdj
      REAL(kind=8) :: rff
      REAL(kind=8) :: qforcei, qforcej
!...... DW
      REAL(kind=8):: SigmaAvg, Tau0Avg_S
      REAL(kind=8):: SgN1_eta, SgN2_eta, SgN3_eta
      REAL(kind=8):: SgN1_mnx, SgN2_mnx, SgN3_mnx
      REAL(kind=8):: SgN1_mny, SgN2_mny, SgN3_mny
      REAL(kind=8):: Tau0SigmaAvg, SigmaQxAvg, SigmaQyAvg
      REAL(kind=8):: AvgSigmaHAbsU, AvgSigmaHAbsV
      REAL(kind=8):: MsFacR_S, MsFacRS_T, MsFacRS_ABC00
      REAL(kind=8):: E0AN1, E0AN2, E0AN3
      REAL(kind=8):: E1AN1, E1AN2, E1AN3
      REAL(kind=8):: E2AN1, E2AN2, E2AN3
      REAL(kind=8):: ESAN1, ESAN2, ESAN3 
      REAL(kind=8):: U1AN1, U1AN2, U1AN3
      REAL(kind=8):: V1AN1, V1AN2, V1AN3
      REAL(kind=8):: ABC00_EAN1, ABC00_EAN2, ABC00_EAN3
      REAL(kind=8):: fwsx1, fwsx2, fwsx3
      REAL(kind=8):: fwsy1, fwsy2, fwsy3
      REAL(kind=8) :: SPM1, SPM2, SPM3 
      LOGICAL :: ETA_SPONGE = .false., GWCE_SPONGE = .false.
      REAL(kind=8) :: Ma2, Cfac = 1.0_SZ, CfacS0 = 1.0_SZ
      
      IF ( LoadAbsLayerSigma ) THEN
         IF (NumNodesAbsLayer(1) > 0 ) ETA_SPONGE = .true.
         IF (SUM(NumNodesAbsLayer) > 0 ) GWCE_SPONGE = .true.
      ENDIF
!..... DW

      call setMessageSource("gwce_new")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
!
!....... DW
!     Switch to the element table to the primary node
      CALL SWITCH_ELTAB_PERBC( IDIREC = 1 ) ;
!.....

      if (subdomainOn.and.enforceBN.eq.2) call checkChange()    ! NCSU Subdomain

!     jgf45.11 Bug fix: calculate the integers OnDiag and OffDiag here
!     instead of inside the GWCE lhs (system matrix) setup, since they
!     are also used in the calculation of the GWCE load vector gwce_lv.

!...  Consistent mass matrix: ILump=0, lumped mass matrix: ILump=1
      OnDiag=(1+ILump)*2        !diagonal coefficient
      OffDiag=(1-ILump)         !off diagonal coefficient
!...
!...  Recompute the GWCE system matrix at the first time step or if any
!...  wetting or drying occurred in the previous time step.
!...
!..... DW, Dev 2019
!.....  if CGWCE_HDP = .true. use H instead of DP in the surface gradient term
      NCIFBLCK: IF(NCChange.GT.0 .or. CGWCE_HDP ) THEN !if any subdomain grid has changed
         NCChange=0
!.....Set up the LHS matrix (for the iterative matrix solver)
         IF ( ILump.eq.0 ) THEN ! default, fully consistent case
            Coef(:,:)=0.0e0
         ELSE ! jgf48.4619: ILump.eq.1, only need the diagonals (Seizo)
            Coefd(:)=0.0e0    ! Only Diagnal
         ENDIF
!
!        jgf47.08 Moved time-varying tau0 subroutine to nodalattr.F
         IF(C2DDI.AND.
     &       ((FullDomainTimeVaryingTau0.or.HighResTimeVaryingTau0)))
     &   THEN
            DO I=1,NP
!              WJP 02.24.2018 Get the infinity norm of the matrix
               TK(I) = max(abs(TKM(1,I))+abs(TKM(3,I)),
     &                     abs(TKM(2,I))+abs(TKM(3,I)))
            ENDDO
            CALL CalculateTimeVaryingTau0(TK, NNeigh, NeiTab, NP)
         ENDIF

         IF(C3D.AND.
     &       ((FullDomainTimeVaryingTau0.or.HighResTimeVaryingTau0)))
     &   THEN
            CALL CalculateTimeVaryingTau0(TK, NNeigh, NeiTab, NP)
         ENDIF

         !jgf48.4619: Added Seizo's handling of Lumping vs non-Lumping
         COEFASSYM : IF ( ILump.eq.0 ) THEN ! default case: fully consistent LHS
            FULLMASS: DO IE=1,NE
               NMI1=nm_mesh(IE,1)
               NMI2=nm_mesh(IE,2)
               NMI3=nm_mesh(IE,3)
               NMJ1=NMI1
               NMJ2=NMI2
               NMJ3=NMI3
               NC1=NodeCode(NMI1)
               NC2=NodeCode(NMI2)
               NC3=NodeCode(NMI3)
               NCEle=NC1*NC2*NC3*NOFF(IE)

               SFacAvg = SFacEle(IE)
!.... Beg DW/WJP:
               SFmxAvg = SFMXEle(IE) ;
               SFmyAvg = SFMYEle(IE) ;
               sfdxfac = (1 - IFSFM)*SFacAvg + IFSFM*SFmxAvg ; 
               sfdyfac = (1 - IFSFM)*1.0_SZ  + IFSFM*SFmyAvg ; 
               SFCXAvg = SFCXEle(IE) ;
               SFCYAvg = SFCYEle(IE) ;
               YCFACAvg = YCSFACEle(IE) ;
               SFCXADJ = SFCXAvg/sfdxfac ; 
               SFCYADJ = YCFACAvg*SFCYAvg/sfdyfac ; 
               FDX1 = FDXE(1,IE)*sfdxfac ; !c FDX1=(Y(NM2)-Y(NM3))*SFacAvg !b1
               FDX2 = FDXE(2,IE)*sfdxfac ; !c FDX2=(Y(NM3)-Y(NM1))*SFacAvg !b2
               FDX3 = FDXE(3,IE)*sfdxfac ; !c FDX3=(Y(NM1)-Y(NM2))*SFacAvg !b3
               FDY1 = FDYE(1,IE)*sfdyfac ; !c  FDY1=X(NM3)-X(NM2) !a1
               FDY2 = FDYE(2,IE)*sfdyfac ; !c  FDY2=X(NM1)-X(NM3) !a2
               FDY3 = FDYE(3,IE)*sfdyfac ; !c  FDY3=X(NM2)-X(NM1) !a3
!.... End WJP/DW
               AreaIE2=Areas(IE)
               AreaIE =AreaIE2*0.5e0
               AreaIE4=AreaIE2*2.0e0

!.... DW, Nov 2019 if IFNL_HDP = 1 use DP + ETA = H in the GWCE
               DPAvg = (DP(NMI1)+DP(NMI2)+DP(NMI3))/3.e0 +   ! BAVG if GWCE_NL_HDP = .FALSE.
     &              IFNLFA*IFNL_HDP*(ETA2(NMI1) + ! HAVG otherwise
     &              ETA2(NMI2)+ETA2(NMI3))/3.e0 ; 
                  
!.....         WJP adjust the dispersion equation through
!              Ali et al. formula
               if (CAliDisp) then
                  Ma2   = G*DPavg/Cs2
                  Cfac = (1.0_SZ - 0.25_SZ*Ma2 - Ad*DPavg**Bd)**2
               endif
               
               GA00DPAvgOAreaIE4 = A00*G*DPAvg*Cfac/AreaIE4
               Tau0Avg = (Tau0Var(NMI1)+Tau0Var(NMI2)+Tau0Var(NMI3))/3.e0

!........... DW, Absorbing layer
!    
               SigmaAvg = 0.E0 ; Tau0SigmaAvg = 0.E0
               Tau0Avg_S = Tau0Avg
               IF ( ETA_SPONGE ) THEN
                  ! \overline( \sigma ) 
                  SigmaAvg = (absorblayer_sigma_eta(NMI1,1) 
     &                     +  absorblayer_sigma_eta(NMI2,1) 
     &                     +  absorblayer_sigma_eta(NMI3,1))/3.E0

                  ! \overline(\tau_{0} \sigma)
                  Tau0SigmaAvg = (absorblayer_sigma_eta(NMI1,1)*Tau0Var(NMI1)
     &              +  absorblayer_sigma_eta(NMI2,1)*Tau0Var(NMI2)
     &              +  absorblayer_sigma_eta(NMI3,1)*Tau0Var(NMI3))/3.E0

                  Tau0Avg_S = Tau0Avg + 
     &             (A00*2.E0*Tau0SigmaAvg*DT + SigmaAvg)*sponge_dis_mthd
               END IF              
!
               MsFacLOnDiag =OnDiag*AreaIE*(1.e0/DT+Tau0Avg_S/2.e0)/DT/12.e0
               MsFacLOffDiag=OffDiag*AreaIE*(1.e0/DT+Tau0Avg_S/2.e0)/DT/12.e0
!
!.......... DW
               DO JN=2,NEIMAX
                  IF(NeiTab(NMI1,JN).EQ.NMJ2) J12=JN
                  IF(NeiTab(NMI1,JN).EQ.NMJ3) J13=JN
                  IF(NeiTab(NMI2,JN).EQ.NMJ1) J21=JN
                  IF(NeiTab(NMI2,JN).EQ.NMJ3) J23=JN
                  IF(NeiTab(NMI3,JN).EQ.NMJ1) J31=JN
                  IF(NeiTab(NMI3,JN).EQ.NMJ2) J32=JN
               END DO
!.........     WJP multiplying by cosine(lat) if doing the 
!              spherical correction 
!              Node-by-Node
               SPM1 = SFCT(NMI1) 
               SPM2 = SFCT(NMI2)
               SPM3 = SFCT(NMI3)
!              Assemble the LHS matrix
               Coef(NMI1,1)  = Coef(NMI1,1)   + (MsFacLOnDiag*SPM1
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX1*FDX1 +
     &                                        SFCYADJ*FDY1*FDY1) )*NCELE
               Coef(NMI1,J12)= Coef(NMI1,J12) + (MsFacLOffDiag*SPM2
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX1*FDX2 + 
     &                                        SFCYADJ*FDY1*FDY2) )*NCELE
               Coef(NMI1,J13)= Coef(NMI1,J13) + (MsFacLOffDiag*SPM3
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX1*FDX3 +
     &                                        SFCYADJ*FDY1*FDY3) )*NCELE

               Coef(NMI2,J21)= Coef(NMI2,J21) + (MsFacLOffDiag*SPM1
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX2*FDX1 + 
     &                                        SFCYADJ*FDY2*FDY1) )*NCELE
               Coef(NMI2,1)  = Coef(NMI2,1)   + (MsFacLOnDiag*SPM2
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX2*FDX2 + 
     &                                        SFCYADJ*FDY2*FDY2) )*NCELE
               Coef(NMI2,J23)= Coef(NMI2,J23) + (MsFacLOffDiag*SPM3
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX2*FDX3 +
     &                                        SFCYADJ*FDY2*FDY3) )*NCELE

               Coef(NMI3,J31)= Coef(NMI3,J31) + (MsFacLOffDiag*SPM1
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX3*FDX1 + 
     &                                        SFCYADJ*FDY3*FDY1) )*NCELE
               Coef(NMI3,J32)= Coef(NMI3,J32) + (MsFacLOffDiag*SPM2
     &                   + GA00DPAvgOAreaIE4*( SFCXADJ*FDX3*FDX2 +
     &                                        SFCYADJ*FDY3*FDY2) )*NCELE
               Coef(NMI3,1)  = Coef(NMI3,1)   + (MsFacLOnDiag*SPM3
     &                  + GA00DPAvgOAreaIE4*( SFCXADJ*FDX3*FDX3 + 
     &                                        SFCYADJ*FDY3*FDY3) )*NCELE
            ENDDO FULLMASS
         ELSE
            !jgf48.4619: Add Seizo's construction of Lumped LHS matrix
            LUMPEDMASS: DO IE=1,NE ! Make LHS Lumped Matrix
               NMI1=nm_mesh(IE,1)
               NMI2=nm_mesh(IE,2)
               NMI3=nm_mesh(IE,3)
               AreaIE2=Areas(IE)
               AreaIE =AreaIE2*0.5e0
               NC1=NodeCode(NMI1)
               NC2=NodeCode(NMI2)
               NC3=NodeCode(NMI3)
               NCEle=NC1*NC2*NC3*NOFF(IE)

               Tau0Avg=(Tau0Var(NMI1)+Tau0Var(NMI2)+Tau0Var(NMI3))/3.e0
!
!.... DW, absorbing layer
               SigmaAvg = 0.E0 ; Tau0sigmaAvg = 0.E0
               Tau0Avg_S = Tau0Avg
               IF ( ETA_SPONGE ) THEN 
                 SigmaAvg = (absorblayer_sigma_eta(NMI1,1) 
     &             +  absorblayer_sigma_eta(NMI2,1) 
     &             +  absorblayer_sigma_eta(NMI3,1))/3.E0
                 Tau0SigmaAvg = (absorblayer_sigma_eta(NMI1,1)*Tau0Var(NMI1)
     &             +  absorblayer_sigma_eta(NMI2,1)*Tau0Var(NMI2)
     &             +  absorblayer_sigma_eta(NMI3,1)*Tau0Var(NMI3))/3.E0

                 Tau0Avg_S = Tau0Avg + 
     &             (A00*2.E0*Tau0SigmaAvg*DT + SigmaAvg)*sponge_dis_mthd
               END IF
               MsFacLOnDiag =OnDiag*AreaIE*(1.e0/DT+Tau0Avg_S/2.e0)/DT/12.e0
!.............DW
!.........     WJP multiplying by cosine(lat) if doing the 
!              spherical correction 
!              Node-by-Node
               SPM1 = SFCT(NMI1) 
               SPM2 = SFCT(NMI2)
               SPM3 = SFCT(NMI3)

               Coefd(NMI1) = Coefd(NMI1) + MsFacLOnDiag*SPM1*NCELE
               Coefd(NMI2) = Coefd(NMI2) + MsFacLOnDiag*SPM2*NCELE
               Coefd(NMI3) = Coefd(NMI3) + MsFacLOnDiag*SPM3*NCELE
            ENDDO LUMPEDMASS
         ENDIF COEFASSYM 

!...  Modify the matrix "COEF" by imposing the elevation specified
!...  boundary conditions while maintaining the symmetry of the system

         IF (ILump.eq.0) THEN
#ifdef CMPI
            EP = PSDOT(NP,Coef(:,1),Coef(:,1))
            EP = SQRT(RNP_GLOBAL*EP)
#else
            EP=0.0E0
            DO I=1,NP
               EP=EP+Coef(I,1)*Coef(I,1)
            ENDDO
            EP=SQRT(EP/NP)
#endif
!...        for each elevation specified boundary node, zero all off diagonal
!...        terms on the row and set diagnoal term to EP
            DO I=1,NETA
               Coef(NBD(I),1)=EP
               DO J=2,NNEIGH(NBD(I))
                  Coef(NBD(I),J)=0.0e0
               ENDDO
            ENDDO
!...        for each elevation specified boundary node, zero all off diagonal
!...        terms on the column but save these to be multiplied by the
!...        boundary value and subtracted from the RHS
            OBCCoef = 0.e0 !jgf53.dev: later code assumes this is initialized to zero
            DO I=1,NETA
               DO J=2,NNeigh(NBD(I))
                  DO IJ=2,NNeigh(NeiTab(NBD(I),J))
                     IF(NBD(I).EQ.NeiTab(NeiTab(NBD(I),J),IJ)) THEN
                        OBCCoef(I,J-1)=Coef(NeiTab(NBD(I),J),IJ)
                        Coef(NeiTab(NBD(I),J),IJ)=0.0e0
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
!.....   DW
            CALL UPDATE_COEF_PERBCSL_CONSISTENT() ;
!.....   DW
!.....      Check that all the diagonal elements in "COEF" are > 0.
            DIE = .FALSE.
            DO I=1,NP
               IF(COEF(I,1).EQ.0.e0) COEF(I,1)=EP
               IF(COEF(I,1).LT.0.e0) THEN                  
                  IF(NSCREEN.NE.0.AND.MYPROC.EQ.0)
     &              WRITE(ScreenUnit,1019) I,COEF(I,1)
                  WRITE(16,1019) I,COEF(I,1)
                  DIE = .TRUE.
               ENDIF
            ENDDO
            IF (DIE) THEN
               ! jgfdebug
               open(899,file='debug.txt',status='replace',action='write')
                  do i=1,np
                     if (coef(i,1).lt.0.e0) then
                        do j=2,nneighele(i)
                           if (neitabele(i,j).ne.0) then
                              write(6,
     & '("Node ",i0," element ",i0," area=",f15.7)') i, j, 
     &  areas(neitabele(i,j))
                           endif
                        enddo
                     endif
                  enddo
               close(899)
               ! end jgfdebug
               CALL gwceTerminate()
            ENDIF

         ELSE  ! jgf48.4619: include Seizo's changes for lumped LHS

            ! Seizo: Explicit scheme can solve localy. (the efect is small?)
            EP=0.0E0
            DO I=1,NP
               EP=EP+Coefd(I)*Coefd(I)
            ENDDO
            EP=SQRT(EP/NP)

            ! set diagonal term to EP
            DO I=1,NETA
               Coefd(NBD(I))=EP
            ENDDO
!......  DW
            CALL UPDATE_COEF_PERBCSL_LUMPED() ;
!......  DW
!.....      Check that all the diagonal elements in "COEFD" are > 0.
            DIE = .FALSE.
            DO I=1,NP
               IF(COEFD(I).EQ.0.e0) COEFD(I)=EP
               IF(COEFD(I).LT.0.e0) THEN
                  IF(NSCREEN.NE.0.AND.MYPROC.EQ.0)
     &              WRITE(ScreenUnit,1019) I,COEFD(I)
                  WRITE(16,1019) I,COEFD(I)
 1019             FORMAT(/,1X,'!!!!!!!!  ERROR  !!!!!!!',
     &              /,1X,'THE DIAGONAL TERM IN THE EQUATION FOR NODE ',I10,
     &                '= ',E15.6,' AND IS <= 0',/)
                  DIE = .TRUE.
               ENDIF
            ENDDO
            IF (DIE) CALL gwceTerminate()
         ENDIF
      ENDIF NCIFBLCK                    !End of GWCE matrix setup
!...
!...  Compute the GWCE load vector GWCE_LV
!...  This is done primarily element by element by forming
!...  temporary vectors and then assembling at the end.
!...  This has been set up to unroll loops to optimize performance
!...  on vector processors.
!...
!...  Elevation and flux boundary conditions are imposed after the
!...  element by element assembly section.
!...

!...  Initialize variables to zero if these forcings are not used
      IF((NWS.NE.0).OR.(NRS.NE.0)) THEN
      ELSE
        WSXAvg=0.e0
        WSYAvg=0.e0
        Pr1N1=0.e0
        Pr1N2=0.e0
        Pr1N3=0.e0
      ENDIF

      ! jgf: Adding support for water level offset as a pseudo
      ! barometric pressure
      if (usingDynamicWaterLevelCorrection.eqv..false.) then
         Off1N1 = 0.e0
         Off1N2 = 0.e0
         Off1N3 = 0.e0
      endif

      IF (CTIP) THEN
      ELSE
        TiPN1=0.e0
        TiPN2=0.e0
        TiPN3=0.e0
      ENDIF

      IF(C3D) THEN
      ELSE
        DispXAvg=0.e0
        DispYAvg=0.e0
      ENDIF

      IF(CBaroclinic) THEN
      ELSE
        BCXAvg=0.e0
        BCYAvg=0.e0
      ENDIF

!...  Compute the Lateral Stress Field using the 2 Part velocity approach (nonsymmetric or symmetric)

      IF ((CGWCE_LS_2PartV).OR.(CGWCE_LS_2PartSV)) THEN

         DO I=1,NP
            LSXX(I)=0.e0
            LSXY(I)=0.e0
            LSYX(I)=0.e0
            LSYY(I)=0.e0
         ENDDO

         DO IE=1,NE
            NM1=nm_mesh(IE,1)
            NM2=nm_mesh(IE,2)
            NM3=nm_mesh(IE,3)
            NC1=NODECODE(NM1)
            NC2=NODECODE(NM2)
            NC3=NODECODE(NM3)
            NCEle=NC1*NC2*NC3*NOFF(IE)
            SFacAvg = SFacEle(IE)
            SFmxAvg=SFMXEle(IE) ;
            SFmyAvg=SFMYEle(IE) ;
            sfdxfac=(1 - IFSFM)*SFacAvg + IFSFM*SFmxAvg ; 
            sfdyfac=(1 - IFSFM)*1.0_SZ  + IFSFM*SFmyAvg ;            
!........ BEG DW:
            FDX1 = FDXE(1,IE)*sfdxfac ; 
            FDX2 = FDXE(2,IE)*sfdxfac ; 
            FDX3 = FDXE(3,IE)*sfdxfac ; 
            FDY1 = FDYE(1,IE)*sfdyfac ; 
            FDY2 = FDYE(2,IE)*sfdyfac ; 
            FDY3 = FDYE(3,IE)*sfdyfac ; 
!........ END DW:
            LSXXGradA=(UU1(NM1)*FDX1+UU1(NM2)*FDX2+UU1(NM3)*FDX3)/2.e0   !A*DUDX
            LSXYGradA=(UU1(NM1)*FDY1+UU1(NM2)*FDY2+UU1(NM3)*FDY3)/2.e0   !A*DUDY
            LSYXGradA=(VV1(NM1)*FDX1+VV1(NM2)*FDX2+VV1(NM3)*FDX3)/2.e0   !A*DVDX
            LSYYGradA=(VV1(NM1)*FDY1+VV1(NM2)*FDY2+VV1(NM3)*FDY3)/2.e0   !A*DVDY
            EVMEle=NCEle*(EVM(NM1)+EVM(NM2)+EVM(NM3))/3.e0
            !
            ! If using Smagorinski vertically-integrated lateral stress coefficient
            IF ((CSmag_Eh).or.(Smag_Comp_Flag)) THEN  
               EVMSmag=EVMEle*
     &             sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA)
     &                 +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
               !tcm v52.30.01 added test for limits
               IF(Smag_Comp_Flag) then
                  IF (EVMSmag .gt. SMAG_UPPER_LIM) EVMSmag = SMAG_UPPER_LIM
                  IF (EVMSmag .lt. SMAG_LOWER_LIM) EVMSmag = SMAG_LOWER_LIM
               ENDIF
               EVMEle=EVMSmag
            ENDIF
            LSXXEle = LSXXGradA*EVMEle
            LSXX(NM1)=LSXX(NM1)+LSXXEle
            LSXX(NM2)=LSXX(NM2)+LSXXEle
            LSXX(NM3)=LSXX(NM3)+LSXXEle
            LSXYEle = LSXYGradA*EVMEle
            LSXY(NM1)=LSXY(NM1)+LSXYEle
            LSXY(NM2)=LSXY(NM2)+LSXYEle
            LSXY(NM3)=LSXY(NM3)+LSXYEle
            LSYXEle = LSYXGradA*EVMEle
            LSYX(NM1)=LSYX(NM1)+LSYXEle
            LSYX(NM2)=LSYX(NM2)+LSYXEle
            LSYX(NM3)=LSYX(NM3)+LSYXEle
            LSYYEle = LSYYGradA*EVMEle
            LSYY(NM1)=LSYY(NM1)+LSYYEle
            LSYY(NM2)=LSYY(NM2)+LSYYEle
            LSYY(NM3)=LSYY(NM3)+LSYYEle
         ENDDO

         DO I=1,NP
            IF(TotalArea(I).NE.0.) THEN
               H2=DP(I)+IFNLFA*ETA2(I)
               H2OTotalArea=H2/TotalArea(I)
               IF (CGWCE_LS_2PartV) THEN          !nonsymmetric
                  LSXX(I)=H2OTotalArea*LSXX(I)
                  LSXY(I)=H2OTotalArea*LSXY(I)
                  LSYX(I)=H2OTotalArea*LSYX(I)
                  LSYY(I)=H2OTotalArea*LSYY(I)
               ENDIF
               IF (CGWCE_LS_2PartSV) THEN         !symmetric
                  LSXX(I)=H2OTotalArea*LSXX(I)
                  LSXY(I)=0.5e0*H2OTotalArea*(LSXY(I)+LSYX(I))
                  LSYX(I)=LSXY(I)
                  LSYY(I)=H2OTotalArea*LSYY(I)
               ENDIF
            ELSE
               LSXX(I)=0.e0
               LSXY(I)=0.e0
               LSYX(I)=0.e0
               LSYY(I)=0.e0
            ENDIF
         ENDDO

      ENDIF

!...  Compute the Lateral Stress Field using the 2 Part flux approach (nonsymmetric or symmetric)

      IF ((CGWCE_LS_2PartQ).OR.(CGWCE_LS_2PartSQ)) THEN

         DO I=1,NP
            LSXX(I)=0.e0
            LSXY(I)=0.e0
            LSYX(I)=0.e0
            LSYY(I)=0.e0
         ENDDO

         DO IE=1,NE
            NM1=nm_mesh(IE,1)
            NM2=nm_mesh(IE,2)
            NM3=nm_mesh(IE,3)
            NC1=NODECODE(NM1)
            NC2=NODECODE(NM2)
            NC3=NODECODE(NM3)
            NCEle=NC1*NC2*NC3*NOFF(IE)
            SFacAvg = SFacEle(IE)
!..... Beg DW:
            SFmxAvg=SFMXEle(IE) ;
            SFmyAvg=SFMYEle(IE) ;
            sfdxfac=(1 - IFSFM)*SFacAvg + IFSFM*SFmxAvg ; 
            sfdyfac=(1 - IFSFM)*1.0_SZ  + IFSFM*SFmyAvg ;            

            FDX1 = FDXE(1,IE)*sfdxfac ; 
            FDX2 = FDXE(2,IE)*sfdxfac ; 
            FDX3 = FDXE(3,IE)*sfdxfac ; 
            FDY1 = FDYE(1,IE)*sfdyfac ; 
            FDY2 = FDYE(2,IE)*sfdyfac ; 
            FDY3 = FDYE(3,IE)*sfdyfac ; 
!..... End DW
            EVMEle = NCEle*(EVM(NM1)+EVM(NM2)+EVM(NM3))/3.e0
            IF ((CSmag_Eh).or.(Smag_Comp_Flag)) THEN  !If using Smagorinski vertically-integrated lateral stress coefficient
               LSXXGradA=(UU1(NM1)*FDX1+UU1(NM2)*FDX2+UU1(NM3)*FDX3)/2.e0
               LSXYGradA=(UU1(NM1)*FDY1+UU1(NM2)*FDY2+UU1(NM3)*FDY3)/2.e0
               LSYXGradA=(VV1(NM1)*FDX1+VV1(NM2)*FDX2+VV1(NM3)*FDX3)/2.e0
               LSYYGradA=(VV1(NM1)*FDY1+VV1(NM2)*FDY2+VV1(NM3)*FDY3)/2.e0
               EVMSmag=EVMEle*
     &            sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA)
     &                +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
               ! tcm v52.30.01 added test for limits
               IF(Smag_Comp_Flag) then
                  IF (EVMSmag .gt. SMAG_UPPER_LIM) EVMSmag = SMAG_UPPER_LIM
                  IF (EVMSmag .lt. SMAG_LOWER_LIM) EVMSmag = SMAG_LOWER_LIM
               ENDIF
               EVMEle=EVMSmag
            ENDIF
            LSXXGradA=(QX1(NM1)*FDX1+QX1(NM2)*FDX2+QX1(NM3)*FDX3)/2.e0
            LSXYGradA=(QX1(NM1)*FDY1+QX1(NM2)*FDY2+QX1(NM3)*FDY3)/2.e0
            LSYXGradA=(QY1(NM1)*FDX1+QY1(NM2)*FDX2+QY1(NM3)*FDX3)/2.e0
            LSYYGradA=(QY1(NM1)*FDY1+QY1(NM2)*FDY2+QY1(NM3)*FDY3)/2.e0
            LSXXEle = LSXXGradA*EVMEle
            LSXX(NM1)=LSXX(NM1)+LSXXEle
            LSXX(NM2)=LSXX(NM2)+LSXXEle
            LSXX(NM3)=LSXX(NM3)+LSXXEle
            LSXYEle = LSXYGradA*EVMEle
            LSXY(NM1)=LSXY(NM1)+LSXYEle
            LSXY(NM2)=LSXY(NM2)+LSXYEle
            LSXY(NM3)=LSXY(NM3)+LSXYEle
            LSYXEle = LSYXGradA*EVMEle
            LSYX(NM1)=LSYX(NM1)+LSYXEle
            LSYX(NM2)=LSYX(NM2)+LSYXEle
            LSYX(NM3)=LSYX(NM3)+LSYXEle
            LSYYEle = LSYYGradA*EVMEle
            LSYY(NM1)=LSYY(NM1)+LSYYEle
            LSYY(NM2)=LSYY(NM2)+LSYYEle
            LSYY(NM3)=LSYY(NM3)+LSYYEle
         ENDDO

         DO I=1,NP
            IF(TotalArea(I).NE.0.) THEN
               IF (CGWCE_LS_2PartQ) THEN          !nonsymmetric
                  LSXX(I)=LSXX(I)/TotalArea(I)
                  LSXY(I)=LSXY(I)/TotalArea(I)
                  LSYX(I)=LSYX(I)/TotalArea(I)
                  LSYY(I)=LSYY(I)/TotalArea(I)
               ENDIF
               IF (CGWCE_LS_2PartSQ) THEN         !symmetric
                  LSXX(I)=LSXX(I)/TotalArea(I)
                  LSXY(I)=0.5e0*(LSXY(I)+LSYX(I))/TotalArea(I)
                  LSYX(I)=LSXY(I)
                  LSYY(I)=LSYY(I)/TotalArea(I)
               ENDIF
            ELSE
               LSXX(I)=0.e0
               LSXY(I)=0.e0
               LSYX(I)=0.e0
               LSYY(I)=0.e0
            ENDIF
         ENDDO

      ENDIF

!...  Assemble the GWCE RHS except for the boundary integral terms

      DO 1037 IE=1,NE

!...     Set nodal values for each element

!orbitt 120322: Localized Advection
         IF (LoadAdvectionState) CALL ADVECTLOCAL(IE)

         NM1=nm_mesh(IE,1)
         NM2=nm_mesh(IE,2)
         NM3=nm_mesh(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
         NCELE=NC1*NC2*NC3*NOFF(IE)
         E0N1=ETA1(NM1)
         E0N2=ETA1(NM2)
         E0N3=ETA1(NM3)
         E1N1=ETA2(NM1)
         E1N2=ETA2(NM2)
         E1N3=ETA2(NM3)
         E1N1SQ=E1N1*E1N1
         E1N2SQ=E1N2*E1N2
         E1N3SQ=E1N3*E1N3
         ESN1=ETAS(NM1)
         ESN2=ETAS(NM2)
         ESN3=ETAS(NM3)
         U1N1=UU1(NM1)
         U1N2=UU1(NM2)
         U1N3=UU1(NM3)
         V1N1=VV1(NM1)
         V1N2=VV1(NM2)
         V1N3=VV1(NM3)
         QX1N1=QX1(NM1)
         QX1N2=QX1(NM2)
         QX1N3=QX1(NM3)
         QY1N1=QY1(NM1)
         QY1N2=QY1(NM2)
         QY1N3=QY1(NM3)
         H1N1=DP(NM1)+IFNLFA*E1N1
         H1N2=DP(NM2)+IFNLFA*E1N2
         H1N3=DP(NM3)+IFNLFA*E1N3
         EVMN1=EVM(NM1)
         EVMN2=EVM(NM2)
         EVMN3=EVM(NM3)
         T0N1=Tau0Var(NM1)
         T0N2=Tau0Var(NM2)
         T0N3=Tau0Var(NM3)
!
!.... DW: begin ! absorbing layer
! 
         E0AN1 = 0.E0 ; E0AN2 = 0.E0 ; E0AN3 = 0.E0 ;
         E1AN1 = 0.E0 ; E1AN2 = 0.E0 ; E1AN3 = 0.E0 ;
         E2AN1 = 0.E0 ; E2AN2 = 0.E0 ; E2AN3 = 0.E0 ;
         ESAN1 = 0.E0 ; ESAN2 = 0.E0 ; ESAN3 = 0.E0 ;
 
         ABC00_EAN1 = 0.E0 ; 
         ABC00_EAN2 = 0.E0 ;
         ABC00_EAN3 = 0.E0 ;

         U1AN1 = 0.E0 ; U1AN2 = 0.E0 ; U1AN3 = 0.E0 ;
         V1AN1 = 0.E0 ; V1AN2 = 0.E0 ; V1AN3 = 0.E0 ;

         SgN1_eta = 0.E0 ; SgN2_eta = 0.E0 ; SgN3_eta = 0.E0 ;
         SgN1_mnx = 0.E0 ; SgN2_mnx = 0.E0 ; SgN3_mnx = 0.E0 ;
         SgN1_mny = 0.E0 ; SgN2_mny = 0.E0 ; SgN3_mny = 0.E0 ;
  
         IF ( GWCE_SPONGE ) THEN 
            SgN1_eta = absorblayer_sigma_eta(NM1,1) ; 
            SgN2_eta = absorblayer_sigma_eta(NM2,1) ; 
            SgN3_eta = absorblayer_sigma_eta(NM3,1) ;

            SgN1_mnx = absorblayer_sigma_mnx(NM1,1) ; 
            SgN2_mnx = absorblayer_sigma_mnx(NM2,1) ; 
            SgN3_mnx = absorblayer_sigma_mnx(NM3,1) ;

            SgN1_mny = absorblayer_sigma_mny(NM1,1) ; 
            SgN2_mny = absorblayer_sigma_mny(NM2,1) ; 
            SgN3_mny = absorblayer_sigma_mny(NM3,1) ;

            E0AN1 = eta0_abslayer(NM1) ;
            E0AN2 = eta0_abslayer(NM2) ; 
            E0AN3 = eta0_abslayer(NM3) ; 

            E1AN1 = eta1_abslayer(NM1) ; 
            E1AN2 = eta1_abslayer(NM2) ; 
            E1AN3 = eta1_abslayer(NM3) ; 

            E2AN1 = eta2_abslayer(NM1) ;
            E2AN2 = eta2_abslayer(NM2) ; 
            E2AN3 = eta2_abslayer(NM3) ;

            ESAN1 = E2AN1 - E0AN1 ;
            ESAN2 = E2AN2 - E0AN2 ;
            ESAN3 = E2AN3 - E0AN3 ;

            ABC00_EAN1 = A00*E2AN1 + B00*E1AN1 + C00*E0AN1 ; 
            ABC00_EAN2 = A00*E2AN2 + B00*E1AN2 + C00*E0AN2 ; 
            ABC00_EAN3 = A00*E2AN3 + B00*E1AN3 + C00*E0AN3 ; 

            U1AN1 = uu1_abslayer(NM1) ; 
            U1AN2 = uu1_abslayer(NM2) ; 
            U1AN3 = uu1_abslayer(NM3) ; 
 
            V1AN1 = vv1_abslayer(NM1) ;
            V1AN2 = vv1_abslayer(NM2) ;
            V1AN3 = vv1_abslayer(NM3) ;          
         END IF
!
         IF ((NWS.NE.0).OR.(NRS.NE.0)) THEN     !wind or radiation stress
            Pr1N1=PR1(NM1)
            Pr1N2=PR1(NM2)
            Pr1N3=PR1(NM3)

            !c Pressure c!
            IF ( NO_MET_IN_SPONGE ) THEN
               IF ( absorblayer_sigma_eta(NM1,1) > 1.0e-9 ) THEN
                  Pr1N1 = PRBCKGRND_MH2O
               END IF
               IF ( absorblayer_sigma_eta(NM2,1) > 1.0e-9 ) THEN
                  Pr1N2 = PRBCKGRND_MH2O
               END IF
               IF ( absorblayer_sigma_eta(NM3,1) > 1.0e-9 ) THEN
                  Pr1N3 = PRBCKGRND_MH2O
               END IF
            END IF
         ENDIF
! ....,... END DW
!
         ! jgf: Adding support for water level offset as a pseudo barometric
         ! pressure
         if (usingDynamicWaterLevelCorrection.eqv..true.) then
            Off1N1 = dynamicWaterLevelCorrection1(NM1)
            Off1N2 = dynamicWaterLevelCorrection1(NM2)
            Off1N3 = dynamicWaterLevelCorrection1(NM3)
         endif
         IF (CTIP) THEN                        !tidal potential
            TiPN1=TiP1(NM1)
            TiPN2=TiP1(NM2)
            TiPN3=TiP1(NM3)
         ENDIF
         IF (C2DDI) THEN                       !2D bottom friction
            !WJP taking into account the symmetrical friction matrix
            BSXN1 = TKM(1,NM1)*QX1N1 + TKM(3,NM1)*QY1N1
            BSYN1 = TKM(3,NM1)*QX1N1 + TKM(2,NM1)*QY1N1
            BSXN2 = TKM(1,NM2)*QX1N2 + TKM(3,NM2)*QY1N2
            BSYN2 = TKM(3,NM2)*QX1N2 + TKM(2,NM2)*QY1N2
            BSXN3 = TKM(1,NM3)*QX1N3 + TKM(3,NM3)*QY1N3
            BSYN3 = TKM(3,NM3)*QX1N3 + TKM(2,NM3)*QY1N3
         ENDIF
         IF (C3D) THEN                         !3D bottom friction
            BSXN1=BSX1(NM1)
            BSXN2=BSX1(NM2)
            BSXN3=BSX1(NM3)
            BSYN1=BSY1(NM1)
            BSYN2=BSY1(NM2)
            BSYN3=BSY1(NM3)
         ENDIF

         AreaIE2=Areas(IE)               !2A
         AreaIE=AreaIE2/2.e0             ! A
         AreaIE4=2.e0*AreaIE2            !4A
         SFacAvg= SFacEle(IE)
!..... DW
         ! derivative: accomodate corrected spherical coordinate
         SFmxAvg=SFMXEle(IE)  ;
         SFmyAvg=SFMYEle(IE) ;
         sfdxfac=(1 - IFSFM)*SFacAvg + IFSFM*SFmxAvg ; 
         sfdyfac=(1 - IFSFM)*1.0_SZ  + IFSFM*SFmyAvg ; 
         ! continuity
         SFCXAvg=SFCXEle(IE) ;
         SFCYAvg=SFCYEle(IE) ;
         YCFACAvg=YCSFACEle(IE) ;
         TANPHIAvg=TANPHIEle(IE) ; 
         SFCXADJ = SFCXAvg/sfdxfac ; 
         SFCYADJ = YCFACAvg*SFCYAvg/sfdyfac ; 

         FDX1 = FDXE(1,IE)*sfdxfac ; !c dx/dlambda*b1*2A
         FDX2 = FDXE(2,IE)*sfdxfac ; !c dx/dlambda*b2*2A
         FDX3 = FDXE(3,IE)*sfdxfac ; !c dx/dlambda*b3*2A
         FDY1 = FDYE(1,IE)*sfdyfac ; !c dy/dphi*a1*2A
         FDY2 = FDYE(2,IE)*sfdyfac ; !c dy/dphi*a1*2A
         FDY3 = FDYE(3,IE)*sfdyfac ; !c dy/dphi*a1*2A
!..... DW

!...     Compute part of several spatial gradients for use below

         E0XGrad2A = 0.0e0
         E0YGrad2A = 0.0e0
         IF (ILump.eq.0) THEN
            E0XGrad2A=E0N1*FDX1+E0N2*FDX2+E0N3*FDX3        !2*Area*deta0/dx
            E0YGrad2A=E0N1*FDY1+E0N2*FDY2+E0N3*FDY3        !2*Area*deta0/dy
         ENDIF
         E1XGrad2A=E1N1*FDX1+E1N2*FDX2+E1N3*FDX3        !2*Area*deta1/dx
         E1YGrad2A=E1N1*FDY1+E1N2*FDY2+E1N3*FDY3        !2*Area*deta1/dy
         Tau0XGrad2A=T0N1*FDX1+T0N2*FDX2+T0N3*FDX3      !2*Area*dTau0/dx
         Tau0YGrad2A=T0N1*FDY1+T0N2*FDY2+T0N3*FDY3      !2*Area*dTau0/dy

!...... BEG: DW, Nov 2019
         Tau0XGrad2A=SFCXADJ*Tau0XGrad2A      !2*Area*dTau0/dx
         Tau0YGrad2A=SFCYADJ*Tau0YGrad2A      !2*Area*dTau0/dy
!...... END DW

!...     Compute the Kolar & Gray lateral stress term extended for spatially varying EVM
!...     ********     Need to think about this more: DW ********
!...     ********     Kolar & Gray does not seem to be suitable for  
!...     ********     the spherical coordinates        
         IF(CGWCE_LS_KGQ) THEN
            EVMXGrad=(EVMN1*FDX1+EVMN2*FDX2+EVMN3*FDX3)/AreaIE2
            EVMYGrad=(EVMN1*FDY1+EVMN2*FDY2+EVMN3*FDY3)/AreaIE2
            EVMAvgODT=((EVMN1+EVMN2+EVMN3)/3.e0)/DT
            MX=(EVMXGrad*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)
     &         +EVMYGrad*(QY1N1*FDX1+QY1N2*FDX2+QY1N3*FDX3)
     &         -EVMAvgODT*(ESN1*FDX1+ESN2*FDX2+ESN3*FDX3))/AreaIE2
            MY=(EVMXGrad*(QX1N1*FDY1+QX1N2*FDY2+QX1N3*FDY3)
     &         +EVMYGrad*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3)
     &         -EVMAvgODT*(ESN1*FDY1+ESN2*FDY2+ESN3*FDY3))/AreaIE2
         ENDIF

!...     Compute the remainder of the 2 Part lateral stress terms

         IF((CGWCE_LS_2PartQ).OR.(CGWCE_LS_2PartV) .OR.
     &        (CGWCE_LS_2PartSQ).OR.(CGWCE_LS_2PartSV)) THEN
!           X-direction
            LSXYEle = (LSXY(NM1)+LSXY(NM2)+LSXY(NM3))/3.e0
            LSYXEle = (LSYX(NM1)+LSYX(NM2)+LSYX(NM3))/3.e0
            MX=(LSXX(NM1)*FDX1+LSXX(NM2)*FDX2+LSXX(NM3)*FDX3
     &         +LSXY(NM1)*FDY1+LSXY(NM2)*FDY2+LSXY(NM3)*FDY3)/AreaIE2
     &         -TANPHIAvg*(LSXYEle + LSYXEle)
!           Y-direction 
            LSXXEle = (LSXX(NM1)+LSXX(NM2)+LSXX(NM3))/3.e0
            LSYYEle = (LSYY(NM1)+LSYY(NM2)+LSYY(NM3))/3.e0
            MY=(LSYX(NM1)*FDX1+LSYX(NM2)*FDX2+LSYX(NM3)*FDX3
     &         +LSYY(NM1)*FDY1+LSYY(NM2)*FDY2+LSYY(NM3)*FDY3)/AreaIE2
     &         +TANPHIAvg*(LSXXEle - LSYYEle)
         ENDIF

!...     Compute the spatial gradients of the velocity dispersion terms if 3D

         IF (C3D) THEN                         !3D bottom friction
            DispX = (DUU1(NM1)*FDX1+DUU1(NM2)*FDX2+DUU1(NM3)*FDX3
     &          +  DUV1(NM1)*FDY1+DUV1(NM2)*FDY2+DUV1(NM3)*FDY3)/AreaIE2
            DispY = (DUV1(NM1)*FDX1+DUV1(NM2)*FDX2+DUV1(NM3)*FDX3
     &          +  DVV1(NM1)*FDY1+DVV1(NM2)*FDY2+DVV1(NM3)*FDY3)/AreaIE2
            DispXAvg = IFNLCT*DispX
            DispYAvg = IFNLCT*DispY
!         ELSEIF (CBaroclinic.and.abs(IDEN).eq.6) THEN
!            ! If 2D and baroclinic and getting information from a 3D BC
!            ! model (the gradient has already been calculated)
!            DispXAvg=(H1N1*VIDISPDXOH(NM1)+H1N2*VIDISPDXOH(NM2)
!     &                                       +H1N3*VIDISPDXOH(NM3))/3.e0
!            DispYAvg=(H1N1*VIDISPDYOH(NM1)+H1N2*VIDISPDYOH(NM2)
!     &                                       +H1N3*VIDISPDYOH(NM3))/3.e0
         ENDIF

!...     Compute elemental averages
         CorifAvg=(Corif(NM1)+Corif(NM2)+Corif(NM3))/3.e0
         Tau0Avg=(T0N1+T0N2+T0N3)/3.e0
         Tau0QXAvg=(T0N1*QX1N1+T0N2*QX1N2+T0N3*QX1N3)/3.e0
         Tau0QYAvg=(T0N1*QY1N1+T0N2*QY1N2+T0N3*QY1N3)/3.e0
         U1Avg=(U1N1+U1N2+U1N3)/3.e0
         V1Avg=(V1N1+V1N2+V1N3)/3.e0
         QX1Avg=(QX1N1+QX1N2+QX1N3)/3.e0
         QY1Avg=(QY1N1+QY1N2+QY1N3)/3.e0
         ESAvg=(ESN1+ESN2+ESN3)/3.e0

!.....   DW, Nov 2019 
         DPAvg = (DP(NM1)+DP(NM2)+DP(NM3))/3.e0 + ! BAVG if GWCE_HDP = .FALSE. 
     &        IFNLFA*IFNL_HDP*(E1N1 + E1N2 + E1N3)/3.e0 ; ! H^{s}_{avg} otherwise
         
         DPAvgS0 = (DP(NM1)+DP(NM2)+DP(NM3))/3.e0 + ! BAVG if GWCE_HDP = .FALSE. 
     &        IFNLFA*IFNL_HDP*(E0N1 + E0N2 + E0N3)/3.e0 ; ! H^{s-1}_{avg} otherwise

!.....   WJP adjust the dispersion equation through
!        Ali et al. formula
         if (CAliDisp) then
            Ma2    = G*DPavg/Cs2
            Cfac   = (1.0_SZ - 0.25_SZ*Ma2 - Ad*DPavg**Bd)**2
            Ma2    = G*DPavg/Cs2
            CfacS0 = (1.0_SZ - 0.25_SZ*Ma2 - Ad*DPavgS0**Bd)**2
         endif

         GDPAvgOAreaIE4 = G*DPAvg*Cfac/AreaIE4
         GDPAvgOAreaIE4_S0 = G*DPAvgS0*CfacS0/AreaIE4
!....  END DW & WJP
         HAvg=(H1N1+H1N2+H1N3)/3.e0
         GHAvg=G*HAvg
         GHAvgOAreaIE2=GHAvg/AreaIE2
         BSXAvg=(BSXN1+BSXN2+BSXN3)/3.e0
         BSYAvg=(BSYN1+BSYN2+BSYN3)/3.e0
         MXAvg=MX           !lateral stresses are constant over an element
         MYAvg=MY           !lateral stresses are constant over an element
!        WJP: Add on the contribution from spherical correction term
!        from the momentum equation
         ! (tan \phi V_{lambda})/R + f
         CorifAvg = CorifAvg + IFNLCT*TANPHIAvg*U1Avg
!
!...... DW absorbing layer
!
         SigmaAvg = (SgN1_eta + SgN2_eta + SgN3_eta)/3.E0
         Tau0SigmaAvg = (SgN1_eta*T0N1 + SgN2_eta*T0N2 + 
     &                   SgN3_eta*T0N3)/3.e0

         ! Use in the momentum equation parts
         SigmaQxAvg = (SgN1_mnx*QX1N1 +
     &         SgN2_mnx*QX1N2 + SgN3_mnx*QX1N3)/3.e0 ;
         SigmaQyAvg = (SgN1_mny*QY1N1 +
     &         SgN2_mny*QY1N2 + SgN3_mny*QY1N3)/3.e0 ;

         AvgSigmaHAbsU = (SgN1_mnx*H1N1*U1AN1 +
     &        SgN2_mnx*H1N2*U1AN2 + SgN3_mnx*H1N3*U1AN3)/3.E0 ; 

         AvgSigmaHAbsV = (SgN1_mny*H1N1*V1AN1 +
     &        SgN2_mny*H1N2*V1AN2 + SgN3_mny*H1N3*V1AN3)/3.E0 ; 

         IF ((NWS.NE.0).OR.(NRS.NE.0)) THEN     !wind or radiation stress
            fwsx1 = 1.E0 ; fwsx2 = 1.E0 ; fwsx3 = 1.E0 ;
            fwsy1 = 1.E0 ; fwsy2 = 1.E0 ; fwsy3 = 1.E0 ;
            IF ( NO_MET_IN_SPONGE ) THEN
               IF ( absorblayer_sigma_mnx(NM1,1) > 1e-9 ) fwsx1 = 0.E0
               IF ( absorblayer_sigma_mnx(NM2,1) > 1e-9 ) fwsx2 = 0.E0
               IF ( absorblayer_sigma_mnx(NM3,1) > 1e-9 ) fwsx3 = 0.E0
           
               IF ( absorblayer_sigma_mny(NM1,1) > 1e-9 ) fwsy1 = 0.E0
               IF ( absorblayer_sigma_mny(NM2,1) > 1e-9 ) fwsy2 = 0.E0
               IF ( absorblayer_sigma_mny(NM3,1) > 1e-9 ) fwsy3 = 0.E0
            END IF
         
            WSXAvg=(fwsx1*WSX1(NM1)+fwsx2*WSX1(NM2)+fwsx3*WSX1(NM3))/3.e0
            WSYAvg=(fwsy1*WSY1(NM1)+fwsy2*WSY1(NM2)+fwsy3*WSY1(NM3))/3.e0
         ENDIF
!------- END DW
! 
         IF (CBaroclinic) THEN
#ifndef NOFSBPG
            BCXAvg = (H1N1*VIDBCPDXOH(NM1)+H1N2*VIDBCPDXOH(NM2)
     &                                    +H1N3*VIDBCPDXOH(NM3))/3.e0
            BCYAvg = (H1N1*VIDBCPDYOH(NM1)+H1N2*VIDBCPDYOH(NM2)
     &                                    +H1N3*VIDBCPDYOH(NM3))/3.e0
#else
            BCXAvg = (DP(NM1)*VIDBCPDXOH(NM1) + DP(NM2)*VIDBCPDXOH(NM2)
     &                + DP(NM3)*VIDBCPDXOH(NM3))/3.e0
            BCYAvg = (DP(NM1)*VIDBCPDYOH(NM1) + DP(NM2)*VIDBCPDYOH(NM2)
     &                + DP(NM3)*VIDBCPDYOH(NM3))/3.e0
#endif
         ENDIF

!...     Compute additional partial factors
!........DW: begin add absorbing layer term
         Tau0Avg_S = Tau0Avg + 
     &     ((A00 + B00)*2.E0*Tau0SigmaAvg*DT + SigmaAvg)*sponge_dis_mthd
   
         MsFacR=AreaIE*(1.e0/DT-Tau0Avg_S/2.e0)/DT/12.e0
         MsFacR_S = (-(A00 + B00 + C00)*AreaIE*
     &          Tau0SigmaAvg/12.E0)*sponge_dis_mthd 

         MsFacRS_ABC00 = (Tau0SigmaAvg*AreaIE/12.E0)*sponge_dis_mthd
         MsFacRS_T = ((SigmaAvg/(2.E0*DT))*AreaIE/12.E0)*sponge_dis_mthd
!......... end DW
!
         GOAreaIE4=G/AreaIE4
         Tau0SpaVar=(QX1Avg*Tau0XGrad2A+QY1Avg*Tau0YGrad2A)/6.e0
         A00pB00=A00+B00

!...     Compute the JX, JY terms less the advection terms
         JXAvg = CorifAvg*QY1Avg
     &          -IFNLFA*GOAreaIE4*(1 - IFNL_HDP)*(
     &               E1N1SQ*FDX1 +E1N2SQ*FDX2 + E1N3SQ*FDX3)
     &          -GHAvgOAreaIE2*((PR1N1-TiPN1-Off1N1)*FDX1
     &          +(PR1N2-TiPN2-Off1N2)*FDX2+(PR1N3-TiPN3-Off1N3)*FDX3)
     &          +WSXAvg-BSXAvg+MXAvg-DispXAvg-BCXAvg+Tau0QXAvg 
!--- begin: DW, absorbing layer
     &          - SigmaQxAvg*sponge_dis_mthd
     &          + AvgSigmaHAbsU*sponge_dis_mthd 
!--- end DW


         JYAvg =-CorifAvg*QX1Avg
     &          -IFNLFA*GOAreaIE4*(1 - IFNL_HDP)*(
     &                E1N1SQ*FDY1+E1N2SQ*FDY2+E1N3SQ*FDY3)
     &          -GHAvgOAreaIE2*((PR1N1-TiPN1-Off1N1)*FDY1
     &          +(PR1N2-TiPN2-Off1N2)*FDY2+(PR1N3-TiPN3-Off1N3)*FDY3)
     &          +WSYAvg-BSYAvg+MYAvg-DispYAvg-BCYAvg+Tau0QYAvg
!--- begin: DW, absorbing layer
     &          -SigmaQYAvg*sponge_dis_mthd 
     &          +AvgSigmaHAbsV*sponge_dis_mthd
!--- end DW

!...     Complete the JX, JY terms depending on the advection formulation
         IF(CGWCE_Advec_NC) THEN        !nonconservative advection
           JXAvg = JXAvg - IFNLCT*(
     &              QX1Avg*(U1N1*FDX1+U1N2*FDX2+U1N3*FDX3)
     &             +QY1Avg*(U1N1*FDY1+U1N2*FDY2+U1N3*FDY3))/AreaIE2
     &             +IFNLCAT*U1Avg*ESAvg/DT
           JYAvg = JYAvg - IFNLCT*(
     &              QX1Avg*(V1N1*FDX1+V1N2*FDX2+V1N3*FDX3)
     &             +QY1Avg*(V1N1*FDY1+V1N2*FDY2+V1N3*FDY3))/AreaIE2
     &             +IFNLCAT*V1Avg*ESAvg/DT
         ENDIF
         IF(CGWCE_Advec_C1) THEN        !conservative v1 advection
           JXAvg = JXAvg - IFNLCT*(
     &              (U1N1*QX1N1*FDX1+U1N2*QX1N2*FDX2+U1N3*QX1N3*FDX3)
     &             +(U1N1*QY1N1*FDY1+U1N2*QY1N2*FDY2+U1N3*QY1N3*FDY3))
     &             /AreaIE2
           JYAvg = JYAvg - IFNLCT*(
     &              (V1N1*QX1N1*FDX1+V1N2*QX1N2*FDX2+V1N3*QX1N3*FDX3)
     &             +(V1N1*QY1N1*FDY1+V1N2*QY1N2*FDY2+V1N3*QY1N3*FDY3))
     &             /AreaIE2
         ENDIF
         IF(CGWCE_Advec_C2) THEN        !conservative v2 advection
           JXAvg = JXAvg - IFNLCT*(
     &              QX1Avg*(U1N1*FDX1+U1N2*FDX2+U1N3*FDX3)
     &             +QY1Avg*(U1N1*FDY1+U1N2*FDY2+U1N3*FDY3)
     &             +U1Avg*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)
     &             +U1Avg*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3))/AreaIE2
           JYAvg = JYAvg - IFNLCT*(
     &              QX1Avg*(V1N1*FDX1+V1N2*FDX2+V1N3*FDX3)
     &             +QY1Avg*(V1N1*FDY1+V1N2*FDY2+V1N3*FDY3)
     &             +V1Avg*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)
     &             +V1Avg*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3))/AreaIE2
         ENDIF

!...... BEG: DW, DEC 2019         
!...     Needs multiplying by cos(lat) in case of the correction 
!        Node by node
         SPM1 = SFCT(NM1) ; 
         SPM2 = SFCT(NM2) ; 
         SPM3 = SFCT(NM3) ;
!...... END DW         

!...     Assemble forcing for node NM1 (local index j=1)

         Temp_LV_A1=

!...     Transient and Tau0 terms from LHS
     &         (OnDiag*ESN1*SPM1 + OffDiag*(ESN2*SPM2+ESN3*SPM3))*MsFacR

!...     Free surface terms from LHS (time levels s-1 & s)
     &        - GDPAvgOAreaIE4_S0*( C00*( SFCXADJ*FDX1*E0XGrad2A+ 
     &                                    SFCYADJ*FDY1*E0YGrad2A ) ) 
     &        - GDPAvgOAreaIE4*(A00pB00*( SFCXADJ*FDX1*E1XGrad2A+
     &                                    SFCYADJ*FDY1*E1YGrad2A ) )

!...     Terms from momentum eqs.
     &        +(SFCXADJ*JXAvg*FDX1 + SFCYADJ*JYAvg*FDY1)*0.5e0

!...     Spatially varying Tau0 terms
     &        +Tau0SpaVar
!
!..... DW: begin
     &        +(OnDiag*E0N1 + OffDiag*(E0N2+E0N3))*MsFacR_S
     &        +(onDiag*ESAN1 + OffDiag*(ESAN2+ESAN3))*MsFacRS_T
     &        +(onDiag*ABC00_EAN1 + 
     &                OffDiag*(ABC00_EAN2 + ABC00_EAN3))*MsFacRS_ABC00
!..... END DW
!

!...     Assemble forcing for node NM2 (local index j=2)
         Temp_LV_A2=

!...     Transient and Tau0 terms from LHS
     &         (OnDiag*ESN2*SPM2 + OffDiag*(ESN1*SPM1+ESN3*SPM3))*MsFacR

!...     Free surface terms from LHS (time levels s-1 & s)
     &        - GDPAvgOAreaIE4_S0*( C00*( SFCXADJ*FDX2*E0XGrad2A +
     &                                    SFCYADJ*FDY2*E0YGrad2A ) ) 
     &        - GDPAvgOAreaIE4*(A00pB00*( SFCXADJ*FDX2*E1XGrad2A +
     &                                    SFCYADJ*FDY2*E1YGrad2A ) )

!...     Terms from momentum eqs.
     &        +(SFCXADJ*JXAvg*FDX2 + SFCYADJ*JYAvg*FDY2)*0.5e0

!...     Spatially varying Tau0 terms
     &        +Tau0SpaVar 
!
!..... DW: begin
     &        +(OnDiag*E0N2 + OffDiag*(E0N1+E0N3))*MsFacR_S
     &        +(onDiag*ESAN2 + OffDiag*(ESAN1+ESAN3))*MsFacRS_T
     &        +(onDiag*ABC00_EAN2 + 
     &                 OffDiag*(ABC00_EAN1 + ABC00_EAN3))*MsFacRS_ABC00
!..... END DW
!


!...     Assemble forcing for node NM3 (local index j=3)

         Temp_LV_A3=

!...     Transient and Tau0 terms from LHS
     &         (OnDiag*ESN3*SPM3 + OffDiag*(ESN1*SPM1+ESN2*SPM2))*MsFacR

!...     Free surface terms from LHS (time levels s-1 & s)
     &        - GDPAvgOAreaIE4_S0*(  C00*( SFCXADJ*FDX3*E0XGrad2A +
     &                                     SFCYADJ*FDY3*E0YGrad2A ) )
     &        - GDPAvgOAreaIE4*( A00pB00*( SFCXADJ*FDX3*E1XGrad2A +
     &                                     SFCYADJ*FDY3*E1YGrad2A ) )

!...     Terms from momentum eqs.
     &        +(SFCXADJ*JXAvg*FDX3 + SFCYADJ*JYAvg*FDY3)*0.5e0

!...     Spatially varying Tau0 terms
     &        +Tau0SpaVar
!
!..... DW: begin, absorbing layer
     &        +(OnDiag*E0N3 + OffDiag*(E0N1+E0N2))*MsFacR_S
     &        +(OnDiag*ESAN3 + OffDiag*(ESAN1+ESAN2))*MsFacRS_T
     &        +(OnDiag*ABC00_EAN3 + 
     &               OffDiag*(ABC00_EAN1 + ABC00_EAN2))*MsFacRS_ABC00
!..... END DW
!

!...     Put these partial products into further elemental storage for a vector computer
!...     These will be put into nodal storage outside of the elemental loop
#ifdef CVEC
         Temp_LV_A(IE,1)=Temp_LV_A1*NCEle
         Temp_LV_A(IE,2)=Temp_LV_A2*NCEle
         Temp_LV_A(IE,3)=Temp_LV_A3*NCEle
#endif

!...     Put these partial products directly into nodal storage for a scalar (non-vector) computer
#ifdef CSCA
         GWCE_LV(NM1)=GWCE_LV(NM1)+Temp_LV_A1*NCEle
         GWCE_LV(NM2)=GWCE_LV(NM2)+Temp_LV_A2*NCEle
         GWCE_LV(NM3)=GWCE_LV(NM3)+Temp_LV_A3*NCEle
#endif

 1037 CONTINUE                  !End of elemental loop

!...  Put load vector elemental values into nodal storage for a vector computer
#ifdef CVEC
      DO IE=1,NE
         NM1=nm_mesh(IE,1)
         NM2=nm_mesh(IE,2)
         NM3=nm_mesh(IE,3)
         GWCE_LV(NM1)=GWCE_LV(NM1)+Temp_LV_A(IE,1)
         GWCE_LV(NM2)=GWCE_LV(NM2)+Temp_LV_A(IE,2)
         GWCE_LV(NM3)=GWCE_LV(NM3)+Temp_LV_A(IE,3)
      END DO
#endif

!...  Save the elevation at the past time step into Eta1 and zero out Eta2
!kmd  Need to save z(s-1) and etas(s-1) for the corrector loop

      DO I=1,NP
         IF(CPRECOR) THEN
            ETAS0(I)=ETAS(I)
            Eta0(I)=Eta1(I)
         END IF
         Eta1(I)=Eta2(I)
         Eta2(I)=0.0e0
      END DO

!...  At elevation boundary condition nodes, determine the elevation at
!...  the s+1 time step
!...
!...  For periodic elevation boundary conditions

      DO J=1,NBFR
         IF(PER(J).EQ.0.) THEN
            NCYC=0
         ELSE
#ifdef IBM
            NCYC=INT(timeh/PER(J),KIND(0.0e0))
#else
            NCYC=INT(timeh/PER(J))
#endif
         ENDIF
         ARGJ=AMIG(J)*(timeh-NCYC*PER(J))+FACE(J)
         RFF=FF(J)*RampElev
         DO I=1,NETA
            ARG=ARGJ-EFA(J,I)
            NBDI=NBD(I)
            Eta2(NBDI)=Eta2(NBDI)+EMO(J,I)*RFF*COS(ARG)
         END DO
      END DO

!...  FOR APERIODIC ELEVATION BOUNDARY CONDITION

      if (subdomainOn) then                     ! NCSU Subdomain
          if(enforceBN.eq.1) call enforceEcb()  ! NCSU Subdomain
      else                                      ! NCSU Subdomain
! 03.27.2018 WJP allowing for periodic plus non-periodic Elevation BCS
          IF (NPEBC) THEN
            IF(TimeLoc.GT.ETIME2) THEN
               ETIME1=ETIME2
               ETIME2=ETIME1+ETIMINC
               DO J=1,NETA
                  ESBIN1(J)=ESBIN2(J)
                  READ(19,*) ESBIN2(J)
               END DO
            ENDIF
            ETRATIO=(TimeLoc-ETIME1)/ETIMINC
            DO I=1,NETA
               NBDI=NBD(I)
               Eta2(NBDI)=Eta2(NBDI) + RampElev
     &           *(ESBIN1(I)+ETRATIO*(ESBIN2(I)-ESBIN1(I)))
            END DO
         ENDIF
      endif                                   ! NCSU Subdomain
      ! 
      ! jgf46.02 Added the ability to include geoid offset on the boundary.
      !
      ! aaltuntas51.48: Deactivated geoidOffset for subdomain boundary 
      ! conditions.
      if (subdomainOn.and.enforceBN.eq.1) then  ! NCSU Subdomain
         continue                                ! NCSU Subdomain
      else 
         IF (LoadGeoidOffset) THEN
            DO I=1,NETA
               ETA2(NBD(I))=ETA2(NBD(I))+GeoidOffset(NBD(I))
            END DO
         ENDIF
      endif                                      ! NCSU Subdomain
!
!
!     jgf48.04 Added an inverted barometer boundary condition so that
!     low pressure systems can cross the boundary without creating an
!     elevation anomaly.
!     jgf52.30.04: Included a parameter from fort.15 to turn this 
!     on and off according to analyst preference. 
      IF (NO_MET_IN_SPONGE.eqv..false.) THEN
         if (invertedBarometerOnElevationBoundary.eqv..true.) then
            DO I=1,NETA
               ETA2(NBD(I)) = ETA2(NBD(I))
     &                  + RampMete*(101300.e0/(RHOWAT0*G) - PR2(NBD(I)))
            END DO
         endif
      ENDIF

      !jgf: Added water level offset to elevation specified boundary
      if (usingDynamicWaterLevelCorrection.eqv..true.) then
         DO I=1,NETA
            ETA2(NBD(I))=ETA2(NBD(I)) + dynamicWaterLevelCorrection1(nbd(i)) 
           END DO
      endif

!   kmd48.33bc add information for the levels of no motion boundary conditions
!              these are considered the steric adjustments.
      IF ((ABS(RES_BC_FLAG).GE.1).AND.(CBaroclinic).AND.(NOPE.GT.0))THEN
!asey 140701: Added the following IF statement.
         IF(BCFLAG_LNM.GT.0)THEN
            DO I=1,NETA
              NBDI=NBD(I)
              ETA2(NBDI) = ETA2(NBDI) + LNM_BC(I)
            END DO
         END IF
      ENDIF

!...  IMPOSE NORMAL FLOW, RADIATION OR GRADIENT BOUNDARY CONDITIONS
!...  ALONG FLOW BOUNDARY TO LOAD VECTOR GWCE_LV(I)

!...  Note 2, Boundary conditions using specified fluxes (LBCODEI < 29)
!...  assume that QN is positive into the domain.  QFORCEJ has a -1
!...  built in and the terms are not explicitly negated. Boundary
!...  conditions using computed fluxes (LBCODEI 30, 40) compute a normal
!...  flux that  is positive out of the domain.  Therefore, to match
!...  the formulation these terms must be explicitly multiplied by -1.

!...Note 3, Eta1 is the latest computed elevation (it was updated above).

      IF((NFLUXF.EQ.1).OR.(NFLUXB.EQ.1).OR.(NFLUXIB.EQ.1)
     &     .OR.(NFLUXGBC.EQ.1).OR.(NFLUXRBC.EQ.1)) THEN
         NBDJ=NBV(1)
         IF (LBCODEI(1).LE.29) QFORCEJ=(QN2(1)-QN0(1))/DT2 +
     &        Tau0VAR(NBDJ)*QN1(1)

         IF (LBCODEI(1).EQ.30) THEN
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(1)
         ENDIF

         IF (LBCODEI(1).EQ.32) THEN
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
!C DW, Below is an original scheme and is correct 
!C            QFORCEJ=(QN1(1)-QN0(1))/DT
!C     &           -CELERITY*(ETAS(NBDJ)-(EN1(1)-EN0(1)))/DT
!C     &           +TAU0VAR(NBDJ)*(QN1(1)-CELERITY*(ETA1(NBDJ)-EN1(1)))
!C DW, a 2nd order extrapolation formula: \eta^{n+1} = 2*\eta^{n} - \eta^{n-1} is used
!C     and second order in time derivative is used
            QFORCEJ=(QN2(1)-QN0(1))/DT2
     &           - CELERITY*(2.0*ETAS(NBDJ) - (EN2(1)-EN0(1)))/DT2
     &           + TAU0VAR(NBDJ)*(QN1(1)-CELERITY*(ETA1(NBDJ)-EN1(1))) 
         ENDIF

         IF((LBCODEI(1).EQ.40).OR.(LBCODEI(1).EQ.41)) QFORCEJ=
     &        -(QN1(1)-QN0(1))/DT - TAU0VAR(NBDJ)*(QN1(1)+QN0(1))/2.e0

!     jgf46.21 Added IBTYPE=52.
         IF(LBCODEI(1).EQ.52) THEN
            QFORCEJ=(QN2(1)-QN0(1))/DT2 + Tau0VAR(NBDJ)*QN1(1)
            IF (IT.GT.FluxSettlingIT) THEN
               HH1=DP(NBDJ)+IFNLFA*Eta1(NBDJ)
               Celerity=SQRT(G*HH1)
               QforceJ=QforceJ - Celerity*(EtaS(NBDJ)/DT
     &              + Tau0Var(NBDJ)*(Eta1(NBDJ)-ElevDisc(1)))
            ENDIF
         ENDIF
!
         DO J=2,NVEL
            NBDI=NBDJ
            NBDJ=NBV(J)
            QFORCEI=QFORCEJ

            IF(LBCODEI(J).LE.29) QFORCEJ=(QN2(J)-QN0(J))/DT2+
     &           Tau0VAR(NBDJ)*QN1(J)
            IF(LBCODEI(J).EQ.30) THEN
               H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
               CELERITY=SQRT(G*H1)
               QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(J)
            ENDIF

            IF(LBCODEI(J).EQ.32) THEN
               H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
               CELERITY=SQRT(G*H1)
!C DW, This is an original formula and is definitely correct (1st order in time)
!C                QFORCEJ=(QN1(J)-QN0(J))/DT
!C     &              -CELERITY*(ETAS(NBDJ)-(EN1(J)-EN0(J)))/DT
!C     &              +TAU0VAR(NBDJ)*(QN1(J)-CELERITY*(ETA1(NBDJ)-EN1(J)))
!C DW, below a 2nd-order extrapolation formula: \eta^{n+1} = 2*\eta^{n} - \eta^{n-1}
!C     and second order in time derivative is used
                QFORCEJ=(QN2(J)-QN0(J))/DT2
     &              -CELERITY*(2.0*ETAS(NBDJ)-(EN2(J)-EN0(J)))/DT2
     &              +TAU0VAR(NBDJ)*(QN1(J)-CELERITY*(ETA1(NBDJ)-EN1(J)))
            ENDIF

            IF((LBCODEI(J).EQ.40).OR.(LBCODEI(J).EQ.41)) QFORCEJ=
     &        -(QN1(J)-QN0(J))/DT - TAU0VAR(NBDJ)*(QN1(J)+QN0(J))/2.e0
!
!     jgf46.21 Added IBTYPE=52
            IF(LBCODEI(J).EQ.52) THEN
               QFORCEJ=(QN2(J)-QN0(J))/DT2 + Tau0VAR(NBDJ)*QN1(J)
               IF (IT.GT.FluxSettlingIT) THEN
                  HH1=DP(NBDJ)+IFNLFA*Eta1(NBDJ)
                  Celerity=SQRT(G*HH1)
                  QforceJ=QforceJ - Celerity*(EtaS(NBDJ)/DT
     &                 + Tau0Var(NBDJ)*(Eta1(NBDJ)-ElevDisc(J)))
               ENDIF
            ENDIF

            NCI=NodeCode(NBDI)
            NCJ=NodeCode(NBDJ)
            BndLenO6NC=NCI*NCJ*BndLen2O3(J-1)/4.e0
            GWCE_LV(NBDI)=GWCE_LV(NBDI)
     &                              + BndLenO6NC*(2.e0*QForceI+QForceJ)
            GWCE_LV(NBDJ)=GWCE_LV(NBDJ)
     &                              + BndLenO6NC*(2.e0*QForceJ+QForceI)
          ENDDO
        ENDIF

!..... DW
      CALL UPDATE_RHS_SLPERBC() ;
!..... DW
!...
!...  IMPOSE ELEVATION BOUNDARY CONDITIONS TO LOAD VECTOR GWCE_LV(I) NOTE; EP
!...  IS THE RMS OF ALL THE DIAGONAL MEMBERS IN THE GWCE.  IT IS USED TO
!...  SCALE THE DIAGONAL ELEMENT FOR THE ELEVATION SPECIFIED BOUNDARY
!...  NODES AND THEREFORE MUST ALSO BE USED TO SCALE THE RHS OF THE
!...  EQUATIONS
!...
      IF ( ILump.eq.0 ) THEN ! default, fully consistent GWCE LHS
         DO I=1,NETA
            NBDI = NBD(I)
            ETAS(NBDI) = ETA2(NBDI)-ETA1(NBDI)
            GWCE_LV(NBDI) = ETAS(NBDI)*NODECODE(NBDI)*EP
            DO J=2,NNEIGH(NBDI)
               GWCE_LV(NEITAB(NBDI,J))=GWCE_LV(NEITAB(NBDI,J))
     &              -ETAS(NBDI)*OBCCOEF(I,J-1)
            END DO
         END DO
      ELSE                   ! ILump.eq.1, lumped GWCE
         DO I=1,NETA
            NBDI=NBD(I)
            ETAS(NBDI)=ETA2(NBDI)-ETA1(NBDI)
            GWCE_LV(NBDI) = ETAS(NBDI)*NODECODE(NBDI)*COEFD(NBDI)
         END DO
      ENDIF
!...
!...  SOLVE GWCE FOR ELEVATION AT NEW TIME LEVEL
!...

!...  UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE

#ifdef CMPI
!...UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE
      IF (ILump.eq.0) THEN ! default, fully consistent LHS
!... DW, periodic bcs
         IF(NPERSEG.GT.0)THEN
             CALL UPDATER_W_PERBC(GWCE_LV,COEF(:,1),DUMY1,2)
         ELSE
             CALL UPDATER(GWCE_LV,COEF(:,1),DUMY1,2)
         ENDIF
      ELSE  ! lumped LHS
         CALL UPDATER(GWCE_LV,COEFD,DUMY1,2)
      ENDIF
#endif

      if (subdomainOn.and.enforceBN.eq.2) call enforceEob() ! NCSU Subdomain
      if (subdomainOn.and.enforceBN.eq.2) call enforceEib() ! NCSU Subdomain
      if (subdomainOn.and.enforceBN.eq.2) call enforceGWCELVob() ! NCSU Subdomain

      IF (ILump.eq.0) THEN ! default, fully consistent LHS
!...  JCG ITERATIVE MATRIX SOLVER
         IPARM(1)=ITMAX
         CALL JCG(NP,MNP,MNEI,NEITAB,COEF,GWCE_LV,ETAS,
     &        IWKSP,NW,WKSP,IPARM,RPARM,IER)
         NUMITR=IPARM(1)
         DO I=1,NP
            ETA2(I) = NODECODE(I)*ETAS(I)+ETA1(I) !COMPUTE NEW ELEVATIONS
         END DO
      ELSE ! lumped LHS
         DO I = 1, NP
            IF (COEFD(I).eq.0.0e0) THEN
               RDIAG = 0.0e0
            ELSE
               RDIAG = 1.0e0 / COEFD(I)
            ENDIF
            ETAS(I) = GWCE_LV(I) * RDIAG
         ENDDO
         NUMITR=0
         DO I=1,NP
            ETA2(I)=NODECODE(I)*ETAS(I)+ETA1(I) !COMPUTE NEW ELEVATIONS
         END DO
      ENDIF

!...... DW, periodic bcs
      CALL UPDATE_PERSLNODES() ;
!     Switch back to the original element table 
      CALL SWITCH_ELTAB_PERBC( IDIREC = 2 ) ;
!......

#ifdef CMPI
      CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif      
     
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN

!**********************************************************************
      CONTAINS

!.... DW, perioidic boundary condition
        SUBROUTINE UPDATE_COEF_PERBCSL_CONSISTENT() 
        IMPLICIT NONE

        INTEGER I, I2

!...  Modify equations associated with the slave nodes. 
!...  For each slave node, zero all off diagonal
!......terms on the row and set diagnoal term to EP
        IF ( NPERSEG > 0 ) THEN
           DO I = 1, NNPERBC
              I2 = IPERCONN(I,2) ;
              Coef(I2,1)=EP ;
              DO J = 2, NNEIGH(I2)
                 Coef(I2,J) = 0.0E0 ;
              END DO
           END DO
        END IF

        RETURN ;
        END SUBROUTINE UPDATE_COEF_PERBCSL_CONSISTENT      

        SUBROUTINE UPDATE_COEF_PERBCSL_LUMPED() 
        IMPLICIT NONE 
      
        INTEGER:: I, I2 
        
        IF ( NPERSEG > 0 ) THEN
           DO I = 1, NNPERBC
              I2 = IPERCONN(I,2) ;
              Coefd(I2)=EP ;  
           END DO
        END IF

        RETURN ;
        END SUBROUTINE  UPDATE_COEF_PERBCSL_LUMPED
      
        SUBROUTINE UPDATE_RHS_SLPERBC() 
        IMPLICIT NONE
        
        INTEGER:: I, I2
        
        IF ( NPERSEG > 0 ) THEN
           DO I = 1, NNPERBC
              I2 = IPERCONN(I,2) ;
              ETAS(I2) = 0.0E0   ; !  ETA1 ==> ETA2  
              GWCE_LV(I2)=ETAS(I2)*NODECODE(I2)*EP
           END DO
        END IF
        
        RETURN ;
        END SUBROUTINE UPDATE_RHS_SLPERBC 
      
        SUBROUTINE UPDATE_PERSLNODES()
        IMPLICIT NONE
       
        INTEGER:: I, I1, I2
        IF ( NPERSEG > 0 ) THEN
!     Update values of slave node ! 
           DO I = 1, NNPERBC
              I1 = IPERCONN(I,1) ; 
              I2 = IPERCONN(I,2) ;
              
              ETAS(I2)=ETAS(I1) ; 
              ETA2(I2)=ETA2(I1) ;
           END DO
        END IF 
        
        RETURN ;
        END SUBROUTINE UPDATE_PERSLNODES
!..... DW
!.....
      END SUBROUTINE GWCE_NEW
!**********************************************************************


!*******************************************************************************
!                                                                              *
!    Subroutine to compute the elevation using the GWCE formluation            *
!    This subroutine is the corrector step for the predictor-corrector         *
!    algorithm and obtains the corrected elevations                            *
!    Re-written to conform to the ADCIRC Theory Report                         *
!                                                                              *
!                            k.d.  06/24/2004                                  *
!                            r.l.  06/22/2005                                  *
!*******************************************************************************

      SUBROUTINE GWCE_new_pc(IT,TimeLoc,TimeH)
!
      USE SIZES
      USE GLOBAL, ONLY : GWCE_LV, COEF, COEFD, OBCCOEF, ETAS, ETA1, ETA2,
     &   ESBIN1, ESBIN2, C2DDI, C3D, CBaroclinic, CGWCE_Advec_C1,
     &   CGWCE_Advec_C2, CGWCE_Advec_NC, CGWCE_LS_2PartQ, CGWCE_LS_2PartSQ,
     &   CGWCE_LS_2PartV, CGWCE_LS_2PartSV, CGWCE_LS_KGQ, CPRECOR,
     &   CSmag_Eh, Smag_Comp_Flag, CTIP, DT, FluxSettlingIT, G, IFNLCAT,
     &   IFNLCT, IFNLFA, NBFR, NODECODE, NRS, NSCREEN, NCCHANGE, IDEN,
     &   NWS, RAMPELEV, RAMPMETE, RES_BC_FLAG, RHOWAT0, screenUnit,
     &   SMAG_LOWER_LIM, SMAG_UPPER_LIM, TK, NOFF, UU1, VV1, QX1, QY1,
     &   PR1, TiP1, BSX1, BSY1, DUU1, DUV1, DVV1, WSX1, WSY1, VIDBCPDXOH,
     &   VIDBCPDYOH, CORIF, PER, FACE, FF, EFA, EMO, AMIG, PR2, LNM_BC,
     &   QN1, EN1, QN0, ElevDisc, QN2, EN0, iLump, BCFLAG_LNM, EN2,
     &   TKM, NPERSEG, NNPERBC, IPERCONN, UU0, VV0, QX0, QY0, TK0, TK2, 
     &   QX2, QY2, UU2, VV2
#ifdef CMPI
     &   , dumy1, dumy2
#endif
      USE MESH, ONLY : NE, NP, NM, DP, X, Y, TotalArea, Areas, NeiTab,
     &                  NNEIGH, SFAC
      USE BOUNDARIES, ONLY : NETA, NFLUXB, NFLUXF, NFLUXGBC, NFLUXRBC,
     &   NFLUXIB, NOPE, NVEL, NBD, NBV, LBCODEI, BndLen2O3, NPEBC
      USE ITPACKV
      USE NodalAttributes, ONLY : FRIC, Tau0Var, HBREAK, FTHETA, FGAMMA,
     &     IFLINBF, IFNLBF, IFHYBF, EVM, LoadAdvectionState, advectlocal
#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE

      INTEGER IE, I, J                           !local loop counters
      INTEGER IT
      INTEGER  NM1, NM2, NM3
      INTEGER NC1, NC2, NC3, NCEle, NCI, NCJ
      INTEGER NCyc
      INTEGER NBDI
      INTEGER OnDiag, OffDiag

      REAL(kind=8) A00pB00
      REAL(kind=8) BCXAvg, BCYAvg
      REAL(kind=8) BndLenO6NC    !BNDLEN2O3NC, NCBND need to be removed from global.f and put in original GWCE subroutine
      REAL(kind=8) BSXN1, BSXN2, BSXN3, BSYN1, BSYN2, BSYN3, BSXAvg, BSYAvg
      REAL(kind=8) CorifAvg
      REAL(kind=8) DPAvg, GDPAvgOAreaIE4
      REAL(kind=8) DispX, DispY, DispXAvg, DispYAvg
      REAL(kind=8) E0N1, E0N2, E0N3, E0XGrad2A, E0YGrad2A
      REAL(kind=8) E1N1, E1N2, E1N3, E1XGrad2A, E1YGrad2A
      REAL(kind=8) E1N1SQ, E1N2SQ, E1N3SQ
      REAL(kind=8) ESN1, ESN2, ESN3, ESAvg
      REAL(kind=8) EVMN1, EVMN2, EVMN3, EVMXGrad, EVMYGrad, EVMAvgODT
      REAL(kind=8) EVMEle, EVMSmag
      REAL(kind=8) GHAvg, GHAvgOAreaIE2, GOAreaIE4
      REAL(kind=8) H1N1, H1N2, H1N3, H2N1, H2N2, H2N3, HAvg, H1, H2
      REAL(kind=8) H2OTotalArea
      REAL(kind=8) LSXXGradA, LSXYGradA, LSYXGradA, LSYYGradA
      REAL(kind=8) LSXXEle, LSXYEle, LSYXEle, LSYYEle
      REAL(kind=8) MsFacR
      REAL(kind=8) MX, MY, MXAvg, MYAvg
      REAL(kind=8) JXAvg, JYAvg
      REAL(kind=8) Pr1N1, Pr1N2, Pr1N3
      REAL(kind=8) QX1N1, QX1N2, QX1N3, QY1N1, QY1N2, QY1N3, QX1Avg, QY1Avg
      REAL(kind=8) SFacAvg
      REAL(kind=8) T0N1,T0N2, T0N3
      REAL(kind=8) Tau0Avg, Tau0QXAvg, Tau0QYAvg
      REAL(kind=8) Tau0XGrad2A, Tau0YGrad2A, Tau0SpaVar
      REAL(kind=8) TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3
      REAL(kind=8) TiPN1, TiPN2, TiPN3
      REAL(kind=8) UV0, UV1, UV2
      REAL(kind=8) U1N1,U1N2,U1N3, U1Avg
      REAL(kind=8) V1N1,V1N2,V1N3, V1Avg
      REAL(kind=8) WSXAvg, WSYAvg

      REAL(8) AreaIE, AreaIE2, AreaIE4
      REAL(8) FDX1, FDX2, FDX3, FDY1, FDY2, FDY3
      REAL(8) TimeLoc, TimeH

!kmd   Added in parameters for the pc algorithm
      REAL(kind=8) BSX0N1, BSX0N2, BSX0N3, BSY0N1
      REAL(kind=8) BSY0N2, BSY0N3, BSX0Avg, BSY0Avg
      REAL(kind=8) BSX2N1, BSX2N2, BSX2N3, BSY2N1
      REAL(kind=8) BSY2N2, BSY2N3, BSX2Avg, BSY2Avg
      REAL(kind=8) E2N1,E2N2,E2N3
      REAL(kind=8) E0N1SQ, E0N2SQ, E0N3SQ
      REAL(kind=8) E2N1SQ, E2N2SQ, E2N3SQ
      REAL(kind=8) H0N1, H0N2, H0N3, H00
      REAL(kind=8) QX0N1, QX0N2, QX0N3, QY0N1, QY0N2, QY0N3, QX0Avg, QY0Avg
      REAL(kind=8) QX2N1, QX2N2, QX2N3, QY2N1, QY2N2, QY2N3, QX2Avg, QY2Avg
      REAL(kind=8) Tau0QX0Avg, Tau0QY0Avg, Tau0QX2Avg, Tau0QY2Avg
      REAL(kind=8) Tau0SpaVar0, Tau0SpaVar2
      REAL(kind=8) U0N1,U0N2,U0N3, U0Avg
      REAL(kind=8) V0N1,V0N2,V0N3, V0Avg
      REAL(kind=8) U2N1,U2N2,U2N3, U2Avg
      REAL(kind=8) V2N1,V2N2,V2N3, V2Avg
      REAL(kind=8) timewtgwce0,timewtgwce1,timewtgwce2,timeagflag
      REAL(kind=8) :: arg, argj
      REAL(kind=8) :: celerity
      REAL(kind=8) :: etratio
      INTEGER :: ier
      INTEGER :: nbdj
      REAL(kind=8) :: rff
      REAL(kind=8) :: qforcei, qforcej

      call setMessageSource("gwce_new_pc")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

!kmd    Must reset the result vector to zero before recomputing
!kmd  the next time level.
      DO I=1,NP
         GWCE_LV(I) =0.E0
      END DO

!     Consistent mass matrix: ILump=0, lumped mass matrix: ILump=1
!     Re-compute these local values

      OnDiag=(1+ILump)*2                         !diagonal coefficient
      OffDiag=(1-ILump)                          !off diagonal coefficient

!...
!...  Compute the GWCE load vector GWCE_LV
!...  This is done primarily element by element by forming
!...  temporary vectors and then assembling at the end.
!...  This has been set up to unroll loops to optimize performance
!...  on vector processors.
!...
!...  Elevation and flux boundary conditions are imposed after the
!...  element by element assembly section.
!...

!...  Initialize variables to zero if these forcings are not used

      IF((NWS.NE.0).OR.(NRS.NE.0)) THEN
      ELSE
        WSXAvg=0.e0
        WSYAvg=0.e0
        Pr1N1=0.e0
        Pr1N2=0.e0
        Pr1N3=0.e0
      ENDIF

      IF (CTIP) THEN
      ELSE
        TiPN1=0.e0
        TiPN2=0.e0
        TiPN3=0.e0
      ENDIF

      IF(C3D) THEN
      ELSE
        DispXAvg=0.e0
        DispYAvg=0.e0
      ENDIF

      IF(CBaroclinic) THEN
      ELSE
        BCXAvg=0.e0
        BCYAvg=0.e0
      ENDIF

!...  Compute the Lateral Stress Field using the 2 Part velocity approach (nonsymmetric or symmetric)

      IF ((CGWCE_LS_2PartV).OR.(CGWCE_LS_2PartSV)) THEN

         DO I=1,NP
            LSXX(I)=0.e0
            LSXY(I)=0.e0
            LSYX(I)=0.e0
            LSYY(I)=0.e0
         ENDDO

         DO IE=1,NE
            NM1=NM(IE,1)
            NM2=NM(IE,2)
            NM3=NM(IE,3)
            NC1=NODECODE(NM1)
            NC2=NODECODE(NM2)
            NC3=NODECODE(NM3)
            NCEle=NC1*NC2*NC3*NOFF(IE)
            SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.e0
            FDX1 = (Y(NM2)-Y(NM3))*SFacAvg               !b1
            FDX2 = (Y(NM3)-Y(NM1))*SFacAvg               !b2
            FDX3 = (Y(NM1)-Y(NM2))*SFacAvg               !b3
            FDY1 = X(NM3)-X(NM2)                         !a1
            FDY2 = X(NM1)-X(NM3)                         !a2
            FDY3 = X(NM2)-X(NM1)                         !a3
            LSXXGradA=(UU1(NM1)*FDX1+UU1(NM2)*FDX2+UU1(NM3)*FDX3)/2.e0   !A*DUDX
            LSXYGradA=(UU1(NM1)*FDY1+UU1(NM2)*FDY2+UU1(NM3)*FDY3)/2.e0   !A*DUDY
            LSYXGradA=(VV1(NM1)*FDX1+VV1(NM2)*FDX2+VV1(NM3)*FDX3)/2.e0   !A*DVDX
            LSYYGradA=(VV1(NM1)*FDY1+VV1(NM2)*FDY2+VV1(NM3)*FDY3)/2.e0   !A*DVDY
            EVMEle=NCEle*(EVM(NM1)+EVM(NM2)+EVM(NM3))/3.e0
            ! If using Smagorinski vertically-integrated lateral stress coefficient
            IF ((CSmag_Eh).or.(Smag_Comp_Flag)) THEN  
               EVMSmag=EVMEle*
     &             sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA)
     &                 +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
               !tcm v52.30.01 added test for limits
               IF(Smag_Comp_Flag) then
                  IF (EVMSmag .gt. SMAG_UPPER_LIM) EVMSmag = SMAG_UPPER_LIM
                  IF (EVMSmag .lt. SMAG_LOWER_LIM) EVMSmag = SMAG_LOWER_LIM
               ENDIF
               EVMEle=EVMSmag
            ENDIF
            LSXXEle = LSXXGradA*EVMEle
            LSXX(NM1)=LSXX(NM1)+LSXXEle
            LSXX(NM2)=LSXX(NM2)+LSXXEle
            LSXX(NM3)=LSXX(NM3)+LSXXEle
            LSXYEle = LSXYGradA*EVMEle
            LSXY(NM1)=LSXY(NM1)+LSXYEle
            LSXY(NM2)=LSXY(NM2)+LSXYEle
            LSXY(NM3)=LSXY(NM3)+LSXYEle
            LSYXEle = LSYXGradA*EVMEle
            LSYX(NM1)=LSYX(NM1)+LSYXEle
            LSYX(NM2)=LSYX(NM2)+LSYXEle
            LSYX(NM3)=LSYX(NM3)+LSYXEle
            LSYYEle = LSYYGradA*EVMEle
            LSYY(NM1)=LSYY(NM1)+LSYYEle
            LSYY(NM2)=LSYY(NM2)+LSYYEle
            LSYY(NM3)=LSYY(NM3)+LSYYEle
         ENDDO

         DO I=1,NP
            IF(TotalArea(I).NE.0.) THEN
               H2=DP(I)+IFNLFA*ETA2(I)
               H2OTotalArea=H2/TotalArea(I)
               IF (CGWCE_LS_2PartV) THEN          !nonsymmetric
                  LSXX(I)=H2OTotalArea*LSXX(I)
                  LSXY(I)=H2OTotalArea*LSXY(I)
                  LSYX(I)=H2OTotalArea*LSYX(I)
                  LSYY(I)=H2OTotalArea*LSYY(I)
               ENDIF
               IF (CGWCE_LS_2PartSV) THEN         !symmetric
                  LSXX(I)=H2OTotalArea*LSXX(I)
                  LSXY(I)=0.5e0*H2OTotalArea*(LSXY(I)+LSYX(I))
                  LSYX(I)=LSXY(I)
                  LSYY(I)=H2OTotalArea*LSYY(I)
               ENDIF
            ELSE
               LSXX(I)=0.e0
               LSXY(I)=0.e0
               LSYX(I)=0.e0
               LSYY(I)=0.e0
            ENDIF
         ENDDO

      ENDIF

!...  Compute the Lateral Stress Field using the 2 Part flux approach (nonsymmetric or symmetric)

      IF ((CGWCE_LS_2PartQ).OR.(CGWCE_LS_2PartSQ)) THEN

         DO I=1,NP
            LSXX(I)=0.e0
            LSXY(I)=0.e0
            LSYX(I)=0.e0
            LSYY(I)=0.e0
         ENDDO

         DO IE=1,NE
            NM1=NM(IE,1)
            NM2=NM(IE,2)
            NM3=NM(IE,3)
            NC1=NODECODE(NM1)
            NC2=NODECODE(NM2)
            NC3=NODECODE(NM3)
            NCEle=NC1*NC2*NC3*NOFF(IE)
            SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.e0
            FDX1 = (Y(NM2)-Y(NM3))*SFacAvg               !b1
            FDX2 = (Y(NM3)-Y(NM1))*SFacAvg               !b2
            FDX3 = (Y(NM1)-Y(NM2))*SFacAvg               !b3
            FDY1 = X(NM3)-X(NM2)                         !a1
            FDY2 = X(NM1)-X(NM3)                         !a2
            FDY3 = X(NM2)-X(NM1)                         !a3
            EVMEle=NCEle*(EVM(NM1)+EVM(NM2)+EVM(NM3))/3.e0
            ! If using Smagorinski vertically-integrated lateral stress coefficient
            IF ((CSmag_Eh).or.(Smag_Comp_Flag)) THEN
               LSXXGradA=(UU1(NM1)*FDX1+UU1(NM2)*FDX2+UU1(NM3)*FDX3)/2.e0
               LSXYGradA=(UU1(NM1)*FDY1+UU1(NM2)*FDY2+UU1(NM3)*FDY3)/2.e0
               LSYXGradA=(VV1(NM1)*FDX1+VV1(NM2)*FDX2+VV1(NM3)*FDX3)/2.e0
               LSYYGradA=(VV1(NM1)*FDY1+VV1(NM2)*FDY2+VV1(NM3)*FDY3)/2.e0
               EVMSmag=EVMEle*
     &            sqrt((LSXXGradA-LSYYGradA)*(LSXXGradA-LSYYGradA)
     &                +(LSYXGradA+LSXYGradA)*(LSYXGradA+LSXYGradA))
               !tcm v52.30.01 added test for limits
               IF(Smag_Comp_Flag) then
                  IF (EVMSmag .gt. SMAG_UPPER_LIM) EVMSmag = SMAG_UPPER_LIM
                  IF (EVMSmag .lt. SMAG_LOWER_LIM) EVMSmag = SMAG_LOWER_LIM
               ENDIF
               EVMEle=EVMSmag
            ENDIF
            LSXXGradA=(QX1(NM1)*FDX1+QX1(NM2)*FDX2+QX1(NM3)*FDX3)/2.e0
            LSXYGradA=(QX1(NM1)*FDY1+QX1(NM2)*FDY2+QX1(NM3)*FDY3)/2.e0
            LSYXGradA=(QY1(NM1)*FDX1+QY1(NM2)*FDX2+QY1(NM3)*FDX3)/2.e0
            LSYYGradA=(QY1(NM1)*FDY1+QY1(NM2)*FDY2+QY1(NM3)*FDY3)/2.e0
            LSXXEle = LSXXGradA*EVMEle
            LSXX(NM1)=LSXX(NM1)+LSXXEle
            LSXX(NM2)=LSXX(NM2)+LSXXEle
            LSXX(NM3)=LSXX(NM3)+LSXXEle
            LSXYEle = LSXYGradA*EVMEle
            LSXY(NM1)=LSXY(NM1)+LSXYEle
            LSXY(NM2)=LSXY(NM2)+LSXYEle
            LSXY(NM3)=LSXY(NM3)+LSXYEle
            LSYXEle = LSYXGradA*EVMEle
            LSYX(NM1)=LSYX(NM1)+LSYXEle
            LSYX(NM2)=LSYX(NM2)+LSYXEle
            LSYX(NM3)=LSYX(NM3)+LSYXEle
            LSYYEle = LSYYGradA*EVMEle
            LSYY(NM1)=LSYY(NM1)+LSYYEle
            LSYY(NM2)=LSYY(NM2)+LSYYEle
            LSYY(NM3)=LSYY(NM3)+LSYYEle
         ENDDO

         DO I=1,NP
            IF(TotalArea(I).NE.0.) THEN
               IF (CGWCE_LS_2PartQ) THEN          !nonsymmetric
                  LSXX(I)=LSXX(I)/TotalArea(I)
                  LSXY(I)=LSXY(I)/TotalArea(I)
                  LSYX(I)=LSYX(I)/TotalArea(I)
                  LSYY(I)=LSYY(I)/TotalArea(I)
               ENDIF
               IF (CGWCE_LS_2PartSQ) THEN         !symmetric
                  LSXX(I)=LSXX(I)/TotalArea(I)
                  LSXY(I)=0.5e0*(LSXY(I)+LSYX(I))/TotalArea(I)
                  LSYX(I)=LSXY(I)
                  LSYY(I)=LSYY(I)/TotalArea(I)
               ENDIF
            ELSE
               LSXX(I)=0.e0
               LSXY(I)=0.e0
               LSYX(I)=0.e0
               LSYY(I)=0.e0
            ENDIF
         ENDDO

      ENDIF

      DO I=1,NP
!kmd
!kmd  Added in the three time levels for the tau term.
!kmd  Every term is updated for the three time levels.
!kmd
            UV0=SQRT(UU0(I)*UU0(I)+VV0(I)*VV0(I))
            UV1=SQRT(UU1(I)*UU1(I)+VV1(I)*VV1(I))
            UV2=SQRT(UU2(I)*UU2(I)+VV2(I)*VV2(I))
            H00=DP(I)+IFNLFA*ETA0(I)
            H1=DP(I)+IFNLFA*ETA1(I)
            H2=DP(I)+IFNLFA*ETA2(I)
            TK0(I)=FRIC(I)*(IFLINBF + (UV0/H00)*(IFNLBF + IFHYBF*
     &           (1+(HBREAK/H00)**FTHETA)**(FGAMMA/FTHETA)))
            TK(I)=FRIC(I)*(IFLINBF + (UV1/H1)*(IFNLBF + IFHYBF*
     &           (1+(HBREAK/H1)**FTHETA)**(FGAMMA/FTHETA)))
            TK2(I)=FRIC(I)*(IFLINBF + (UV2/H2)*(IFNLBF + IFHYBF*
     &           (1+(HBREAK/H2)**FTHETA)**(FGAMMA/FTHETA)))
       END DO

!kmd      Added in the time weights
!...     Time weights for the nonlinear terms in the GWCE for
!...       the corrector step
          timewtgwce0=0.33e0
          timewtgwce1=0.34e0
          timewtgwce2=0.33e0
          timeagflag=1.0e0

!...  Assemble the GWCE RHS except for the boundary integral terms
!kmd  Renumber the GWCE loop for the corrector step

      DO 1038 IE=1,NE

!...     Set nodal values for each element
!kmd
!kmd  Define the needed product terms at three time levels
!kmd


!orbitt 120322: Localized Advection
         IF(LoadAdvectionState)CALL ADVECTLOCAL(IE)

         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
         NCELE=NC1*NC2*NC3*NOFF(IE)
         E0N1=ETA0(NM1)
         E0N2=ETA0(NM2)
         E0N3=ETA0(NM3)
         E1N1=ETA1(NM1)
         E1N2=ETA1(NM2)
         E1N3=ETA1(NM3)
         E2N1=ETA2(NM1)
         E2N2=ETA2(NM2)
         E2N3=ETA2(NM3)
         E0N1SQ=E0N1*E0N1
         E0N2SQ=E0N2*E0N2
         E0N3SQ=E0N3*E0N3
         E1N1SQ=E1N1*E1N1
         E1N2SQ=E1N2*E1N2
         E1N3SQ=E1N3*E1N3
         E2N1SQ=E2N1*E2N1
         E2N2SQ=E2N2*E2N2
         E2N3SQ=E2N3*E2N3
         ESN1=ETAS0(NM1)
         ESN2=ETAS0(NM2)
         ESN3=ETAS0(NM3)
         U0N1=UU0(NM1)
         U0N2=UU0(NM2)
         U0N3=UU0(NM3)
         U1N1=UU1(NM1)
         U1N2=UU1(NM2)
         U1N3=UU1(NM3)
         U2N1=UU2(NM1)
         U2N2=UU2(NM2)
         U2N3=UU2(NM3)
         V0N1=VV0(NM1)
         V0N2=VV0(NM2)
         V0N3=VV0(NM3)
         V1N1=VV1(NM1)
         V1N2=VV1(NM2)
         V1N3=VV1(NM3)
         V2N1=VV2(NM1)
         V2N2=VV2(NM2)
         V2N3=VV2(NM3)
         QX0N1=QX0(NM1)
         QX0N2=QX0(NM2)
         QX0N3=QX0(NM3)
         QX1N1=QX1(NM1)
         QX1N2=QX1(NM2)
         QX1N3=QX1(NM3)
         QX2N1=QX2(NM1)
         QX2N2=QX2(NM2)
         QX2N3=QX2(NM3)
         QY0N1=QY0(NM1)
         QY0N2=QY0(NM2)
         QY0N3=QY0(NM3)
         QY1N1=QY1(NM1)
         QY1N2=QY1(NM2)
         QY1N3=QY1(NM3)
         QY2N1=QY2(NM1)
         QY2N2=QY2(NM2)
         QY2N3=QY2(NM3)
         H0N1=DP(NM1)+IFNLFA*E0N1
         H0N2=DP(NM2)+IFNLFA*E0N2
         H0N3=DP(NM3)+IFNLFA*E0N3
         H1N1=DP(NM1)+IFNLFA*E1N1
         H1N2=DP(NM2)+IFNLFA*E1N2
         H1N3=DP(NM3)+IFNLFA*E1N3
         H2N1=DP(NM1)+IFNLFA*E2N1
         H2N2=DP(NM2)+IFNLFA*E2N2
         H2N3=DP(NM3)+IFNLFA*E2N3
         EVMN1=EVM(NM1)
         EVMN2=EVM(NM2)
         EVMN3=EVM(NM3)
         T0N1=Tau0Var(NM1)
         T0N2=Tau0Var(NM2)
         T0N3=Tau0Var(NM3)
         IF((NWS.NE.0).OR.(NRS.NE.0)) THEN     !wind or radiation stress
            Pr1N1=PR1(NM1)
            Pr1N2=PR1(NM2)
            Pr1N3=PR1(NM3)
         ENDIF
         IF (CTIP) THEN                        !tidal potential
            TiPN1=TiP1(NM1)
            TiPN2=TiP1(NM2)
            TiPN3=TiP1(NM3)
         ENDIF
         IF (C2DDI) THEN                       !2D bottom friction
            BSX0N1=TK0(NM1)*QX0N1
            BSY0N1=TK0(NM1)*QY0N1
            BSX0N2=TK0(NM2)*QX0N2
            BSY0N2=TK0(NM2)*QY0N2
            BSX0N3=TK0(NM3)*QX0N3
            BSY0N3=TK0(NM3)*QY0N3
            BSXN1=TK(NM1)*QX1N1
            BSYN1=TK(NM1)*QY1N1
            BSXN2=TK(NM2)*QX1N2
            BSYN2=TK(NM2)*QY1N2
            BSXN3=TK(NM3)*QX1N3
            BSYN3=TK(NM3)*QY1N3
            BSX2N1=TK2(NM1)*QX2N1
            BSY2N1=TK2(NM1)*QY2N1
            BSX2N2=TK2(NM2)*QX2N2
            BSY2N2=TK2(NM2)*QY2N2
            BSX2N3=TK2(NM3)*QX2N3
            BSY2N3=TK2(NM3)*QY2N3

         ENDIF
         IF (C3D) THEN                         !3D bottom friction
            BSXN1=BSX1(NM1)
            BSXN2=BSX1(NM2)
            BSXN3=BSX1(NM3)
            BSYN1=BSY1(NM1)
            BSYN2=BSY1(NM2)
            BSYN3=BSY1(NM3)
         ENDIF

         AreaIE2=Areas(IE)               !2A
         AreaIE=AreaIE2/2.e0             ! A
         AreaIE4=2.e0*AreaIE2            !4A

         SFacAvg=(SFAC(NM1)+SFAC(NM2)+SFAC(NM3))/3.e0

         FDX1 = (Y(NM2)-Y(NM3))*SFacAvg !b1 = 2*Area*dphi1/dx
         FDX2 = (Y(NM3)-Y(NM1))*SFacAvg !b2 = 2*Area*dphi2/dx
         FDX3 = (Y(NM1)-Y(NM2))*SFacAvg !b3 = 2*Area*dphi3/dx
         FDY1 = X(NM3)-X(NM2)           !a1 = 2*Area*dphi1/dy
         FDY2 = X(NM1)-X(NM3)           !a2 = 2*Area*dphi2/dy
         FDY3 = X(NM2)-X(NM1)           !a3 = 2*Area*dphi3/dy

!...     Compute part of several spatial gradients for use below

         E0XGrad2A=E0N1*FDX1+E0N2*FDX2+E0N3*FDX3        !2*Area*deta0/dx
         E0YGrad2A=E0N1*FDY1+E0N2*FDY2+E0N3*FDY3        !2*Area*deta0/dy
         E1XGrad2A=E1N1*FDX1+E1N2*FDX2+E1N3*FDX3        !2*Area*deta1/dx
         E1YGrad2A=E1N1*FDY1+E1N2*FDY2+E1N3*FDY3        !2*Area*deta1/dy
         Tau0XGrad2A=T0N1*FDX1+T0N2*FDX2+T0N3*FDX3      !2*Area*dTau0/dx
         Tau0YGrad2A=T0N1*FDY1+T0N2*FDY2+T0N3*FDY3      !2*Area*dTau0/dy

!...     Compute the Kolar & Gray lateral stress term extended for spatially varying EVM

         IF(CGWCE_LS_KGQ) THEN
            EVMXGrad=(EVMN1*FDX1+EVMN2*FDX2+EVMN3*FDX3)/AreaIE2
            EVMYGrad=(EVMN1*FDY1+EVMN2*FDY2+EVMN3*FDY3)/AreaIE2
            EVMAvgODT=((EVMN1+EVMN2+EVMN3)/3.e0)/DT
            MX=(EVMXGrad*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)
     &         +EVMYGrad*(QY1N1*FDX1+QY1N2*FDX2+QY1N3*FDX3)
     &         -EVMAvgODT*(ESN1*FDX1+ESN2*FDX2+ESN3*FDX3))/AreaIE2
            MY=(EVMXGrad*(QX1N1*FDY1+QX1N2*FDY2+QX1N3*FDY3)
     &         +EVMYGrad*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3)
     &         -EVMAvgODT*(ESN1*FDY1+ESN2*FDY2+ESN3*FDY3))/AreaIE2
         ENDIF

!...     Compute the remainder of the 2 Part lateral stress terms

         IF((CGWCE_LS_2PartQ).OR.(CGWCE_LS_2PartV)) THEN
            MX=(LSXX(NM1)*FDX1+LSXX(NM2)*FDX2+LSXX(NM3)*FDX3
     &         +LSXY(NM1)*FDY1+LSXY(NM2)*FDY2+LSXY(NM3)*FDY3)/AreaIE2
            MY=(LSYX(NM1)*FDX1+LSYX(NM2)*FDX2+LSYX(NM3)*FDX3
     &         +LSYY(NM1)*FDY1+LSYY(NM2)*FDY2+LSYY(NM3)*FDY3)/AreaIE2
         ENDIF

!...     Compute the spatial gradients of the velocity dispersion terms if 3D
         IF (C3D) THEN !3D velocity dispersion
            DispX=(DUU1(NM1)*FDX1+DUU1(NM2)*FDX2+DUU1(NM3)*FDX3
     &            +DUV1(NM1)*FDY1+DUV1(NM2)*FDY2+DUV1(NM3)*FDY3)/AreaIE2
            DispY=(DUV1(NM1)*FDX1+DUV1(NM2)*FDX2+DUV1(NM3)*FDX3
     &            +DVV1(NM1)*FDY1+DVV1(NM2)*FDY2+DVV1(NM3)*FDY3)/AreaIE2
         ENDIF

!...     Compute elemental averages

         CorifAvg=(Corif(NM1)+Corif(NM2)+Corif(NM3))/3.e0
         Tau0Avg=(T0N1+T0N2+T0N3)/3.e0
         Tau0QX0Avg=(T0N1*QX0N1+T0N2*QX0N2+T0N3*QX0N3)/3.e0
         Tau0QY0Avg=(T0N1*QY0N1+T0N2*QY0N2+T0N3*QY0N3)/3.e0
         Tau0QXAvg=(T0N1*QX1N1+T0N2*QX1N2+T0N3*QX1N3)/3.e0
         Tau0QYAvg=(T0N1*QY1N1+T0N2*QY1N2+T0N3*QY1N3)/3.e0
         Tau0QX2Avg=(T0N1*QX2N1+T0N2*QX2N2+T0N3*QX2N3)/3.e0
         Tau0QY2Avg=(T0N1*QY2N1+T0N2*QY2N2+T0N3*QY2N3)/3.e0
         U0Avg=(U0N1+U0N2+U0N3)/3.e0
         V0Avg=(V0N1+V0N2+V0N3)/3.e0
         U1Avg=(U1N1+U1N2+U1N3)/3.e0
         V1Avg=(V1N1+V1N2+V1N3)/3.e0
         U2Avg=(U2N1+U2N2+U2N3)/3.e0
         V2Avg=(V2N1+V2N2+V2N3)/3.e0
         QX0Avg=(QX0N1+QX0N2+QX0N3)/3.e0
         QY0Avg=(QY0N1+QY0N2+QY0N3)/3.e0
         QX1Avg=(QX1N1+QX1N2+QX1N3)/3.e0
         QY1Avg=(QY1N1+QY1N2+QY1N3)/3.e0
         QX2Avg=(QX2N1+QX2N2+QX2N3)/3.e0
         QY2Avg=(QY2N1+QY2N2+QY2N3)/3.e0
         ESAvg=(ESN1+ESN2+ESN3)/3.e0
         DPAvg=(DP(NM1)+DP(NM2)+DP(NM3))/3.e0
         GDPAvgOAreaIE4=G*DPAvg/AreaIE4
         HAvg=(H1N1+H1N2+H1N3)/3.e0
         GHAvg=G*HAvg
         GHAvgOAreaIE2=GHAvg/AreaIE2
         BSX0Avg=(BSX0N1+BSX0N2+BSX0N3)/3.e0
         BSY0Avg=(BSY0N1+BSY0N2+BSY0N3)/3.e0
         BSXAvg=(BSXN1+BSXN2+BSXN3)/3.e0
         BSYAvg=(BSYN1+BSYN2+BSYN3)/3.e0
         BSX2Avg=(BSX2N1+BSX2N2+BSX2N3)/3.e0
         BSY2Avg=(BSY2N1+BSY2N2+BSY2N3)/3.e0
         MXAvg=MX           !lateral stresses are constant over an element
         MYAvg=MY           !lateral stresses are constant over an element
         IF((NWS.NE.0).OR.(NRS.NE.0)) THEN     !wind or radiation stress
            WSXAvg=(WSX1(NM1)+WSX1(NM2)+WSX1(NM3))/3.e0
            WSYAvg=(WSY1(NM1)+WSY1(NM2)+WSY1(NM3))/3.e0
         ENDIF
         IF (C3D.or.(CBaroclinic.and.abs(IDEN).eq.5)) THEN !3D velocity dispersion
            DispXAvg=IFNLCT*DispX
            DispYAvg=IFNLCT*DispY
         ENDIF
         IF(CBaroclinic) THEN
            BCXAvg=(H1N1*VIDBCPDXOH(NM1)+H1N2*VIDBCPDXOH(NM2)
     &                                       +H1N3*VIDBCPDXOH(NM3))/3.e0
            BCYAvg=(H1N1*VIDBCPDYOH(NM1)+H1N2*VIDBCPDYOH(NM2)
     &                                       +H1N3*VIDBCPDYOH(NM3))/3.e0
         ENDIF

!...     Compute additional partial factors

         MsFacR=AreaIE*(1.e0/DT-Tau0Avg/2.e0)/DT/12.e0
         GOAreaIE4=G/AreaIE4
         Tau0SpaVar0=(QX0Avg*Tau0XGrad2A+QY0Avg*Tau0YGrad2A)/6.e0
         Tau0SpaVar=(QX1Avg*Tau0XGrad2A+QY1Avg*Tau0YGrad2A)/6.e0
         Tau0SpaVar2=(QX2Avg*Tau0XGrad2A+QY2Avg*Tau0YGrad2A)/6.e0
         A00pB00=A00+B00

!...     Compute the JX, JY terms less the advection terms
         JXAvg = timewtgwce0*CorifAvg*QY0Avg+timewtgwce1*CorifAvg
     &                          *QY1Avg+timewtgwce2*CorifAvg*QY2Avg
     &          -IFNLFA*GOAreaIE4*timewtgwce0*(E0N1SQ*FDX1+E0N2SQ*FDX2
     &                                        +E0N3SQ*FDX3)
     &          -IFNLFA*GOAreaIE4*timewtgwce1*(E1N1SQ*FDX1+E1N2SQ*FDX2
     &                                        +E1N3SQ*FDX3)
     &          -IFNLFA*GOAreaIE4*timewtgwce2*(E2N1SQ*FDX1+E2N2SQ*FDX2
     &                                        +E2N3SQ*FDX3)
     &          -GHAvgOAreaIE2*((PR1N1-TiPN1)*FDX1
     &                     +(PR1N2-TiPN2)*FDX2+(PR1N3-TiPN3)*FDX3)
     &          +WSXAvg-timewtgwce0*BSX0Avg-timewtgwce1*BSXAvg
     &                                        -timewtgwce2*BSX2Avg
     &          +MXAvg-DispXAvg-BCXAvg
     &          +timewtgwce0*Tau0QX0Avg+timewtgwce1*Tau0QXAvg
     &                                 +timewtgwce2*Tau0QX2Avg
         JYAvg =-timewtgwce0*CorifAvg*QX0Avg-timewtgwce1*CorifAvg
     &                          *QX1Avg-timewtgwce2*CorifAvg*QX2Avg
     &          -IFNLFA*GOAreaIE4*timewtgwce0*(E0N1SQ*FDY1+E0N2SQ*FDY2
     &                                        +E0N3SQ*FDY3)
     &          -IFNLFA*GOAreaIE4*timewtgwce1*(E1N1SQ*FDY1+E1N2SQ*FDY2
     &                                        +E1N3SQ*FDY3)
     &          -IFNLFA*GOAreaIE4*timewtgwce2*(E2N1SQ*FDY1+E2N2SQ*FDY2
     &                                        +E2N3SQ*FDY3)
     &          -GHAvgOAreaIE2*((PR1N1-TiPN1)*FDY1
     &                    +(PR1N2-TiPN2)*FDY2+(PR1N3-TiPN3)*FDY3)
     &          +WSYAvg-timewtgwce0*BSY0Avg-timewtgwce1*BSYAvg
     &                                      -timewtgwce2*BSY2Avg
     &          +MYAvg-DispYAvg-BCYAvg
     &          +timewtgwce0*Tau0QY0Avg+timewtgwce1*Tau0QYAvg
     &                                    +timewtgwce2*Tau0QY2Avg

!...     Complete the JX, JY terms depending on the advection formulation
         IF(CGWCE_Advec_NC) THEN        !nonconservative advection
           JXAvg = JXAvg - IFNLCT*timewtgwce0*(
     &              QX0Avg*(U0N1*FDX1+U0N2*FDX2+U0N3*FDX3)
     &             +QY0Avg*(U0N1*FDY1+U0N2*FDY2+U0N3*FDY3))/AreaIE2
     &             - IFNLCT*timewtgwce1
     &             *(QX1Avg*(U1N1*FDX1+U1N2*FDX2+U1N3*FDX3)
     &             +QY1Avg*(U1N1*FDY1+U1N2*FDY2+U1N3*FDY3))/AreaIE2
     &             - IFNLCT*timewtgwce2
     &             *(QX2Avg*(U2N1*FDX1+U2N2*FDX2+U2N3*FDX3)
     &             +QY2Avg*(U2N1*FDY1+U2N2*FDY2+U2N3*FDY3))/AreaIE2
     &             +IFNLCAT*(timewtgwce0*U0Avg+timewtgwce1*U1Avg
     &             +timewtgwce2*U2Avg)
     &             *(timeagflag*0.5e0*((E2N1-E0N1+E2N2-E0N2
     &             +E2N3-E0N3)/DT)
     &             +(1.E0-timeagflag)*ESAvg/DT)
           JYAvg = JYAvg - IFNLCT*timewtgwce0*(
     &              QX0Avg*(V0N1*FDX1+V0N2*FDX2+V0N3*FDX3)
     &             +QY0Avg*(V0N1*FDY1+V0N2*FDY2+V0N3*FDY3))/AreaIE2
     &             - IFNLCT*timewtgwce1
     &             *(QX1Avg*(V1N1*FDX1+V1N2*FDX2+V1N3*FDX3)
     &             +QY1Avg*(V1N1*FDY1+V1N2*FDY2+V1N3*FDY3))/AreaIE2
     &             - IFNLCT*timewtgwce2
     &             *(QX2Avg*(V2N1*FDX1+V2N2*FDX2+V2N3*FDX3)
     &             +QY2Avg*(V2N1*FDY1+V2N2*FDY2+V2N3*FDY3))/AreaIE2
     &             +IFNLCAT*(timewtgwce0*V0Avg+timewtgwce1*V1Avg
     &             +timewtgwce2*V2Avg)
     &             *(timeagflag*0.5e0*((E2N1-E0N1+E2N2-E0N2
     &             +E2N3-E0N3)/DT)
     &             +(1.E0-timeagflag)*ESAvg/DT)
         ENDIF
         IF(CGWCE_Advec_C1) THEN        !conservative v1 advection
           JXAvg = JXAvg - IFNLCT*(
     &              (U1N1*QX1N1*FDX1+U1N2*QX1N2*FDX2
     &                              +U1N3*QX1N3*FDX3)
     &             +(U1N1*QY1N1*FDY1+U1N2*QY1N2*FDY2
     &                              +U1N3*QY1N3*FDY3))/AreaIE2
           JYAvg = JYAvg - IFNLCT*(
     &              (V1N1*QX1N1*FDX1+V1N2*QX1N2*FDX2
     &                              +V1N3*QX1N3*FDX3)
     &             +(V1N1*QY1N1*FDY1+V1N2*QY1N2*FDY2
     &                              +V1N3*QY1N3*FDY3))/AreaIE2
         ENDIF
         IF(CGWCE_Advec_C2) THEN        !conservative v2 advection
           JXAvg = JXAvg - IFNLCT*(
     &              QX1Avg*(U1N1*FDX1+U1N2*FDX2+U1N3*FDX3)
     &             +QY1Avg*(U1N1*FDY1+U1N2*FDY2+U1N3*FDY3)
     &             +U1Avg*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)
     &             +U1Avg*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3))/AreaIE2
           JYAvg = JYAvg - IFNLCT*(
     &              QX1Avg*(V1N1*FDX1+V1N2*FDX2+V1N3*FDX3)
     &             +QY1Avg*(V1N1*FDY1+V1N2*FDY2+V1N3*FDY3)
     &             +V1Avg*(QX1N1*FDX1+QX1N2*FDX2+QX1N3*FDX3)
     &             +V1Avg*(QY1N1*FDY1+QY1N2*FDY2+QY1N3*FDY3))/AreaIE2
         ENDIF


!...     Assemble forcing for node NM1 (local index j=1)

         Temp_LV_A1=

!...     Transient and Tau0 terms from LHS
     &         (OnDiag*ESN1 + OffDiag*(ESN2+ESN3))*MsFacR

!...     Free surface terms from LHS (time levels s-1 & s)
     &        -GDPAvgOAreaIE4*(  C00  *(FDX1*E0XGrad2A+FDY1*E0YGrad2A)
     &                        +A00pB00*(FDX1*E1XGrad2A+FDY1*E1YGrad2A))

!...     Terms from momentum eqs.
     &        +(JXAvg*FDX1 + JYAvg*FDY1)/2.e0

!...     Spatially varying Tau0 terms
     &        +timewtgwce0*Tau0SpaVar0+timewtgwce1*Tau0SpaVar
     &                    +timewtgwce2*Tau0SpaVar2


!...     Assemble forcing for node NM2 (local index j=2)

         Temp_LV_A2=

!...     Transient and Tau0 terms from LHS
     &         (OnDiag*ESN2 + OffDiag*(ESN1+ESN3))*MsFacR

!...     Free surface terms from LHS (time levels s-1 & s)
     &        -GDPAvgOAreaIE4*(  C00  *(FDX2*E0XGrad2A+FDY2*E0YGrad2A)
     &                        +A00pB00*(FDX2*E1XGrad2A+FDY2*E1YGrad2A))

!...     Terms from momentum eqs.
     &        +(JXAvg*FDX2 + JYAvg*FDY2)/2.e0

!...     Spatially varying Tau0 terms
     &        +timewtgwce0*Tau0SpaVar0+timewtgwce1*Tau0SpaVar
     &                   +timewtgwce2*Tau0SpaVar2


!...     Assemble forcing for node NM3 (local index j=3)

         Temp_LV_A3=

!...     Transient and Tau0 terms from LHS
!...    (consistent mass matrix: ILump=0, lumped mass matrix: ILump=1)
     &         (OnDiag*ESN3 + OffDiag*(ESN1+ESN2))*MsFacR

!...     Free surface terms from LHS (time levels s-1 & s)
     &        -GDPAvgOAreaIE4*(  C00  *(FDX3*E0XGrad2A+FDY3*E0YGrad2A)
     &                        +A00pB00*(FDX3*E1XGrad2A+FDY3*E1YGrad2A))

!...     Terms from momentum eqs.
     &        +(JXAvg*FDX3 + JYAvg*FDY3)/2.e0

!...     Spatially varying Tau0 terms
     &        +timewtgwce0*Tau0SpaVar0+timewtgwce1*Tau0SpaVar
     &                    +timewtgwce2*Tau0SpaVar2


!...     Put these partial products into further elemental storage for a vector computer
!...     These will be put into nodal storage outside of the elemental loop
#ifdef CVEC
         Temp_LV_A(IE,1)=Temp_LV_A1*NCEle
         Temp_LV_A(IE,2)=Temp_LV_A2*NCEle
         Temp_LV_A(IE,3)=Temp_LV_A3*NCEle
#endif

!...     Put these partial products directly into nodal storage for a scalar (non-vector) computer
#ifdef CSCA
         GWCE_LV(NM1)=GWCE_LV(NM1)+Temp_LV_A1*NCEle
         GWCE_LV(NM2)=GWCE_LV(NM2)+Temp_LV_A2*NCEle
         GWCE_LV(NM3)=GWCE_LV(NM3)+Temp_LV_A3*NCEle
#endif

!        IF(IE.EQ.1) THEN
!           WRITE(101,*) ' '
!           WRITE(101,*) '  ************* GWCE Load Vector ************'
!           WRITE(101,*) '  Time Step = ',IT
!           ENDIF
!        WRITE(101,*) IE, ESN1, ESN2, ESN3
!        WRITE(101,*) IE, TEMP_LV_A1, TEMP_LV_A2, TEMP_LV_A3

!kmd  Change the number of the loop
 1038 CONTINUE      !End of elemental loop


!...  Put load vector elemental values into nodal storage for a vector computer
#ifdef CVEC
      DO IE=1,NE
         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         GWCE_LV(NM1)=GWCE_LV(NM1)+Temp_LV_A(IE,1)
         GWCE_LV(NM2)=GWCE_LV(NM2)+Temp_LV_A(IE,2)
         GWCE_LV(NM3)=GWCE_LV(NM3)+Temp_LV_A(IE,3)
      END DO
#endif


!...  Save the elevation at the past time step into Eta1 and zero out Eta2
!kmd
!kmd  Already did this and haven't advanced in time yet, so just
!kmd  comment out the save elevation line. But we do need to zero
!kmd  out eta2 because it's involved in the summation below.
!kmd
      DO I=1,NP
!         Eta1(I)=Eta2(I)
         Eta2(I)=0.0e0
      END DO

!...  At elevation boundary condition nodes, determine the elevation at
!...  the s+1 time step
!...
!...  For periodic elevation boundary conditions

      DO J=1,NBFR
         IF(PER(J).EQ.0.) THEN
            NCYC=0
         ELSE
#ifdef IBM
            NCYC=INT(timeh/PER(J),KIND(0.0e0))
#else
            NCYC=INT(timeh/PER(J))
#endif
         ENDIF
         ARGJ=AMIG(J)*(timeh-NCYC*PER(J))+FACE(J)
         RFF=FF(J)*RampElev
         DO I=1,NETA
            ARG=ARGJ-EFA(J,I)
            NBDI=NBD(I)
            Eta2(NBDI)=Eta2(NBDI)+EMO(J,I)*RFF*COS(ARG)
         END DO
      END DO

!...  FOR APERIODIC ELEVATION BOUNDARY CONDITION

       IF (NPEBC) THEN
         IF(TimeLoc.GT.ETIME2) THEN
            ETIME1=ETIME2
            ETIME2=ETIME1+ETIMINC
            DO J=1,NETA
               ESBIN1(J)=ESBIN2(J)
               READ(19,*) ESBIN2(J)
            END DO
         ENDIF
         ETRATIO=(TimeLoc-ETIME1)/ETIMINC
         DO I=1,NETA
            NBDI=NBD(I)
            Eta2(NBDI)=RampElev
     &           *(ESBIN1(I)+ETRATIO*(ESBIN2(I)-ESBIN1(I)))
         END DO
      ENDIF


!...  IMPOSE NORMAL FLOW, RADIATION OR GRADIENT BOUNDARY CONDITIONS
!...  ALONG FLOW BOUNDARY TO LOAD VECTOR GWCE_LV(I)

!...  Note 2, Boundary conditions using specified fluxes (LBCODEI < 29)
!...  assume that QN is positive into the domain.  QFORCEJ has a -1
!...  built in and the terms are not explicitly negated. Boundary
!...  conditions using computed fluxes (LBCODEI 30, 40) compute a normal
!...  flux that  is positive out of the domain.  Therefore, to match
!...  the formulation these terms must be explicitly multiplied by -1.

!...Note 3, Eta1 is the latest computed elevation (it was updated above).

      IF((NFLUXF.EQ.1).OR.(NFLUXB.EQ.1).OR.(NFLUXIB.EQ.1)
     &     .OR.(NFLUXGBC.EQ.1).OR.(NFLUXRBC.EQ.1)) THEN
         NBDJ=NBV(1)
         IF(LBCODEI(1).LE.29) QFORCEJ=(QN2(1)-QN0(1))/DT2 +
     &        Tau0VAR(NBDJ)*QN1(1)

         IF(LBCODEI(1).EQ.30) THEN
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(1)
            ENDIF

         IF(LBCODEI(1).EQ.32) THEN
            H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
            CELERITY=SQRT(G*H1)
            QFORCEJ=(QN1(1)-QN0(1))/DT
     &           -CELERITY*(ETAS(NBDJ)-(EN1(1)-EN0(1)))/DT
     &           +TAU0VAR(NBDJ)*(QN1(1)-CELERITY*(ETA1(NBDJ)-EN1(1)))
            ENDIF

         IF((LBCODEI(1).EQ.40).OR.(LBCODEI(1).EQ.41)) QFORCEJ=
     &        -(QN1(1)-QN0(1))/DT - TAU0VAR(NBDJ)*(QN1(1)+QN0(1))/2.e0

         DO J=2,NVEL
            NBDI=NBDJ
            NBDJ=NBV(J)
            QFORCEI=QFORCEJ

            IF(LBCODEI(J).LE.29) QFORCEJ=(QN2(J)-QN0(J))/DT2+
     &           Tau0VAR(NBDJ)*QN1(J)

            IF(LBCODEI(J).EQ.30) THEN
               H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
               CELERITY=SQRT(G*H1)
               QFORCEJ=-CELERITY*ETAS(NBDJ)/DT - Tau0VAR(NBDJ)*QN1(J)
               ENDIF

            IF(LBCODEI(J).EQ.32) THEN
               H1=DP(NBDJ)+IFNLFA*ETA1(NBDJ)
               CELERITY=SQRT(G*H1)
               QFORCEJ=(QN1(J)-QN0(J))/DT
     &              -CELERITY*(ETAS(NBDJ)-(EN1(J)-EN0(J)))/DT
     &              +TAU0VAR(NBDJ)*(QN1(J)-CELERITY*(ETA1(NBDJ)-EN1(J)))
               ENDIF

            IF((LBCODEI(J).EQ.40).OR.(LBCODEI(J).EQ.41)) QFORCEJ=
     &        -(QN1(J)-QN0(J))/DT - TAU0VAR(NBDJ)*(QN1(J)+QN0(J))/2.e0

            NCI=NodeCode(NBDI)
            NCJ=NodeCode(NBDJ)
            BndLenO6NC=NCI*NCJ*BndLen2O3(J-1)/4.e0
            GWCE_LV(NBDI)=GWCE_LV(NBDI)
     &                              + BndLenO6NC*(2.e0*QForceI+QForceJ)
            GWCE_LV(NBDJ)=GWCE_LV(NBDJ)
     &                              + BndLenO6NC*(2.e0*QForceJ+QForceI)
          ENDDO
        ENDIF

!...
!...  IMPOSE ELEVATION BOUNDARY CONDITIONS TO LOAD VECTOR GWCE_LV(I) NOTE; EP
!...  IS THE RMS OF ALL THE DIAGONAL MEMBERS IN THE GWCE.  IT IS USED TO
!...  SCALE THE DIAGONAL ELEMENT FOR THE ELEVATION SPECIFIED BOUNDARY
!...  NODES AND THEREFORE MUST ALSO BE USED TO SCALE THE RHS OF THE
!...  EQUATIONS
!...
      DO I=1,NETA
         NBDI=NBD(I)
         ETAS(NBDI)=ETA2(NBDI)-ETA1(NBDI)
         GWCE_LV(NBDI)=ETAS(NBDI)*NODECODE(NBDI)*EP
         DO J=2,NNEIGH(NBDI)
            GWCE_LV(NEITAB(NBDI,J))=GWCE_LV(NEITAB(NBDI,J))
     &           -ETAS(NBDI)*OBCCOEF(I,J-1)
         END DO
      END DO

!...
!...  SOLVE GWCE FOR ELEVATION AT NEW TIME LEVEL
!...

!...  UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE

#ifdef CMPI
!...UPDATE LOAD VECTOR INITIAL GUESS and DIAGONAL FOR GWCE SOLVE
      CALL UPDATER(GWCE_LV,COEF(:,1),DUMY1,2)
#endif

!...  JCG ITERATIVE MATRIX SOLVER
      IPARM(1)=ITMAX
      CALL JCG(NP,MNP,MNEI,NEITAB,COEF,GWCE_LV,ETAS,
     &     IWKSP,NW,WKSP,IPARM,RPARM,IER)

      NUMITR=IPARM(1)
      DO I=1,NP
         ETA2(I)=NODECODE(I)*ETAS(I)+ETA1(I) !COMPUTE NEW ELEVATIONS
      END DO

!     UPDATE ELEVATIONS

#ifdef CMPI
      CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif

!kmd
!kmd  Eta2 values are now corrected elevations at time level s+1.
!kmd
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
!**********************************************************************
      END SUBROUTINE GWCE_NEW_PC
!**********************************************************************

!*************************************************************************
!   Subroutine to check if elemental slope limiting is needed
!       Chris Massey, USACE-ERDC-CHL Dec. 8, 2014 -- Removed this
!         section from timestep subroutine and put it in its own
!         subroutine.
!     jgf53.dev: Moved to GWCE module. 
!*************************************************************************
      subroutine check_slopes(it,TimeLoc)
      USE SIZES, ONLY : SZ, mnproc, myproc
      use global, only : eta2,nodecode,NOFF,ESLCOUNT,screenUnit,
     &      setMessageSource, unsetMessageSource, scratchMessage,
     &      allMessage,logMessage, DEBUG, ECHO, INFO, WARNING, ERROR,
     &      nodes_lg, IFSFM
      use mesh, only : nm, ne, areas, FDXE, FDYE, 
     &      SFacEle, SFMYEle, SFMXEle
      use NodalAttributes, ONLY : LoadEleSlopeLim,
     &     elemental_slope_limiter_active,
     &     elemental_slope_limiter_grad_max,
     &     elemental_slope_limiter_max_exceeded
      implicit none
      INTEGER, intent(in) :: IT
      Real(8), intent(in) :: TimeLoc
      INTEGER IE,I
      INTEGER NM1, NM2, NM3, NM123
      INTEGER NC1, NC2, NC3, NCEle, NCI
      REAL(kind=8) DEta2DX, DEta2DY, DEta2Mag
      REAL(kind=8) SFacAvg, SFmxAvg, SFmyAvg, sfdxfac, sfdyfac 
      REAL(kind=8) FDX1,FDX2,FDX3,FDY1,FDY2,FDY3
      REAL(8) :: AreaIE2
      integer nodeNumber ! fulldomain node number where wse slope is exceeded
      

      call setMessageSource("check_slopes")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ! jgf51.47: Per Rick's specification, start off by deactivating
      ! any slope limitation triggered in the previous time step;
      ! we have a fresh start each time step in determining whether
      ! to activate slope limiting at each node. 
      ! jgf51.50: Added this back to the subroutine version of the
      ! slope limiter. 
      elemental_slope_limiter_active(:) = .false.

!obell  CHECK IF THE LOCAL GRADIENT FOR ELEMENTS IS EXCEEDED AND TRIGGER THE
!...  ELEMENTAL SLOPE LIMITER ACCORDINGLY (original routine by Crystal Fulcher)

      DO IE=1,NE
         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NC1=NODECODE(NM1)
         NC2=NODECODE(NM2)
         NC3=NODECODE(NM3)
         NCEle=NC1*NC2*NC3*NOFF(IE)
         IF(NCEle.EQ.0)THEN
            CYCLE  ! this element is dry, go to the next one
         ENDIF
         AreaIE2 = Areas(IE)
!.... Beg WP:
         SFacAvg = SFacEle(IE)
!..... BEG DW/WJP
         SFmxAvg = SFMXEle(IE) ; 
         SFmyAvg = SFMYEle(IE) ; 
         sfdxfac = (1 - IFSFM)*SFacAvg + IFSFM*SFmxAvg ;
         sfdyfac = (1 - IFSFM)*1.0_SZ + IFSFM*SFmyAvg ;
!..... END DW/WJP
         FDX1 = FDXE(1,IE)*sfdxfac ; !c FDX1=(Y(NM2)-Y(NM3))*MX !b1*mx
         FDX2 = FDXE(2,IE)*sfdxfac ; !c FDX2=(Y(NM3)-Y(NM1))*MX !b2*mx
         FDX3 = FDXE(3,IE)*sfdxfac ; !c FDX3=(Y(NM1)-Y(NM2))*MX !b3*mx
         FDY1 = FDYE(1,IE)*sfdyfac ; !c FDY1=(X(NM3)-X(NM2))*MY !a1*my
         FDY2 = FDYE(2,IE)*sfdyfac ; !c FDY2=(X(NM1)-X(NM3))*MY !a2*my
         FDY3 = FDYE(3,IE)*sfdyfac ; !c FDY3=(X(NM2)-X(NM1))*MY !a3*my         
!.... end WP
         dEta2Dx  = (Eta2(NM1)*FDX1+Eta2(NM2)*FDX2+Eta2(NM3)*FDX3)
     &                 /AreaIE2
         dEta2Dy  = (Eta2(NM1)*FDY1+Eta2(NM2)*FDY2+Eta2(NM3)*FDY3)
     &                 /AreaIE2
         dEta2Mag = sqrt(dEta2Dx*dEta2Dx + dEta2Dy*dEta2Dy)
         !
         ! jgf51.51: Now that the slope limiter gets reset at 
         ! every time step, I had to rewrite the logging so that
         ! a log message is only written the first time the slope
         ! limiter is activated at a node during a particular run.
         DO I=1,3
            ! If the limiter is on already, go to the next node. 
            IF (elemental_slope_limiter_active(NM(IE,I))) CYCLE
            ! Compare the elemental slope to the maximum elemental gradient.
            grad_check: IF (dEta2Mag.GE.
     &           ABS(elemental_slope_limiter_grad_max(NM(IE,I)))) THEN
               ! jgf51.51: Log the fulldomain node number.
               nodeNumber = nm(ie,i)
               if (mnproc.gt.1) then
                  nodeNumber = nodes_lg(nm(ie,i))
               endif
               ! zc - If gradmax is positive or zero, activate slope
               ! limiting. 
               grad_case: if (elemental_slope_limiter_grad_max(nm(ie,i))
     &               .ge.0.0e0) then
                  ! If it is the first time that the slope limiter
                  ! has been activated at this node, write a log 
                  ! message. 
                  eslcount(nm(ie,i)) = eslcount(nm(ie,i))+1 ! for output file
                  if (eslcount(nm(ie,i)).eq.1) then
                     write(scratchMessage,1983) nodeNumber,dEta2Mag, 
     &                     elemental_slope_limiter_grad_max(NM(IE,I)),
     &                     it, timeLoc
                     call allMessage(INFO,scratchMessage)
                  endif
                  elemental_slope_limiter_active(NM(IE,I)) = .true.
               else
                  ! Just print log message the first time the 
                  ! gradient is exceeded. 
                  IF (elemental_slope_limiter_max_exceeded(NM(IE,I))
     &                  .eqv..false.) THEN
                     write(scratchMessage,1984) nodeNumber,dEta2Mag, 
     &                     elemental_slope_limiter_grad_max(NM(IE,I)),
     &                     it, timeLoc
                     call allMessage(INFO,scratchMessage)
                     elemental_slope_limiter_max_exceeded(NM(IE,I)) =
     &                      .true.
                  endif
               endif grad_case
            endif grad_check
         enddo ! loop around nodes of an element
      enddo ! loop over the elements


1983  format('Elemental slope limiter turned on at fulldomain node ',i0,
     &  ' where the elemental slope is ',1pE12.4E3,
     &  ' and the maximum elemental slope is ',1pE12.4E3,
     &  ' on time step ',i0,' and time = ',e15.8,'.')

1984  format('Maximum elemental slope exceeded at fulldomain node ',i0,
     &  ' where the elemental slope is ',1pE12.4E3,
     &  ' and the maximum elemental slope is ',1pE12.4E3,
     &  ' on time step ',i0,' and time = ',e15.8,'.')
     
!...  END CHECKING ELEMENT GRADIENTS

#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !-----------------------------------------------------------------
      END SUBROUTINE CHECK_SLOPES
      !-----------------------------------------------------------------   
!
!***********************************************************************

!***********************************************************************
!  Apply Elemental Slope Limiter 
!    Chris Massey, USACE-ERDC-CHL, Dec. 8, 2014 
!       Made into a subroutine
!     jgf53.dev: Moved to GWCE module. 
!***********************************************************************
!
      SUBROUTINE APPLY_SLOPE_LIMITS(ETA2Lim,LocNP)
      USE SIZES, ONLY : SZ
      use global, only : nodecode,NOFF,IFNLFA,
     &      setMessageSource, unsetMessageSource, allMessage,
     &      logMessage, DEBUG, ECHO, INFO, WARNING, ERROR
      use mesh, only : ne, nm, areas, totalArea
      use NodalAttributes, ONLY : LoadEleSlopeLim,
     &     elemental_slope_limiter_active
      implicit none
      integer :: IE,NM1,NM2,NM3,NC1,NC2,NC3,NCEle
      integer, intent(in) :: LocNP
      REAL(kind=8) :: EtaN1,EtaN2,EtaN3,EtaN123
      real(kind=8), intent(inout) :: Eta2lim(LocNP)
      REAL(8) :: AreaEle
      REAL(kind=8), ALLOCATABLE, SAVE :: elevSum(:) ! used if elemental slope limiter is active
      LOGICAL, SAVE :: firstCall = .true.


      call setMessageSource("check_slopes")
#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      IF (LoadEleSlopeLim.eqv..true.) THEN
         IF (firstCall.eqv..true.) THEN
            allocate(elevSum(LocNP))
            firstCall = .false.
            elevSum(:) = 0.e0
         ENDIF
      ENDIF


!       ELEMENTAL SLOPE LIMITER
!
!        CHECK TO SEE IF A HIGH GRADIENT HAS BEEN DETECTED IN THE PREVIOUS
!        TIME STEP AND APPLY SMOOTHING IF DESIRED. -ZC

      IF (LoadEleSlopeLim.eqv..true.) THEN
!
         elevSum(:) = 0.e0
         DO IE=1,NE
            NM1=NM(IE,1)
            NM2=NM(IE,2)
            NM3=NM(IE,3)
            NC1=NODECODE(NM1)
            NC2=NODECODE(NM2)
            NC3=NODECODE(NM3)
            NCEle=NC1*NC2*NC3*NOFF(IE)
            EtaN1=IFNLFA*Eta2Lim(NM1)
            EtaN2=IFNLFA*Eta2Lim(NM2)
            EtaN3=IFNLFA*Eta2Lim(NM3)
            AreaEle=NCEle*Areas(IE)/2.e0
            EtaN123=(EtaN1+EtaN2+EtaN3)/3.e0
            elevSum(NM1)=elevSum(NM1)+AreaEle*EtaN123
            elevSum(NM2)=elevSum(NM2)+AreaEle*EtaN123
            elevSum(NM3)=elevSum(NM3)+AreaEle*EtaN123
         ENDDO

!obell   CHECK TO SEE IF A HIGH GRADIENT HAS BEEN DETECTED IN THE PREVIOUS
!        TIME STEP AND APPLY SMOOTHING IF DESIRED. THIS ROUTINE PARALLELS THE
!        ABOVE ROUTINE.
         WHERE ((elemental_slope_limiter_active).and.
     &         (TotalArea.ne.0.e0))
               Eta2Lim = elevSum / TotalArea
         END WHERE
      ENDIF

!... Will apply the updating outside the subroutine
!#ifdef CMPI
!      CALL UPDATER(ETA2Lim,DUMY1,DUMY2,1)
!#endif


#if defined(TIMESTEP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()


      RETURN

      END SUBROUTINE APPLY_SLOPE_LIMITS

!***********************************************************************


!-----------------------------------------------------------------------
!      S U B R O U T I N E   G  W  C  E   T E R M I N A T E 
!-----------------------------------------------------------------------
!     jgf: Placed this subroutine here to eliminate dependency of this
!     module on adcirc_mod module. Bombs out of mpi gracefully when
!     there is an issue. 
!-----------------------------------------------------------------------
      subroutine gwceTerminate(NO_MPI_FINALIZE)
#ifdef CMPI
      use messenger
#endif
      implicit none
      logical, optional :: no_mpi_finalize
!
      call setMessageSource("gwceTerminate")
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      call allMessage(INFO,"ADCIRC Terminating.")

#ifdef CMPI
      subdomainFatalError = .true.
      if (present(no_mpi_finalize)) then
        call msg_fini(no_mpi_finalize)
      else
        call msg_fini()
      endif
#endif
      CALL EXIT(1)
!
#if defined(ADCIRC_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
!-----------------------------------------------------------------------
      end subroutine gwceTerminate
!-----------------------------------------------------------------------



!----------------------------------------------------------------------
!----------------------------------------------------------------------
      end module gwce
!----------------------------------------------------------------------
!----------------------------------------------------------------------
