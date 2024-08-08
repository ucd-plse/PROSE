!******************************************************************************
! PADCIRC VERSION 45.12 03/17/2006                                            *
!  last changes in this file VERSION 45.08                                    *
!... TCM V50.66.01 --ADDING TIME DEPENDENT BATHYMETRY                         *
!                 ADDED NDDT,BTIMINC,DP1,DP2,BTIME1,BTIME2                    *
!                                                                             *
! This module declares all global variables that are not exclusive to the 3D  *
! routines.  The 3D exclusive variables are declared in global_3DVS           *
!                                                                             *
!******************************************************************************
!
      MODULE GLOBAL
      USE SIZES, ONLY : SZ 
      USE KDTREE2_MODULE   !v49.48.01 tcm fast search algorithm
#ifdef DATETIME
      use datetime_module, only: datetime
#endif
      IMPLICIT NONE
      SAVE
!...
!...SET GLOBAL PHYSICAL CONSTANTS
!...
      REAL(kind=8) :: g = 9.80665 !default gravitational acceleration
!.....nominal density of water RHOWAT0
      REAL(kind=8), PARAMETER ::  RhoWat0=1000.E0
!.....Sigma T value of reference density
      REAL(kind=8), PARAMETER ::  SigT0=RHOWAT0-1000.E0
!.....Background Atmospheric Pressure in mb
      REAL(kind=8), PARAMETER :: PRBCKGRND = 1013.0e0   !tcm 20100617 v49.16 Added
      REAL(kind=8), PARAMETER :: TBCKGRND = 288.15e0    !wjp degrees kelvin 
      REAL(kind=8), PARAMETER :: RAir = 287.058e0        !wjp specific gas
!onstant 
      REAL(kind=8) :: rhoAir = 1.293e0 !PRBCKGRND*100e0/(RAir*TBCKGRND)  !wjp Ideal Gas Law kg/m^3
      REAL(8) :: windReduction = 0.9e0
      REAL(8) :: one2ten = 0.8928e0
      REAL(8) :: ten2one = 1.e0/0.8928e0
!asey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
      REAL(kind=8) :: WaveWindMultiplier = 1.E0
!.....PI etc
      REAL(8), PARAMETER ::  PI=3.141592653589793E0
      REAL(8), PARAMETER ::  DEG2RAD = PI/180.E0  ! degrees to radians
      REAL(8), PARAMETER ::  RAD2DEG = 180.E0/PI  ! radians to degrees
      REAL(8), PARAMETER ::  e=2.718281828459045e0
!.....Days to seconds conversion
      REAL(8), PARAMETER :: day2sec=24.e0*3600.e0
      REAL(8), PARAMETER :: sec2day=1.e0/day2sec
      ! jgf50.38.03: Additional useful constants
      REAL(8), PARAMETER :: Rearth = 6378206.4e0 ! radius of earth (m)
      REAL(8), PARAMETER :: nm2m = 1852.e0 ! nautical miles to meters
      REAL(8), PARAMETER :: m2nm = 1.e0/nm2m ! meters to nautical miles
      REAL(8), PARAMETER :: kt2ms = nm2m / 3600.0 ! knots to m/s
      REAL(8), PARAMETER :: ms2kt = 1.e0 / kt2ms ! m/s to knots
      REAL(8), PARAMETER :: omega = 7.29212d-5
      REAL(8), PARAMETER :: mb2pa = 100.e0
!..   Ali's full dispersion coefficient for M2
      REAL(kind=8) :: Bd = 0.23394_SZ  ! exponent of depth (default)
      REAL(kind=8) :: Ad = 0.0050189_SZ ! coefficient of depth (default)
      REAL(kind=8) :: Cs = 1500.0_SZ ! Speed of sound in water
      REAL(kind=8) :: Cs2 ! Speed of sound in water^2
!     jgf46.21 Added support for IBTYPE=52.
      REAL(kind=8),ALLOCATABLE, TARGET ::   ElevDisc(:)
      REAL(kind=8) FluxSettlingTime
      INTEGER  FluxSettlingIT
!... v49.48 tcm -- added for coupling with STWAVE (NRS=4)
      REAL(8)  :: STARTWAVE   ! Start of Wave Service
      REAL(8)  :: ENDWAVE     ! End of Wave Service
      LOGICAL  :: CPL2STWAVE = .FALSE.  !Coupled to STWAVE (NRS=4)
      LOGICAL  :: Flag_ElevError = .FALSE.  
!      real(kind=8) rampriver,drampriver,rampriver1,rampriver2
!     jgf46.08 Fine grained ramp functions (jgf46.21 split flux b.c.s)
      REAL(kind=8) RampExtFlux,DRampExtFlux ! Ramp for external flux b.c.s
      REAL(kind=8) RampIntFlux,DRampIntFlux ! Ramp for internal flux b.c.s
      REAL(kind=8) RampElev,DRampElev    ! Ramp for elevation boundary conditions.
      REAL(kind=8) RampTip,DRampTip      ! Ramp for tidal potential
      REAL(kind=8) RampMete,DRampMete    ! Ramp for wind and atmospheric pressure
      REAL(kind=8) RampWRad,DRampWRad    ! Ramp for wave radiation stress
!orbitt 1203022: Added Zachs Fix for Assigning a Start Time to Mete Ramping
      REAL(kind=8) DUnRampMete           ! Time to not have ramped Mete
!
!     jgf46.10 user-controlled warning, dump and error levels on elevations
      REAL(kind=8) WarnElev       ! elevation at which a warning is issued
      INTEGER iWarnElevDump   ! 1 to dump warning elevs to fort.69
      LOGICAL WarnElevDump    ! .True. if iWarnElevDump.eq.1
      INTEGER WarnElevDumpLimit ! max dumps; ADCIRC terminates if exceeded
      INTEGER WarnElevDumpCounter ! counter for number of elev dumps to fort.69
      REAL(kind=8) ErrorElev      ! ADCIRC terminates if this elev is exceeded
      INTEGER :: screenUnit = 6 ! jgf46.19 I/O unit where screen output is sent
      INTEGER :: Terminate_LocalProc = 0 ! zc - allow a single processor to initiate code shutdown
!
!     Variables related to hotstarting.
      type HOTSTART_t
        character(20)     :: filename
        integer           :: lun
      end type HOTSTART_t

      type(HOTSTART_t):: hss
      !       
      ! iLump=0 for consistent mass matrix
      ! iLump=1 for lumped mass matrix 
      integer :: ilump 

      real(kind=8),allocatable,target ::   mom_lv_x(:)
      real(kind=8),allocatable ::   mom_lv_y(:)

      !-----------------------------------------------------------------
      !        F U L L   D O M A I N   V A R I A B L E S   
      !                  O N L Y   U S E D   
      !        I N   P A R A L L E L   E X E C U T I O N
      !-----------------------------------------------------------------
      ! jgf51.21.25: The variables in this section are only used in
      ! parallel execution but they are declared here (rather than in
      ! the messenger module) because they are used in many different 
      ! modules, including harm.F, write_output.F, globalio.F et al.
      !-----------------------------------------------------------------
      integer :: comm        ! MPI communicator.
      integer :: np_g, ne_g  ! global (full domain) number of nodes, elements

!     jgf48.03 Arrays for mapping subdomain nodes and elements to full
!     domain ... used in globalio module but allocated and populated
!     in messenger module
      integer,allocatable,target :: nodes_lg(:)
      integer,allocatable,target :: imap_el_lg(:) 
!
!     global number of elevation, velocity, and meteorological stations
      integer :: nstae_g, nstav_g, nstam_g
      integer :: nstac_g     ! number of concentration stations
!
!     jgf49.44: These variables represent fulldomain arrays and are
!     used in hstart.F for reading fulldomain hotstart files as well as in
!     write_output.F/writeHotstart and
!     writer.F/writeHotstart_through_HSwriter
!     for writing fulldomain hotstart files.
      REAL(kind=8), ALLOCATABLE, TARGET :: ETA1_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET :: ETA2_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET :: EtaDisc_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET :: UU2_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET :: VV2_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET :: CH1_g(:)
      INTEGER, ALLOCATABLE, TARGET :: labels_g(:)
      INTEGER, ALLOCATABLE, TARGET :: NodeCode_g(:)
      INTEGER, ALLOCATABLE :: NNodeCode_g(:)
      INTEGER, ALLOCATABLE, TARGET :: NOFF_g(:)
      ! full domain min and max files in parallel
      !tcm v51.20.01 added variable for time of max (time is referenced off TimeLoc)
      REAL(kind=8),ALLOCATABLE, TARGET :: HOT_ETAMAX_G(:)
      REAL(kind=8),ALLOCATABLE, TARGET :: ETAMAX_Time(:) 
      REAL(kind=8),ALLOCATABLE, TARGET :: ETAMAX_Time_G(:) 
      REAL(kind=8),ALLOCATABLE, TARGET :: HOT_UMAX_G(:)
      REAL(kind=8),ALLOCATABLE, TARGET :: UMAX_Time(:) 
      REAL(kind=8),ALLOCATABLE, TARGET :: UMAX_Time_G(:) 
      REAL(kind=8),ALLOCATABLE, TARGET :: HOT_PrMin_G(:)
      REAL(kind=8),ALLOCATABLE, TARGET :: PrMin_Time(:) 
      REAL(kind=8),ALLOCATABLE, TARGET :: PrMin_Time_G(:) 
      REAL(kind=8),ALLOCATABLE, TARGET :: HOT_WVNOutMax_G(:)
      REAL(kind=8),ALLOCATABLE, TARGET :: WVNOutMax_Time(:) 
      REAL(kind=8),ALLOCATABLE, TARGET :: WVNOutMax_Time_G(:) 
      REAL(kind=8),ALLOCATABLE, TARGET :: HOT_RSNMax_g(:)
      REAL(kind=8),ALLOCATABLE, TARGET :: RSNMax_Time(:) 
      REAL(kind=8),ALLOCATABLE, TARGET :: RSNMax_Time_G(:) 
     
#if defined CSWAN || defined ADCSWAN
      !...Global SWAN Arrays
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_HSOut_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_DIROut_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TM01Out_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TPSOut_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_WINDXOut_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_WINDYOut_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TM02Out_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TMM10Out_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_HSMaxOut_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_DIRMaxOut_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TM01MaxOut_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TPSMaxOut_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_WINDMaxOut_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TM02MaxOut_g(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TMM10MaxOut_g(:)
#endif
      
      ! mappings from subdomain to fulldomain for stations
      integer,allocatable,target :: imap_stae_lg(:)
      integer,allocatable,target :: imap_stav_lg(:)
      integer,allocatable,target :: imap_stam_lg(:)
      integer,allocatable,target :: imap_stac_lg(:)
!
!     jgf48.03 Variable, arrays, and datastructure that are used by both
!     globalio and writer modules
      integer, parameter :: BUFSIZE_MAX = 131072
      real(kind=8) :: buf(BUFSIZE_MAX)
      real(kind=8) :: resultBuf(BUFSIZE_MAX)
      integer :: integerBuffer(BUFSIZE_MAX)
      integer :: integerResultBuffer(BUFSIZE_MAX)
      integer                    :: float_type
!
!     jgf48.03 Variables shared between messenger and writer modules
      integer ::  realtype, dbletype
      !
      ! there is a separate mpi communicator for each dedicated
      ! writer processor; their communicator ids are stored in the 
      ! following array
      integer, allocatable ::  comm_writer(:)
      ! communicators for hotstart file writing !st3
      integer, allocatable ::  comm_writeh(:), comm_hsleep(:)
      !
      ! = 0 if this is a compute proc; > 0 if this is a writer node
      integer :: writer_id
      ! following signals are used in messenger.f and writer.f
      integer, parameter :: sig_val   = 1
      integer, parameter :: sig_ini   = 0
      integer, parameter :: sig_write = 100
      integer, parameter :: sig_mesh  = 14   ! signal to read mesh for xdmf output
      integer, parameter :: sig_term  = 999  ! signal to terminate adc_writer
      integer, parameter :: sig_dl    = 333  ! signal to download new BC3D data
      integer, parameter :: sig_bcini = 330  ! signal to initialize 3DBC
!  tcm v49.48.01
      integer, parameter :: sig_pause = 99 ! signal to pause adc_writer
      !----------------------------------------------------------------
!     End of variables shared between messenger and other modules
      !----------------------------------------------------------------
!     
      INTEGER             :: FileFmtVersion      ! File format version.

!..RJW merged 09/02/2008 Casey 071219: Added the following variable declaration for 3D mass balance.
!             RESELEM is the number of processors on which each element appears.
      INTEGER, ALLOCATABLE :: RESELEM(:)
!.
#ifdef CMPI
      INTEGER, ALLOCATABLE ::   IDUMY(:)
      REAL(kind=8),ALLOCATABLE ::   DUMY1(:),DUMY2(:)
#endif
      !jgf50.44: Added variable to record the time that the min or max
      ! has occurred
      REAL(kind=8),ALLOCATABLE, TARGET ::  ETA1(:), ETA2(:), UU2(:), VV2(:)
      REAL(kind=8),ALLOCATABLE, TARGET ::  ETAMAX(:) ! v46.50 sb 11/11/2006
      REAL(kind=8),ALLOCATABLE, TARGET ::  UMAX(:) ! v46.50 sb 11/11/2006
      REAL(kind=8),ALLOCATABLE, TARGET ::  ET00(:), UU00(:), VV00(:)
      REAL(kind=8),ALLOCATABLE, TARGET ::  RMP00(:), RMU00(:), RMV00(:)
      REAL(kind=8),ALLOCATABLE, TARGET ::  RMICE00(:), CICEOUT(:) !v49.64.01 tcm -added for ice concentration

      REAL(kind=8),ALLOCATABLE, TARGET ::  Pr2(:)
#ifdef CSWAN
      REAL(kind=8),ALLOCATABLE, TARGET ::  SWAN_WDragCo(:) !MVB
#endif
      REAL(kind=8),ALLOCATABLE, TARGET ::  windDragOut(:) !MVB 

      REAL(kind=8),ALLOCATABLE, TARGET ::  PrMin(:) ! v46.50 sb 11/11/2006
      REAL(kind=8),ALLOCATABLE, TARGET ::  WVNOutMax(:) ! v46.50 sb 11/11/2006
      
      real(kind=8), allocatable, target :: wvnxout(:) ! output wind velocity in x-dir at each node 
      real(kind=8), allocatable, target :: wvnyout(:) ! output wind velocity in y-dir at each node 



      
#if defined CSWAN || defined ADCSWAN     
      !...Local SWAN Arrays
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_HSOut(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_DIROut(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TM01Out(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TPSOut(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_WINDXOut(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_WINDYOut(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TM02Out(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TMM10Out(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_HSMaxOut(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_DIRMaxOut(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TM01MaxOut(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TPSMaxOut(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_WINDMaxOut(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TM02MaxOut(:)
      REAL(kind=8), ALLOCATABLE, TARGET      :: SWAN_TMM10MaxOut(:)
#endif      

!  tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
!#ifdef CSWAN
!asey 090302: Added the following arrays for output of radiation stress gradients.
      REAL(kind=8),ALLOCATABLE, TARGET :: RSNXOUT(:), RSNYOUT(:)
!#endif
      LOGICAL :: outputWindDrag = .false. ! true if wind drag data should be written to output file

#if defined CSWAN || defined ADCSWAN      
!obell 20120510: Added logicals for turning on/off swan output files   
      LOGICAL                      :: SWAN_OutputHS = .FALSE.
      LOGICAL                      :: SWAN_OutputDIR = .FALSE.
      LOGICAL                      :: SWAN_OutputTM01 = .FALSE.
      LOGICAL                      :: SWAN_OutputTPS = .FALSE.
      LOGICAL                      :: SWAN_OutputWIND = .FALSE.
      LOGICAL                      :: SWAN_OutputTM02 = .FALSE.
      LOGICAL                      :: SWAN_OutputTMM10 = .FALSE.
      LOGICAL                      :: SWAN_OutputAgg(7)
#endif

!obell - Added for Time Varying Weirs, Set to default values
      LOGICAL                      :: USE_TVW = .FALSE.     !..From namelist, used time varying weirs
      CHARACTER(80)                :: TVW_FILE = 'fort.tvw' !..From namelist, file for time varying weirs
      REAL(kind=8),ALLOCATABLE,TARGET  :: TVW(:)                !..Current elevation of weirs, for output


!obell -  Added the following for output of ESL at end of run
      REAL(kind=8),ALLOCATABLE, TARGET :: ESLCOUNT(:)

      REAL(kind=8),ALLOCATABLE, TARGET ::  CC00(:)
      REAL(kind=8),ALLOCATABLE, TARGET ::  EtaDisc(:) ! vjp 9/2006 added for hotstart file
      REAL(kind=8),ALLOCATABLE, TARGET ::   CH1(:)
      REAL(kind=8),ALLOCATABLE ::   ETAS(:), UU0(:),VV0(:)
      real(kind=8),allocatable ::   GWCE_LV(:)
      REAL(kind=8),ALLOCATABLE ::   UU1(:),VV1(:)
      REAL(kind=8),ALLOCATABLE ::   QX1(:)
      REAL(kind=8),ALLOCATABLE,TARGET ::   QX2(:)
      REAL(kind=8),ALLOCATABLE ::   QY1(:)
      REAL(kind=8),ALLOCATABLE,TARGET ::   QY2(:)
      REAL(kind=8),ALLOCATABLE ::   QX0(:),QY0(:)
      REAL(kind=8),ALLOCATABLE ::   TK0(:),TK2(:)
      REAL(kind=8),ALLOCATABLE ::  DP1(:),DP2(:) !TCM V50.66.01 FOR TIME VARYING BATHYMETRY

      REAL(kind=8),ALLOCATABLE ::   CORIF(:)
      REAL(kind=8),ALLOCATABLE ::   TPK(:),FFT(:)
      REAL(kind=8),ALLOCATABLE ::   FACET(:),ETRF(:)
      REAL(kind=8),ALLOCATABLE ::   FFF(:),FFACE(:)
      REAL(kind=8),ALLOCATABLE ::   ESBIN1(:),ESBIN2(:)
      REAL(kind=8),ALLOCATABLE ::   TEMP_LV_A(:,:),TEMP_LV_B(:,:)
      REAL(kind=8),ALLOCATABLE ::   QN0(:),QN1(:)
      REAL(kind=8),ALLOCATABLE, TARGET :: QN2(:)
      REAL(kind=8),ALLOCATABLE ::   EN0(:),EN1(:),EN2(:)
      REAL(kind=8),ALLOCATABLE ::   QNAM(:,:),QNPH(:,:)
      REAL(kind=8),ALLOCATABLE ::   ENAM(:,:),ENPH(:,:)
      REAL(kind=8),ALLOCATABLE ::   QNIN1(:),QNIN2(:)
      REAL(kind=8),ALLOCATABLE ::   ENIN1(:),ENIN2(:)
      !
      ! spatial interpolating factors for elevation stations
      real(kind=8),allocatable, target :: staie1(:),staie2(:),staie3(:)
      real(8), allocatable, target :: xev(:) ! x coords of current vel stations
      real(8), allocatable, target :: yev(:) ! y coords of current vel stations
      real(8), allocatable, target :: slev(:)! longitudes of current vel sta
      real(8), allocatable, target :: sfev(:)! latitudes of current vel sta

      ! spatial interpolating factors for velocity stations
      real(kind=8),allocatable, target ::   staiv1(:),staiv2(:),staiv3(:)
      real(8), allocatable, target ::   xec(:),yec(:),slec(:),sfec(:)
      ! spatial interpolating factors for concentration stations
      real(kind=8),allocatable ::   staic1(:),staic2(:),staic3(:)

      ! spatial interpolating factors for meteorological stations
      real(kind=8),allocatable, target ::   staim1(:),staim2(:),staim3(:)
      real(8), allocatable, target :: xem(:) ! x coords of met stations
      real(8), allocatable, target :: yem(:) ! y coords of met stations
      real(8), allocatable, target :: slem(:) ! longitudes of met stations
      real(8), allocatable, target :: sfem(:) ! latitudes of met stations

      REAL(kind=8),ALLOCATABLE ::   TRANS_LV_B(:),TRANS_LV_A(:)
      REAL(kind=8),ALLOCATABLE ::   SOURSIN(:)
      REAL(kind=8),ALLOCATABLE ::   WSX1(:),WSY1(:),PR1(:)
      REAL(kind=8),ALLOCATABLE ::   WSX2(:),WSY2(:)
      REAL(kind=8),ALLOCATABLE ::   PRN1(:)
      REAL(kind=8),ALLOCATABLE ::   PRN2(:)
      REAL(kind=8),ALLOCATABLE ::   RSNX1(:),RSNY1(:),RSNX2(:),RSNY2(:)
      REAL(kind=8),ALLOCATABLE ::   CICE1(:),CICE2(:)   !v49.64.01 tcm -- added
                                                    !for ice concentration
      REAL(kind=8),ALLOCATABLE, TARGET :: RSNMax(:) ! v46.50 sb 11/11/2006

      REAL(kind=8),ALLOCATABLE :: TK(:)     ! WJP bottom frictiom
      REAL(kind=8),ALLOCATABLE :: TKM(:,:)  ! symmetrical friction matrix
      REAL(8),ALLOCATABLE :: EMO(:,:),EFA(:,:)
      REAL(8),ALLOCATABLE,TARGET :: XEL(:)  ! x-coords of elevation stations
      REAL(8),ALLOCATABLE,TARGET :: YEL(:)  ! y-coords of elevation stations
      REAL(8),ALLOCATABLE,TARGET :: SLEL(:) ! latitudes of elevation stations
      REAL(8),ALLOCATABLE,TARGET :: SFEL(:) ! longitudes of elevation stations


      REAL(kind=8),ALLOCATABLE,TARGET ::   DUU1(:),DUV1(:),DVV1(:)
      REAL(kind=8),ALLOCATABLE,TARGET ::   BSX1(:),BSY1(:)
      REAL(kind=8),ALLOCATABLE ::   TIP1(:),TIP2(:)
      REAL(kind=8),ALLOCATABLE ::   L_N(:,:)    ! WJP The array for TIP multiplier
      REAL(kind=8),ALLOCATABLE ::   SALTAMP(:,:),SALTPHA(:,:)
      REAL(kind=8),ALLOCATABLE ::   OBCCOEF(:,:),COEF(:,:)
      REAL(kind=8),ALLOCATABLE ::   COEFD(:) !jgf48.4619 from Seizo (Lumped GWCE)
      INTEGER :: NBFR, NCOR
!    kmd48.33bc - added variables for the levels of no motion and top temperature boundary condition
      REAL(kind=8),ALLOCATABLE ::   LNM_BC1(:), LNM_BC2(:), LNM_BC(:) ! level of no motion
      REAL(kind=8),ALLOCATABLE ::   q_heat1(:), q_heat2(:), HFLUX(:)  ! top temperature
      REAL(kind=8),ALLOCATABLE, TARGET :: Sponge(:)  ! sponge information for boundary conditions

      INTEGER,ALLOCATABLE, TARGET ::    NOFF(:), NODECODE(:)
      INTEGER,ALLOCATABLE ::    NNODECODE(:)
      INTEGER,ALLOCATABLE ::    NIBNODECODE(:)
      !
      ! lists of elements containing each station 
      integer, allocatable, target :: nnc(:) ! concentration station elements
      integer, allocatable, target :: nne(:) ! elevation station elements
      integer, allocatable, target :: nnv(:) ! velocity station elements
      integer, allocatable, target :: nnm(:) ! meteorological station elements
      !
      INTEGER,ALLOCATABLE ::    IPV(:)
!
!.....for buoyancy forcing
      REAL(kind=8),ALLOCATABLE ::   DASigT(:), DASalt(:), DATemp(:)
      REAL(kind=8),ALLOCATABLE ::   VIDBCPDXOH(:), VIDBCPDYOH(:)
      REAL(kind=8),ALLOCATABLE ::   SigTS(:), NM(:), NB(:), MLD(:)
!.....for scalar transport
      REAL(kind=8),ALLOCATABLE ::   DAConc(:)
      
!     The vertically-integrated momentum dispersion terms     
      REAL(kind=8),ALLOCATABLE ::   VIDISPDXOH(:), VIDISPDYOH(:)

      REAL(kind=8) DT
!                                
      INTEGER NOLICA,NOLIFA,IHOT,IDEN
      INTEGER :: NSCREEN = 1
!...
!...DECLARE REAL(8) AND CHAR VARIABLES, EQUIVALENCES
!...

      REAL(8) STATIM,REFTIM
      REAL(8) QTIME1,QTIME2
      REAL(8), TARGET :: DTDP      

      REAL(8) FTIMINC,RSTIME1,RSTIME2,RSTIMINC

      REAL(8) DIST
      REAL(8),ALLOCATABLE :: AMIG(:),AMIGT(:),FAMIG(:)
      REAL(8),ALLOCATABLE :: PER(:),PERT(:),FPER(:)
      REAL(8),ALLOCATABLE :: FREQ(:),FF(:),FACE(:)


      REAL(8) DTDPHS  !kmd48.33 - added time step for hot start
      CHARACTER(80)             :: RUNDES, RUNID, AFRIC
      CHARACTER(4)              :: RDES4(20),RID4(20)
      CHARACTER(8)              :: RDES8(10),RID8(10),AID8(10)
      CHARACTER(10)             :: ALPHA ! periodic constituent name
      CHARACTER(10),ALLOCATABLE :: ELEVALPHA(:) ! periodic constituent name
      CHARACTER(5),ALLOCATABLE  :: TIPOTAG(:),BOUNTAG(:),FBOUNTAG(:)

!.ral v45.02 to fix compiler problem, see also change in read_input.F
!      EQUIVALENCE (RDES4(1),RDES8(1),RUNDES), (RID4(1),RID8(1),RUNID),
!     &            (AID4(1),AID8(1),AGRID)

!...
!...DECLARE AND INITIALIZE LOGICAL VARIABLES
!...  jgf45.06 Default choices are set to .TRUE. in read_input.F
!     jgf45.07 The definitions of many of these are in timestep.F
      LOGICAL :: C2DDI               = .FALSE.
      LOGICAL :: C3D                 = .FALSE.
      LOGICAL :: C3DDSS              = .FALSE.
      LOGICAL :: C3DVS               = .FALSE.
      LOGICAL :: C2D_BTrans          = .FALSE.
      LOGICAL :: C2D_PTrans          = .FALSE.
      LOGICAL :: CBaroclinic         = .FALSE.
      LOGICAL :: CSmag_Eh            = .FALSE.
      LOGICAL :: CGWCE_New           = .FALSE.
      LOGICAL :: CGWCE_Lump          = .FALSE.
      LOGICAL :: CGWCE_LS_KGQ        = .FALSE.
      LOGICAL :: CGWCE_LS_2PartQ     = .FALSE.
      LOGICAL :: CGWCE_LS_2PartV     = .FALSE.
      LOGICAL :: CGWCE_LS_2PartSQ    = .FALSE.
      LOGICAL :: CGWCE_LS_2PartSV    = .FALSE.
      LOGICAL :: CGWCE_Advec_NC      = .FALSE.
      LOGICAL :: CGWCE_Advec_C1      = .FALSE.
      LOGICAL :: CGWCE_Advec_C2      = .FALSE.
      LOGICAL :: CGWCE_HDP           = .FALSE.
      LOGICAL :: CME_Orig            = .FALSE.
      LOGICAL :: CME_New_NC          = .FALSE.
      LOGICAL :: CME_New_C1          = .FALSE.
      LOGICAL :: CME_New_C2          = .FALSE.
      LOGICAL :: CPRECOR             = .FALSE.
      LOGICAL :: CME_LS_IBPQ         = .FALSE.
      LOGICAL :: CME_LS_IBPV         = .FALSE.
      LOGICAL :: CME_LS_IBPSQ        = .FALSE.
      LOGICAL :: CME_LS_IBPSV        = .FALSE.
      LOGICAL :: CME_LS_2PartQ       = .FALSE.
      LOGICAL :: CME_LS_2PartV       = .FALSE.
      LOGICAL :: CME_AreaInt_Orig    = .FALSE.
      LOGICAL :: CME_AreaInt_Corr    = .FALSE.
      LOGICAL :: CTIP                = .FALSE.
      LOGICAL :: CAliDisp            = .FALSE.
!    kmd48.33bc  - logical statement for hot start using different time step
      ! This is .true. if the time step size of the current run is
      ! different from the time step size in the hotstart file; also
      ! true for IDEN=1, 3D baroclinic runs 
      LOGICAL :: CHOTHS              = .FALSE.

! kmd - Evan's updates for rivers above MSL
      LOGICAL :: River_above_MSL     = .FALSE.

!     jgf48.4627 When only meteorological output is requested, skip past
!     the GWCE and Momentum equations, so that ADCIRC runs faster and
!     only calculates the requested quantities.
      LOGICAL :: METONLY             = .FALSE.
!
!     jgf51.14: Moved the parameters that control output timing 
!     from read_input.F to here so that they can be used to recalculate
!     starting and ending time step for output when the time step
!     size has changed after a hotstart. 
      REAL(kind=8) TOUTSC
      REAL(kind=8) TOUTSE,TOUTFE
      REAL(kind=8) TOUTSV,TOUTFV
      REAL(kind=8) TOUTFGC, TOUTFGE, TOUTFGV, TOUTFGW, TOUTF_TVW
      REAL(kind=8) TOUTFM, TOUTSGC, TOUTSGE, TOUTSGV, TOUTSGW, TOUTSM, TOUTS_TVW

      REAL(kind=8),ALLOCATABLE ::   RPIPEWL1AVG(:),RPIPEWL2AVG(:)
      INTEGER  NODEDRYMIN,NODEWETMIN
      INTEGER  IM
      INTEGER :: NABOUT = 0
      INTEGER  NFFR
      INTEGER  NFOVER,NOLICAT,NOUTC
      INTEGER  NOUTE,NSPOOLE,NOUTV,NSPOOLV
      INTEGER  NRAMP,NRS,NSTAE,NSTAV,NT,NTCYFE
      INTEGER  NTCYFV,NTCYSE,NTCYSV,NTIF,NTIP,NTRSPE,NTRSPV
      INTEGER  NTCYS_TVW,NTCYF_TVW
      INTEGER  NWS
      INTEGER  ICSTP

!orbitt 120322: Create Elemental Versions of IFNLCAT and IFNLCT
      INTEGER  IFNLCAT, IFNLCT, IFNLFA, IFNLCTE, IFNLCATE
      INTEGER :: IFNL_HDP = 0 ! use H (1) or D (0) in LHS surface gradient term
      INTEGER  IFWIND
      integer :: IFSFM = 0 ! corrected spherical coordinates (1) or not (0)
      integer :: IFSPROTS = 0 ! rotated spherical coordinates (1) or not (0)
      ! 
      ! fulldomain (a.k.a. global) file position counters
      integer, target :: igcp ! concentration
      integer, target :: igep ! elevation
      integer, target :: igpp ! barometric pressure
      integer, target :: igvp ! water velocity
      integer, target :: igwp ! wind velocity
      integer, target :: itau0p ! position in tau0 output file for variable tau0
      integer, target :: igsp ! sponge layer output
      integer, target :: igdp ! debug elevation
      ! station file position counters
      integer, target :: ivstp ! file position in velocity station file
      integer, target :: iwstp ! file position in wind velocity station file
      integer, target :: ipstp ! position in barometric pressure station file
      integer, target :: iestp ! file position in station elevation file 
      integer, target :: ig_tvw ! file position in weir elevation file
      integer, target :: nodeCodeOutFilePos
      integer, target :: noffOutFilePos
      integer, target :: windDragFilePosition
      ! jgf52.08.03: Some files (e.g., initiallydry.63 and endrisinginun.63)
      ! are only written once, but are written using the same subroutines
      ! that write time varying output files; the output subroutine increments
      ! the file position counter; so we have to have a fake file position 
      ! counter even for output files that are only written once.
      integer, target :: fauxFilePosition 
      ! jgf: added support for water level offset
      integer, target :: dynamicWaterLevelCorrectionFilePos
      integer, target :: dynamicWaterLevelCorrectionStaFilePos

      INTEGER  IHOTSTP, IMHS
      INTEGER  ITIME, ITEMPSTP
      INTEGER  ITHS, ITITER    
      
      INTEGER  LUMPT ! =1/0; LUMP/DO NOT LUMP THE TRANSPORT EQN
      INTEGER  NDSETSC
      INTEGER  NCChange !jgf45.06 flag indicating wetting / drying has occurred
      INTEGER  NDSETSE, NDSETSV, NDSETSW, NDSETS_TVW
      INTEGER  NHSINC, NHSTAR
      LOGICAL :: NHOUTONCE = .false.
!  tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
!#ifdef CSWAN 
!asey 090302: Added a variable for output of radiation stress gradients.
      INTEGER, TARGET :: IGRadS
!#endif
      INTEGER  NOUTGC, NOUTGE, NOUTGV, NOUTGW, NOUTM, NOUT_TVW
      INTEGER, TARGET :: NSCOUC, NSCOUE, NSCOUGC, NSCOUGE, NSCOUGV ! output spool counters
      INTEGER, TARGET :: NSCOUGW, NSCOURS, NSCOUGP, NSCOUM, NSCOUV ! output spool counters
      INTEGER, TARGET :: NSCOUGI, NSCOUGT, NSCOUGS  ! output spool counters
      INTEGER, TARGET :: NSCOU_TVW ! output spool counters
      INTEGER, TARGET :: NSCOUP ! output spool counter for fort.71 (bar. press. sta.)
      INTEGER, TARGET :: NSCOUI ! output spool counter for fort.91 (ice sta.)     
      integer, target :: nodeCodeOutSpoolCounter
      integer, target :: noffOutSpoolCounter
      integer, target :: initiallyDrySpoolCounter
      integer, target :: windDragSpoolCounter
      integer, target :: dynamicWaterLevelCorrectionSpoolCounter
      integer, target :: dynamicWaterLevelCorrectionStaSpoolCounter
      
      INTEGER :: NSPOOLC, NSPOOLGC, NSPOOLGE ! output increments
      INTEGER :: NSPOOLGV, NSPOOLGW, NSPOOLM, NSPOOL_TVW ! output increments
      INTEGER  NSTAC, NSTAM, NTCYFC, NTCYFGC, NTCYFGE
      INTEGER  NTCYFGV, NTCYFGW
      INTEGER  NTCYFM, NTCYSC, NTCYSGC, NTCYSGE
      INTEGER  NTCYSGV, NTCYSGW, NTCYSM

      INTEGER  NTRSPC, NTRSPM, NUMITR, NW, NWSEGWI

      INTEGER  NP_GLOBAL
!     kmd48.33bc - added variables for 3D boundary conditions
      INTEGER  RES_BC_FLAG  ! boundary condition flag for salinity and temperature
      INTEGER  BCFLAG_LNM   ! boundary condition flag for levels of no motion
      INTEGER  BCFLAG_TEMP  ! boundary condition flag for heat flux (top temperature)
      REAL(kind=8) RBCTIMEINC, RBCTIME1, RBCTIME2  ! time information for level of no motion
      REAL(kind=8) SBCTIMEINC, SBCTIME1, SBCTIME2  ! time information for salinity
      REAL(kind=8) TBCTIMEINC, TBCTIME1, TBCTIME2  ! time information for temperature
      REAL(kind=8) TTBCTIMEINC, TTBCTIME1, TTBCTIME2  ! time information for heat flux (top temperature)
      REAL(kind=8) BCSTATIM, RBCRATIO
      REAL(kind=8) SBCSTATIM, SBCRATIO
      REAL(kind=8) TBCSTATIM, TBCRATIO
      REAL(kind=8) TTBCSTATIM, TTBCRATIO
      REAL(kind=8) SPONGEDIST  ! distance information for the sponge layer
      LOGICAL :: OUTPUTSPONGE  =  .FALSE.
      REAL(kind=8) RIVBCTIMINC, RIVBCTIME1, RIVBCTIME2  ! time information for river - baroclinic
      REAL(kind=8) BCRivRATIO, RIVBCSTATIM
      REAL(kind=8) BEDSTR
      REAL(kind=8) COND
      REAL(kind=8) RAMP, RBARWL
      REAL(kind=8) RBARWL1,RBARWL1F, RBARWL2, RBARWL2F
      REAL(kind=8) TIMEIT, TOUTFC
      REAL(kind=8) CORI
      REAL(kind=8) DAY,DRAMP
      REAL(kind=8) H0
      REAL(kind=8) RNDAY
      REAL(8)  RNP_GLOBAL
      REAL(8)  REFSEC   ! required to run in either 32-bit or 64-bit

      LOGICAL :: useNetCDF = .false. ! .true. if netcdf is used at all
      LOGICAL :: useNetCDFOutput = .false. ! only for output, not hotstart
      LOGICAL :: useXDMF = .false.
! CF Start addition for NETCDF I/O 6/21/06
!     -----------------
!     Attribute general
!     -----------------
      CHARACTER(80) :: title       ! description of the dataset
      CHARACTER(80) :: institution ! where the original data was produced
      CHARACTER(80) :: source      ! method of production of the original data
      CHARACTER(80) :: history     ! provides an audit trail for modifications to the original data
      CHARACTER(80) :: base_date
#ifdef DATETIME
      type(datetime) :: basedatetime
#endif
      CHARACTER(80) :: comments    ! miscellaneous information about the data
      CHARACTER(80) :: host
      CHARACTER(80) :: convention
      CHARACTER(80) :: contact     ! user's institution/ e-mail address
      CHARACTER(80) :: references  ! published/web-based references that describe the data
      CHARACTER(7),  ALLOCATABLE :: STATNUMB(:,:)
      CHARACTER(50), ALLOCATABLE, TARGET  :: STATNAME(:)
      CHARACTER(7),  ALLOCATABLE :: STATNUMBV(:,:)
      CHARACTER(50), ALLOCATABLE, TARGET  :: STATNAMEV(:)
      INTEGER, PARAMETER :: SNUMLEN = 7
      INTEGER, PARAMETER :: SNAMLEN = 50
      CHARACTER(7),  ALLOCATABLE  :: STATNUMBM(:,:)
      CHARACTER(50), ALLOCATABLE, TARGET  :: STATNAMEM(:)
      CHARACTER(50), ALLOCATABLE, TARGET  :: STATNAMEC(:)

      INTEGER NCIESTP, NCIPSTP,NCIVSTP,NCIWSTP
      INTEGER NCIGEP,NCIGPP, NCIGVP,NCIGWP
      INTEGER NCCOUE,NCCOUM,NCCOUV,NCCOUW
      INTEGER NCCOUGE,NCCOUGW,NCCOUGV

#ifdef ADCNETCDF
      LOGICAL, PARAMETER :: netcdf_avail = .true. 
#else
      LOGICAL, PARAMETER :: netcdf_avail = .false. ! .true. if NetCDF support was compiled in
#endif
      
! CF Finish addition for NETCDF I/O

! TCM v49.64.01 Addition for Ice Concentration Fields
      INTEGER :: NCICE
      INTEGER, TARGET :: IGIP ! position in fulldomain Ice output file
      INTEGER, TARGET :: IICESTP ! position in station Ice output file
      REAL(8) :: CICE_TIMINC,CICE_TIME1,CICE_TIME2

! TCM V50.66.01 ADDITIONS FOR TIME VARYING BATHYMETRY
      INTEGER :: NDDT  !FLAG TO USE TIME VARYING BATHYMETRY
      REAL(8) :: BTIME1,BTIME2  !TIMES FOR TIME VARYING BATHYMETRY Records
      REAL(8) :: BTIME_END  !time that 
      !variables related to global bathymetry output file
      INTEGER, TARGET :: IGBP ! position in global bathy output file
      INTEGER :: NDSETSB,NTCYFGB,NTCYSGB,NTRSPB
      INTEGER, TARGET :: NSCOUGB
      !variables related to station bathymetry output file
      INTEGER :: NTCYSB,NTCYFB
      INTEGER, TARGET :: NSCOUB
      INTEGER, TARGET :: IBSTP  ! position in station bathy output file      
      REAL(8) :: BTIMINC ! TIME INCREMENT (SECONDS) FOR BATHYMETRY CHANGES USED FOR TIME VARYING BATHYMETRY
      REAL(8) :: BCHGTIMINC ! time increment (seconds) over which bathymetry changes during a btimeinc interval
      REAL(kind=8),ALLOCATABLE, TARGET ::  DP00(:)

! TCM v52.30.01 ADDITIONS FOR USING SMAGORINSKY eddy MODEL OPTIONS
      REAL(kind=8) :: SMAG_UPPER_LIM = 100.e0  ! m^2/s (default upper limit)
      REAL(kind=8) :: SMAG_LOWER_LIM = 1.0d-8  ! m^2/s (default lower limit)
      LOGICAL :: SMAG_COMP_FLAG = .FALSE.   !BY DEFAULT DO NOT USE SMAG
! TCM v53.01.02 Addition for special nearest neighbor station interpolation when
!           station location occurs in a partially wet / partially dry element
      LOGICAL :: StatPartWetFix = .FALSE.   !.true. if one wants to handle partially wet/dry elements in a special way
      INTEGER :: How2FixStatPartWet = 0     ! other numbers could be used later for different ways to interpolate,
                                            ! instead of nearest neighbor
                                            ! 0 is use nearest neighbor if wet and H > 0.8H0
                                            ! 1 is use nearest neighbor if wet regardless if H > 0.8H0

      ! 
      ! all info needed for self describing dataset in XDMF
      type xdmfMetaData_t
         integer*8 :: xdmfFortranObj  ! object associated with the file
         integer :: depthID    ! index of the depth attribute
         integer :: geometryID ! index of the Geometry data in the output file
         real(8), allocatable :: tempCoord(:,:) ! holds Geometry data
         integer :: topologyID ! index of the Topology data in the output file
         integer, allocatable :: xdmf_nm(:,:)  ! 0-offset connectivity array
         logical :: createdIDs        ! .true. if infoIDs have been created 
         character(80) :: variable_name 
         integer :: variable_name_id
         character(80) :: long_name 
         integer :: long_name_id
         character(80) :: standard_name
         integer :: standard_name_id
         character(80) :: coordinates
         integer :: coordinates_id
         character(80) :: units
         integer :: units_id
         character(80) :: positive
         integer :: positive_id
         real(kind=8), allocatable :: data_array3(:,:)
      end type xdmfMetaData_t
      !
      ! info required for dedicated writers
      type writerData_t
         logical :: bufferInitialized  !.true. if memory has been allocated
         logical :: bufferLoaded       !.true. if buffered data must be written
         real(kind=8) :: myTime
         integer :: myTimeStep
      end type writerData_t
!
!     jgf48.03 the following data structure is used to describe
!     the size and shape of the output data.
!     TCM48.4618 Added ConsiderWetDry, alternate_value and field_name 
!     elements to the structure
      type OutputDataDescript_t
        integer           :: lun
        CHARACTER(1024)   :: file_name ! name of output file
        character(20)     :: field_name !Name of the field associated with the data structure
        integer           :: specifier  ! output format from fort.15 (NOUTE,etc)
        REAL(kind=8)          :: initial_value
        logical           :: ConsiderWetDry ! Flag for considering nodes wet/dry status
        real(kind=8)          :: alternate_value  ! Value that can be used as an alternate value
        REAL(kind=8), pointer :: x_coord(:) ! x-coord or longitude (station or mesh)
        REAL(kind=8), pointer :: y_coord(:) ! y-coord or lat (station or mesh)
        ! num_items_per_record is number of "tuples", i.e.
        ! it is equal to 1 if scalar, 2 if 2D vector, 3 if 3D vector, etc
        integer           :: num_items_per_record
        ! array length for the full domain
        integer           :: num_fd_records
        ! array length for this domain (in serial, refers to full domain,
        ! in parallel, refers to the subdomain processed by this PE)
        integer           :: num_records_this
        integer,  pointer :: imap(:)   ! local-to-global mapping
        integer,  pointer :: iarray(:)
        REAL(kind=8), pointer :: array(:)  ! scalar data or 1st component of vector
        REAL(kind=8), pointer :: array2(:) ! 2nd component of vector
        REAL(kind=8), pointer :: array3(:) ! 3rd component of vector
        REAL(kind=8), pointer :: array2D(:,:) ! subdomain matrix
        logical :: isInteger ! .true. if this is an integer array
        INTEGER, pointer  :: iarray_g(:) ! full domain integer
        REAL(kind=8), pointer :: array_g(:)  ! full domain data scalar etc
        REAL(kind=8), pointer :: array2_g(:) ! full domain 2nd vec component
        REAL(kind=8), pointer :: array3_g(:) ! full domain 3rd vec component
        REAL(kind=8), pointer :: array2D_g(:,:) ! fulldomain matrix
        REAL(kind=8), pointer :: hotstart(:) ! from hotstart (used by min max)
        REAL(kind=8), pointer :: hotstart_g(:) ! from hotstart (used by min max)
        logical :: writeFlag ! .true. if the data should be written during this run
        logical :: isStation ! .true. if the data represent a recording station
        logical :: divideByDepth ! .true. if the solution must be divided by depth b/f output
        integer, allocatable :: writerFormats(:) ! formats that dedicated writers can write
        logical :: useWriter ! .true. if a dedicated writer should be used to write these data
        integer :: file_extension  ! end of the file name       
        character(len=32) :: file_basename  ! usually 'fort' but sometimes 'rads' etc
        real(kind=8), pointer :: interped_array(:) ! array that stations spatially interpolate
        real(kind=8), pointer :: interped_array2(:) ! array2 that stations spatially interpolate
        integer :: startTimeStep ! time step after which output should begin
        integer :: endTimeStep   ! last time step when output should occur    
        integer :: outputTimeStepIncrement ! num time steps between outputs          
        integer, pointer :: spoolCounter ! counts time steps between outputs
        integer, pointer :: filepos ! keeps track of position in output files
        integer, pointer :: elements(:) ! array of elements that contain the stations
        real(kind=8), pointer :: interp_fac1(:) ! array of interp. factors for stations
        real(kind=8), pointer :: interp_fac2(:) ! array of interp. factors for stations        
        real(kind=8), pointer :: interp_fac3(:) ! array of interp. factors for stations
        logical :: initialized ! .true. if initialization has occurred (netcdf, xdmf)
        logical :: minmax_timestamp ! .true. if a time stamp will be written to a max file
        logical :: readMaxMin ! .true. if max/min file should be read upon hotstart
        type(xdmfMetaData_t) :: xdmfMD ! used if this data will be written to XDMF
        type(writerData_t) :: wMD ! used if this data will be written by dedicated writer
      end type OutputDataDescript_t
      !
      ! inunThresh is used to set a minimum inundation depth for 
      ! accumulating inundationTime for detailed inundation output
      real(kind=8) :: inunThresh = 0.6e0
      !
      ! inundationOutput can be set to .true. in read_input()
      logical :: inundationOutput = .false. 
!
!     jgf49.44: Log levels, in order from largest amount of log messages
!     written (DEBUG) to fewest log messages written (ERROR). Compared
!     to the value of NABOUT when determining which messages to write
!     to the screen or to log files.
      INTEGER, PARAMETER :: DEBUG = -1  ! write all messages and echo input
      INTEGER, PARAMETER :: ECHO = 0    ! echo input, plus write all non-debug
      INTEGER, PARAMETER :: INFO = 1    ! don't echo input; write all non-debug
      INTEGER, PARAMETER :: WARNING = 2 ! don't echo input; write only warn/err
      INTEGER, PARAMETER :: ERROR = 3   ! don't echo input; only fatal msgs
!
      CHARACTER(len=10), dimension(5) :: logLevelNames
      CHARACTER(len=50), dimension(50) :: messageSources ! subroutine names
      CHARACTER(len=1024) :: scratchMessage ! used for formatted messages
      CHARACTER(len=1024) :: scratchFormat  ! used for Fortran format strings
      INTEGER :: sourceNumber ! index into messageSources for current sub

!     v49.48.01 tcm 20110126 Variables related to kdtree searches
      INTEGER :: SRCHDP=100
      TYPE(KDTREE2), POINTER :: TREE
      TYPE(KDTREE2_RESULT), ALLOCATABLE :: KDRESULTS(:)

      !jgf: water level offset at two time levels
      real(kind=8), allocatable :: dynamicWaterLevelCorrection1(:)
      real(kind=8), allocatable, target :: dynamicWaterLevelCorrection2(:)
      logical :: usingDynamicWaterLevelCorrection = .false.
      real(kind=8), allocatable, target :: dynamicWaterLevelCorrectionSta(:)

!----- DW --- periodic bc 
      INTEGER:: NPERSEG, NNPERBC
      INTEGER, ALLOCATABLE:: IPERCONN(:,:), PERBC_IDN_MAP(:)
!-------------------end of data declarations----------------------------------C

      integer, parameter :: max_clk = 40 ! max allowable number of clocks to register
      integer, private :: n_clk = 0 ! number of registered clocks
      integer, private :: clk_count_begin(max_clk) = -99999
      integer, private :: clk_count_end(max_clk) = -99999
      integer, private :: clk_count_total(max_clk) = 0
      integer, private :: clk_count_rate
      character(len=50), private :: clock_name(max_clk) = "null"

      CONTAINS

      function register_clock(name)
        character(len=*),  intent(in) :: name
        integer :: register_clock
        n_clk = n_clk+1
        if (n_clk == max_clk) then
            stop "ERROR: max number of clocks exceeded"
        endif
        clock_name(n_clk) = name
        register_clock = n_clk
      end function register_clock

      subroutine begin_clock(id)
        integer, intent(in) :: id
        if (clk_count_begin(id) == -99999) then
            call system_clock(clk_count_begin(id))
        else
            print *, "ERROR: begin_clock called twice before calling end_clock for ", trim(clock_name(id))
            stop
        endif
      end subroutine begin_clock

      subroutine end_clock(id)
        integer, intent(in) :: id
        if (clk_count_begin(id) /= -99999) then
            call system_clock(clk_count_end(id), clk_count_rate)
            clk_count_total(id) = clk_count_total(id) + clk_count_end(id) - clk_count_begin(id)
            clk_count_begin(id) = -99999
        else
            print *, "ERROR: end_clock called before begin clock for ", trim(clock_name(id))
            stop
        endif
      end subroutine end_clock

      subroutine log_clocks()
        integer :: id
        do id=1, max_clk
          if (trim(clock_name(id)) == "null") exit
          if (clk_count_begin(id) /= -99999) then
              print *, "ERROR: log_clock called before end_clock for ", trim(clock_name(id))
              stop
          else
              print *, "clock ", trim(clock_name(id)), ": ", real(clk_count_total(id), kind=8)/clk_count_rate
          endif
        enddo
      end subroutine log_clocks

!--------------------------------------------------------------------
!     S U B R O U T I N E    I N I T   L O G G I N G
!--------------------------------------------------------------------
!     jgf49.44: Initialize the names for the logging levels and the counter
!     for the current subroutine.
!--------------------------------------------------------------------
      SUBROUTINE initLogging()
      IMPLICIT NONE
!...
!...  OPEN STATEMENT FOR UNIT 16 OUTPUT FILE (ADCIRC LOG FILE)
!...
      sourceNumber = 0
      logLevelNames(1) = "DEBUG"
      logLevelNames(2) = "ECHO"
      logLevelNames(3) = "INFO"
      logLevelNames(4) = "WARNING"
      logLevelNames(5) = "ERROR"

!--------------------------------------------------------------------
      END SUBROUTINE initLogging
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!     S U B R O U T I N E    O P E N   L O G   F I L E
!--------------------------------------------------------------------
!     jgf50.65: Open the log file; this must be called after make dirname
!     so that we know where to put this log file.
!--------------------------------------------------------------------
      SUBROUTINE openLogFile()
      use sizes, only : localDir
      IMPLICIT NONE
!...
!...  OPEN STATEMENT FOR UNIT 16 OUTPUT FILE (ADCIRC LOG FILE)
!...
!     cms51.06: moved fort.33 file from read_input.F to here
#ifdef CMPI
      OPEN(16,FILE=trim(localdir)//'/'//'fort.16', ACTION='WRITE',
     &    STATUS='REPLACE')
      OPEN(33,FILE=trim(localdir)//'/'//'fort.33', ACTION='WRITE',
     &    STATUS='REPLACE')
#else
      OPEN(16,FILE='fort.16', ACTION='WRITE', STATUS='REPLACE')
      OPEN(33,FILE='fort.33', ACTION='WRITE', STATUS='REPLACE')
#endif
!--------------------------------------------------------------------
      END SUBROUTINE openLogFile
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!     S U B R O U T I N E    S C R E E N
!--------------------------------------------------------------------
!     jgf49.44: General purpose subroutine to write a message to
!     the screen with a certain "logging level", and subject to the
!     user's selection of where to write screen output. The logging
!     level is controlled by NABOUT from the fort.15 file. The actual
!     destination of messages written to the screen is controlled by
!     NSCREEN from the fort.15 file.
!
!     In parallel, only the processor with rank 0 actually writes
!     the message.
!
!     This subroutine assumes that the global variable "caller" has
!     been set to the name of the subroutine calling it. Therefore,
!     the setMessageSource subroutine must be called at the beginning
!     of the subroutine that calls this one, and unsetMessageSource
!     must be called at the end.
!--------------------------------------------------------------------
      SUBROUTINE screenMessage(level, message)
      USE SIZES, ONLY : myProc
      IMPLICIT NONE
      INTEGER, intent(in) :: level
      CHARACTER(*), intent(in) :: message
      INTEGER j  ! loop counter for stack
!
      IF (myProc.eq.0) THEN
         IF (NSCREEN.ne.0) THEN
            IF (level.ge.NABOUT) THEN
#ifdef FULL_STACK
               write(screenUnit,331,advance="no")
     &             trim(logLevelNames(level+2)),
     &             (trim(messageSources(j)),j=1,sourceNumber)
               write(screenUnit,332) trim(message)
#else
               write(screenUnit,333) trim(logLevelNames(level+2)),
     &             trim(messageSources(sourceNumber)), trim(message)
#endif
!#ifdef FLUSH_MESSAGES
               flush(screenUnit)
!#endif
            ENDIF
         ENDIF
      ENDIF
331   FORMAT(A,": ",A,50(:,"->",A))
332   FORMAT(": ",A)
333   FORMAT(A,": ",A,": ",A)
!--------------------------------------------------------------------
      END SUBROUTINE screenMessage
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!     S U B R O U T I N E    L O G   M E S S A G E
!--------------------------------------------------------------------
!     jgf49.44: General purpose subroutine to write a message to
!     the fort.16 file. In parallel, processors of all ranks will
!     write the message to their own subdomain fort.16 files.
!
!     This subroutine assumes that the global variable "caller" has
!     been set to the name of the subroutine calling it. Therefore,
!     the setMessageSource subroutine must be called at the beginning
!     of the subroutine that calls this one, and unsetMessageSource
!     must be called at the end.
!--------------------------------------------------------------------
      SUBROUTINE logMessage(level, message)
      IMPLICIT NONE
      INTEGER, intent(in) :: level
      CHARACTER(*), intent(in) :: message
      INTEGER j  ! loop counter for stack
!
      IF (level.ge.NABOUT) THEN
#ifdef FULL_STACK
         write(16,331,advance="no") trim(logLevelNames(level+2)),
     &       (trim(messageSources(j)),j=1,sourceNumber)
         write(16,332) trim(message)
#else
         write(16,333) trim(logLevelNames(level+2)),
     &       trim(messageSources(sourceNumber)), trim(message)
#endif
!#ifdef FLUSH_MESSAGES
               flush(16)
!#endif
      ENDIF
331   FORMAT(A,": ",A,50(:,"->",A))
332   FORMAT(": ",A)
333   FORMAT(A,": ",A,": ",A)
!--------------------------------------------------------------------
      END SUBROUTINE logMessage
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!     S U B R O U T I N E   A L L    M E S S A G E
!--------------------------------------------------------------------
!     jgf49.44: General purpose subroutine to write a message to
!     both the screen and to the fort.16 log file.
!--------------------------------------------------------------------
      SUBROUTINE allMessage(level, message)
      IMPLICIT NONE
      INTEGER, intent(in) :: level
      CHARACTER(*), intent(in) :: message
!
      call screenMessage(level, message)
      call logMessage(level, message)

!--------------------------------------------------------------------
      END SUBROUTINE allMessage
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!     S U B R O U T I N E   S E T   M E S S A G E   S O U R C E
!--------------------------------------------------------------------
!     jgf49.44: Sets the name of the subroutine that is writing
!     log and/or screen messages. Must use at the start of any subroutine
!     that calls screen, logMessage, or allMessage.
!--------------------------------------------------------------------
      SUBROUTINE setMessageSource(source)
      IMPLICIT NONE
      CHARACTER(*), intent(in) :: source
!
      sourceNumber = sourceNumber + 1
      messageSources(sourceNumber) = source
!--------------------------------------------------------------------
      END SUBROUTINE setMessageSource
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!     S U B R O U T I N E   U N S E T   M E S S A G E   S O U R C E
!--------------------------------------------------------------------
!     jgf49.44: Removes the name of the subroutine that is no longer
!     writing log and/or screen messages. Must use at the end of
!     any subroutine that calls screen, logMessage, or allMessage.
!--------------------------------------------------------------------
      SUBROUTINE unsetMessageSource()
      IMPLICIT NONE
!
      sourceNumber = sourceNumber - 1
!--------------------------------------------------------------------
      END SUBROUTINE unsetMessageSource
!--------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E   O P E N  F I L E  F O R  R E A D
!-----------------------------------------------------------------------
!     jgf50.16 Added general subroutine for opening an existing input
!     file for reading. Includes error checking.
!-----------------------------------------------------------------------
      SUBROUTINE openFileForRead(lun, filename, errorIO)
      IMPLICIT NONE
      INTEGER, intent(in) :: lun   ! fortran logical unit number
      CHARACTER(*), intent(in) :: filename ! full pathname of file
      INTEGER, intent(out) :: errorIO  ! zero if the file opened successfully
      LOGICAL :: fileFound    ! .true. if the file is present
!
      call setMessageSource("openFileForRead")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
!
      errorIO = 0
!
!     Check to see if file exists
      write(scratchMessage,21) lun
  21  format("Searching for file to open on unit ",i0,"...")
      call logMessage(INFO,trim(scratchMessage))
      inquire(FILE=trim(filename),EXIST=fileFound)
      if (fileFound.eqv..false.) then
         write(scratchMessage,23) trim(filename)
  23     format("The file '",A,"' was not found.")
         call allMessage(INFO,scratchMessage)
         errorIO = 1
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         return  ! file not found
      else
         write(scratchMessage,24) trim(filename)
  24     format("The file '",A,"' was found. The file will be opened.")
         call logMessage(INFO,trim(scratchMessage))
      endif
!
!     Open existing file
      OPEN(lun,FILE=trim(filename),STATUS='OLD',
     & ACTION='READ',IOSTAT=errorIO)
      if (errorIO.ne.0) then
          write(scratchMessage,25) trim(filename)
  25      format("Could not open the file '",A,"'.")
          call allMessage(ERROR,trim(scratchMessage))
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
          call allMessage(DEBUG,"Return.")
#endif
          call unsetMessageSource()
          return  ! file found but could not be opened
      else
         write(scratchMessage,26) trim(filename)
  26     format("The file '",A,"' was opened successfully.")
         call logMessage(INFO,trim(scratchMessage))
      endif

#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
!-----------------------------------------------------------------------
      END SUBROUTINE openFileForRead
!-----------------------------------------------------------------------

!     ----------------------------------------------------------------
!      F U N C T I O N   S P H E R I C A L   D I S T A N C E
!     ----------------------------------------------------------------
!     jgf49.1001 Function to get the distance along the surface of
!     a sphere (the earth's surface in this case).
!     ----------------------------------------------------------------
      REAL(kind=8) FUNCTION sphericalDistance(dx, dy, y1, y2)
      IMPLICIT NONE
      REAL(kind=8), intent(in) :: DX    ! longitude distance in radians
      REAL(kind=8), intent(in) :: DY    ! latitude distance in radians
      REAL(kind=8), intent(in) :: y1    ! degrees latitude of starting point
      REAL(kind=8), intent(in) :: y2    ! degrees latitude of ending point
!
!     compute the distances based on haversine formula for
!     distance along a sphere
      sphericalDistance = Rearth *
     &  ( 2.0e0*ASIN(
     &     sqrt(sin(DY/2.0e0)**2.0e0
     &      + cos(y1*DEG2RAD)*cos(y2*DEG2RAD)*sin(DX/2.0e0)**2.0e0
     &   )))

      RETURN
!     ----------------------------------------------------------------
      END FUNCTION sphericalDistance
!     ----------------------------------------------------------------

!...
!...Allocate space for Arrays dimensioned by MNE and MNP
!...
      SUBROUTINE ALLOC_MAIN1a()
      USE SIZES, ONLY : MNP, MNE
      IMPLICIT NONE
!  
      !jgfdebug TODO: rationalize initialization of these variables
      ALLOCATE ( ETAS(MNP))
      etas = -99999.e0
      allocate ( gwce_lv(mnp))
      gwce_lv = -99999.e0
      allocate ( mom_lv_x(mnp))
      mom_lv_x = -99999.e0
      allocate ( mom_lv_y(mnp))
      mom_lv_y = -99999.e0
      ALLOCATE ( UU1(MNP),VV1(MNP))
      uu1(:) = -99999.e0
      vv1(:) = -99999.e0
      ALLOCATE ( QX1(MNP),QY1(MNP))
      qx1(:) = -99999.e0
      qy1(:) = -99999.e0
      ALLOCATE ( NNODECODE(MNP))
      nnodecode(:) = -99999
      ALLOCATE ( NOFF(MNE) )
      noff(:) = -99999
      ALLOCATE ( NIBNODECODE(MNP))
      nibnodecode = -99999
      ALLOCATE ( CH1(MNP),TRANS_LV_B(MNP),TRANS_LV_A(MNP))
      ch1 = -99999.e0
      trans_lv_a = -99999.e0
      trans_lv_b = -99999.e0
      ALLOCATE ( SOURSIN(MNP))
      soursin = -99999.e0
      ALLOCATE ( TK(MNP), TKM(3,MNP))
      tk  = -99999.e0
      tkm = -99999.e0
      ALLOCATE ( UU2(MNP),VV2(MNP))
      uu2(:) = -99999.e0
      vv2(:) = -99999.e0
      ALLOCATE ( QX2(MNP),QY2(MNP))
      qx2(:) = -99999.e0
      qy2(:) = -99999.e0
      ALLOCATE ( ETA1(MNP),ETA2(MNP))
      eta1(:) = -99999.e0
      eta2(:) = -99999.e0
      ALLOCATE ( CORIF(MNP))
      corif(:) = -99999.e0      
      ALLOCATE ( NODECODE(MNP))
      nodecode = -99999
      ALLOCATE ( DASigT(MNP),DATemp(MNP),DASalt(MNP))
      DASigT = -99999.e0
      DATemp = -99999.e0
      DASalt = -99999.e0
      ALLOCATE ( VIDBCPDXOH(MNP), VIDBCPDYOH(MNP))
      VIDBCPDXOH = -99999.e0
      VIDBCPDYOH = -99999.e0
      ALLOCATE ( NB(MNP), NM(MNP), SigTS(MNP), MLD(MNP) )
      NB    = -99999.e0
      NM    = -99999.e0
      SigTS = -99999.e0
      MLD   = -99999.e0
      ALLOCATE ( VIDISPDXOH(MNP), VIDISPDYOH(MNP))
      VIDISPDXOH = -99999.e0
      VIDISPDYOH = -99999.e0
      ALLOCATE ( DAConc(MNP))
      DAConc = -99999.e0
!kmd   Added in parameters for the pc algorithm
      ALLOCATE ( UU0(MNP), VV0(MNP))
      ALLOCATE ( QX0(MNP), QY0(MNP) )
      ALLOCATE ( TK0(MNP), TK2(MNP))

!    kmd48.33bc - arrays for baroclinic
      ALLOCATE (q_heat1(MNP), q_heat2(MNP), HFLUX(MNP) )
      ALLOCATE (sponge(MNP) )

!    Allocate array for ESL i/o - only needed if
!    loadEleSlopeLim.eqv..true
      ALLOCATE (ESLCOUNT(MNP))
      ESLCOUNT = 0.e0

#ifdef CMPI
      ALLOCATE ( IDUMY(MNP), DUMY1(MNP), DUMY2(MNP) )
#endif

#ifdef CVEC
      ALLOCATE ( TEMP_LV_A(MNE,3),TEMP_LV_B(MNE,3) )
#endif
      ALLOCATE(EtaDisc(MNP))
      EtaDisc(:) = -99999.e0
      ALLOCATE ( ETAMAX(MNP)) ! v46.50 sb 11/11/2006
      ALLOCATE ( ETAMAX_TIME(MNP)) ! tcm v51.20.01 
      ALLOCATE ( UMAX(MNP)) ! v46.50 sb 11/11/2006
      ALLOCATE ( UMAX_TIME(MNP)) ! tcm v51.20.01 
      ETAMAX(:) = -99999.e0 ! v46.50 sb 11/11/2006
      EtaMax_Time(:) = 0.e0
      UMAX(:) = 0.e0 ! v46.50 sb 11/11/2006
      UMax_Time(:) = 0.e0
!-----------------------------------------------------------------------
      END SUBROUTINE ALLOC_MAIN1a
!-----------------------------------------------------------------------

!...
!...Allocate space for Arrays dimensioned by MNOPE and MNETA
!...
      SUBROUTINE ALLOC_MAIN2()
      USE SIZES, ONLY : MNETA
      IMPLICIT NONE
!
      ALLOCATE ( ESBIN1(MNETA),ESBIN2(MNETA))

!    kmd48.33bc - arrays for boundary conditions
      IF (CBAROCLINIC) THEN
         ALLOCATE ( LNM_BC1(MNETA), LNM_BC2(MNETA) )
         ALLOCATE ( LNM_BC(MNETA) )
      END IF
!
      RETURN
      END SUBROUTINE ALLOC_MAIN2


!...
!...Allocate space for nonperiodic zero and nonzero normal flow
!...boundary arrays including barriers
!...
      SUBROUTINE ALLOC_MAIN3()
      USE SIZES, ONLY : MNVEL
      IMPLICIT NONE
!
      ALLOCATE ( QN0(MNVEL),QN1(MNVEL),QN2(MNVEL))
      ALLOCATE ( EN0(MNVEL),EN1(MNVEL),EN2(MNVEL))
      ALLOCATE ( ENIN1(MNVEL),ENIN2(MNVEL))
      ALLOCATE ( QNIN1(MNVEL),QNIN2(MNVEL))
!
!     jgf46.21 Added support for IBTYPE=52.
      ALLOCATE ( ElevDisc(MNVEL))
!
!     jgf49.58 Explicitly initialize variables that would otherwise be
!     used without initialization.
      EN2(:) = -99999.e0

      END SUBROUTINE ALLOC_MAIN3

!...
!...Allocate space for tidal potential terms
!...
      SUBROUTINE ALLOC_MAIN4a()
      USE SIZES, ONLY : MNTIF, MNP
      IMPLICIT NONE
!
      ALLOCATE ( TPK(MNTIF),AMIGT(MNTIF),FFT(MNTIF) )
      ALLOCATE ( FACET(MNTIF),PERT(MNTIF),ETRF(MNTIF) )
      ALLOCATE ( TIPOTAG(MNTIF) )

      IF ( CTIP ) THEN
         ALLOCATE( TIP1(MNP),TIP2(MNP),L_N(0:2,MNP))
      ENDIF
!
      RETURN
      END SUBROUTINE ALLOC_MAIN4a


!...
!...Allocate space for Earth load/self-attraction tide
!...
      SUBROUTINE ALLOC_MAIN4b()
      USE SIZES, ONLY : MNTIF, MNP
      IMPLICIT NONE
!
      ALLOCATE ( SALTAMP(MNTIF,MNP),SALTPHA(MNTIF,MNP) )
!
      RETURN
      END SUBROUTINE ALLOC_MAIN4b


!...
!...Allocate space for Arrays dimensioned by MNBFR
!...
      SUBROUTINE ALLOC_MAIN5()
      USE SIZES, ONLY : MNBFR, MNETA
      IMPLICIT NONE
!
      ALLOCATE ( AMIG(MNBFR),PER(MNBFR))
      ALLOCATE ( FF(MNBFR),FACE(MNBFR))
      ALLOCATE ( EMO(MNBFR,MNETA),EFA(MNBFR,MNETA))
      ALLOCATE ( BOUNTAG(MNBFR) )
!
      RETURN
      END SUBROUTINE ALLOC_MAIN5


!...
!...Allocate space for periodic normal flow boundary conditions
!...
      SUBROUTINE ALLOC_MAIN6()
      USE SIZES, ONLY : MNFFR, MNVEL
      IMPLICIT NONE
!
      ALLOCATE ( QNAM(MNFFR,MNVEL),QNPH(MNFFR,MNVEL))
      ALLOCATE ( ENAM(MNFFR,MNVEL),ENPH(MNFFR,MNVEL))
      ALLOCATE ( FBOUNTAG(MNFFR))
      ALLOCATE ( FAMIG(MNFFR), FFF(MNFFR), FFACE(MNFFR), FPER(MNFFR) )
!
      RETURN
      END SUBROUTINE ALLOC_MAIN6


!...
!...Allocate space for arrays used for station elevation output
!...
      SUBROUTINE ALLOC_MAIN7()
      USE SIZES, ONLY : MNSTAE
      IMPLICIT NONE
!
      ALLOCATE ( NNE(MNSTAE),ET00(MNSTAE))
      ALLOCATE ( STAIE1(MNSTAE),STAIE2(MNSTAE),STAIE3(MNSTAE))
      ALLOCATE ( XEL(MNSTAE),YEL(MNSTAE),SLEL(MNSTAE),SFEL(MNSTAE))

!...  tcm v50.66.01 added for time varying bathymetry
      IF (NDDT.NE.0) ALLOCATE (DP00(MNSTAE) )
!
      RETURN
      END SUBROUTINE ALLOC_MAIN7


!...
!...Allocate space for arrays used for station velocity output
!...
      SUBROUTINE ALLOC_MAIN8()
      USE SIZES, ONLY : MNSTAV
      IMPLICIT NONE
!
      ALLOCATE ( XEV(MNSTAV),YEV(MNSTAV),SLEV(MNSTAV),SFEV(MNSTAV))
      ALLOCATE ( NNV(MNSTAV))
      ALLOCATE ( UU00(MNSTAV),VV00(MNSTAV))
      ALLOCATE ( STAIV1(MNSTAV),STAIV2(MNSTAV),STAIV3(MNSTAV))
!
      RETURN
      END SUBROUTINE ALLOC_MAIN8


!...
!...Allocate space for arrays used for station concentration output
!...
      SUBROUTINE ALLOC_MAIN9()
      USE SIZES, ONLY : MNSTAC
      IMPLICIT NONE
!
      ALLOCATE ( XEC(MNSTAC),YEC(MNSTAC),SLEC(MNSTAC),SFEC(MNSTAC))
      ALLOCATE ( NNC(MNSTAC))
      ALLOCATE ( CC00(MNSTAC))
      ALLOCATE ( STAIC1(MNSTAC),STAIC2(MNSTAC),STAIC3(MNSTAC))
!
      RETURN
      END SUBROUTINE ALLOC_MAIN9


!...
!...Allocate space for arrays used for station meteorological output
!...
      SUBROUTINE ALLOC_MAIN10()
      USE SIZES, ONLY : MNSTAM
      IMPLICIT NONE
!
      ALLOCATE ( XEM(MNSTAM),YEM(MNSTAM),SLEM(MNSTAM),SFEM(MNSTAM))
      ALLOCATE ( NNM(MNSTAM))
      ALLOCATE ( RMU00(MNSTAM),RMV00(MNSTAM),RMP00(MNSTAM))
      if (NCICE.NE.0) then
         ALLOCATE ( RMICE00(MNSTAM) )    ! v49.64.01 tcm added for ice concentration
      endif
      ALLOCATE ( STAIM1(MNSTAM),STAIM2(MNSTAM),STAIM3(MNSTAM))
!
      RETURN
      END SUBROUTINE ALLOC_MAIN10


!...
!...Allocate space for arrays needed by GWCE matrix and iterative solver
!...

      SUBROUTINE ALLOC_MAIN11_LUMPED()
      USE SIZES, ONLY : MNP
      IMPLICIT NONE
!     jgf48.4619: Array used in lumped GWCE (from Seizo)
      ALLOCATE( COEFD(MNP) )
      coefd = -99999.e0
      RETURN
      END SUBROUTINE ALLOC_MAIN11_LUMPED




      SUBROUTINE ALLOC_MAIN11()
      USE SIZES, ONLY : MNETA,MNEI,MNP
      IMPLICIT NONE

!  Arrays used by JCG iterative solver

      ALLOCATE( OBCCOEF(MNETA,MNEI-1),COEF(MNP,MNEI))
      obccoef = -99999.e0
      coef = -99999.e0
!
      RETURN
      END SUBROUTINE ALLOC_MAIN11


!...
!...Allocate space for wind forcing
!...
      SUBROUTINE ALLOC_MAIN12()
      USE SIZES, ONLY : MNP
      IMPLICIT NONE
!
      ALLOCATE ( WSX1(MNP),WSY1(MNP),PR1(MNP) )
      ALLOCATE ( WSX2(MNP),WSY2(MNP),PR2(MNP) )

      WSX1(:) = -99999.0e0
      WSY1(:) = -99999.0e0
      PR1(:) = -99999.0e0
      WSX2(:) = -99999.0e0
      WSY2(:) = -99999.0e0
      PR2(:) = -99999.0e0

      allocate ( wvnxout(mnp),wvnyout(mnp) )
      wvnxout(:) = -99999.e0
      wvnyout(:) = -99999.e0

      ALLOCATE ( RSNX1(MNP),RSNY1(MNP),RSNX2(MNP),RSNY2(MNP) )
      RSNX1(:) = 0.e0 !jgf51.46: Added initialization of RSNX1, RSNY1
      RSNY1(:) = 0.e0 ! must be initialized to zero for SWAN
      RSNX2(:) = 0.e0
      RSNY2(:) = 0.e0

      if (NCICE.NE.0) then
         ALLOCATE ( CICE1(MNP),CICE2(MNP),CICEOUT(MNP) ) !tcm v49.64.01 added
                                                        ! for ice concentration
         CICE1(:) = 0.0e0 !tcm v49.64.01 added for ice concentration
         CICE2(:) = 0.0e0 !tcm v49.64.01 added for ice concentration
         CICEOUT(:) = 0.0e0 !tcm v49.64.01 added for ice concentration
      ENDIF
      
#ifdef CSWAN
      ALLOCATE ( SWAN_WDragCo(MNP) ) !MVB
      SWAN_WDragCo(:) = 0.0e0 !MVB
#endif
      allocate (windDragOut(mnp))
      windDragOut(:) = 0.e0

!  tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
      IF ((ABS(NRS).EQ.3).or.(ABS(NRS).EQ.4)) then
!#ifdef CSWAN
!asey 090302: Added the next line for output of radiation stress gradients.
         ALLOCATE ( RSNXOUT(MNP),RSNYOUT(MNP) )
         RSNXOUT(:) = 0.e0  !tcm v50.75 added initialization to zero
         RSNYOUT(:) = 0.e0
      ENDIF
!#endif
      ALLOCATE ( PrMin(MNP)) ! v46.50 sb 11/11/2006
      PrMin(:) = 99999.e0
      ALLOCATE (PrMin_Time(MNP) ) !tcm v51.20.01
      PrMin_Time(:) = 0.e0
      ALLOCATE ( WVNOUTMax(MNP)) ! v46.50 sb 11/11/2006
      WVNOUTMAX(:) = 0.e0
      ALLOCATE (WVNOutMax_Time(MNP) ) !tcm v51.20.01
      WVNOutMax_Time(:) = 0.e0
      ALLOCATE ( RSNMax(MNP))               ! v46.50 sb 11/11/2006
      RSNMax(:) = 0.e0
      ALLOCATE (RSNMax_Time(MNP) ) !tcm v51.20.01
      RSNMax_Time(:) = 0.e0

!
      RETURN
      END SUBROUTINE ALLOC_MAIN12

! TCM V50.66.01 ADDITIONS FOR TIME VARYING BATHYMETRY (ALLOC_MAIN13)
!...
!...Allocate space for TIME VARYING BATHYMETRY
!...
      SUBROUTINE ALLOC_MAIN13()
      USE SIZES, ONLY : MNP
      IMPLICIT NONE
      
      ALLOCATE ( DP1(MNP), DP2(MNP) )
      DP1(:) = 0.E0
      DP2(:) = 0.E0

      END SUBROUTINE ALLOC_MAIN13

      !-----------------------------------------------------------------
      !                      S U B R O U T I N E   
      !                G E T    F O R M A T    I N T E G E R
      !-----------------------------------------------------------------
      ! Sets file format integer that corresponds to the file 
      ! format string.
      !-----------------------------------------------------------------
      subroutine getFormatInteger(typeString, typeHolder)
      implicit none
      character(len=20), intent(in) :: typeString
      integer, intent(out) :: typeHolder
      !
      select case(trim(typeString))
      case('OFF')
         typeHolder = 0 
      case('ASCII')
         typeHolder = 1
      case('BINARY')
         typeHolder = 2 
      case('NETCDF3')
         typeHolder = 3 
      case('SPARSE_ASCII')
         typeHolder = 4 
      case('NETCDF4')
         typeHolder = 5 
      case('XDMF')
         typeHolder = 7
      case default
         write(6,'("WARNING: Unrecognized format ",i0,".")') trim(typeString)
      end select
      !----------------------------------------------------------------------
      end subroutine getFormatInteger
      !----------------------------------------------------------------------


!-------------------------------------------------------------------

      function VERSION_NUMBER(Major, Minor, Rev) result(vn)

      implicit none
      integer   :: vn, major, minor, rev

      vn = ior(ior(ishft(major,20),ishft(minor,10)),rev)

      return
      end function VERSION_NUMBER

!-------------------------------------------------------------------

      function CMP_VERSION_NUMBERS(a,b) result(match)

      implicit none
      integer  a,b
      logical  :: match

      match = ishft(a,-10) == ishft(b,-10)

      return
      end function CMP_VERSION_NUMBERS
!-------------------------------------------------------------------

      !=================================================================
      !=================================================================
      !=================================================================
      !      =====                                            =====
      !      =====            SUBROUTINE PARSE                =====
      !      =====                                            =====
      !=================================================================
      !=================================================================
      !=================================================================

      !=================================================================
      ! This subroutines parse station names so they can passed on to
      ! the NetCDF I/O module of ADCIRC and be written out in the netCDF
      ! output files.
      !
      ! Revision history:
      !
      ! Date      Programmer                       Description of change
      ! ----      ----------                       ---------------------
      ! 12/16/07  Cristina Forbes, PSGS @ UNC-IMS  Wrote parser
      !                                            based on character
      !                                            manipulation subroutines
      !                                            developed at the
      !                                            University of Oslo, Norway
      ! v49.29    Jason Fleming  Fixes to avoid placing leading exclamation
      !                          points in the station names, and to make use
      !                          of trailing exclamation points as station name
      !                          terminators or comment separators.
      ! v50.36    Chris Massey   Added logic to parse to handle tabs and commas
      !                          and not just spaces.
      !=================================================================
      subroutine parse (line, lvar)
      IMPLICIT NONE
      CHARACTER(132) LINE,L
      CHARACTER(50) LVAR(3)
      INTEGER LS, L1, LEN1, I
      INTEGER L2,L3,L4,PB
      INTEGER NAMSTART  ! position of beginning of station name
      INTEGER NAMEND    ! position of end of station name
      INTEGER BANGEND   ! location of terminating "!", if any
      CHARACTER(len=1) :: tab_character
!
      call setMessageSource("parse")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      tab_character = achar(9)
      LS=1
      L=TRIM(ADJUSTL(LINE))
      CALL COMPACT(L)
      LEN1=LEN(L)
!     EXTRACT STATION DATA FROM LINE
      L2 = 0
      L3 = 0
      L4 = 0
      DO I=1,2
!        DO I=1,3
         if (index(L," ").ne.0) then
            L2=LEN(TRIM(ADJUSTL(L(1:INDEX(L," ")-1))))          !space
         endif
         if (index(L,tab_character).ne.0) then
            L3=LEN(Trim(adjustl(L(1:index(L,tab_character)-1))))     !tab
         endif
         if (index(L,",").ne.0) then
            L4=LEN(trim(adjustl(L(1:index(L,",")-1))))          !comma
         endif
         !Determine starting position based on which delimiter was found first
         L1 = L2
         PB = 1
         IF(L1.GT.L3.AND.L3.NE.0) THEN
            L1 = L3
            PB = 2
         ENDIF
         IF(L1.GT.L4.AND.L4.NE.0) THEN
           L1=L4
           PB = 2
         ENDIF

         LVAR(I)=ADJUSTL(L(LS:L1))
         L=ADJUSTL(L(L1+1:LEN(L)))
         CALL COMPACT(L)

      END DO
!
!     Find start of station name, ignoring leading "!" if it is there
      NAMSTART=1
      IF ( L(1:1).EQ."!" ) THEN
         NAMSTART=2
      ENDIF
!
!     Find end of station name, ignoring the first occurrence of an
!     embedded "!" as well as any characters to the right of that, since
!     they are assumed to be a comment.
      NAMEND=LEN(L)
      BANGEND=INDEX(L(NAMSTART:NAMEND),"!")
      IF ( BANGEND.NE.0 ) THEN
         NAMEND = BANGEND + (NAMSTART-2)
      ENDIF
      LVAR(3)=ADJUSTL(L(NAMSTART:NAMEND)) !jgf49.29: trim leading "!" and spaces
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !=================================================================
      END SUBROUTINE parse
      !=================================================================

      !=================================================================
       SUBROUTINE COMPACT(STR)
!     tcm50.36 SUBROUTINE TO CONVERT MULTIPLE SPACES AND TABS TO SINGLE
!              SPACES; DELETES CONTROL CHARACTERS;REMOVES INITIAL SPACES.
!              Adapted from stringmod.f90 from the website:
!              http://www.gbenthien.net/strings/index.html
      !=================================================================
       IMPLICIT NONE

       INTEGER :: I,K,ISP,LENSTR,ICH
       CHARACTER(LEN=*):: STR
       CHARACTER(LEN=1):: CH
       CHARACTER(LEN=LEN_TRIM(STR)):: OUTSTR

      call setMessageSource("compact")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
       STR=ADJUSTL(STR)
       LENSTR=LEN_TRIM(STR)
       OUTSTR=' '
       ISP=0
       K=0

       DO I=1,LENSTR
          CH=STR(I:I)
          ICH=IACHAR(CH)

          SELECT CASE(ICH)

          CASE(9,32)     ! SPACE OR TAB CHARACTER
             IF(ISP==0) THEN
                K=K+1
                OUTSTR(K:K)=' '
             END IF
             ISP=1

          CASE(33:)      ! NOT A SPACE, QUOTE, OR CONTROL CHARACTER
             K=K+1
             OUTSTR(K:K)=CH
             ISP=0

          END SELECT

       END DO

       STR=ADJUSTL(OUTSTR)
!
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !=================================================================
       END SUBROUTINE COMPACT
      !=================================================================


      !=================================================================
        SUBROUTINE i2a(in2t, res)
!     jgf49.31 Subroutine to convert an integer number to a string. Moved
!     here from read_input.F to make it more generally available in different
!     modules.
      !=================================================================
        IMPLICIT NONE
!// Arguments
        INTEGER, INTENT(IN)            :: in2t
        CHARACTER(LEN=*), INTENT(OUT)  :: res
!// Local variables
        INTEGER                        :: i,j,k
        CHARACTER(LEN=80)              :: sbuf
        sbuf = ' '
        res = ' '
        i = in2t
        k = 1
        DO
        j = MOD(i,10)
        sbuf(k:k) = ACHAR(j+48)
        i = i / 10
        IF (i <= 0) THEN
        EXIT
        END IF
        k = k + 1
        END DO
        k = LEN_TRIM(sbuf)
        IF (k>1) THEN
        j = k
        i = 1
        DO i = 1, k
         res(i:i) = sbuf(j:j)
         j = j - 1
        END DO
        ELSE
        res(1:1) = sbuf(1:1)
        END IF
      !=================================================================
        END SUBROUTINE i2a
      !=================================================================
!
      !=================================================================
      FUNCTION front_trim(buf) RESULT(res)
!     jgf49.31 Subroutine to return the index of the first non-space
!     character of a string.
      !=================================================================
      IMPLICIT NONE
      !// Arguments
      CHARACTER(LEN=*), INTENT(IN)   :: buf
      INTEGER                        :: res
      !// Local variables
      INTEGER                        :: i, lng
      res = 1
      lng = LEN_TRIM(buf)
      DO i = 1, lng
         IF (buf(i:i) .NE. ' ') THEN
            res = i
            RETURN
         END IF
      END DO
      !=================================================================
      END FUNCTION front_trim
      !=================================================================
!
      !=================================================================
      FUNCTION a2d(buf) RESULT(res)
!     jgf49.31 Subroutine to convert a double precision real number to
!     a string. Moved here from read_input.F to make it more generally
!     available to different parts of ADCIRC.
      !=================================================================
      IMPLICIT NONE
      !// Arguments
        CHARACTER(LEN=*)               :: buf
        DOUBLE PRECISION               :: res
!// Local variables
        INTEGER                        :: i,j,k,itmp, foffs, lng
        INTEGER                        :: p, q
        DOUBLE PRECISION               :: dtmp
        LOGICAL                        :: neg
        LOGICAL                        :: exponential
        neg = .FALSE.
        exponential = .FALSE.
        lng = LEN_TRIM(buf)
        foffs = front_trim(buf)
        k = INDEX(buf,'.')
        p = INDEX(buf,'e')
        q = INDEX(buf,'E')
        IF(p /= -1 .AND. p /= 0) THEN
        exponential = .TRUE.
        END IF
        IF(q /= -1 .AND. q /= 0) THEN
        exponential = .TRUE.
        p = q
        END IF
        IF (k /= -1 .AND. k /= 0) THEN
!// We have a floating point number
        itmp = 0
!// Get the integer part of the number
        DO i = foffs, k - 1
        IF(buf(i:i) .EQ. '-') THEN
        neg = .TRUE.
        CONTINUE
        END IF
        IF(buf(i:i) .GE. '0' .AND. buf(i:i) .LE. '9' .AND.
     &  buf(i:i) .NE. ' ') THEN
        itmp = itmp * 10
        itmp = itmp + (IACHAR(buf(i:i)) - 48)
        END IF
        END DO
        res = DBLE(itmp)
        dtmp = 0.
        q = 0
        IF(.NOT. exponential) THEN
!// We do not have an exponential number
        DO i = LEN_TRIM(buf), k+1, -1
        itmp = (IACHAR(buf(i:i)) - 48)
        dtmp = dtmp + FLOAT(itmp)
        dtmp = dtmp / 10.
        END DO
        ELSE
!// We have an exponential number
        DO i = p-1, k+1, -1
        itmp = (IACHAR(buf(i:i)) - 48)
        dtmp = dtmp + FLOAT(itmp)
        dtmp = dtmp / 10.
        END DO
        q = a2i(buf(p+1:LEN_TRIM(buf)))
        END IF
        res = res + dtmp
        IF(exponential) THEN
        res = res * 10.**q
        END IF
        IF(neg) THEN
        res = res * (-1)
        END IF
        ELSE
!// We have an integer
        res = DBLE(a2i(buf))
        END IF
      !=================================================================
        END FUNCTION a2d
      !=================================================================
!
      !=================================================================
      FUNCTION a2f(buf) RESULT(res)
!     jgf49.31 Subroutine to convert a single precision real number to
!     a string. Moved here from read_input.F to make it more generally
!     available to different parts of ADCIRC.
      !=================================================================
      IMPLICIT NONE
      !// Arguments
      CHARACTER(LEN=*)               :: buf
      REAL                           :: res
      !// Local variables
      INTEGER                        :: i,j,k,itmp, foffs, lng
      INTEGER                        :: p, q
      REAL                           :: rtmp
      LOGICAL                        :: neg
      LOGICAL                        :: exponential
      neg = .FALSE.
      exponential = .FALSE.
      lng = LEN_TRIM(buf)
      foffs = front_trim(buf)
      k = INDEX(buf,'.')
      p = INDEX(buf,'e')
      q = INDEX(buf,'E')
      IF(p /= -1 .AND. p /= 0) THEN
         exponential = .TRUE.
      END IF
      IF(q /= -1 .AND. q /= 0) THEN
         exponential = .TRUE.
         p = q
      END IF
      IF (k /= -1 .AND. k /= 0) THEN
         !// We have a floating point number
         itmp = 0
         !// Get the integer part of the number
         DO i = foffs, k - 1
            IF(buf(i:i) .EQ. '-') THEN
               neg = .TRUE.
               CONTINUE
            END IF
            IF(buf(i:i) .GE. '0' .AND. buf(i:i) .LE. '9' .AND.
     &           buf(i:i) .NE. ' ') THEN
               itmp = itmp * 10
               itmp = itmp + (IACHAR(buf(i:i)) - 48)
            END IF
         END DO
         res = DBLE(itmp)
         rtmp = 0.
         q = 0
         IF(.NOT. exponential) THEN
!// We do not have an exponential number
            DO i = LEN_TRIM(buf), k+1, -1
               itmp = (IACHAR(buf(i:i)) - 48)
               rtmp = rtmp + FLOAT(itmp)
               rtmp = rtmp / 10.
            END DO
         ELSE
!// We have an exponential number
            DO i = p-1, k+1, -1
               itmp = (IACHAR(buf(i:i)) - 48)
               rtmp = rtmp + FLOAT(itmp)
               rtmp = rtmp / 10.
            END DO
            q = a2i(buf(p+1:LEN_TRIM(buf)))
         END IF
         res = res + rtmp
         IF(exponential) THEN
            res = res * 10.**q
         END IF
         IF(neg) THEN
            res = res * (-1)
         END IF
      ELSE
         !// We have an integer
         res = REAL(a2i(buf))
      END IF
      !=================================================================
      END FUNCTION a2f
      !=================================================================
!
      !=================================================================
      FUNCTION a2i(buf) RESULT(res)
!     jgf49.31 Subroutine to convert a string to an integer.
!     Moved here from read_input.F to make it more generally
!     available to different parts of ADCIRC.
      !=================================================================
      IMPLICIT NONE
      !// Arguments
      CHARACTER(LEN=*), INTENT(IN)   :: buf
      INTEGER                        :: res
      !// Local variables
      INTEGER                        :: i, foffs, lng
      LOGICAL                        :: neg
      neg = .FALSE.
      lng = LEN_TRIM(buf)
      foffs = front_trim(buf)
      res = 0
      DO i = foffs, lng
         IF(buf(i:i) .EQ. '-') THEN
            neg = .TRUE.
           CONTINUE
         END IF
         IF(buf(i:i) .GE. '0' .AND. buf(i:i) .LE. '9' .AND.
     &        buf(i:i) .NE. ' ') THEN
            res = res * 10
            res = res + (IACHAR(buf(i:i)) - 48)
         END IF
      END DO
      IF(neg) THEN
         res = res * (-1)
      END IF
      !=================================================================
      END FUNCTION a2i
      !=================================================================

!-----------------------------------------------------------------------
!     @brief Internal norm2 function
!     If the compiler supports Fortran 2008 then the intrinsic is used. 
!     Otherwise, the implementation is provided here. Use the compiler 
!     flag -DNOF2008 to enable the internal implementation. norm2
!     calculates the Euclidean vector norm (L2 norm) of an array
!     @param[in] x array to calculate the L2 norm of 
!     @result euclidean vector norm of x
!-----------------------------------------------------------------------
      pure function adcirc_norm2(x)
!-----------------------------------------------------------------------
        implicit none
        REAL(kind=8), intent(in) :: x(:)
        real(kind=8)             :: adcirc_norm2
#ifdef NOF2008
        intrinsic            :: dot_product, sqrt
        adcirc_norm2 = sqrt(dot_product(x,x))
#else 
        intrinsic            :: norm2
        adcirc_norm2 = norm2(x)
#endif
!-----------------------------------------------------------------------
      end function adcirc_norm2
!-----------------------------------------------------------------------

!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
      FUNCTION toLowercase(string) RESULT(lowercase)
          IMPLICIT NONE
          CHARACTER(*),INTENT(IN)  :: string
          CHARACTER(LEN(string))   :: lowercase
          CHARACTER(26), PARAMETER :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
          CHARACTER(26), PARAMETER :: low = 'abcdefghijklmnopqrstuvwxyz'
          INTEGER                  :: I
          INTEGER                  :: IDX
          lowercase = string
          DO I = 1,LEN_TRIM(lowercase)
            idx = INDEX(cap,lowercase(I:I))
            IF(idx.GT.0)lowercase(I:I)=low(IDX:IDX)
          ENDDO
      END FUNCTION toLowercase
      
      FUNCTION toUppercase(string) RESULT(uppercase)
          IMPLICIT NONE
          CHARACTER(*),INTENT(IN)  :: string
          CHARACTER(LEN(string))   :: uppercase
          CHARACTER(26), PARAMETER :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
          CHARACTER(26), PARAMETER :: low = 'abcdefghijklmnopqrstuvwxyz'
          INTEGER                  :: I
          INTEGER                  :: IDX
          uppercase = string
          DO I = 1,LEN_TRIM(uppercase)
            idx = INDEX(low,uppercase(I:I))
            IF(idx.GT.0)uppercase(I:I)=cap(IDX:IDX)
          ENDDO
      END FUNCTION toUppercase


!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
      END MODULE GLOBAL
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
