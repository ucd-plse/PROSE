PROGRAM MOM6_ODA_offline_driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A parallel offline driver for MOM6 data assimilation
!
! Authors: Matthew.Harrison@noaa.gov
!          Feiyu.Liu@noaa.gov and
!          Tony.Rosati@noaa.gov
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use MOM_file_parser, only : get_param, log_param, param_file_type, close_param_file
  use MOM_get_input, only : get_MOM_input, directories
  use MOM_domains, only : MOM_domains_init, MOM_infra_init, MOM_infra_end, clone_MOM_domain, create_group_pass
  use MOM_grid, only : ocean_grid_type, MOM_grid_init
  use MOM_io, only : slasher, MOM_io_init, ASCII_FILE, READONLY_FILE
  use MOM_io, only : check_nml_error, file_exists, open_file, close_file, io_infra_init
  use MOM_io, only : vardesc, var_desc
  use MOM_error_handler, only : FATAL, WARNING, MOM_error
  use MOM_hor_index, only : hor_index_type, hor_index_init
  use MOM_grid, only : ocean_grid_type
  use MOM_grid_initialize, only : set_grid_metrics
  use MOM_fixed_initialization, only : MOM_initialize_fixed, MOM_initialize_topography
  use MOM_coord_initialization, only : MOM_initialize_coord
  use MOM_transcribe_grid,       only : copy_dyngrid_to_MOM_grid, copy_MOM_grid_to_dyngrid
  use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
  use MOM_verticalGrid, only : verticalGrid_type, verticalGridInit
  use MOM_diag_mediator, only : diag_ctrl
  !use MOM_diagnostics, only : transport_diag_IDS, surface_diag_IDs
  use MOM_variables, only : thermo_var_ptrs, vertvisc_type, accel_diag_ptrs, cont_diag_ptrs
  use MOM_restart, only : register_restart_field, query_initialized, restart_init, MOM_restart_CS, save_restart
  use MOM_state_initialization, only : MOM_initialize_state
  use MOM_tracer_registry,       only : tracer_registry_type
  use MOM_sponge,                only : sponge_CS
  use MOM_ALE_sponge,            only : ALE_sponge_CS
  use MOM_ALE,                   only : ALE_CS
  use MOM_open_boundary, only         : ocean_OBC_type
  use MOM_checksum_packages, only : MOM_state_chksum, MOM_thermo_chksum, MOM_state_stats
  use time_manager_mod, only: time_type, set_time, set_date, JULIAN, NOLEAP, NO_CALENDAR, set_calendar_type
  use time_manager_mod, only: time_manager_init, get_time
  use ensemble_manager_mod, only : get_ensemble_size, ensemble_manager_init, ensemble_pelist_setup
  use ensemble_manager_mod, only : get_ensemble_id, get_ensemble_pelist
  use mpp_mod, only : set_current_pelist => mpp_set_current_pelist
  use mpp_mod, only : get_current_pelist => mpp_get_current_pelist
  use mpp_mod, only : set_root_pe => mpp_set_root_pe, sync_self => mpp_sync_self
  use mpp_mod, only : declare_pelist => mpp_declare_pelist, pe => mpp_pe, npes => mpp_npes
  use fms_mod, only : fms_end
  use diag_manager_mod, only : diag_manager_init, diag_manager_end, get_base_date
  use MOM_oda_driver_mod, only : ODA_CS, oda, init_oda, oda_end
  use MOM_oda_driver_mod, only : set_prior_tracer, set_analysis_time, get_posterior_tracer, save_obs_diff
  use oda_types_mod, only : ocean_profile_type
  implicit none

#include <MOM_memory.h>

!! A scaled-down version of the MOM control structure, including the variables that describe
!! the state of the ocean.
type :: MOM_control_struct
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: &
    h, &      !< layer thickness (m or kg/m2 (H))
    T, &      !< potential temperature (degrees C)
    S         !< salinity (ppt)
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    u         !< zonal velocity component (m/s)
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    v         !< meridional velocity (m/s)
  type(ocean_grid_type) :: G  !< structure containing metrics and grid info
  type(verticalGrid_type), pointer :: &
    GV => NULL()              !< structure containing vertical grid info
  type(thermo_var_ptrs) :: tv !< structure containing pointers to available
                              !! thermodynamic fields
  type(time_type), pointer :: Time   !< pointer to ocean clock
  type(MOM_restart_CS), pointer :: restart_CSp !<pointer set in this routine to the restart control structure
  type(tracer_registry_type),    pointer :: tracer_Reg             => NULL()
  type(ocean_OBC_type),          pointer :: OBC                    => NULL()
  type(sponge_CS),               pointer :: sponge_CSp             => NULL()
  type(ALE_sponge_CS),           pointer :: ALE_sponge_CSp         => NULL()
  type(ALE_CS),                  pointer :: ALE_CSp                => NULL()
end type MOM_control_struct

  type(ocean_grid_type), pointer :: Grid=> NULL() !<temporary pointer to ensemble member horizontal grid
  type(verticalGrid_type), pointer :: GV=> NULL() !< pointer to the ocean vertical grid type
  type(dyn_horgrid_type), pointer :: dG=> NULL() !<temporary pointer to ensemble member horizontal grid
  type(MOM_control_struct), pointer :: CS=> NULL() !<pointer to an element in CSp
  type(ODA_CS), pointer :: odaCS => NULL() !<pointer to ocean data assimilation control structure
  type(hor_index_type), pointer :: HI=> NULL() !<pointer to ensemble member indexing object
  type(time_type) :: time, time_start, time_start_segment !<time types
  type(param_file_type), pointer, dimension(:) :: PFp => NULL() !<pointer to parameter file list for ensemble members
  type(param_file_type), pointer :: PF=> NULL() !<pointer to parameter file
  type(thermo_var_ptrs), pointer :: tv => NULL() !<pointer to an ocean thermo_var type
  real, dimension(:,:,:), pointer :: h => NULL()
  type(directories) :: dirs !<directories object
  type(vardesc) :: vd !<variable descriptor
  type(ocean_profile_type), pointer :: Prof=> NULL() !<pointer to an ocean profile
  integer, dimension(:,:), pointer :: ens_pelist=> NULL() !<pointer to ocean pelist for each ensemble member
  integer :: ensemble_id
  integer :: ens_size, npes_ocn, npes_atm, npes_lnd, npeS_ice, ens_info(6)
  integer, dimension(:), allocatable :: ocn_pelist, atm_pelist, lnd_pelist, ice_pelist, global_pelist
  integer :: calendar_type=-1
  character(len=16) :: month='jan',calendar='julian'
  !<it is currently hard-coded since we are only running an ocean
  !<model for this test and do not need to reserve separate PEs for the atm/land/ice
  logical :: get_increment = .true.
  integer :: date_init(6), date(6), years=0, months=0, days=0, hours=0, minutes=0, seconds=0
  integer :: yr, mon, day, hr, min, sec
  integer :: is, ie, js, je, isd, ied, jsd, jed, IscB, IecB, JscB, JecB, IsdB, IedB, JsdB, JedB
  integer :: ierr, io_status, n, nz, unit, i
  integer :: X_FLAGS, Y_FLAGS
  integer :: profile_fid
  character(len=128) :: mesg
  !logical, dimension(:), pointer :: ens_pe(:)=>NULL() !< True if current PE is associated with ensemble member(n)

  namelist /test_driver_nml/ date_init, date, calendar, get_increment

  call MOM_infra_init()
  call io_infra_init()
  call time_manager_init()
  call ensemble_manager_init()
! Read namelist
  date(:)=0
  calendar_type=-1
  call open_file(unit,'input.nml',form=ASCII_FILE,action=READONLY_FILE)
  read(unit,test_driver_nml,iostat=io_status)
  call close_file(unit)
  ierr = check_nml_error(io_status,'test_driver_nml')
  ens_info = get_ensemble_size()
  ens_size = ens_info(1); npes_ocn = ens_info(2)
  npes_atm=npes_ocn; npes_lnd=npes_ocn; npes_ice=npes_ocn
  ! All components for an individual ensemble member will run on an identical PElist
  allocate(ocn_pelist(npes_ocn))
  allocate(atm_pelist(npes_ocn))
  allocate(lnd_pelist(npes_ocn))
  allocate(ice_pelist(npes_ocn))
  allocate(ens_pelist(ens_size,npes_ocn))
  call ensemble_pelist_setup(.false.,npes_atm,npes_ocn,npes_lnd,npes_ice,&
          atm_pelist,ocn_pelist,lnd_pelist,ice_pelist)
  call get_ensemble_pelist(ens_pelist)
  ensemble_id = get_ensemble_id()
  write(mesg,'(a,i2.2,a,i2)') 'Number of ensemble members= ',ens_size,': Number of PEs per member= ',npes_ocn

  if (file_exists(trim(dirs%restart_input_dir)//'coupler.res')) then
    call open_file(unit,trim(dirs%restart_input_dir)//'coupler.res',form=ASCII_FILE,action=READONLY_FILE)
    read(unit,*) calendar_type; read(unit,*) date_init; read(unit,*) date
    call close_file(unit)
  else
    ! Otherwise, if a coupler.res file does not exist, then use
    ! information provided by test_driver_nml
    if (calendar(1:6)=='julian') calendar_type=JULIAN
    if (calendar(1:11)=='no_calendar') calendar_type=NO_CALENDAR
    if (calendar(1:11)=='noleap') calendar_type=NOLEAP
    if (calendar_type<0) call MOM_error(FATAL,'paricles_driver: Invalid namelist value for calendar')
  endif
  ! Set the calendar type and initial and current dates
  call set_calendar_type(calendar_type)
  if (sum(date).eq.0) date=date_init
  call diag_manager_init(TIME_INIT=date) ! initialize diag_manager
  !----- always override initial/base date with diag_manager value -----
  call get_base_date(date_init(1),date_init(2),date_init(3),date_init(4),&
       date_init(5),date_init(6))
!----- set initial and current time types ------
  Time_start = set_date (date_init(1), date_init(2), date_init(3), &
       date_init(4), date_init(5), date_init(6))
  Time=set_date(date(1),date(2),date(3),date(4),date(5),date(6))
  Time_start_segment=Time

! Allocate storage for ensemble members.
  call set_current_pelist(ens_pelist(ensemble_id,:))
  call set_root_pe(ens_pelist(ensemble_id,1))
  allocate(PF)
  call get_MOM_input(PF,dirs,ensemble_num=ensemble_id) !< read PF for this member
  allocate(Grid)
  call MOM_domains_init(Grid%domain,PF,symmetric=.false.) !< initialize storage for horizontal and vertical grids and decompisition
  allocate(HI)
  call hor_index_init(Grid%Domain, HI, PF, &
       local_indexing=.false.)
  allocate(CS)
  call verticalGridInit( PF, CS%GV )
  call create_dyn_horgrid(dG,HI)
  call clone_MOM_domain(Grid%Domain, dG%Domain,symmetric=.false.)
  is   = dG%isc   ; ie   = dG%iec  ; js   = dG%jsc  ; je   = dG%jec
  isd  = dG%isd   ; ied  = dG%ied  ; jsd  = dG%jsd  ; jed  = dG%jed
  IsdB = dG%IsdB  ; IedB = dG%IedB ; JsdB = dG%JsdB ; JedB = dG%JedB
  nz = CS%GV%ke
! Allocate and initialize storage fo MOM state vector
  allocate(CS%u(IsdB:IedB,jsd:jed,nz))   ; CS%u(:,:,:) = 0.0
  allocate(CS%v(isd:ied,JsdB:JedB,nz))   ; CS%v(:,:,:) = 0.0
  allocate(CS%h(isd:ied,jsd:jed,nz))     ; CS%h(:,:,:) = CS%GV%Angstrom
  allocate(CS%T(isd:ied,jsd:jed,nz))   ; CS%T(:,:,:) = 0.0
  allocate(CS%S(isd:ied,jsd:jed,nz))   ; CS%S(:,:,:) = 0.0
  CS%tv%T => CS%T ; CS%tv%S => CS%S
! Initialize grids and vertical coordinate
!    call MOM_initialize_fixed(dG,CS%OBC,PF,.false.,dirs%output_directory)
  call set_grid_metrics(dG,PF)
  call MOM_initialize_topography(dg%bathyT,dG%max_depth,dG,PF)
  call MOM_initialize_coord(CS%GV, PF, .false., &
       dirs%output_directory, CS%tv, dG%max_depth)
  call MOM_grid_init(Grid, PF, HI,global_indexing=.true.)
  call copy_dyngrid_to_MOM_grid(dG, Grid)
  Grid%ke = CS%GV%ke; Grid%g_Earth = CS%GV%g_Earth
! Initialize restart info (h,T and S are required)
  call restart_init(PF,CS%restart_CSp)
  vd = var_desc("h",'m',"Layer Thickness")
  call register_restart_field(CS%h, vd, .true., CS%restart_CSp)
  vd = var_desc("u","m s-1","Zonal velocity",'u','L')
  call register_restart_field(CS%u, vd, .false., CS%restart_CSp)
  vd = var_desc("v","m s-1","Meridional velocity",'v','L')
  call register_restart_field(CS%v, vd, .false., CS%restart_CSp)
  vd = var_desc("Temp","degC","Potential Temperature")
  call register_restart_field(CS%tv%T, vd, .true., CS%restart_CSp)
  vd = var_desc("Salt","PPT","Salinity")
  call register_restart_field(CS%tv%S, vd, .true., CS%restart_CSp)
!! Initialize MOM state vector override changes to intent(inout) Time
!! argument.
  call MOM_initialize_state(CS%u, CS%v, CS%h, CS%tv, Time, Grid, CS%GV, PF, &
       dirs, CS%restart_CSp, CS%ALE_CSp, CS%tracer_Reg, &
       CS%sponge_CSp, CS%ALE_sponge_CSp, CS%OBC,Time_in=Time_start_segment)
!! Checksum the model state vector
  write(mesg,*) 'Initial thermo state for ensemble member ',ensemble_id
  call MOM_thermo_chksum(mesg,CS%tv,Grid)
! Close PF for this member
  call close_param_file(PF)
!  call set_current_pelist()
  call get_time(Time,sec,days)
  call init_oda(Time, Grid, CS%GV, odaCS )
  call set_analysis_time(Time,odaCS) !< Set analysis time
  call set_prior_tracer(Time, Grid, CS%GV, CS%h, CS%tv, odaCS)
  call oda(Time,odaCS) !<read observations and calculate ensemble increments or posterior
  call save_obs_diff('first_guess_profiles.nc' , odaCS)
  call get_posterior_tracer(Time, odaCS, Grid, GV, h, tv, increment=get_increment)
  !print *,'00001x',ens_pelist(ensemble_id,:)
!  call set_current_pelist(ens_pelist(ensemble_id,:))
!  call set_root_pe(ens_pelist(ensemble_id,1))
  write(mesg,*) 'Posterior thermo state for ensemble member after remapping back to model domain', ensemble_id
  call MOM_thermo_chksum(mesg,tv,Grid)
  call diag_manager_end(Time) ! close diag_manager
  deallocate(odaCS)
  deallocate(tv,atm_pelist,lnd_pelist,ice_pelist,ens_pelist)

  call fms_end()

  stop
END PROGRAM MOM6_ODA_offline_driver
