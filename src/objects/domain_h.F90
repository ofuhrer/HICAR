module domain_interface
  use, intrinsic :: iso_c_binding, only : C_DOUBLE_COMPLEX
  use mpi
  use options_interface,        only : options_t
  use boundary_interface,       only : boundary_t
  use grid_interface,           only : grid_t
  use variable_interface,       only : variable_t
  use data_structures,          only : interpolable_type, tendencies_type, index_type
  use halo_interface,           only : halo_t
  use timer_interface,          only : timer_t
  use flow_object_interface,    only : flow_obj_t
  use sparse_lbc_reader,        only : sparse_lbc_reader_t
  use icar_constants,               only : kMAX_STORAGE_VARS, kVARS, kMAX_NAME_LENGTH, MAXLEVELS
  implicit none

  private
  public :: domain_t, auto_dz

  type , extends(flow_obj_t) :: domain_t
    type(grid_t)         :: grid,   grid8w,  u_grid,   v_grid
    type(grid_t)         :: column_grid
    type(grid_t)         :: global_grid_2d, global_grid, global_grid8w
    type(grid_t)         :: neighbor_grid_2d, neighbor_grid, neighbor_grid8w
    type(grid_t)         :: grid2d, u_grid2d, v_grid2d
    type(grid_t)         :: grid_monthly, grid_soil
    type(grid_t)         :: grid_snow, grid_snow_i, grid_snowsoil, grid_fm
    type(grid_t)         :: grid_wind_height
    type(grid_t)         :: grid_soilcomp, grid_gecros, grid_croptype
    type(grid_t)         :: grid_hlm, grid_Sx !! MJ added
    type(grid_t)         :: grid_lake , grid_lake_soisno, grid_lake_soi, grid_lake_soisno_1
    type(halo_t)         :: halo

    ! note that not all variables are allocated at runtime, physics packages must request a variable be created
    ! though variables considered "required" are requested by the domain object itself (e.g. terrain)
    ! core model species to be advected

    ! these data are stored on the domain wide grid even if this process is only looking at a subgrid
    ! these variables are necessary with linear winds, especially with spatially variable dz, to compute the LUT

    type(tendencies_type) :: tend

    type(index_type) :: vars_to_out(kMAX_STORAGE_VARS)
    
    ! Array listing variables to advect with pointers to local data
    type(index_type), allocatable :: adv_vars(:), exch_vars(:)
    integer :: n_adv_2d, n_adv_3d, n_exch_2d, n_exch_3d
    
    type(interpolable_type) :: geo
    type(interpolable_type) :: geo_agl
    type(interpolable_type) :: geo_u
    type(interpolable_type) :: geo_v

    real :: smooth_height, dx
    integer :: nsmooth

    ! complex(C_DOUBLE_COMPLEX),  allocatable :: terrain_frequency(:,:) ! FFT(terrain)

    type(variable_t), allocatable :: forcing_hi(:)
    ! forcing_hi%data_*d and forcing_hi%dqdt_*d retain the exact
    ! interpolated left and right forcing records, respectively.  This flag
    ! distinguishes the first right-endpoint load from later interval
    ! advances, when the previous right endpoint must be promoted to the
    ! next interval's left endpoint.
    logical :: forcing_interval_ready = .False.
    type(sparse_lbc_reader_t) :: sparse_lbc

    ! Compact, restart-persistent definition of the static theta_bar(z)
    ! reference used by split-potential-temperature advection.  The expanded
    ! 3-D field is reconstructed from these two one-dimensional tables after
    ! every process start.  Without this state, a restart rebuilds theta_bar
    ! from the later atmospheric state and changes the first advection step.
    integer :: adv_theta_ref_n = 0
    real :: adv_theta_ref_z(MAXLEVELS) = 0.0
    real :: adv_theta_ref_theta(MAXLEVELS) = 0.0

    type(variable_t), allocatable :: vars_1d(:)
    type(variable_t), allocatable :: vars_2d(:)
    type(variable_t), allocatable :: vars_3d(:)
    type(variable_t), allocatable :: vars_4d(:)

    type(index_type) :: var_indx(kMAX_STORAGE_VARS), forcing_var_indx(kMAX_STORAGE_VARS)
    
    ! MPI communicator object for doing parallel communications among domain objects
    integer, public :: compute_comms

    ! timers used to track the time spent doing various operations
    type(timer_t) :: initialization_timer, total_timer, input_timer, &
                        output_timer, physics_timer, wind_timer, mp_timer, &
                        adv_timer, rad_timer, lsm_timer, pbl_timer, exch_timer, &
                        send_timer, ret_timer, wait_timer, forcing_timer, diagnostic_timer, wind_bal_timer, &
                        flux_timer, flux_corr_timer, sum_timer, adv_wind_timer, cpu_gpu_timer, nest_timer

    ! contains the size of the domain (or the local tile?)
    integer :: nx, ny, nz, nx_global, ny_global
    integer :: ximg, ximages, yimg, yimages
    logical :: north_boundary = .True.
    logical :: south_boundary = .True.
    logical :: east_boundary = .True.
    logical :: west_boundary = .True.
    integer :: FILTER_WIDTH = 7

    ! Map-scale factors m = (nominal dx) / (true ground distance), computed
    ! at init from the hi-res lat/lon fields (init_map_factors). All 1.0
    ! when use_map_factors is off or the lat/lon fields are degenerate
    ! (idealized grids) — kernels multiply unconditionally; x*1.0 is exact.
    ! _u arrays live on the u grid, _v on the v grid; mapfac_mxy = m_x*m_y
    ! at mass points (the cell-area factor). max_mapfac feeds a
    ! conservative CFL correction in compute_dt.
    real, allocatable :: mapfac_mx_u(:,:), mapfac_my_u(:,:)
    real, allocatable :: mapfac_mx_v(:,:), mapfac_my_v(:,:)
    real, allocatable :: mapfac_mxy(:,:)
    real :: max_mapfac = 1.0

    ! The hourly climate fields themselves hold running time integrals between
    ! output events. Only the current ten-minute scalar-speed integrals need
    ! separate scratch storage. Checkpoints are constrained to hourly output
    ! boundaries, so no partial accumulator state crosses a restart.
    logical :: wind_climatology_enabled = .False.
    real :: wind_climatology_hour_seconds = 0.0
    real :: wind_climatology_tenminute_seconds = 0.0
    real, allocatable :: wind_speed_tenminute_sum_agl(:,:,:)
    real, allocatable :: wind_speed_tenminute_sum_10m(:,:)

    ! store the start (s) and end (e) for the i,j,k dimensions
    integer ::  ids,ide, jds,jde, kds,kde, & ! for the entire model domain    (d)
                ims,ime, jms,jme, kms,kme, & ! for the memory in these arrays (m)
                its,ite, jts,jte, kts,kte, & ! for the data tile to process   (t)
                ihs,ihe, jhs,jhe, khs,khe    ! for the neighborhood arrays for non-local calculations (h)

    integer :: neighborhood_max ! The maximum neighborhood radius in indices
    

  contains
    procedure :: init => init_domain
    procedure :: release
    procedure :: enforce_limits

    procedure :: batch_exch
    procedure :: halo_3d_send
    procedure :: halo_3d_retrieve
    procedure :: halo_2d_send
    procedure :: halo_2d_retrieve

    procedure :: get_initial_conditions
    procedure :: diagnostic_update
    procedure :: update_wind_height_diagnostics
    procedure :: accumulate_wind_climatology
    procedure :: finalize_wind_climatology
    procedure :: reset_wind_climatology
    procedure :: interpolate_forcing
    procedure :: update_delta_fields
    procedure :: apply_forcing
    procedure :: initialize_sparse_lbc
    procedure :: synchronize_sparse_lbc
    procedure :: apply_sparse_lbc
    procedure :: forcing_phase_at
    procedure :: read_land_variables

    procedure :: update_device
    procedure :: update_host

  end type

  integer, parameter :: space_dimension=3

  interface

    ! Set default component values
    module subroutine init_domain(this, options, nest_indx)
        implicit none
        class(domain_t), intent(inout) :: this
        type(options_t), intent(inout) :: options
        integer,         intent(in)    :: nest_indx
    end subroutine init_domain

    ! finalize domain object, freeing halo mpi windows
    module subroutine release(this)
        implicit none
        class(domain_t), intent(inout) :: this
    end subroutine

    module subroutine batch_exch(this, two_d, exch_only)
        implicit none
        class(domain_t), intent(inout) :: this
        logical, optional,   intent(in) :: two_d, exch_only
  end subroutine

    module subroutine halo_3d_send(this, exch_only)
        implicit none
        class(domain_t), intent(inout) :: this
        logical, optional,   intent(in) :: exch_only
    end subroutine
    
    module subroutine halo_3d_retrieve(this, exch_only)
      implicit none
      class(domain_t), intent(inout) :: this
      logical, optional,   intent(in) :: exch_only
  end subroutine

  module subroutine halo_2d_send(this)
    implicit none
    class(domain_t), intent(inout) :: this
  end subroutine

  module subroutine halo_2d_retrieve(this)
    implicit none
    class(domain_t), intent(inout) :: this
  end subroutine

    ! read initial atmospheric conditions from forcing data
    module subroutine get_initial_conditions(this, forcing, options)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(boundary_t), intent(inout) :: forcing
        type(options_t),  intent(in)    :: options
    end subroutine

    module subroutine diagnostic_update(this, forcing_update, thermo_only)
      implicit none
      class(domain_t),  intent(inout)   :: this
      logical, intent(in), optional    :: forcing_update
      logical, intent(in), optional    :: thermo_only
    end subroutine

    module subroutine update_wind_height_diagnostics(this)
      implicit none
      class(domain_t), intent(inout) :: this
    end subroutine

    module subroutine accumulate_wind_climatology(this, dt_seconds)
      implicit none
      class(domain_t), intent(inout) :: this
      real, intent(in) :: dt_seconds
    end subroutine

    module subroutine finalize_wind_climatology(this)
      implicit none
      class(domain_t), intent(inout) :: this
    end subroutine

    module subroutine reset_wind_climatology(this)
      implicit none
      class(domain_t), intent(inout) :: this
    end subroutine

    module subroutine interpolate_forcing(this, forcing, update)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(boundary_t), intent(inout) :: forcing
        logical,          intent(in),   optional :: update
    end subroutine

    ! Exchange subdomain boundary information
    !module subroutine halo_exchange(this, two_d, exch_only)
    !    implicit none
    !    class(domain_t), intent(inout) :: this
    !    logical, optional,   intent(in) :: two_d, exch_only
    !end subroutine

    ! Make sure no hydrometeors are getting below 0
    module subroutine enforce_limits(this,update_in)
        implicit none
        class(domain_t), intent(inout) :: this
        logical, optional, intent(in)  :: update_in
    end subroutine

    module subroutine update_delta_fields(this)
        implicit none
        class(domain_t),    intent(inout) :: this
    end subroutine

    module subroutine apply_forcing(this, options, dt)
        implicit none
        class(domain_t),    intent(inout) :: this
        type(options_t), intent(in)       :: options
        real, intent(in)                  :: dt
    end subroutine

    module subroutine initialize_sparse_lbc(this, options)
        implicit none
        class(domain_t), intent(inout) :: this
        type(options_t), intent(in) :: options
    end subroutine initialize_sparse_lbc

    module subroutine synchronize_sparse_lbc(this)
        implicit none
        class(domain_t), intent(inout) :: this
    end subroutine synchronize_sparse_lbc

    module subroutine apply_sparse_lbc(this, dt)
        implicit none
        class(domain_t), intent(inout) :: this
        real, intent(in) :: dt
    end subroutine apply_sparse_lbc

    module function forcing_phase_at(this, offset_seconds) result(phase)
        implicit none
        class(domain_t), intent(in) :: this
        real, intent(in)            :: offset_seconds
        real                        :: phase
    end function

    module subroutine read_land_variables(this, options)
        implicit none
        class(domain_t), intent(inout) :: this
        type(options_t), intent(in)    :: options
    end subroutine

    module subroutine update_device(this)
        implicit none
        class(domain_t), intent(inout) :: this
    end subroutine

    module subroutine update_host(this)
        implicit none
        class(domain_t), intent(inout) :: this
    end subroutine

    module subroutine auto_dz(options)
        implicit none
        type(options_t), intent(inout) :: options
    end subroutine

  end interface

end module
