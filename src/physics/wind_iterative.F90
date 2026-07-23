!>------------------------------------------------------------
!! Native HICAR iterative wind solver.
!!
!! Right-preconditioned BiCGStab + vertical-line block-Jacobi smoothing
!! over a 15-point stencil with terrain-following boundary
!! conditions. Pure OpenACC + MPI — vendor-portable, no external
!! solver library.
!!
!! Phase 1 (this file): classic BiCGStab, MPI_Allreduce on dot
!! products, MPI Sendrecv halos, pure-stencil SpMV (no CSR).
!! Multi-nest cache and pipelined/NCCL paths land in later phases.
!!
!! Matrix structure mirrors wind_iterative_amgx.F90 exactly so the
!! two solvers are interchangeable for validation.
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!------------------------------------------------------------

module wind_iterative
    use iso_c_binding
    use domain_interface,  only : domain_t
    use icar_constants,    only : STD_OUT_PE, kVARS, kITERATIVE_WINDS
    use options_interface, only : options_t
    use iso_fortran_env
    use mpi
    use openacc
    use string,        only : str
    use debug_module,  only : domain_check_winds
    use wind_hypre,    only : wind_hypre_available, wind_hypre_invalidate, wind_hypre_solve, wind_hypre_apply
    use wind_multilevel, only : horizontal_tile_transfer_t, galerkin_tile_stencil_t, &
                                vertical_line_factor_t, assemble_colored_tile_galerkin, &
                                owned_coarse_interval
    use wind_multilevel_mpi, only : horizontal_halo_exchange_t
#ifdef USE_NCCL
    use nccl_interface, only : nccl_comm_init, nccl_comm_destroy, &
                               nccl_group_start, nccl_group_end, &
                               nccl_send_double, nccl_recv_double, &
                               nccl_allreduce_double_sum
#endif

    implicit none
    private
    public :: init_iter_winds, calc_iter_winds, finalize_iter_winds
    public :: probe_lambda_pattern, probe_zero_corrections, probe_apply_corrections, &
              probe_record, probe_finalize, probe_random_pattern, probe_compare_operator
    public :: multilevel_preconditioner_is_ready
    public :: multilevel_preconditioner_smoke

    logical :: initialized_iter_winds = .false.
    logical :: structure_uploaded = .false.

    ! the stencil coefficients have
    ! been replaced by the exact composition of the model's divergence
    ! operator with the velocity-correction operator (obtained by lattice
    ! probing — see wind.F90::calibrate_projection_operator). spmv/extract_diagonal then
    ! use the probed 15-point operator on all physical rows and identity
    ! rows on the ghost planes (the correction operator no longer reads
    ! ghost lambda under w_to_grid).
    logical :: operator_probed = .false.
    ! Emit one post-calibration backend-selection report.  This makes a
    ! fallback before HYPRE is entered visible in the model stdout without
    ! adding a line for every subsequent wind solve.
    logical :: hypre_selection_reported = .false.
    ! An assembled external matrix is immutable after the calibration probe.
    ! Check its action once per probe, then avoid an extra distributed SpMV
    ! and a line of stdout for every physical time step.
    logical :: hypre_operator_verified = .false.
    ! Opt-in structural diagnostics for the calibrated projection operator.
    ! The audit is enabled only when HICAR_WIND_OPERATOR_AUDIT is set to a
    ! non-zero value, and runs once per calibrated nest.  It never changes
    ! the solver tolerance or acceptance path.
    logical :: operator_audit_enabled = .false.
    logical :: operator_audit_done = .false.
    logical :: krylov_audit_written = .false.
    logical :: bootstrap_rhs_saved = .false.
    logical :: multilevel_requested = .false.
    logical :: multilevel_setup_attempted = .false.
    logical :: multilevel_ready = .false.
    logical :: multilevel_arrays_uploaded = .false.
    integer :: operator_audit_max_iters = 0
    character(len=512) :: operator_audit_file = 'hicar_wind_operator_audit.csv'


    real, parameter :: deg2rad = 0.017453293
    real, parameter :: rad2deg = 57.2957779371

    integer :: wind_solver_max_iters = 1500
    real(c_double), parameter :: bicg_tol_abs = 1.0e-10_c_double
    real(c_double), parameter :: bicg_tol_rel = 1.0e-5_c_double
    real(c_double), parameter :: breakdown_eps = 1.0e-30_c_double

    ! Number of Richardson sweeps per vertical-line block-Jacobi apply.  Each
    ! block is an exact Thomas solve through one terrain-following column, so
    ! the dominant fine-vertical-grid coupling is removed before Krylov sees it.
    !
    ! Adaptive retry: if a solve diverges/stagnates, calc_iter_winds bumps
    ! this up to MAX_PREC_SWEEPS, retrying after each bump. On successful
    ! convergence we reset to BASE_PREC_SWEEPS for the next call. Mirrors the
    ! AMGX module's prec_max_iters retry. Per-nest state is cached.
    integer, parameter :: BASE_PREC_SWEEPS = 2
    integer, parameter :: MAX_PREC_SWEEPS  = 4
    integer :: precond_n_sweeps = BASE_PREC_SWEEPS

    ! Per-solve status/residual/timing printing. Off by default under
    ! RANS (one solve per physics step makes it far too verbose); on in
    ! debug mode or under the diagnostic iterative solver (solves only at
    ! wind updates). Convergence-failure/retry warnings always print.
    logical :: verbose_solver = .False.
    ! Emit a rank-0 heartbeat while a verbose solve is running.  Swiss-scale
    ! solves can legitimately take minutes before the final summary, so the
    ! heartbeat makes a live, converging run distinguishable from a hang.
    integer, parameter :: SOLVER_PROGRESS_INTERVAL = 100
    ! Keep the native flexible-GMRES workspace bounded.  Unlike the external
    ! host FGMRES path this is an explicit GPU allocation, so the national
    ! decomposition has a predictable (2*(restart+1)) vector footprint.
    ! Before the exact recursive hierarchy, 12- and 50-vector restarts lost too
    ! much Krylov information, while restart 100 still failed and exhausted
    ! national-scale memory.  The verified hierarchy now converges the 250 m
    ! regression and the hard 701x701 regional bridge in 5--6 iterations, so a
    ! bounded 20-vector space leaves ample margin without allocating roughly
    ! 60 GiB of basis vectors per Swiss-domain compute rank.
    integer, parameter :: FGMRES_RESTART = 20
    ! Coarse global sketch used only by the opt-in Arnoldi audit.  Signed
    ! basis sums on this fixed grid let the offline analyzer reconstruct the
    ! spatial envelope of harmonic Ritz vectors without exporting full 3-D
    ! Krylov fields.
    integer, parameter :: AUDIT_X_BINS = 16
    integer, parameter :: AUDIT_Y_BINS = 12
    integer, parameter :: AUDIT_Z_BINS = 8

    ! 15-point stencil coefficients (same names as AMGX module)
    real, allocatable, dimension(:,:,:) :: A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, &
                                           H_coef, I_coef, J_coef, K_coef, L_coef, M_coef, N_coef, O_coef
    real    :: dx
    real, allocatable, dimension(:,:,:) :: div, dz_if, jaco, dzdx, dzdy, sigma, alpha
    real, allocatable, dimension(:,:)   :: dzdx_surf, dzdy_surf
    real, allocatable, dimension(:,:,:) :: jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag
    real, allocatable, dimension(:)     :: adv_dz_col

    ! Domain decomposition / boundary-cell scalars (PETSc/AMGX convention)
    integer :: hs, i_s, i_e, k_s, k_e, j_s, j_e
    integer :: ims, ime, jms, jme, ids, ide, jds, jde
    integer :: xs, ys, zs, xm, ym, zm, mx, my, mz
    integer :: n_rows, n_rows_global

    ! Krylov / RHS / preconditioner state — 3D with halos so SpMV reads neighbours directly.
    ! Allocation is (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1) which covers:
    !   - the owned cell range xs..xs+xm-1 (incl. BC layer at boundary ranks),
    !   - one lateral halo cell in i and j on non-boundary edges.
    !
    ! Classic right-preconditioned BiCGStab vectors. Earlier Phase-2 versions tracked
    ! u = M^{-1} r and q = M^{-1} v via recurrences to skip fresh preconditioner applies,
    ! but the recurrences accumulated roundoff and could stall convergence at high iter
    ! counts. Reverted to fresh M^{-1} on p and s. Merged 5-value 2nd allreduce is kept
    ! so we still pay only 2 allreduces per iter.
    real(c_double), allocatable, dimension(:,:,:) :: x_sol      ! solution
    real(c_double), allocatable, dimension(:,:,:) :: r_vec      ! residual
    real(c_double), allocatable, dimension(:,:,:) :: r_hat      ! shadow residual r̂_0
    real(c_double), allocatable, dimension(:,:,:) :: p_vec, p_hat
    real(c_double), allocatable, dimension(:,:,:) :: v_vec
    real(c_double), allocatable, dimension(:,:,:) :: s_vec, s_hat
    real(c_double), allocatable, dimension(:,:,:) :: t_vec
    real(c_double), allocatable, dimension(:,:,:) :: rhs        ! right-hand side b
    real(c_double), allocatable, dimension(:,:,:) :: D_inv      ! inverse Thomas pivots for vertical-line blocks
    real(c_double), allocatable, dimension(:,:,:) :: line_cprime ! Thomas upper factors for vertical-line blocks
    real(c_double), allocatable, dimension(:,:,:) :: prec_res   ! scratch for multi-sweep Richardson residual

    ! Persistent halo face buffers — allocated once at first solve, reused every iter.
    ! 8 buffers (4 directions x send/recv) so all 4 sends and all 4 recvs can be in flight
    ! simultaneously via Isend/Irecv + Waitall. Cuts halo latency from 4*alpha (sequential
    ! Sendrecvs) to ~1*alpha — major win at high MPI rank counts.
    real(c_double), allocatable, target, dimension(:,:) :: east_send,  east_recv     ! shape (k_s-1:k_e+1, j_s-1:j_e+1)
    real(c_double), allocatable, target, dimension(:,:) :: west_send,  west_recv     ! shape (k_s-1:k_e+1, j_s-1:j_e+1)
    real(c_double), allocatable, target, dimension(:,:) :: north_send, north_recv    ! shape (i_s-1:i_e+1, k_s-1:k_e+1)
    real(c_double), allocatable, target, dimension(:,:) :: south_send, south_recv    ! shape (i_s-1:i_e+1, k_s-1:k_e+1)

    integer :: solver_comm
    integer :: solver_rank = -1
    integer :: east_neighbor = -1, west_neighbor = -1, north_neighbor = -1, south_neighbor = -1

    type :: multilevel_deep_level_t
        type(horizontal_tile_transfer_t) :: transfer_from_parent
        type(galerkin_tile_stencil_t) :: stencil
        type(vertical_line_factor_t) :: line_factor
        type(horizontal_halo_exchange_t) :: halo
        logical :: arrays_uploaded = .false.
        real(c_double), allocatable :: weight(:,:,:), weight_halo(:,:,:)
        real(c_double), allocatable :: halo_x(:,:,:)
        real(c_double), allocatable :: b(:,:,:), x(:,:,:), ax(:,:,:)
        real(c_double), allocatable :: r(:,:,:), correction(:,:,:)
    end type multilevel_deep_level_t

    integer, parameter :: MAX_ML_DEEP_LEVELS = 12
    integer :: ml_deep_count = 0
    integer :: ml_assembly_child = 0
    type(multilevel_deep_level_t), allocatable :: ml_deep(:)

    ! Feature-gated, distributed Petrov-Galerkin V-cycle.  These arrays use
    ! uniform local indexing rather than aliasing the solver vectors, whose
    ! physical-edge bounds differ by rank.  The separation is deliberate:
    ! it keeps transfer and coarse-grid halo ownership decomposition-exact.
    type(horizontal_tile_transfer_t) :: ml_transfer
    type(galerkin_tile_stencil_t) :: ml_stencil
    type(vertical_line_factor_t) :: ml_line_factor
    type(horizontal_halo_exchange_t) :: ml_fine_halo, ml_coarse_halo
    real(c_double), allocatable :: ml_fine_owned(:,:,:)
    real(c_double), allocatable :: ml_fine_residual(:,:,:), ml_fine_weight(:,:,:)
    real(c_double), allocatable :: ml_coarse_halo_x(:,:,:)
    real(c_double), allocatable :: ml_coarse_weight(:,:,:), ml_coarse_weight_halo(:,:,:)
    real(c_double), allocatable :: ml_coarse_b(:,:,:)
    real(c_double), allocatable :: ml_coarse_x(:,:,:), ml_coarse_ax(:,:,:)
    real(c_double), allocatable :: ml_coarse_r(:,:,:), ml_coarse_correction(:,:,:)

    ! Multi-nest cache: each nest's full state lives in its own slot. On nest
    ! context switch (nest_manager.F90:switch_nest_context), init_iter_winds
    ! saves the active state and restores the target. O(1) switching once each
    ! nest has been visited once. Mirrors the AMGX module's cache (which is
    ! itself based on the long-standing physics module pattern in HICAR).
    integer, parameter :: MAX_NESTS = 4
    integer :: active_nest_indx = -1

    type :: hicar_cache_t
        logical :: valid = .false.
        ! Grid scalars
        integer :: i_s, i_e, k_s, k_e, j_s, j_e
        integer :: ims, ime, jms, jme, ids, ide, jds, jde
        integer :: xs, ys, zs, xm, ym, zm, mx, my, mz
        integer :: hs, n_rows, n_rows_global
        real    :: dx
        integer :: solver_rank, east_neighbor, west_neighbor, north_neighbor, south_neighbor
        integer :: solver_comm
        ! Flags
        logical :: structure_uploaded = .false.   ! default-init: read at restore_from_cache
        logical :: operator_probed = .false.
        logical :: operator_audit_enabled = .false.
        logical :: operator_audit_done = .false.
        logical :: krylov_audit_written = .false.
        logical :: bootstrap_rhs_saved = .false.
        logical :: multilevel_requested = .false.
        integer :: wind_solver_max_iters
        integer :: precond_n_sweeps = BASE_PREC_SWEEPS
        ! Geometry (single precision, owned-range allocations)
        real, allocatable, dimension(:,:,:) :: A_coef, B_coef, C_coef, D_coef, E_coef
        real, allocatable, dimension(:,:,:) :: F_coef, G_coef, H_coef, I_coef, J_coef
        real, allocatable, dimension(:,:,:) :: K_coef, L_coef, M_coef, N_coef, O_coef
        real, allocatable, dimension(:,:,:) :: div, dz_if, jaco, dzdx, dzdy, sigma, alpha
        real, allocatable, dimension(:,:)   :: dzdx_surf, dzdy_surf
        real, allocatable, dimension(:,:,:) :: jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag
        real, allocatable, dimension(:)     :: adv_dz_col
        ! Krylov / RHS / precond (double precision, halo'd allocations)
        real(c_double), allocatable, dimension(:,:,:) :: x_sol, r_vec, r_hat
        real(c_double), allocatable, dimension(:,:,:) :: p_vec, p_hat, v_vec
        real(c_double), allocatable, dimension(:,:,:) :: s_vec, s_hat, t_vec
        real(c_double), allocatable, dimension(:,:,:) :: rhs, D_inv, line_cprime, prec_res
        ! Halo face buffers (target attribute only needed on the live module-level
        ! variables for c_loc — not on cache slots, which are pure host storage)
        real(c_double), allocatable, dimension(:,:) :: east_send,  east_recv
        real(c_double), allocatable, dimension(:,:) :: west_send,  west_recv
        real(c_double), allocatable, dimension(:,:) :: north_send, north_recv
        real(c_double), allocatable, dimension(:,:) :: south_send, south_recv
#ifdef USE_NCCL
        real(c_double), allocatable, dimension(:) :: sigma_dev, red5_dev, rho0_dev
#endif
    end type hicar_cache_t

    type(hicar_cache_t) :: domain_cache(MAX_NESTS)

#ifdef USE_NCCL
    ! NCCL communicator + stream. Stream is aliased to OpenACC's sync queue so
    ! all NCCL ops are stream-ordered with the OpenACC pack/unpack/reduction kernels
    ! that produce inputs and consume outputs — no explicit cross-stream sync needed.
    type(c_ptr) :: nccl_comm   = c_null_ptr
    type(c_ptr) :: nccl_stream = c_null_ptr
    logical     :: nccl_initialized = .false.
    ! Device-resident allreduce buffers (one for sigma, one for the 5-value reduction).
    ! Allocated once at first solve so NCCL operates entirely on device pointers.
    real(c_double), allocatable, target, dimension(:) :: sigma_dev   ! length 1
    real(c_double), allocatable, target, dimension(:) :: red5_dev    ! length 5
    real(c_double), allocatable, target, dimension(:) :: rho0_dev    ! length 2 (rho_0 + ||r0||²)
#endif

    ! Per-solve timing accumulators (set to 0 at start of each bicgstab_solve, summary
    ! printed at end). MPI_Wtime is portable across CPU and GPU paths. With synchronous
    ! OpenACC kernels, MPI_Wtime around a kernel call measures the kernel's wall time.
    real(c_double) :: t_spmv_acc      = 0.0_c_double
    real(c_double) :: t_halo_acc      = 0.0_c_double
    real(c_double) :: t_precond_acc   = 0.0_c_double
    real(c_double) :: t_allreduce_acc = 0.0_c_double
    real(c_double) :: t_vecops_acc    = 0.0_c_double
    real(c_double) :: t_total_acc     = 0.0_c_double

contains

    logical function multilevel_preconditioner_is_ready()
        multilevel_preconditioner_is_ready = multilevel_ready
    end function multilevel_preconditioner_is_ready


    logical function multilevel_preconditioner_smoke(domain)
        type(domain_t), intent(in) :: domain
        real(c_double) :: local_norm, global_norm, local_rhs_norm, global_rhs_norm
        real(c_double) :: local_residual_norm, global_residual_norm
        integer :: i, j, k, ierr

        multilevel_preconditioner_smoke = .false.
        if (.not. multilevel_requested .or. .not. operator_probed .or. .not. structure_uploaded) return
        call build_line_preconditioner()
        if (.not. multilevel_setup_attempted) call setup_multilevel_preconditioner(domain)
        if (.not. multilevel_ready) return
        !$acc parallel loop gang vector collapse(3) present(rhs)
        do j = ys, ys+ym-1
            do k = zs, zs+zm-1
                do i = xs, xs+xm-1
                    if (i <= 0 .or. i >= mx-1 .or. j <= 0 .or. j >= my-1 .or. &
                        k <= 0 .or. k >= mz-1) then
                        rhs(i,k,j) = 0.0_c_double
                    else
                        rhs(i,k,j) = sin(0.071_c_double*real(3*i+5*k+7*j,c_double))
                    endif
                enddo
            enddo
        enddo
        call apply_multilevel_preconditioner(rhs, x_sol, domain)
        call vec_norm2_local(x_sol, local_norm)
        call MPI_Allreduce(local_norm, global_norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call exchange_krylov_halos(x_sol, domain)
        call spmv(x_sol, t_vec)
        !$acc parallel loop gang vector collapse(3) present(rhs,t_vec,prec_res)
        do j = ys, ys+ym-1
            do k = zs, zs+zm-1
                do i = xs, xs+xm-1
                    prec_res(i,k,j) = rhs(i,k,j)-t_vec(i,k,j)
                enddo
            enddo
        enddo
        call vec_norm2_local(rhs, local_rhs_norm)
        call vec_norm2_local(prec_res, local_residual_norm)
        call MPI_Allreduce(local_rhs_norm, global_rhs_norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(local_residual_norm, global_residual_norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                           solver_comm, ierr)
        multilevel_preconditioner_smoke = global_norm > 0.0_c_double .and. &
            global_norm < huge(global_norm) .and. global_residual_norm < global_rhs_norm
    end function multilevel_preconditioner_smoke

    !>------------------------------------------------------------
    !! Per-nest initialisation — called once per nest at startup.
    !! Phase 1: no multi-nest cache, so this just allocates state
    !! the first time and otherwise returns. Calling for a second
    !! nest will currently overwrite state for the first.
    !!------------------------------------------------------------
    subroutine init_iter_winds(domain, options)
        implicit none
        type(domain_t),  intent(in) :: domain
        type(options_t), intent(in) :: options
        integer :: ierr
        integer :: nprocs, device_num
        integer :: target_nest
        integer :: my_rank_in_comm
        integer :: env_status, env_length
        character(len=32) :: audit_env, audit_max_iters_env, multilevel_env
        character(len=512) :: audit_file_env
        integer :: audit_read_status
#ifdef USE_NCCL
        integer(c_int) :: nccl_rc
#endif

        ! NCCL communicator is shared across nests (created once on first call).
        ! It's tied to compute_comms which we assume is the same across nests.
#ifdef USE_NCCL
        if (.not. nccl_initialized) then
            call MPI_Comm_size(domain%compute_comms, nprocs,         ierr)
            call MPI_Comm_rank(domain%compute_comms, my_rank_in_comm, ierr)
            device_num = acc_get_device_num(acc_device_nvidia)
            nccl_rc = nccl_comm_init(nccl_comm, nprocs, my_rank_in_comm, &
                                      domain%compute_comms, device_num)
            if (nccl_rc /= 0 .and. STD_OUT_PE) then
                print*, "WARNING: nccl_comm_init failed in wind_iterative, rc=", nccl_rc
            endif
            nccl_stream = transfer(acc_get_cuda_stream(acc_async_sync), nccl_stream)
            nccl_initialized = .true.
        endif
#endif

        target_nest = domain%nest_indx

        ! Already on this nest — nothing to do
        if (target_nest == active_nest_indx) return

        ! Save current state to its cache slot (if any)
        if (active_nest_indx > 0 .and. active_nest_indx <= MAX_NESTS) then
            call save_to_cache(active_nest_indx)
        endif

        ! Restore from cache if we've seen this nest before — O(1) context switch
        if (target_nest > 0 .and. target_nest <= MAX_NESTS .and. &
            domain_cache(target_nest)%valid) then
            call restore_from_cache(target_nest)
            active_nest_indx = target_nest
            return
        endif

        ! --- Fresh init for a nest seen for the first time ---
        solver_comm = domain%compute_comms

        ! Cache MPI rank + neighbour ranks (exchange_krylov_halos hits these ~200x/solve)
        call MPI_Comm_rank(solver_comm, solver_rank, ierr)
        east_neighbor  = solver_rank + 1
        west_neighbor  = solver_rank - 1
        north_neighbor = solver_rank + domain%grid%ximages
        south_neighbor = solver_rank - domain%grid%ximages

        wind_solver_max_iters = options%wind%wind_solver_iterations

        audit_env = ''
        call get_environment_variable('HICAR_WIND_OPERATOR_AUDIT', audit_env, &
                                      length=env_length, status=env_status)
        operator_audit_enabled = env_status == 0 .and. env_length > 0 .and. &
                                 trim(adjustl(audit_env)) /= '0'
        operator_audit_done = .false.
        krylov_audit_written = .false.
        bootstrap_rhs_saved = .false.
        operator_audit_max_iters = 0
        audit_max_iters_env = ''
        call get_environment_variable('HICAR_WIND_OPERATOR_AUDIT_MAX_ITERS', &
                                      audit_max_iters_env, length=env_length, status=env_status)
        if (operator_audit_enabled .and. env_status == 0 .and. env_length > 0) then
            read(audit_max_iters_env(:min(env_length, len(audit_max_iters_env))), *, &
                 iostat=audit_read_status) operator_audit_max_iters
            if (audit_read_status /= 0 .or. operator_audit_max_iters < 0) then
                if (STD_OUT_PE) write(output_unit,'(A,A)') &
                    ' Invalid HICAR_WIND_OPERATOR_AUDIT_MAX_ITERS: ', &
                    trim(audit_max_iters_env)
                error stop
            endif
        endif
        operator_audit_file = 'hicar_wind_operator_audit.csv'
        audit_file_env = ''
        call get_environment_variable('HICAR_WIND_OPERATOR_AUDIT_FILE', audit_file_env, &
                                      length=env_length, status=env_status)
        if (env_status == 0 .and. env_length > 0) then
            operator_audit_file = trim(audit_file_env(:min(env_length, len(audit_file_env))))
        endif

        verbose_solver = options%general%debug .or. &
                         (options%physics%windtype == kITERATIVE_WINDS)

        multilevel_env = ''
        call get_environment_variable('HICAR_WIND_MULTILEVEL', multilevel_env, &
                                      length=env_length, status=env_status)
        multilevel_requested = env_status == 0 .and. env_length > 0 .and. &
                               trim(adjustl(multilevel_env)) /= '0'
        multilevel_setup_attempted = .false.
        multilevel_ready = .false.

        call init_module_vars(domain)

        n_rows_global = mx * my * mz
        n_rows        = xm * ym * zm

        initialized_iter_winds  = .true.
        structure_uploaded = .false.
        active_nest_indx   = target_nest

        if (STD_OUT_PE) print*, "HICAR native wind solver initialised for nest ", target_nest
    end subroutine init_iter_winds


    !>------------------------------------------------------------
    !! Main entry: allocate/update the solver state, then solve and apply the
    !! wind correction.  The optional setup-only path lets wind.F90 allocate
    !! the probe workspace before the first exact D o G calibration, avoiding
    !! a redundant solve with the approximate analytic bootstrap operator.
    !!------------------------------------------------------------
    subroutine calc_iter_winds(domain, alpha_in, div_in, adv_den, setup_only)
        implicit none
        type(domain_t), intent(inout) :: domain
        real, dimension(ims:ime, domain%kms:domain%kme, jms:jme), intent(in) :: alpha_in, div_in
        logical, intent(in) :: adv_den
        logical, intent(in), optional :: setup_only

        integer :: i, j, k
        real    :: alpha_min, alpha_max
        logical :: varying_alpha
        integer :: status, n_iters, apply_status, calibrated_max_iters
        integer :: nan_count
        real(c_double) :: res0, res_final, local_norm2, global_norm2, target_norm, max_x_global
        real(c_double) :: local_apply_stats(2), global_apply_stats(2), operator_error
        integer :: ierr

        ! Copy alpha and div into module-resident arrays on GPU
        !$acc parallel loop gang vector collapse(3) present(alpha, div, alpha_in, div_in)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_s, i_e
                    div(i,k,j)   = div_in(i,k,j)
                    alpha(i,k,j) = alpha_in(i,k,j)
                enddo
            enddo
        enddo

        ! Debug: if the input divergence contains any NaN, dump the min/max of the
        ! state fields that feed the divergence so the source can be traced.
        nan_count = 0
        !$acc parallel loop gang vector collapse(3) reduction(+:nan_count) present(div)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_s, i_e
                    if (div(i,k,j) /= div(i,k,j)) nan_count = nan_count + 1
                enddo
            enddo
        enddo
        if (nan_count > 0) then
            associate(density     => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d,     &
                      pressure    => domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d,    &
                      temperature => domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d, &
                      qv          => domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,  &
                      u           => domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d,            &
                      v           => domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d,            &
                      w_grid      => domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d)
                !$acc update host(density, pressure, temperature, qv, u, v, w_grid)
                write(*,*) '--------------- HICAR wind solver: NaN in input divergence (rank=', &
                           solver_rank, ', count=', nan_count, ') ---------------'
                write(*,*) "  density:     min=", minval(density),     " max=", maxval(density)
                write(*,*) "  pressure:    min=", minval(pressure),    " max=", maxval(pressure)
                write(*,*) "  temperature: min=", minval(temperature), " max=", maxval(temperature)
                write(*,*) "  qv:          min=", minval(qv),          " max=", maxval(qv)
                write(*,*) "  u:           min=", minval(u),           " max=", maxval(u)
                write(*,*) "  v:           min=", minval(v),           " max=", maxval(v)
                write(*,*) "  w:           min=", minval(w_grid),      " max=", maxval(w_grid)
            end associate
        endif

        if (.not. structure_uploaded) then
            ! First call: build coefficient arrays (CPU host needs alpha)
            !$acc update host(alpha)

            call initialize_coefs(domain)

            ! Allocate Krylov / RHS / vertical-line preconditioner state on host then push to device.
            ! Range matches AMGX's lambda_3d: (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1).
            allocate(x_sol  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(r_vec  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(r_hat  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(p_vec  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(p_hat  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(v_vec  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(s_vec  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(s_hat  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(t_vec  (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(rhs     (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(D_inv   (i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(line_cprime(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(prec_res(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))
#ifdef USE_NCCL
            ! Device-resident scalar buffers for NCCL allreduces.
            allocate(sigma_dev(1)); sigma_dev = 0.0_c_double
            allocate(red5_dev (5)); red5_dev  = 0.0_c_double
            allocate(rho0_dev (2)); rho0_dev  = 0.0_c_double
            !$acc enter data copyin(sigma_dev, red5_dev, rho0_dev)
#endif

            ! Persistent halo face buffers — 8 separate buffers so all 4 directions can be
            ! in flight at once (Isend/Irecv + Waitall pattern).
            allocate(east_send (k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(east_recv (k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(west_send (k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(west_recv (k_s-1:k_e+1, j_s-1:j_e+1))
            allocate(north_send(i_s-1:i_e+1, k_s-1:k_e+1))
            allocate(north_recv(i_s-1:i_e+1, k_s-1:k_e+1))
            allocate(south_send(i_s-1:i_e+1, k_s-1:k_e+1))
            allocate(south_recv(i_s-1:i_e+1, k_s-1:k_e+1))

            ! Zero everything to make halo cells well-defined
            x_sol = 0.0_c_double; r_vec = 0.0_c_double; r_hat = 0.0_c_double
            p_vec = 0.0_c_double; p_hat = 0.0_c_double; v_vec = 0.0_c_double
            s_vec = 0.0_c_double; s_hat = 0.0_c_double; t_vec = 0.0_c_double
            rhs   = 0.0_c_double; D_inv = 1.0_c_double; line_cprime = 0.0_c_double; prec_res = 0.0_c_double
            east_send  = 0.0_c_double; east_recv  = 0.0_c_double
            west_send  = 0.0_c_double; west_recv  = 0.0_c_double
            north_send = 0.0_c_double; north_recv = 0.0_c_double
            south_send = 0.0_c_double; south_recv = 0.0_c_double

            !$acc enter data copyin(x_sol, r_vec, r_hat, p_vec, p_hat, v_vec, &
            !$acc                   s_vec, s_hat, t_vec, rhs, D_inv, line_cprime, prec_res, &
            !$acc                   east_send, east_recv, west_send, west_recv, &
            !$acc                   north_send, north_recv, south_send, south_recv)
            !$acc enter data copyin(A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, &
            !$acc                   H_coef, I_coef, J_coef, K_coef, L_coef, M_coef, N_coef, O_coef)

            structure_uploaded = .true.
        else
            ! Subsequent call: detect varying alpha on GPU; refresh coefs only when needed
            alpha_min =  HUGE(1.0)
            alpha_max = -HUGE(1.0)
            !$acc parallel loop gang vector collapse(3) reduction(min:alpha_min) reduction(max:alpha_max) present(alpha)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_s, i_e
                        alpha_min = min(alpha_min, alpha(i,k,j))
                        alpha_max = max(alpha_max, alpha(i,k,j))
                    enddo
                enddo
            enddo
            varying_alpha = (alpha_max > alpha_min)
            ! Under the probed operator the coefficients ARE the exact
            ! alpha-dependent composition; the analytic refresh would
            ! clobber them. Alpha changes are handled by re-probing
            ! (wind.F90::calibrate_projection_operator).
            if (varying_alpha .and. .not. operator_probed) call update_coefs_gpu()
        endif

        ! Calibration only needs the allocated vectors and uploaded stencil
        ! storage.  Return before constructing a preconditioner or solving the
        ! approximate analytic operator; probe_finalize will replace its
        ! coefficients with the exact discrete D o G operator.
        if (present(setup_only)) then
            if (setup_only) return
        endif

        ! Build / refresh the vertical-line block-Jacobi preconditioner from
        ! the current coefficients.
        call build_line_preconditioner()
        if (multilevel_requested .and. operator_probed .and. .not. multilevel_setup_attempted) &
            call setup_multilevel_preconditioner(domain)

        ! Build RHS on GPU (3D layout: rhs(i,k,j) = -2*div for interior, 0 at BCs)
        call compute_rhs_3d()

        if (operator_probed .and. operator_audit_enabled .and. .not. operator_audit_done) then
            call audit_operator_structure(domain)
            operator_audit_done = .true.
        endif

        ! Initial guess: warm-start from the previous solve's lambda.
        ! bicgstab_solve computes r0 = b - A*x0 properly and the
        ! convergence target is normalized by ||b||, so a correlated
        ! guess (RANS per-step solves) cuts iterations directly, while an
        ! uncorrelated one (diagnostic hourly solves) just behaves like a
        ! cold start. x_sol is reset after operator probing.

        !$acc wait

        ! The calibrated operator can use HYPRE's algebraic GPU hierarchy when
        ! the executable was built with HICAR_HYPRE.  The analytic bootstrap
        ! stays on the established native path; it is only the exact probed
        ! matrix that is handed to the external backend.
        if (STD_OUT_PE .and. .not. hypre_selection_reported) then
            write(*,'(A,L1,A,L1)') ' HICAR wind backend selection: operator_probed=', operator_probed, &
                ' hypre_available=', wind_hypre_available()
            flush(output_unit)
            hypre_selection_reported = .true.
        endif
        if (operator_probed) then
            calibrated_max_iters = wind_solver_max_iters
            if (operator_audit_enabled .and. operator_audit_max_iters > 0) then
                calibrated_max_iters = min(calibrated_max_iters, operator_audit_max_iters)
                if (STD_OUT_PE) then
                    write(output_unit,'(A,I0)') &
                        ' HICAR operator audit calibrated-solve iteration cap=', calibrated_max_iters
                    flush(output_unit)
                endif
            endif
            call fgmres_line_solve(domain, calibrated_max_iters, status, n_iters, res0, res_final)
            if (status == 0) then
                call verify_true_residual(domain, res_final, max_x_global)
                ! FGMRES, like the existing BiCGStab implementation, uses
                ! ||b|| as its relative scale.  Do not tighten that criterion
                ! after a warm start by substituting the smaller ||r_0||.
                call vec_norm2_local(rhs, local_norm2)
                call MPI_Allreduce(local_norm2, global_norm2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
                target_norm = max(bicg_tol_abs, bicg_tol_rel * sqrt(global_norm2))
                if (STD_OUT_PE) then
                    write(*,'(A,I0,A,ES12.4,A,ES12.4)') ' HICAR native FGMRES+line: iterations=', n_iters, &
                        ' true_residual=', res_final, ' target=', target_norm
                    flush(output_unit)
                endif
                if (res_final > target_norm) status = 3
            else if (STD_OUT_PE) then
                write(*,'(A,I0,A,I0)') ' HICAR native FGMRES+line failed: status=', status, &
                    ' iterations=', n_iters
                flush(output_unit)
            endif
        else if (operator_probed .and. wind_hypre_available()) then
            !$acc update host(rhs, x_sol)
            call vec_norm2_local(rhs, local_norm2)
            call MPI_Allreduce(local_norm2, global_norm2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
            res0 = sqrt(global_norm2)
            target_norm = max(bicg_tol_abs, bicg_tol_rel * res0)
            call wind_hypre_solve(solver_comm, xs, ys, zs, xm, ym, zm, mx, my, mz, &
                                  i_s, i_e, k_s, k_e, j_s, j_e, A_coef, B_coef, C_coef, D_coef, E_coef, &
                                  F_coef, G_coef, H_coef, I_coef, J_coef, K_coef, L_coef, M_coef, N_coef, O_coef, &
                                  rhs, x_sol, wind_solver_max_iters, bicg_tol_rel, status, n_iters, res_final)
            !$acc update device(x_sol)
            if (status == 0 .and. .not. hypre_operator_verified) then
                call wind_hypre_apply(xs, ys, zs, xm, ym, zm, i_s, i_e, k_s, k_e, j_s, j_e, x_sol, r_vec, apply_status)
                if (apply_status == 0) then
                    call exchange_krylov_halos(x_sol, domain)
                    call spmv(x_sol, t_vec)
                    !$acc update host(t_vec)
                    local_apply_stats = 0.0_c_double
                    do j = j_s, j_e
                        do k = k_s, k_e
                            do i = i_s, i_e
                                local_apply_stats(1) = local_apply_stats(1) + (t_vec(i,k,j)-r_vec(i,k,j))**2
                                local_apply_stats(2) = local_apply_stats(2) + t_vec(i,k,j)**2
                            enddo
                        enddo
                    enddo
                    call MPI_Allreduce(local_apply_stats, global_apply_stats, 2, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
                    operator_error = sqrt(global_apply_stats(1) / max(global_apply_stats(2), tiny(1.0_c_double)))
                    if (STD_OUT_PE) then
                        write(*,'(A,ES12.4)') ' HICAR HYPRE/native operator relative L2=', operator_error
                        flush(output_unit)
                    endif
                    if (operator_error <= 1.0e-10_c_double) then
                        hypre_operator_verified = .true.
                    else
                        if (STD_OUT_PE) write(*,'(A,ES12.4)') ' HICAR HYPRE matrix rejected: operator mismatch=', operator_error
                        status = 4
                    endif
                else
                    if (STD_OUT_PE) write(*,'(A,I0)') ' HICAR HYPRE matrix apply failed: status=', apply_status
                    status = apply_status
                endif
            endif
            if (status == 0) then
                call verify_true_residual(domain, res_final, max_x_global)
                if (STD_OUT_PE) then
                    write(*,'(A,I0,A,ES12.4,A,ES12.4)') ' HICAR HYPRE FGMRES+AMG: iterations=', n_iters, &
                        ' true_residual=', res_final, ' target=', target_norm
                    flush(output_unit)
                endif
                if (res_final > target_norm) then
                    if (STD_OUT_PE) write(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') &
                        ' HICAR HYPRE FGMRES+AMG rejected by true residual: ', res_final, ' target=', target_norm, &
                        ' max |x|=', max_x_global
                    status = 3
                endif
            else if (STD_OUT_PE) then
                write(*,'(A,I0,A,I0)') ' HICAR HYPRE FGMRES+AMG failed: status=', status, &
                    ' iterations=', n_iters
                flush(output_unit)
            endif
        else
            call bicgstab_solve(domain, wind_solver_max_iters, status, n_iters, res0, res_final)
        endif

        if (STD_OUT_PE .and. verbose_solver) then
            write(*,*) ' HICAR BiCGStab status=', status, ' iterations=', n_iters, &
                       ' precond_n_sweeps=', precond_n_sweeps
            write(*,*) '   Residual at iter 0:    ', res0
            write(*,*) '   Residual at final iter:', res_final
            if (t_total_acc > 0.0_c_double) then
                write(*,'(A)')          '   --- per-solve timing breakdown (rank 0, seconds) ---'
                write(*,'(A,F10.4,A)')  '     total           : ', t_total_acc, ''
                write(*,'(A,F10.4,A,F6.1,A)') '     SpMV            : ', t_spmv_acc, &
                    '   (', 100.0_c_double * t_spmv_acc / t_total_acc, ' %)'
                write(*,'(A,F10.4,A,F6.1,A)') '     halo exchange   : ', t_halo_acc, &
                    '   (', 100.0_c_double * t_halo_acc / t_total_acc, ' %)'
                write(*,'(A,F10.4,A,F6.1,A)') '     preconditioner  : ', t_precond_acc, &
                    '   (', 100.0_c_double * t_precond_acc / t_total_acc, ' %)'
                write(*,'(A,F10.4,A,F6.1,A)') '     allreduce       : ', t_allreduce_acc, &
                    '   (', 100.0_c_double * t_allreduce_acc / t_total_acc, ' %)'
                write(*,'(A,F10.4,A,F6.1,A)') '     vector ops      : ', t_vecops_acc, &
                    '   (', 100.0_c_double * t_vecops_acc / t_total_acc, ' %)'
                write(*,'(A,F10.4,A,F6.1,A)') '     other           : ', &
                    t_total_acc - (t_spmv_acc + t_halo_acc + t_precond_acc + t_allreduce_acc + t_vecops_acc), &
                    '   (', 100.0_c_double * (t_total_acc - (t_spmv_acc + t_halo_acc + t_precond_acc + t_allreduce_acc + t_vecops_acc)) / t_total_acc, ' %)'
                if (n_iters > 0) then
                    write(*,'(A,F10.6,A)')  '     per-iter total  : ', t_total_acc / real(n_iters, c_double), ' s'
                endif
            endif
        endif

        ! A capped calibrated solve is an audit-only execution: the Krylov
        ! artifact has been written, but the iterate is not authorized for
        ! physics.  Abort the full MPI application so CPU I/O ranks cannot
        ! remain blocked after compute ranks stop, and so no output can be
        ! mistaken for a scientifically accepted run.
        if (operator_probed .and. operator_audit_enabled .and. operator_audit_max_iters > 0) then
            if (STD_OUT_PE) then
                write(output_unit,'(A,I0,A,I0,A,ES12.4)') &
                    ' HICAR operator audit-only exit: status=', status, &
                    ' iterations=', n_iters, ' residual=', res_final
                flush(output_unit)
            endif
            call MPI_Abort(MPI_COMM_WORLD, 86, ierr)
            error stop
        endif

        ! Adaptive preconditioner retry: any non-converged solve must be
        ! retried while a stronger configured sweep count remains.  Residual
        ! reduction alone is not an acceptance criterion: a solve can shrink
        ! by much more than 100x and still miss the requested tolerance.
        ! Each retry restarts from x_sol = 0 (bicgstab_solve assumes x0=0 —
        ! see r_vec = rhs init).
        ! Mirrors wind_iterative_amgx.F90's prec_max_iters retry.
        if (status /= 0 .and. n_iters > 0) then
            do while (status /= 0 .and. precond_n_sweeps < MAX_PREC_SWEEPS)
                precond_n_sweeps = precond_n_sweeps + 1
                if (STD_OUT_PE) write(*,*) ' Convergence unsatisfactory, retrying solve with precond_n_sweeps=', &
                                           precond_n_sweeps

                ! Restart from x = 0 so bicgstab_solve's r0 = b - A*x0 = b assumption holds.
                call vec_zero(x_sol)
                !$acc wait

                call bicgstab_solve(domain, wind_solver_max_iters, status, n_iters, res0, res_final)

                if (STD_OUT_PE) then
                    write(*,*) '  Retry status=', status, ' iterations=', n_iters
                    write(*,*) '  Residual at iter 0:     ', res0
                    write(*,*) '  Residual at final iter: ', res_final
                endif
            end do

            ! Only the Krylov solver's success status proves convergence.  A
            ! residual reduction alone may still be far above its requested
            ! tolerance, and must never be used to update the winds.
            if (status == 0) then
                if (STD_OUT_PE) write(*,*) ' Retry converged, resetting precond_n_sweeps=', BASE_PREC_SWEEPS
            endif
            precond_n_sweeps = BASE_PREC_SWEEPS

            ! A non-zero status after every retry is fatal.  Continuing with
            ! the failed iterate silently contaminates the following physics
            ! step, even when the residual happened to decrease somewhat.
            if (status /= 0) then
                if (STD_OUT_PE) then
                    associate(density => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                              u       => domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d, &
                              v       => domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d, &
                              w_grid  => domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d, &
                              w_real  => domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d, &
                              alpha_d => domain%vars_3d(domain%var_indx(kVARS%wind_alpha)%v)%data_3d)
                        write(*,*) '---------------HICAR BiCGStab DIVERGED (precond_n_sweeps=', &
                                   precond_n_sweeps, ')-----------'
                        !$acc update host(density, u, v, w_grid, w_real, alpha_d, div)
                        write(*,*) "max abs val of density: ", maxval(abs(density))
                        write(*,*) "max abs val of u:       ", maxval(abs(u))
                        write(*,*) "max abs val of v:       ", maxval(abs(v))
                        write(*,*) "max abs val of w_grid:  ", maxval(abs(w_grid))
                        write(*,*) "max abs val of w_real:  ", maxval(abs(w_real))
                        write(*,*) "max abs val of alpha:   ", maxval(abs(alpha_d))
                        write(*,*) "max abs val of div:     ", maxval(abs(div))
                        call domain_check_winds(domain, "Solver Diverged: ", dqdt=.True.)
                    end associate
                endif
                stop
            endif
        endif

        ! A status that did not meet the retry trigger (for example a failed
        ! independent true-residual check after substantial reduction) is
        ! still not an acceptable wind correction.  Never fall through to
        ! calc_updated_winds with it.
        if (status /= 0) then
            if (STD_OUT_PE) write(*,'(A,I0)') ' HICAR wind solve rejected after acceptance gate: status=', status
            stop
        endif

        ! Preserve the accepted analytic-bootstrap RHS in an existing Krylov
        ! work vector.  The calibrated FGMRES path does not otherwise use
        ! r_hat, so this costs no additional national-scale allocation and
        ! lets the audit compare the bootstrap and failing calibrated RHSs.
        if (.not. operator_probed .and. operator_audit_enabled) then
            call vec_copy(r_hat, rhs)
            bootstrap_rhs_saved = .true.
        endif

        call calc_updated_winds(domain, adv_den)
    end subroutine calc_iter_winds


    !>------------------------------------------------------------
    !! Right-preconditioned BiCGStab with merged reductions.
    !!
    !! Solves A x = b with Block-Jacobi (diagonal) preconditioner M.
    !!
    !! Per iteration: 2 SpMVs, 2 halo exchanges, 2 MPI_Iallreduce
    !! calls, 2 preconditioner applies (fresh M^{-1} on p and s — no
    !! auxiliary recurrences, so no roundoff accumulation that could
    !! stall convergence at high iter count).
    !!
    !! Convergence test mirrors AMGX's COMBINED_REL_INI_ABS:
    !!   converged when ||r|| <= max(tol_abs, tol_rel * ||r0||).
    !! ||r_new||^2 is computed locally from the packed reduction:
    !!   ||r_new||^2 = <s,s> - 2*omega*<t,s> + omega^2*<t,t>
    !! so no separate convergence allreduce is needed.
    !!
    !! 2 allreduces/iter:
    !!   1) sigma = <r̂, v>                                  (1 double)
    !!   2) <t,s>, <t,t>, <r̂,s>, <r̂,t>, <s,s>               (5 doubles)
    !!------------------------------------------------------------
    subroutine bicgstab_solve(domain, max_iters, status_out, n_iters_out, res0_out, res_final_out)
        implicit none
        type(domain_t), intent(in)  :: domain
        integer,        intent(in)  :: max_iters
        integer,        intent(out) :: status_out      ! 0 = converged, 1 = iter cap, 2 = breakdown
        integer,        intent(out) :: n_iters_out
        real(c_double), intent(out) :: res0_out, res_final_out

        integer :: it, ierr
        real(c_double) :: rho, rho_new, alpha_s, beta_s, omega
        real(c_double) :: sigma_local, sigma_global
        real(c_double) :: red5_local(5), red5_global(5)
        real(c_double) :: ts, tt, rs, rt, ss
        real(c_double) :: rnorm_global, rnorm_squared, target_norm, max_x_global
        real(c_double) :: rho0_pack(2), b_norm2
        real(c_double) :: t0_solve, t0_region

        status_out    = 1
        n_iters_out   = 0
        res0_out      = 0.0_c_double
        res_final_out = 0.0_c_double

        ! Reset timing accumulators for this solve
        t_spmv_acc      = 0.0_c_double
        t_halo_acc      = 0.0_c_double
        t_precond_acc   = 0.0_c_double
        t_allreduce_acc = 0.0_c_double
        t_vecops_acc    = 0.0_c_double
        t0_solve        = MPI_Wtime()

        ! --- Initial setup ---
        t0_region = MPI_Wtime()
        ! r0 = b - A*x0. Computed properly (one SpMV) so the caller may
        ! warm-start from the previous solve's lambda; with a zero x0 the
        ! SpMV result is zero and r0 = b as before. x_sol halos may be one
        ! iteration stale from the previous solve — refresh before SpMV.
        call exchange_krylov_halos(x_sol, domain)
        call spmv(x_sol, t_vec)
        call vec_axpby_into(r_vec, 1.0_c_double, rhs, -1.0_c_double, t_vec)
        ! r̂0 = r0 (shadow residual; held fixed across the solve)
        call vec_copy(r_hat, r_vec)
        ! p_0 = r_0
        call vec_copy(p_vec, r_vec)

        ! Initial reductions: rho_0 = <r̂, r_0>, ||r_0||², plus ||b||² for
        ! the convergence target. The target is normalized by ||b|| (the
        ! physical scale of the divergence to remove), NOT by ||r_0||:
        ! otherwise a warm start that shrinks r_0 would tighten its own
        ! target and gain nothing.
        call vec_dot_local(r_hat, r_vec, rho0_pack(1))
        call vec_norm2_local(r_vec, rho0_pack(2))
        call vec_norm2_local(rhs, b_norm2)
        t_vecops_acc = t_vecops_acc + (MPI_Wtime() - t0_region)

        t0_region = MPI_Wtime()
#ifdef USE_NCCL
        rho0_dev(1) = rho0_pack(1); rho0_dev(2) = rho0_pack(2)
        !$acc update device(rho0_dev)
        !$acc host_data use_device(rho0_dev)
        ierr = nccl_allreduce_double_sum(c_loc(rho0_dev), c_loc(rho0_dev), 2, nccl_comm, nccl_stream)
        !$acc end host_data
        !$acc update host(rho0_dev)
        rho0_pack(1) = rho0_dev(1); rho0_pack(2) = rho0_dev(2)
        ! ||b||² reduced separately (host MPI; once per solve)
        call MPI_Allreduce(MPI_IN_PLACE, b_norm2, 1, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, solver_comm, ierr)
#else
        call MPI_Allreduce(MPI_IN_PLACE, rho0_pack, 2, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, b_norm2, 1, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, solver_comm, ierr)
#endif
        t_allreduce_acc = t_allreduce_acc + (MPI_Wtime() - t0_region)

        rho           = rho0_pack(1)
        rnorm_global  = sqrt(rho0_pack(2))
        res0_out      = rnorm_global
        res_final_out = rnorm_global
        target_norm   = max(bicg_tol_abs, bicg_tol_rel * sqrt(b_norm2))

        if (STD_OUT_PE .and. verbose_solver) then
            write(output_unit,'(A,ES12.4,A,ES12.4,A,I0,A)') &
                ' HICAR BiCGStab start: residual=', rnorm_global, &
                ' target=', target_norm, ' max_iterations=', max_iters, '.'
            flush(output_unit)
        endif

        if (rnorm_global <= target_norm) then
            status_out  = 0
            n_iters_out = 0
            t_total_acc = MPI_Wtime() - t0_solve
            return
        endif

        omega   = 1.0_c_double
        alpha_s = 1.0_c_double

        do it = 1, max_iters

            ! ----- First half-step: p_hat = M^{-1} p, v = A p_hat -----
            t0_region = MPI_Wtime()
            call apply_precond(p_vec, p_hat, domain)
            t_precond_acc = t_precond_acc + (MPI_Wtime() - t0_region)

            t0_region = MPI_Wtime()
            call exchange_krylov_halos(p_hat, domain)
            t_halo_acc = t_halo_acc + (MPI_Wtime() - t0_region)

            t0_region = MPI_Wtime()
            call spmv(p_hat, v_vec)
            t_spmv_acc = t_spmv_acc + (MPI_Wtime() - t0_region)

            ! ============== ALLREDUCE 1: sigma = <r̂, v> ==================
            t0_region = MPI_Wtime()
            call vec_dot_local(r_hat, v_vec, sigma_local)
            t_vecops_acc = t_vecops_acc + (MPI_Wtime() - t0_region)

            t0_region = MPI_Wtime()
#ifdef USE_NCCL
            sigma_dev(1) = sigma_local
            !$acc update device(sigma_dev)
            !$acc host_data use_device(sigma_dev)
            ierr = nccl_allreduce_double_sum(c_loc(sigma_dev), c_loc(sigma_dev), 1, nccl_comm, nccl_stream)
            !$acc end host_data
            !$acc update host(sigma_dev)
            sigma_global = sigma_dev(1)
#else
            call MPI_Allreduce(sigma_local, sigma_global, 1, MPI_DOUBLE_PRECISION, &
                               MPI_SUM, solver_comm, ierr)
#endif
            t_allreduce_acc = t_allreduce_acc + (MPI_Wtime() - t0_region)

            if (abs(sigma_global) < breakdown_eps) then
                if (STD_OUT_PE) write(*,*) ' BiCGStab breakdown: <r̂, v> ~ 0 at iter ', it
                status_out  = 2
                n_iters_out = it
                exit
            endif
            alpha_s = rho / sigma_global

            ! s = r - alpha v
            t0_region = MPI_Wtime()
            call vec_axpby_into(s_vec, 1.0_c_double, r_vec, -alpha_s, v_vec)
            t_vecops_acc = t_vecops_acc + (MPI_Wtime() - t0_region)

            ! ----- Second half-step: s_hat = M^{-1} s, t = A s_hat -----
            t0_region = MPI_Wtime()
            call apply_precond(s_vec, s_hat, domain)
            t_precond_acc = t_precond_acc + (MPI_Wtime() - t0_region)

            t0_region = MPI_Wtime()
            call exchange_krylov_halos(s_hat, domain)
            t_halo_acc = t_halo_acc + (MPI_Wtime() - t0_region)

            t0_region = MPI_Wtime()
            call spmv(s_hat, t_vec)
            t_spmv_acc = t_spmv_acc + (MPI_Wtime() - t0_region)

            ! ============== ALLREDUCE 2: 5-value packed =================
            ! <t,s>, <t,t>, <r̂,s>, <r̂,t>, <s,s>
            ! From these we get omega, rho_new (via <r̂,r_new> = <r̂,s> - omega <r̂,t>),
            ! and ||r_new||² (via ss - 2 omega ts + omega² tt) — convergence rides for free.
            t0_region = MPI_Wtime()
            call vec_dots_fused(t_vec, s_vec, r_hat, &
                                red5_local(1), red5_local(2), red5_local(3), &
                                red5_local(4), red5_local(5))
            t_vecops_acc = t_vecops_acc + (MPI_Wtime() - t0_region)

            t0_region = MPI_Wtime()
#ifdef USE_NCCL
            red5_dev(1:5) = red5_local(1:5)
            !$acc update device(red5_dev)
            !$acc host_data use_device(red5_dev)
            ierr = nccl_allreduce_double_sum(c_loc(red5_dev), c_loc(red5_dev), 5, nccl_comm, nccl_stream)
            !$acc end host_data
            !$acc update host(red5_dev)
            red5_global(1:5) = red5_dev(1:5)
#else
            call MPI_Allreduce(red5_local, red5_global, 5, MPI_DOUBLE_PRECISION, &
                               MPI_SUM, solver_comm, ierr)
#endif
            t_allreduce_acc = t_allreduce_acc + (MPI_Wtime() - t0_region)

            ts = red5_global(1)
            tt = red5_global(2)
            rs = red5_global(3)
            rt = red5_global(4)
            ss = red5_global(5)

            if (tt < breakdown_eps) then
                if (STD_OUT_PE) write(*,*) ' BiCGStab breakdown: <t,t> ~ 0 at iter ', it
                status_out  = 2
                n_iters_out = it
                exit
            endif
            omega   = ts / tt
            rho_new = rs - omega * rt

            ! ||r_new||² = <s - omega t, s - omega t> = ss - 2*omega*ts + omega²*tt
            rnorm_squared = ss - 2.0_c_double*omega*ts + omega*omega*tt
            if (rnorm_squared < 0.0_c_double) rnorm_squared = 0.0_c_double  ! roundoff guard
            rnorm_global  = sqrt(rnorm_squared)
            res_final_out = rnorm_global
            n_iters_out   = it

            if (STD_OUT_PE .and. verbose_solver .and. mod(it, SOLVER_PROGRESS_INTERVAL) == 0) then
                write(output_unit,'(A,I0,A,ES12.4,A,ES12.4,A,F10.2,A)') &
                    ' HICAR BiCGStab progress: iteration=', it, &
                    ' residual=', rnorm_global, ' target=', target_norm, &
                    ' elapsed_s=', MPI_Wtime() - t0_solve, '.'
                flush(output_unit)
            endif

            ! x = x + alpha*p_hat + omega*s_hat   (single fused kernel)
            t0_region = MPI_Wtime()
            call vec_axpy2(x_sol, alpha_s, p_hat, omega, s_hat)
            t_vecops_acc = t_vecops_acc + (MPI_Wtime() - t0_region)

            if (rnorm_global <= target_norm) then
                status_out = 0
                exit
            endif

            if (abs(omega) < breakdown_eps) then
                if (STD_OUT_PE) write(*,*) ' BiCGStab breakdown: omega ~ 0 at iter ', it
                status_out  = 2
                exit
            endif

            ! r_new = s - omega*t
            t0_region = MPI_Wtime()
            call vec_axpby_into(r_vec, 1.0_c_double, s_vec, -omega, t_vec)
            t_vecops_acc = t_vecops_acc + (MPI_Wtime() - t0_region)

            if (abs(rho) < breakdown_eps) then
                if (STD_OUT_PE) write(*,*) ' BiCGStab breakdown: rho ~ 0 at iter ', it
                status_out  = 2
                exit
            endif
            beta_s = (rho_new / rho) * (alpha_s / omega)

            ! p_new = r + beta * (p - omega * v)
            t0_region = MPI_Wtime()
            call vec_p_update(p_vec, r_vec, beta_s, omega, v_vec)
            t_vecops_acc = t_vecops_acc + (MPI_Wtime() - t0_region)

            rho = rho_new
        enddo

        ! Recurrence residuals can be spuriously small for a non-normal
        ! operator after a long Krylov run.  Verify the physical algebraic
        ! residual independently before accepting or applying this solution.
        if (status_out == 0) then
            call verify_true_residual(domain, rnorm_global, max_x_global)
            res_final_out = rnorm_global
            if (rnorm_global > target_norm) then
                if (STD_OUT_PE) write(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') &
                    ' HICAR BiCGStab rejected false recurrence convergence: true residual=', rnorm_global, &
                    ' target=', target_norm, ' max |x|=', max_x_global
                status_out = 3
            endif
        endif

        t_total_acc = MPI_Wtime() - t0_solve

    end subroutine bicgstab_solve


    !> Recompute ||b-Ax|| from the final solution.  This must stay separate
    !! from the BiCGStab recurrence because residual drift can otherwise admit
    !! an enormous, physically invalid correction as "converged".
    subroutine verify_true_residual(domain, residual_norm, max_x_global)
        implicit none
        type(domain_t), intent(in) :: domain
        real(c_double), intent(out) :: residual_norm, max_x_global
        real(c_double) :: local_norm2, global_norm2, max_x_local
        integer :: i, j, k, ierr

        call exchange_krylov_halos(x_sol, domain)
        call spmv(x_sol, t_vec)
        call vec_axpby_into(r_vec, 1.0_c_double, rhs, -1.0_c_double, t_vec)
        call vec_norm2_local(r_vec, local_norm2)

        max_x_local = 0.0_c_double
        !$acc parallel loop gang vector collapse(3) reduction(max:max_x_local) present(x_sol)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    max_x_local = max(max_x_local, abs(x_sol(i,k,j)))
                enddo
            enddo
        enddo
        call MPI_Allreduce(local_norm2, global_norm2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(max_x_local, max_x_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, solver_comm, ierr)
        residual_norm = sqrt(global_norm2)
    end subroutine verify_true_residual


    !>------------------------------------------------------------
    !! 15-point stencil SpMV: y = A x.
    !!
    !! Split into 4 kernels by row category to remove branches from
    !! the dominant interior loop:
    !!   1) Lateral identity rows (i=0, i=mx-1, j=0, j=my-1)
    !!   2) Top BC k=mz-1 (interior i,j only)
    !!   3) Bottom BC k=0 (interior i,j only)
    !!   4) Interior 15-pt stencil — branchless, vectorizable
    !!
    !! Interior optimisations:
    !!   - D/E/F/G coefs are all 1/dx² (constants) — use a scalar,
    !!     skip the 4 array reads.
    !!   - A coef is derivable as -4/dx² - B - C — compute on the
    !!     fly, skip the A_coef array read.
    !! Net: 5 fewer single-precision reads per interior cell out of
    !! ~140 bytes/cell traffic (~14% reduction), plus the branchless
    !! inner loop auto-vectorises cleanly.
    !!------------------------------------------------------------
    !> Restarted right-preconditioned flexible GMRES for the calibrated wind
    !! operator.  The Arnoldi basis and its preconditioned images live on the
    !! GPU only for this solve and are released afterwards.  This deliberately
    !! avoids the unbounded host FGMRES allocation observed on Switzerland.
    subroutine fgmres_line_solve(domain, max_iters, status_out, n_iters_out, res0_out, res_final_out)
        implicit none
        type(domain_t), intent(in) :: domain
        integer, intent(in) :: max_iters
        integer, intent(out) :: status_out, n_iters_out
        real(c_double), intent(out) :: res0_out, res_final_out
        real(c_double), allocatable :: v_basis(:,:,:,:), z_basis(:,:,:,:)
        real(c_double) :: h(FGMRES_RESTART+1,FGMRES_RESTART), h_local(FGMRES_RESTART+1)
        real(c_double) :: h_raw(FGMRES_RESTART+1,FGMRES_RESTART)
        real(c_double) :: cs(FGMRES_RESTART), sn(FGMRES_RESTART), g(FGMRES_RESTART+1)
        real(c_double) :: ycoef(FGMRES_RESTART), dots_local(FGMRES_RESTART)
        real(c_double) :: beta, bnorm2, target_norm, tmp, denom
        integer :: ierr, cycle, j, i, used, total

        status_out = 1; n_iters_out = 0; res0_out = 0.0_c_double; res_final_out = 0.0_c_double
        t_total_acc = 0.0_c_double
        call exchange_krylov_halos(x_sol, domain)
        call spmv(x_sol, t_vec)
        call vec_axpby_into(r_vec, 1.0_c_double, rhs, -1.0_c_double, t_vec)
        call vec_norm2_local(r_vec, beta)
        call vec_norm2_local(rhs, bnorm2)
        call MPI_Allreduce(MPI_IN_PLACE, beta, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, bnorm2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        beta = sqrt(max(beta, 0.0_c_double))
        target_norm = max(bicg_tol_abs, bicg_tol_rel * sqrt(max(bnorm2, 0.0_c_double)))
        res0_out = beta; res_final_out = beta
        if (STD_OUT_PE .and. verbose_solver) then
            write(output_unit,'(A,ES12.4,A,ES12.4,A,I0,A)') ' HICAR native FGMRES+line start: residual=', beta, &
                ' target=', target_norm, ' restart=', FGMRES_RESTART, '.'
            flush(output_unit)
        endif
        if (beta <= target_norm) then; status_out = 0; return; endif

        allocate(v_basis(i_s-1:i_e+1,k_s-1:k_e+1,j_s-1:j_e+1,FGMRES_RESTART+1))
        allocate(z_basis(i_s-1:i_e+1,k_s-1:k_e+1,j_s-1:j_e+1,FGMRES_RESTART))
        v_basis = 0.0_c_double; z_basis = 0.0_c_double
        !$acc enter data copyin(v_basis, z_basis)
        total = 0
        do cycle = 1, max(1, (max_iters + FGMRES_RESTART - 1) / FGMRES_RESTART)
            if (total >= max_iters) exit
            call vec_axpby_into(v_basis(:,:,:,1), 1.0_c_double/beta, r_vec, 0.0_c_double, r_vec)
            h = 0.0_c_double; h_raw = 0.0_c_double
            cs = 0.0_c_double; sn = 0.0_c_double; g = 0.0_c_double; g(1) = beta
            used = min(FGMRES_RESTART, max_iters-total)
            do j = 1, used
                call apply_precond(v_basis(:,:,:,j), z_basis(:,:,:,j), domain)
                call exchange_krylov_halos(z_basis(:,:,:,j), domain)
                call spmv(z_basis(:,:,:,j), t_vec)
                dots_local = 0.0_c_double
                do i = 1, j
                    call vec_dot_local(v_basis(:,:,:,i), t_vec, dots_local(i))
                enddo
                call MPI_Allreduce(dots_local, h_local, j, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
                do i = 1, j
                    h(i,j) = h_local(i)
                    call vec_axpby_into(t_vec, -h(i,j), v_basis(:,:,:,i), 1.0_c_double, t_vec)
                enddo
                ! Reorthogonalize.  On the long, non-normal calibrated
                ! operator a single classical Gram-Schmidt pass loses enough
                ! orthogonality that a restarted solve can report no useful
                ! progress after the first few cycles.  Accumulate the second
                ! projection in the same Hessenberg column rather than
                ! discarding it.
                dots_local = 0.0_c_double
                do i = 1, j
                    call vec_dot_local(v_basis(:,:,:,i), t_vec, dots_local(i))
                enddo
                call MPI_Allreduce(dots_local, h_local, j, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
                do i = 1, j
                    h(i,j) = h(i,j) + h_local(i)
                    call vec_axpby_into(t_vec, -h_local(i), v_basis(:,:,:,i), 1.0_c_double, t_vec)
                enddo
                call vec_norm2_local(t_vec, h(j+1,j))
                call MPI_Allreduce(MPI_IN_PLACE, h(j+1,j), 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
                h(j+1,j) = sqrt(max(h(j+1,j), 0.0_c_double))
                h_raw(1:j+1,j) = h(1:j+1,j)
                if (h(j+1,j) > breakdown_eps .and. j < FGMRES_RESTART+1) then
                    call vec_axpby_into(v_basis(:,:,:,j+1), 1.0_c_double/h(j+1,j), t_vec, 0.0_c_double, t_vec)
                endif
                do i = 1, j-1
                    tmp = cs(i)*h(i,j) + sn(i)*h(i+1,j)
                    h(i+1,j) = -sn(i)*h(i,j) + cs(i)*h(i+1,j)
                    h(i,j) = tmp
                enddo
                denom = sqrt(h(j,j)*h(j,j) + h(j+1,j)*h(j+1,j))
                if (denom <= breakdown_eps) then; used = j; exit; endif
                cs(j) = h(j,j)/denom; sn(j) = h(j+1,j)/denom
                h(j,j) = denom; h(j+1,j) = 0.0_c_double
                tmp = cs(j)*g(j) + sn(j)*g(j+1)
                g(j+1) = -sn(j)*g(j) + cs(j)*g(j+1); g(j) = tmp
                total = total + 1; n_iters_out = total; res_final_out = abs(g(j+1))
                if (STD_OUT_PE .and. verbose_solver .and. mod(total,SOLVER_PROGRESS_INTERVAL) == 0) then
                    write(output_unit,'(A,I0,A,ES12.4,A,ES12.4,A)') ' HICAR native FGMRES+line progress: iteration=', total, &
                        ' residual=', res_final_out, ' target=', target_norm, '.'
                    flush(output_unit)
                endif
                if (res_final_out <= target_norm) then; used = j; exit; endif
            enddo
            if (operator_audit_enabled .and. .not. krylov_audit_written) then
                call write_krylov_audit(h_raw, used, v_basis, r_vec)
                krylov_audit_written = .true.
            endif
            ycoef = 0.0_c_double
            do i = used, 1, -1
                ycoef(i) = g(i)
                do j = i+1, used; ycoef(i) = ycoef(i) - h(i,j)*ycoef(j); enddo
                ycoef(i) = ycoef(i) / h(i,i)
            enddo
            do i = 1, used
                ! Use the in-place AXPY primitive. Passing x_sol both as an
                ! output and an INTENT(IN) argument to axpby is forbidden
                ! aliasing in Fortran and let the accelerator compiler discard
                ! the restart correction on the 250 m run.
                call vec_axpy(x_sol, ycoef(i), z_basis(:,:,:,i))
            enddo
            call exchange_krylov_halos(x_sol, domain)
            call spmv(x_sol, t_vec)
            call vec_axpby_into(r_vec, 1.0_c_double, rhs, -1.0_c_double, t_vec)
            call vec_norm2_local(r_vec, beta)
            call MPI_Allreduce(MPI_IN_PLACE, beta, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
            beta = sqrt(max(beta, 0.0_c_double)); res_final_out = beta
            if (beta <= target_norm) then; status_out = 0; exit; endif
            if (used < min(FGMRES_RESTART, max_iters-total)) exit
        enddo
        !$acc exit data delete(v_basis, z_basis)
        deallocate(v_basis, z_basis)
    end subroutine fgmres_line_solve

    subroutine spmv(x, y)
        implicit none
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: x
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: y
        integer :: i, j, k
        integer :: i_int_lo, i_int_hi, j_int_lo, j_int_hi
        real(c_double) :: denom, dzdx_s, dzdy_s, dz_kp1, two_dx
        real(c_double) :: inv_dx2, neg_4_inv_dx2
        real(c_double) :: b_val, c_val, a_val

        two_dx        = 2.0_c_double * real(dx, c_double)
        inv_dx2       = 1.0_c_double / (real(dx, c_double) * real(dx, c_double))
        neg_4_inv_dx2 = -4.0_c_double * inv_dx2

        ! Interior (i,j) ranges (avoiding lateral BC cells)
        i_int_lo = max(xs,         1)
        i_int_hi = min(xs + xm - 1, mx - 2)
        j_int_lo = max(ys,         1)
        j_int_hi = min(ys + ym - 1, my - 2)

        ! All physical rows use the full probed 15-point coefficients
        ! (no hardcoded D/E/F/G = 1/dx^2 shortcuts, no analytic vertical
        ! BC rows); the ghost planes k=0 and k=mz-1 are identity rows
        ! (the w_to_grid correction operator never reads ghost lambda).
        if (operator_probed) then
            ! Lateral identity rows. NOTE: duplicated in the analytic
            ! branch below — keep the two copies in sync if editing.
            if (xs == 0) then
                !$acc parallel loop gang vector collapse(2) present(x, y)
                do j = ys, ys + ym - 1
                    do k = zs, zs + zm - 1
                        y(0, k, j) = x(0, k, j)
                    enddo
                enddo
            endif
            if (xs + xm - 1 == mx - 1) then
                !$acc parallel loop gang vector collapse(2) present(x, y)
                do j = ys, ys + ym - 1
                    do k = zs, zs + zm - 1
                        y(mx-1, k, j) = x(mx-1, k, j)
                    enddo
                enddo
            endif
            if (ys == 0) then
                !$acc parallel loop gang vector collapse(2) present(x, y)
                do k = zs, zs + zm - 1
                    do i = xs, xs + xm - 1
                        y(i, k, 0) = x(i, k, 0)
                    enddo
                enddo
            endif
            if (ys + ym - 1 == my - 1) then
                !$acc parallel loop gang vector collapse(2) present(x, y)
                do k = zs, zs + zm - 1
                    do i = xs, xs + xm - 1
                        y(i, k, my-1) = x(i, k, my-1)
                    enddo
                enddo
            endif

            if (i_int_lo <= i_int_hi .and. j_int_lo <= j_int_hi) then
                ! Ghost-plane identity rows
                !$acc parallel loop gang vector collapse(2) present(x, y)
                do j = j_int_lo, j_int_hi
                    do i = i_int_lo, i_int_hi
                        y(i, 0,    j) = x(i, 0,    j)
                        y(i, mz-1, j) = x(i, mz-1, j)
                    enddo
                enddo

                !$acc parallel loop gang vector collapse(3) &
                !$acc present(x, y, A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, &
                !$acc         H_coef, I_coef, J_coef, K_coef, L_coef, M_coef, N_coef, O_coef)
                do j = j_int_lo, j_int_hi
                    do k = 1, mz - 2
                        do i = i_int_lo, i_int_hi
                            y(i,k,j) =                                                       &
                                  real(O_coef(i,k,j), c_double) * x(i,   k-1, j-1)           &
                                + real(G_coef(i,k,j), c_double) * x(i,   k,   j-1)           &
                                + real(M_coef(i,k,j), c_double) * x(i,   k+1, j-1)           &
                                + real(K_coef(i,k,j), c_double) * x(i-1, k-1, j  )           &
                                + real(C_coef(i,k,j), c_double) * x(i,   k-1, j  )           &
                                + real(J_coef(i,k,j), c_double) * x(i+1, k-1, j  )           &
                                + real(E_coef(i,k,j), c_double) * x(i-1, k,   j  )           &
                                + real(A_coef(i,k,j), c_double) * x(i,   k,   j  )           &
                                + real(D_coef(i,k,j), c_double) * x(i+1, k,   j  )           &
                                + real(I_coef(i,k,j), c_double) * x(i-1, k+1, j  )           &
                                + real(B_coef(i,k,j), c_double) * x(i,   k+1, j  )           &
                                + real(H_coef(i,k,j), c_double) * x(i+1, k+1, j  )           &
                                + real(N_coef(i,k,j), c_double) * x(i,   k-1, j+1)           &
                                + real(F_coef(i,k,j), c_double) * x(i,   k,   j+1)           &
                                + real(L_coef(i,k,j), c_double) * x(i,   k+1, j+1)
                        enddo
                    enddo
                enddo
            endif
            return
        endif

        ! ============== 1. Lateral identity rows ==============
        ! West face (i=0): only owned at west-boundary ranks
        if (xs == 0) then
            !$acc parallel loop gang vector collapse(2) present(x, y)
            do j = ys, ys + ym - 1
                do k = zs, zs + zm - 1
                    y(0, k, j) = x(0, k, j)
                enddo
            enddo
        endif
        ! East face (i=mx-1)
        if (xs + xm - 1 == mx - 1) then
            !$acc parallel loop gang vector collapse(2) present(x, y)
            do j = ys, ys + ym - 1
                do k = zs, zs + zm - 1
                    y(mx-1, k, j) = x(mx-1, k, j)
                enddo
            enddo
        endif
        ! South face (j=0)
        if (ys == 0) then
            !$acc parallel loop gang vector collapse(2) present(x, y)
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    y(i, k, 0) = x(i, k, 0)
                enddo
            enddo
        endif
        ! North face (j=my-1)
        if (ys + ym - 1 == my - 1) then
            !$acc parallel loop gang vector collapse(2) present(x, y)
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    y(i, k, my-1) = x(i, k, my-1)
                enddo
            enddo
        endif

        ! ============== 2. Top BC k=mz-1 ==============
        ! 2-pt stencil: y = (x(i,k,j) - x(i,k-1,j)) / dz_if(i,k,j)
        ! Restrict to interior (i,j) so we don't double-write lateral cells.
        if (i_int_lo <= i_int_hi .and. j_int_lo <= j_int_hi) then
            !$acc parallel loop gang vector collapse(2) present(x, y, dz_if)
            do j = j_int_lo, j_int_hi
                do i = i_int_lo, i_int_hi
                    y(i, mz-1, j) = (x(i, mz-1, j) - x(i, mz-2, j)) / real(dz_if(i, mz-1, j), c_double)
                enddo
            enddo
        endif

        ! ============== 3. Bottom BC k=0 ==============
        ! 10-pt terrain-following stencil
        if (i_int_lo <= i_int_hi .and. j_int_lo <= j_int_hi) then
            !$acc parallel loop gang vector collapse(2) &
            !$acc present(x, y, dzdx_surf, dzdy_surf, alpha, jaco, dz_if) &
            !$acc private(denom, dzdx_s, dzdy_s, dz_kp1)
            do j = j_int_lo, j_int_hi
                do i = i_int_lo, i_int_hi
                    dzdx_s = real(dzdx_surf(i,j), c_double)
                    dzdy_s = real(dzdy_surf(i,j), c_double)
                    dz_kp1 = real(dz_if(i,1,j), c_double)
                    denom  = 2.0_c_double * (dzdx_s**2 + dzdy_s**2 + real(alpha(i,k_s,j),c_double)**2) &
                             / real(jaco(i,k_s,j), c_double)
                    y(i,0,j) =                                                                       &
                          ( dzdy_s/(denom*two_dx)) * x(i,   0, j-1)                                  &
                        + ( dzdy_s/(denom*two_dx)) * x(i,   1, j-1)                                  &
                        + ( dzdx_s/(denom*two_dx)) * x(i-1, 0, j  )                                  &
                        + (-1.0_c_double / dz_kp1) * x(i,   0, j  )                                  &
                        + (-dzdx_s/(denom*two_dx)) * x(i+1, 0, j  )                                  &
                        + ( dzdx_s/(denom*two_dx)) * x(i-1, 1, j  )                                  &
                        + ( 1.0_c_double / dz_kp1) * x(i,   1, j  )                                  &
                        + (-dzdx_s/(denom*two_dx)) * x(i+1, 1, j  )                                  &
                        + (-dzdy_s/(denom*two_dx)) * x(i,   0, j+1)                                  &
                        + (-dzdy_s/(denom*two_dx)) * x(i,   1, j+1)
                enddo
            enddo
        endif

        ! ============== 4. Interior 15-pt stencil (branchless, vectorisable) ==============
        ! Skip A/D/E/F/G array reads:
        !   D = E = F = G = 1/dx²  (use scalar)
        !   A = -4/dx² - B - C     (compute on the fly)
        if (i_int_lo <= i_int_hi .and. j_int_lo <= j_int_hi) then
            !$acc parallel loop gang vector collapse(3) &
            !$acc present(x, y, B_coef, C_coef, H_coef, I_coef, J_coef, K_coef, &
            !$acc         L_coef, M_coef, N_coef, O_coef) &
            !$acc private(b_val, c_val, a_val)
            do j = j_int_lo, j_int_hi
                do k = 1, mz - 2
                    do i = i_int_lo, i_int_hi
                        b_val = real(B_coef(i,k,j), c_double)
                        c_val = real(C_coef(i,k,j), c_double)
                        a_val = neg_4_inv_dx2 - b_val - c_val

                        y(i,k,j) =                                                       &
                              real(O_coef(i,k,j), c_double) * x(i,   k-1, j-1)           &
                            + inv_dx2                       * x(i,   k,   j-1)           &
                            + real(M_coef(i,k,j), c_double) * x(i,   k+1, j-1)           &
                            + real(K_coef(i,k,j), c_double) * x(i-1, k-1, j  )           &
                            + c_val                         * x(i,   k-1, j  )           &
                            + real(J_coef(i,k,j), c_double) * x(i+1, k-1, j  )           &
                            + inv_dx2                       * x(i-1, k,   j  )           &
                            + a_val                         * x(i,   k,   j  )           &
                            + inv_dx2                       * x(i+1, k,   j  )           &
                            + real(I_coef(i,k,j), c_double) * x(i-1, k+1, j  )           &
                            + b_val                         * x(i,   k+1, j  )           &
                            + real(H_coef(i,k,j), c_double) * x(i+1, k+1, j  )           &
                            + real(N_coef(i,k,j), c_double) * x(i,   k-1, j+1)           &
                            + inv_dx2                       * x(i,   k,   j+1)           &
                            + real(L_coef(i,k,j), c_double) * x(i,   k+1, j+1)
                    enddo
                enddo
            enddo
        endif

    end subroutine spmv


    !>------------------------------------------------------------
    !! Build one vertical tridiagonal block per horizontal grid column.
    !!
    !! Horizontal and diagonal-vertical stencil entries are deliberately left
    !! to the outer Krylov iteration; B/C plus each column's boundary rows are
    !! factored exactly with Thomas elimination.  This is the natural stronger
    !! local preconditioner for the thin, vertically stretched Alpine grid.
    !! D_inv stores inverse pivots and line_cprime the upper factors.
    !!------------------------------------------------------------
    subroutine build_line_preconditioner()
        implicit none
        integer :: i, j, k
        real(c_double) :: diag, lower, upper, pivot

        !$acc parallel loop gang vector collapse(2) &
        !$acc present(D_inv, line_cprime, A_coef, B_coef, C_coef, dz_if) &
        !$acc private(diag, lower, upper, pivot)
        do j = ys, ys + ym - 1
            do i = xs, xs + xm - 1
                if (i <= 0 .or. j <= 0 .or. i >= mx-1 .or. j >= my-1) then
                    do k = zs, zs + zm - 1
                        D_inv(i,k,j) = 1.0_c_double
                        line_cprime(i,k,j) = 0.0_c_double
                    enddo
                else
                    ! Bottom row: identity for the probed operator; otherwise
                    ! retain its two-point vertical boundary condition.
                    if (operator_probed) then
                        diag = 1.0_c_double; upper = 0.0_c_double
                    else
                        diag = -1.0_c_double / real(dz_if(i,1,j), c_double)
                        upper =  1.0_c_double / real(dz_if(i,1,j), c_double)
                    endif
                    if (abs(diag) < breakdown_eps) diag = sign(breakdown_eps, diag)
                    D_inv(i,0,j) = 1.0_c_double / diag
                    line_cprime(i,0,j) = upper * D_inv(i,0,j)

                    do k = 1, mz - 2
                        lower = real(C_coef(i,k,j), c_double)
                        diag  = real(A_coef(i,k,j), c_double)
                        upper = real(B_coef(i,k,j), c_double)
                        pivot = diag - lower * line_cprime(i,k-1,j)
                        if (abs(pivot) < breakdown_eps) pivot = sign(breakdown_eps, pivot)
                        D_inv(i,k,j) = 1.0_c_double / pivot
                        line_cprime(i,k,j) = upper * D_inv(i,k,j)
                    enddo

                    ! Top row mirrors the bottom treatment.
                    if (operator_probed) then
                        diag = 1.0_c_double; lower = 0.0_c_double
                    else
                        diag  =  1.0_c_double / real(dz_if(i,mz-1,j), c_double)
                        lower = -1.0_c_double / real(dz_if(i,mz-1,j), c_double)
                    endif
                    pivot = diag - lower * line_cprime(i,mz-2,j)
                    if (abs(pivot) < breakdown_eps) pivot = sign(breakdown_eps, pivot)
                    D_inv(i,mz-1,j) = 1.0_c_double / pivot
                    line_cprime(i,mz-1,j) = 0.0_c_double
                endif
            enddo
        enddo
    end subroutine build_line_preconditioner
    !>------------------------------------------------------------
    !! operator probing (orchestrated by wind.F90::calibrate_projection_operator).
    !!
    !! The exact discrete projection requires the solver matrix to be
    !! the composition A = 2 * D o G of the model's divergence operator
    !! D (calc_divergence) with the velocity-correction operator G
    !! (calc_updated_winds, w_to_grid form). Rather than deriving the
    !! composed coefficients by hand, they are extracted numerically:
    !! apply G then D to 27 lattice-colored indicator fields (3x3x3
    !! coloring, so no two stencil supports overlap) and read off all
    !! 15 coefficients per cell. By construction the result is
    !! bit-consistent with whatever D and G implement.
    !!------------------------------------------------------------


    !> Fill x_sol with the (ca,cb,cc) coloring on physical cells
    !! (global-index lattice), zero on ghost planes and the lateral
    !! identity ring — matching solve-time lambda there.
    subroutine probe_lambda_pattern(ca, cb, cc)
        implicit none
        integer, intent(in) :: ca, cb, cc
        integer :: i, j, k

        do j = j_s-1, j_e+1
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    if (i >= 1 .and. i <= mx-2 .and. k >= 1 .and. k <= mz-2 .and. &
                        j >= 1 .and. j <= my-2 .and. &
                        mod(i,3) == ca .and. mod(k,3) == cb .and. mod(j,3) == cc) then
                        x_sol(i,k,j) = 1.0_c_double
                    else
                        x_sol(i,k,j) = 0.0_c_double
                    endif
                enddo
            enddo
        enddo
        !$acc update device(x_sol)
    end subroutine probe_lambda_pattern

    !> Zero the u/v/w dqdt workspace so calc_updated_winds leaves the
    !! pure correction G(lambda) in it.
    subroutine probe_zero_corrections(domain)
        implicit none
        type(domain_t), intent(inout) :: domain
        integer :: i, j, k

        associate(u_q => domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d, &
                  v_q => domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d, &
                  w_q => domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d, &
                  kms => domain%kms, kme => domain%kme)
        !$acc parallel default(present)
        !$acc loop gang vector collapse(3)
        do j = jms, jme
            do k = kms, kme
                do i = ims, ime+1
                    u_q(i,k,j) = 0.0
                enddo
            enddo
        enddo
        !$acc loop gang vector collapse(3)
        do j = jms, jme+1
            do k = kms, kme
                do i = ims, ime
                    v_q(i,k,j) = 0.0
                enddo
            enddo
        enddo
        !$acc loop gang vector collapse(3)
        do j = jms, jme
            do k = kms, kme
                do i = ims, ime
                    w_q(i,k,j) = 0.0
                enddo
            enddo
        enddo
        !$acc end parallel
        end associate
    end subroutine probe_zero_corrections

    !> Apply the correction operator G from the current x_sol.
    subroutine probe_apply_corrections(domain, adv_den)
        implicit none
        type(domain_t), intent(inout) :: domain
        logical, intent(in) :: adv_den
        call calc_updated_winds(domain, adv_den)
    end subroutine probe_apply_corrections

    !> Record the probed response T = 2*div(G(pattern)) into the stencil
    !! coefficient arrays (host side; finalize pushes to device). Any
    !! response at an offset outside the 15-point footprint is leakage
    !! (a discretization wider than the stencil) and is reported.
    subroutine probe_record(domain, Tfield, ca, cb, cc, max_leak)
        implicit none
        type(domain_t), intent(in) :: domain
        real, intent(in) :: Tfield(domain%ims:domain%ime, domain%kms:domain%kme, domain%jms:domain%jme)
        integer, intent(in) :: ca, cb, cc
        real, intent(inout) :: max_leak
        integer :: i, j, k, di, dk, dj
        real :: tv

        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_s, i_e
                    di = mod(ca - mod(i,3) + 3, 3); if (di == 2) di = -1
                    dk = mod(cb - mod(k,3) + 3, 3); if (dk == 2) dk = -1
                    dj = mod(cc - mod(j,3) + 3, 3); if (dj == 2) dj = -1
                    tv = 2.0 * Tfield(i,k,j)
                    if (dk == 0) then
                        if (di == 0 .and. dj == 0) then
                            A_coef(i,k,j) = tv
                        else if (di == 1 .and. dj == 0) then
                            D_coef(i,k,j) = tv
                        else if (di == -1 .and. dj == 0) then
                            E_coef(i,k,j) = tv
                        else if (di == 0 .and. dj == 1) then
                            F_coef(i,k,j) = tv
                        else if (di == 0 .and. dj == -1) then
                            G_coef(i,k,j) = tv
                        else
                            max_leak = max(max_leak, abs(tv))
                        endif
                    else if (dk == 1) then
                        if (di == 0 .and. dj == 0) then
                            B_coef(i,k,j) = tv
                        else if (di == 1 .and. dj == 0) then
                            H_coef(i,k,j) = tv
                        else if (di == -1 .and. dj == 0) then
                            I_coef(i,k,j) = tv
                        else if (di == 0 .and. dj == 1) then
                            L_coef(i,k,j) = tv
                        else if (di == 0 .and. dj == -1) then
                            M_coef(i,k,j) = tv
                        else
                            max_leak = max(max_leak, abs(tv))
                        endif
                    else
                        if (di == 0 .and. dj == 0) then
                            C_coef(i,k,j) = tv
                        else if (di == 1 .and. dj == 0) then
                            J_coef(i,k,j) = tv
                        else if (di == -1 .and. dj == 0) then
                            K_coef(i,k,j) = tv
                        else if (di == 0 .and. dj == 1) then
                            N_coef(i,k,j) = tv
                        else if (di == 0 .and. dj == -1) then
                            O_coef(i,k,j) = tv
                        else
                            max_leak = max(max_leak, abs(tv))
                        endif
                    endif
                enddo
            enddo
        enddo
    end subroutine probe_record

    !> Switch the solver to the probed operator and push coefficients
    !! to the device.
    subroutine probe_finalize(max_leak)
        implicit none
        real, intent(in) :: max_leak
        operator_probed = .true.
        operator_audit_done = .false.
        krylov_audit_written = .false.
        hypre_selection_reported = .false.
        hypre_operator_verified = .false.
        call wind_hypre_invalidate()
        !$acc update device(A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, &
        !$acc               H_coef, I_coef, J_coef, K_coef, L_coef, M_coef, N_coef, O_coef)
        ! x_sol holds the last probe pattern — useless as a warm-start
        ! guess for the next solve; reset it.
        call vec_zero(x_sol)
        if (STD_OUT_PE) then
            write(*,*) "Poison projection operator calibrated by probing (A = 2*D o G)."
            write(*,*) "  off-stencil leakage max |T| = ", max_leak, " (should be ~0)"
        endif
    end subroutine probe_finalize


    !> Fill the solve vector with a deterministic, non-coloured pseudo-random
    !! field.  This deliberately crosses every rank interface and is used only
    !! to verify that the matrix reconstructed from coloured probes is the
    !! same operator as the direct G then D application.
    subroutine probe_random_pattern()
        implicit none
        integer :: i, j, k

        do j = j_s-1, j_e+1
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    if (i >= 1 .and. i <= mx-2 .and. k >= 1 .and. k <= mz-2 .and. &
                        j >= 1 .and. j <= my-2) then
                        x_sol(i,k,j) = sin(12.9898_c_double * real(i, c_double) + &
                                               78.233_c_double * real(k, c_double) + &
                                               37.719_c_double * real(j, c_double))
                    else
                        x_sol(i,k,j) = 0.0_c_double
                    endif
                enddo
            enddo
        enddo
        !$acc update device(x_sol)
    end subroutine probe_random_pattern


    !> Compare the calibrated matrix with a direct distributed application of
    !! A = 2 D o G to the current probe vector.  The relative L2, absolute
    !! maximum, and rank-interface maximum expose a bad colouring/halo
    !! reconstruction before a Krylov failure is attributed to conditioning.
    subroutine probe_compare_operator(domain, direct_div)
        implicit none
        type(domain_t), intent(in) :: domain
        real, intent(in) :: direct_div(domain%ims:domain%ime, domain%kms:domain%kme, domain%jms:domain%jme)
        integer :: i, j, k, ierr
        real(c_double) :: local_stats(4), global_stats(4), matrix_value, direct_value, error_value
        logical :: rank_interface

        call exchange_krylov_halos(x_sol, domain)
        call spmv(x_sol, t_vec)
        !$acc update host(t_vec)

        local_stats = 0.0_c_double
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_s, i_e
                    if (i < 1 .or. i > mx-2 .or. k < 1 .or. k > mz-2 .or. j < 1 .or. j > my-2) cycle
                    matrix_value = t_vec(i,k,j)
                    direct_value = 2.0_c_double * real(direct_div(i,k,j), c_double)
                    error_value = matrix_value - direct_value
                    local_stats(1) = local_stats(1) + error_value * error_value
                    local_stats(2) = local_stats(2) + direct_value * direct_value
                    local_stats(3) = max(local_stats(3), abs(error_value))
                    rank_interface = (.not. domain%west_boundary  .and. i == i_s) .or. &
                                     (.not. domain%east_boundary  .and. i == i_e) .or. &
                                     (.not. domain%south_boundary .and. j == j_s) .or. &
                                     (.not. domain%north_boundary .and. j == j_e)
                    if (rank_interface) local_stats(4) = max(local_stats(4), abs(error_value))
                enddo
            enddo
        enddo
        call MPI_Allreduce(local_stats, global_stats, 2, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(local_stats(3), global_stats(3), 1, MPI_DOUBLE_PRECISION, MPI_MAX, solver_comm, ierr)
        call MPI_Allreduce(local_stats(4), global_stats(4), 1, MPI_DOUBLE_PRECISION, MPI_MAX, solver_comm, ierr)
        if (STD_OUT_PE) then
            write(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') ' Projection operator equivalence: relative L2=', &
                sqrt(global_stats(1) / max(global_stats(2), tiny(1.0_c_double))), ' max abs=', global_stats(3), &
                ' rank-interface max abs=', global_stats(4)
        endif
        call vec_zero(x_sol)
    end subroutine probe_compare_operator


    !> Measure structural properties of the calibrated A = 2 D o G operator.
    !! Three deterministic vector pairs are tested in the Euclidean,
    !! Jacobian-weighted, and inverse-Jacobian-weighted inner products.  All
    !! trial vectors vanish on the solver's identity boundary rows, so the
    !! reported Rayleigh quotients describe the physical projection block.
    subroutine audit_operator_structure(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        integer, parameter :: N_WEIGHTS = 3, N_SAMPLES = 3
        integer :: sample, weight_kind, i, j, k, ierr
        real(c_double) :: defect, rayleigh_x, rayleigh_y, cosine_xy
        real(c_double) :: symmetry_defect(N_WEIGHTS)
        real(c_double) :: rayleigh_min(N_WEIGHTS), rayleigh_max(N_WEIGHTS)
        real(c_double) :: rhs_defect(N_WEIGHTS), rhs_rayleigh_bootstrap(N_WEIGHTS)
        real(c_double) :: rhs_rayleigh_current(N_WEIGHTS), rhs_cosine(N_WEIGHTS)
        real(c_double) :: local_constant_response, global_constant_response
        real(c_double) :: local_constant_diagonal, global_constant_diagonal
        real(c_double) :: local_constant_count, global_constant_count
        real(c_double) :: local_boundary_error, global_boundary_error
        logical :: is_boundary

        symmetry_defect = 0.0_c_double
        rayleigh_min = huge(1.0_c_double)
        rayleigh_max = -huge(1.0_c_double)
        rhs_defect = 0.0_c_double
        rhs_rayleigh_bootstrap = 0.0_c_double
        rhs_rayleigh_current = 0.0_c_double
        rhs_cosine = 0.0_c_double

        do sample = 1, N_SAMPLES
            call fill_audit_vector(p_vec, sample, 0)
            call fill_audit_vector(v_vec, sample, 1)
            call exchange_krylov_halos(p_vec, domain)
            call spmv(p_vec, p_hat)
            call exchange_krylov_halos(v_vec, domain)
            call spmv(v_vec, s_hat)

            do weight_kind = 1, N_WEIGHTS
                call audit_weighted_pair(p_vec, p_hat, v_vec, s_hat, weight_kind, &
                                         defect, rayleigh_x, rayleigh_y, cosine_xy)
                symmetry_defect(weight_kind) = max(symmetry_defect(weight_kind), defect)
                rayleigh_min(weight_kind) = min(rayleigh_min(weight_kind), rayleigh_x, rayleigh_y)
                rayleigh_max(weight_kind) = max(rayleigh_max(weight_kind), rayleigh_x, rayleigh_y)
            enddo
        enddo

        ! A constant-vector response diagnoses row sums and the boundary
        ! identity constraints.  The physical response is normalized by the
        ! diagonal Frobenius scale to remain meaningful across resolutions.
        !$acc parallel loop gang vector collapse(3) present(p_vec)
        do j = j_s-1, j_e+1
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    if (i >= xs .and. i <= xs+xm-1 .and. &
                        k >= zs .and. k <= zs+zm-1 .and. &
                        j >= ys .and. j <= ys+ym-1) then
                        p_vec(i,k,j) = 1.0_c_double
                    else
                        p_vec(i,k,j) = 0.0_c_double
                    endif
                enddo
            enddo
        enddo
        call exchange_krylov_halos(p_vec, domain)
        call spmv(p_vec, p_hat)
        local_constant_response = 0.0_c_double
        local_constant_diagonal = 0.0_c_double
        local_constant_count = 0.0_c_double
        local_boundary_error = 0.0_c_double
        !$acc parallel loop gang vector collapse(3) &
        !$acc reduction(+:local_constant_response,local_constant_diagonal,local_constant_count) &
        !$acc reduction(max:local_boundary_error) present(p_hat, A_coef) private(is_boundary)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    is_boundary = i <= 0 .or. i >= mx-1 .or. j <= 0 .or. j >= my-1 .or. &
                                  k <= 0 .or. k >= mz-1
                    if (is_boundary) then
                        local_boundary_error = max(local_boundary_error, abs(p_hat(i,k,j)-1.0_c_double))
                    else
                        local_constant_response = local_constant_response + p_hat(i,k,j)*p_hat(i,k,j)
                        local_constant_diagonal = local_constant_diagonal + real(A_coef(i,k,j),c_double)**2
                        local_constant_count = local_constant_count + 1.0_c_double
                    endif
                enddo
            enddo
        enddo
        call MPI_Allreduce(local_constant_response, global_constant_response, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(local_constant_diagonal, global_constant_diagonal, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(local_constant_count, global_constant_count, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(local_boundary_error, global_boundary_error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, solver_comm, ierr)

        if (bootstrap_rhs_saved) then
            call exchange_krylov_halos(r_hat, domain)
            call spmv(r_hat, p_hat)
            call exchange_krylov_halos(rhs, domain)
            call spmv(rhs, s_hat)
            do weight_kind = 1, N_WEIGHTS
                call audit_weighted_pair(r_hat, p_hat, rhs, s_hat, weight_kind, &
                                         rhs_defect(weight_kind), &
                                         rhs_rayleigh_bootstrap(weight_kind), &
                                         rhs_rayleigh_current(weight_kind), &
                                         rhs_cosine(weight_kind))
            enddo
        endif

        if (STD_OUT_PE) then
            write(output_unit,'(A)') ' HICAR calibrated wind operator structural audit:'
            do weight_kind = 1, N_WEIGHTS
                write(output_unit,'(A,A,A,ES12.4,A,ES12.4,A,ES12.4)') &
                    '   inner_product=', trim(audit_weight_name(weight_kind)), &
                    ' symmetry_defect=', symmetry_defect(weight_kind), &
                    ' sampled_rayleigh_min=', rayleigh_min(weight_kind), &
                    ' sampled_rayleigh_max=', rayleigh_max(weight_kind)
            enddo
            write(output_unit,'(A,ES12.4,A,ES12.4,A,ES12.4)') &
                '   constant_response_relative_diagonal=', &
                sqrt(global_constant_response/max(global_constant_diagonal,tiny(1.0_c_double))), &
                ' constant_response_rms=', &
                sqrt(global_constant_response/max(global_constant_count,1.0_c_double)), &
                ' boundary_identity_max_error=', global_boundary_error
            if (bootstrap_rhs_saved) then
                do weight_kind = 1, N_WEIGHTS
                    write(output_unit,'(A,A,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
                        '   rhs_inner_product=', trim(audit_weight_name(weight_kind)), &
                        ' adjoint_defect=', rhs_defect(weight_kind), &
                        ' bootstrap_rayleigh=', rhs_rayleigh_bootstrap(weight_kind), &
                        ' calibrated_rayleigh=', rhs_rayleigh_current(weight_kind), &
                        ' rhs_cosine=', rhs_cosine(weight_kind)
                enddo
            else
                write(output_unit,'(A)') '   bootstrap RHS unavailable; RHS comparison skipped.'
            endif
            flush(output_unit)
        endif
    end subroutine audit_operator_structure


    !> Deterministic smooth-plus-oscillatory trial field for distributed
    !! bilinear-form tests.  Boundary rows and halos are zero by construction.
    subroutine fill_audit_vector(v, sample, family)
        implicit none
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: v
        integer, intent(in) :: sample, family
        integer :: i, j, k
        real(c_double) :: phase

        phase = real(7*sample + 13*family, c_double)
        !$acc parallel loop gang vector collapse(3) present(v)
        do j = j_s-1, j_e+1
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    if (i >= 1 .and. i <= mx-2 .and. k >= 1 .and. k <= mz-2 .and. &
                        j >= 1 .and. j <= my-2) then
                        v(i,k,j) = sin((0.011_c_double + 0.001_c_double*sample) * real(i,c_double) + &
                                       (0.037_c_double + 0.002_c_double*family) * real(k,c_double) + &
                                       (0.017_c_double + 0.001_c_double*sample) * real(j,c_double) + phase) + &
                                     0.25_c_double*cos(0.173_c_double*real(i,c_double) - &
                                                           0.119_c_double*real(k,c_double) + &
                                                           0.071_c_double*real(j,c_double) + 0.5_c_double*phase)
                    else
                        v(i,k,j) = 0.0_c_double
                    endif
                enddo
            enddo
        enddo
    end subroutine fill_audit_vector


    !> Bilinear symmetry, Rayleigh, and correlation statistics under a chosen
    !! diagonal scalar-space weight.  weight_kind 1/2/3 selects I/J/J^{-1}.
    subroutine audit_weighted_pair(x, ax, y, ay, weight_kind, defect, rayleigh_x, rayleigh_y, cosine_xy)
        implicit none
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in) :: x, ax, y, ay
        integer, intent(in) :: weight_kind
        real(c_double), intent(out) :: defect, rayleigh_x, rayleigh_y, cosine_xy
        real(c_double) :: local_stats(9), global_stats(9), weight
        real(c_double) :: denom
        integer :: i, j, k, ierr

        local_stats = 0.0_c_double
        !$acc parallel loop gang vector collapse(3) &
        !$acc reduction(+:local_stats(1),local_stats(2),local_stats(3),local_stats(4), &
        !$acc             local_stats(5),local_stats(6),local_stats(7),local_stats(8),local_stats(9)) &
        !$acc present(x, ax, y, ay, jaco) private(weight)
        do j = max(ys,1), min(ys+ym-1,my-2)
            do k = 1, mz - 2
                do i = max(xs,1), min(xs+xm-1,mx-2)
                    select case (weight_kind)
                    case (2)
                        weight = max(abs(real(jaco(i,k,j),c_double)), tiny(1.0_c_double))
                    case (3)
                        weight = 1.0_c_double / max(abs(real(jaco(i,k,j),c_double)), tiny(1.0_c_double))
                    case default
                        weight = 1.0_c_double
                    end select
                    local_stats(1) = local_stats(1) + weight*x(i,k,j)*ay(i,k,j)
                    local_stats(2) = local_stats(2) + weight*ax(i,k,j)*y(i,k,j)
                    local_stats(3) = local_stats(3) + weight*x(i,k,j)*x(i,k,j)
                    local_stats(4) = local_stats(4) + weight*y(i,k,j)*y(i,k,j)
                    local_stats(5) = local_stats(5) + weight*ax(i,k,j)*ax(i,k,j)
                    local_stats(6) = local_stats(6) + weight*ay(i,k,j)*ay(i,k,j)
                    local_stats(7) = local_stats(7) + weight*x(i,k,j)*ax(i,k,j)
                    local_stats(8) = local_stats(8) + weight*y(i,k,j)*ay(i,k,j)
                    local_stats(9) = local_stats(9) + weight*x(i,k,j)*y(i,k,j)
                enddo
            enddo
        enddo
        call MPI_Allreduce(local_stats, global_stats, 9, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        denom = sqrt(max(global_stats(3)*global_stats(6),0.0_c_double)) + &
                sqrt(max(global_stats(5)*global_stats(4),0.0_c_double))
        defect = abs(global_stats(1)-global_stats(2)) / max(denom,tiny(1.0_c_double))
        rayleigh_x = global_stats(7) / max(global_stats(3),tiny(1.0_c_double))
        rayleigh_y = global_stats(8) / max(global_stats(4),tiny(1.0_c_double))
        cosine_xy = global_stats(9) / &
                    max(sqrt(max(global_stats(3)*global_stats(4),0.0_c_double)),tiny(1.0_c_double))
    end subroutine audit_weighted_pair


    !> Export the first calibrated Arnoldi relation and RHS projections.  The
    !! small CSV is written by rank zero only; all vector inner products are
    !! globally reduced first.  H_raw represents the right-preconditioned
    !! operator A M^{-1}, which is the operator whose restart behaviour the
    !! recycling solver must address.  Coarse signed spatial sums for each
    !! basis vector are also exported.  Combining them with a Ritz coefficient
    !! vector reconstructs that mode's global 16x12x8 mean-field envelope.
    subroutine write_krylov_audit(h_raw, used, v_basis, current_residual)
        implicit none
        real(c_double), intent(in) :: h_raw(FGMRES_RESTART+1,FGMRES_RESTART)
        integer, intent(in) :: used
        real(c_double), dimension(i_s-1:, k_s-1:, j_s-1:, :), &
            intent(in) :: v_basis
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), &
            intent(in) :: current_residual
        real(c_double) :: bootstrap_local(FGMRES_RESTART), bootstrap_global(FGMRES_RESTART)
        real(c_double) :: current_local(FGMRES_RESTART), current_global(FGMRES_RESTART)
        real(c_double), allocatable :: spatial_local(:,:), spatial_global(:,:)
        real(c_double), allocatable :: cell_count_local(:), cell_count_global(:)
        integer :: i, j, k, basis, bin_index, bx, by, bz, n_spatial_bins
        integer :: ierr, audit_unit, io_status

        bootstrap_local = 0.0_c_double
        bootstrap_global = 0.0_c_double
        current_local = 0.0_c_double
        current_global = 0.0_c_double
        do i = 1, used
            if (bootstrap_rhs_saved) then
                call vec_dot_local(v_basis(:,:,:,i), r_hat, bootstrap_local(i))
            endif
            call vec_dot_local(v_basis(:,:,:,i), current_residual, current_local(i))
        enddo
        call MPI_Allreduce(bootstrap_local, bootstrap_global, used, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(current_local, current_global, used, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)

        n_spatial_bins = AUDIT_X_BINS * AUDIT_Y_BINS * AUDIT_Z_BINS
        allocate(spatial_local(used,n_spatial_bins), spatial_global(used,n_spatial_bins))
        allocate(cell_count_local(n_spatial_bins), cell_count_global(n_spatial_bins))
        spatial_local = 0.0_c_double
        spatial_global = 0.0_c_double
        cell_count_local = 0.0_c_double
        cell_count_global = 0.0_c_double

        ! Arnoldi vectors are device-resident during the solve.  This one-time
        ! audit transfer is bounded by the existing restart allocation and
        ! happens only after the first cycle has completed.
        !$acc update host(v_basis)
        do j = max(ys,1), min(ys+ym-1,my-2)
            by = min(AUDIT_Y_BINS-1, (j*AUDIT_Y_BINS)/max(my,1))
            do k = 1, mz-2
                bz = min(AUDIT_Z_BINS-1, (k*AUDIT_Z_BINS)/max(mz,1))
                do i = max(xs,1), min(xs+xm-1,mx-2)
                    bx = min(AUDIT_X_BINS-1, (i*AUDIT_X_BINS)/max(mx,1))
                    bin_index = 1 + bx + AUDIT_X_BINS*(by + AUDIT_Y_BINS*bz)
                    cell_count_local(bin_index) = cell_count_local(bin_index) + 1.0_c_double
                enddo
            enddo
        enddo
        do basis = 1, used
            do j = max(ys,1), min(ys+ym-1,my-2)
                by = min(AUDIT_Y_BINS-1, (j*AUDIT_Y_BINS)/max(my,1))
                do k = 1, mz-2
                    bz = min(AUDIT_Z_BINS-1, (k*AUDIT_Z_BINS)/max(mz,1))
                    do i = max(xs,1), min(xs+xm-1,mx-2)
                        bx = min(AUDIT_X_BINS-1, (i*AUDIT_X_BINS)/max(mx,1))
                        bin_index = 1 + bx + AUDIT_X_BINS*(by + AUDIT_Y_BINS*bz)
                        spatial_local(basis,bin_index) = spatial_local(basis,bin_index) + v_basis(i,k,j,basis)
                    enddo
                enddo
            enddo
        enddo
        call MPI_Allreduce(cell_count_local, cell_count_global, n_spatial_bins, &
                           MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(spatial_local, spatial_global, used*n_spatial_bins, &
                           MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)

        if (solver_rank == 0) then
            open(newunit=audit_unit, file=trim(operator_audit_file), status='replace', &
                 action='write', iostat=io_status)
            if (io_status /= 0) then
                write(output_unit,'(A,A,A,I0)') ' HICAR wind audit could not open ', &
                    trim(operator_audit_file), ' iostat=', io_status
                flush(output_unit)
                deallocate(spatial_local, spatial_global, cell_count_local, cell_count_global)
                return
            endif
            write(audit_unit,'(A)') 'record,name,index_i,index_j,value'
            write(audit_unit,'(A,I0)') 'metadata,arnoldi_dimension,,,', used
            write(audit_unit,'(A,I0)') 'metadata,restart,,,', FGMRES_RESTART
            write(audit_unit,'(A,I0)') 'metadata,bootstrap_rhs_saved,,,', merge(1,0,bootstrap_rhs_saved)
            write(audit_unit,'(A,I0)') 'metadata,spatial_x_bins,,,', AUDIT_X_BINS
            write(audit_unit,'(A,I0)') 'metadata,spatial_y_bins,,,', AUDIT_Y_BINS
            write(audit_unit,'(A,I0)') 'metadata,spatial_z_bins,,,', AUDIT_Z_BINS
            do j = 1, used
                do i = 1, j+1
                    write(audit_unit,'(A,I0,A,I0,A,ES24.16)') 'hessenberg,H,', i, ',', j, ',', h_raw(i,j)
                enddo
            enddo
            do i = 1, used
                write(audit_unit,'(A,I0,A,ES24.16)') 'projection,bootstrap,', i, ',,', bootstrap_global(i)
                write(audit_unit,'(A,I0,A,ES24.16)') 'projection,current,', i, ',,', current_global(i)
            enddo
            do bin_index = 1, n_spatial_bins
                write(audit_unit,'(A,I0,A,ES24.16)') 'spatial,cell_count,0,', bin_index, ',', &
                    cell_count_global(bin_index)
            enddo
            do basis = 1, used
                do bin_index = 1, n_spatial_bins
                    write(audit_unit,'(A,I0,A,I0,A,ES24.16)') 'spatial,basis_sum,', basis, ',', &
                        bin_index, ',', spatial_global(basis,bin_index)
                enddo
            enddo
            close(audit_unit)
            write(output_unit,'(A,A,A,I0)') ' HICAR wind Krylov audit written to ', &
                trim(operator_audit_file), ' with Arnoldi dimension ', used
            flush(output_unit)
        endif
        deallocate(spatial_local, spatial_global, cell_count_local, cell_count_global)
    end subroutine write_krylov_audit


    function audit_weight_name(weight_kind) result(name)
        implicit none
        integer, intent(in) :: weight_kind
        character(len=16) :: name
        select case (weight_kind)
        case (2)
            name = 'jacobian'
        case (3)
            name = 'inverse_jacobian'
        case default
            name = 'euclidean'
        end select
    end function audit_weight_name


    subroutine setup_multilevel_preconditioner(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        real(c_double), allocatable :: test_x(:,:,:), direct_ax(:,:,:), host_stencil_ax(:,:,:)
        real(c_double), allocatable :: host_halo_reference(:,:,:)
        real(c_double) :: local_error, global_error, local_reference, global_reference
        real(c_double) :: relative_error, host_relative_error, halo_relative_error, minimum_pivot
        integer :: i, j, k, gi, gj, gk, ierr, line_status, deep_status
        integer :: west, east, south, north

        call release_multilevel_preconditioner()
        multilevel_setup_attempted = .true.
        if (.not. multilevel_requested .or. .not. operator_probed) return
        if (solver_rank == 0) then
            write(output_unit,'(A,I0,A,I0,A,I0)') &
                ' HICAR multilevel setup start: fine=', mx, 'x', my, 'x', mz
            flush(output_unit)
        endif

        call ml_transfer%init(mx, my, xs, xm, ys, ym, &
                              fix_lateral_boundaries=.true., fix_vertical_boundaries=.true.)
        west  = merge(west_neighbor,  MPI_PROC_NULL, xs > 0)
        east  = merge(east_neighbor,  MPI_PROC_NULL, xs+xm < mx)
        south = merge(south_neighbor, MPI_PROC_NULL, ys > 0)
        north = merge(north_neighbor, MPI_PROC_NULL, ys+ym < my)
        call ml_fine_halo%init(xm, ym, mz, solver_comm, west, east, south, north)
        call ml_coarse_halo%init(ml_transfer%nx_c_local, ml_transfer%ny_c_local, mz, &
                                 solver_comm, west, east, south, north)

        allocate(ml_fine_owned(xm,mz,ym), &
                 ml_fine_residual(0:xm+1,mz,0:ym+1), &
                 ml_fine_weight(0:xm+1,mz,0:ym+1), &
                 ml_coarse_halo_x(0:ml_transfer%nx_c_local+1,mz,0:ml_transfer%ny_c_local+1), &
                 ml_coarse_weight_halo(0:ml_transfer%nx_c_local+1,mz,0:ml_transfer%ny_c_local+1), &
                 ml_coarse_weight(ml_transfer%nx_c_local,mz,ml_transfer%ny_c_local), &
                 ml_coarse_b(ml_transfer%nx_c_local,mz,ml_transfer%ny_c_local), &
                 ml_coarse_x(ml_transfer%nx_c_local,mz,ml_transfer%ny_c_local), &
                 ml_coarse_ax(ml_transfer%nx_c_local,mz,ml_transfer%ny_c_local), &
                 ml_coarse_r(ml_transfer%nx_c_local,mz,ml_transfer%ny_c_local), &
                 ml_coarse_correction(ml_transfer%nx_c_local,mz,ml_transfer%ny_c_local))
        ml_fine_owned = 0.0_c_double
        ml_fine_residual = 0.0_c_double
        ml_fine_weight = 0.0_c_double
        do j = 1, ym
            gj = ys+j-1
            do k = 1, mz
                gk = zs+k-1
                do i = 1, xm
                    gi = xs+i-1
                    if (gi <= 0 .or. gi >= mx-1 .or. gj <= 0 .or. gj >= my-1 .or. &
                        gk <= 0 .or. gk >= mz-1) then
                        ml_fine_weight(i,k,j) = 1.0_c_double
                    else
                        ml_fine_weight(i,k,j) = 1.0_c_double / &
                            max(abs(real(jaco(gi,gk,gj),c_double)), 1.0e-12_c_double)
                    endif
                enddo
            enddo
        enddo
        call ml_fine_halo%exchange(ml_fine_weight)
        ml_coarse_halo_x = 0.0_c_double
        ml_coarse_weight_halo = 0.0_c_double
        ml_coarse_weight = 0.0_c_double
        ml_coarse_b = 0.0_c_double
        ml_coarse_x = 0.0_c_double
        ml_coarse_ax = 0.0_c_double
        ml_coarse_r = 0.0_c_double
        ml_coarse_correction = 0.0_c_double

        !$acc enter data copyin(ml_fine_weight) &
        !$acc            create(ml_fine_owned, ml_fine_residual, &
        !$acc                   ml_coarse_halo_x, ml_coarse_weight_halo, ml_coarse_weight, &
        !$acc                   ml_coarse_b, ml_coarse_x, &
        !$acc                   ml_coarse_ax, ml_coarse_r, ml_coarse_correction)
        multilevel_arrays_uploaded = .true.
        call ml_transfer%upload_device()
        call ml_coarse_halo%upload_device()
        call ml_transfer%build_owned_coarse_weights_device(ml_fine_weight, ml_coarse_weight)
        !$acc update self(ml_coarse_weight)
        ml_coarse_weight_halo = 0.0_c_double
        ml_coarse_weight_halo(1:ml_transfer%nx_c_local,:,1:ml_transfer%ny_c_local) = ml_coarse_weight
        call ml_coarse_halo%exchange(ml_coarse_weight_halo)
        !$acc update device(ml_coarse_weight_halo)

        if (solver_rank == 0) then
            write(output_unit,'(A,I0,A,I0,A)') ' HICAR multilevel assembling exact ', &
                ml_transfer%nx_c_global, 'x', ml_transfer%ny_c_global, ' coarse R A P stencil'
            flush(output_unit)
        endif
        call assemble_colored_tile_galerkin(ml_transfer%nx_c_global, ml_transfer%ny_c_global, &
            ml_transfer%x_c_first, ml_transfer%y_c_first, ml_transfer%nx_c_local, &
            ml_transfer%ny_c_local, mz, .true., .true., apply_coarse_rap, ml_stencil)
        call ml_line_factor%factorize(ml_stencil, line_status, minimum_pivot)
        if (line_status /= 0) then
            if (solver_rank == 0) write(output_unit,'(A,ES12.4)') &
                ' HICAR multilevel rejected: singular coarse vertical line, pivot=', minimum_pivot
            call release_multilevel_preconditioner()
            multilevel_setup_attempted = .true.
            return
        endif
        call ml_stencil%upload_device()
        call ml_line_factor%upload_device()

        allocate(test_x(ml_transfer%nx_c_local,mz,ml_transfer%ny_c_local), &
                 direct_ax(ml_transfer%nx_c_local,mz,ml_transfer%ny_c_local))
        do j = 1, ml_transfer%ny_c_local
            do k = 1, mz
                do i = 1, ml_transfer%nx_c_local
                    test_x(i,k,j) = sin(0.173_c_double*real( &
                        3*(ml_transfer%x_c_first+i-1)+5*(k-1)+7*(ml_transfer%y_c_first+j-1),c_double))
                    if (ml_transfer%x_c_first+i-1 == 0 .or. &
                        ml_transfer%x_c_first+i-1 == ml_transfer%nx_c_global-1 .or. &
                        ml_transfer%y_c_first+j-1 == 0 .or. &
                        ml_transfer%y_c_first+j-1 == ml_transfer%ny_c_global-1 .or. &
                        k == 1 .or. k == mz) test_x(i,k,j) = 0.0_c_double
                enddo
            enddo
        enddo
        call apply_coarse_rap(test_x, direct_ax)
        ml_coarse_halo_x = 0.0_c_double
        ml_coarse_halo_x(1:ml_transfer%nx_c_local,:,1:ml_transfer%ny_c_local) = test_x
        allocate(host_stencil_ax(ml_transfer%nx_c_local,mz,ml_transfer%ny_c_local), &
                 host_halo_reference(0:ml_transfer%nx_c_local+1,mz,0:ml_transfer%ny_c_local+1))
        call ml_coarse_halo%exchange(ml_coarse_halo_x)
        host_halo_reference = ml_coarse_halo_x
        call ml_stencil%apply_owned(ml_coarse_halo_x, host_stencil_ax)
        local_error = sum((host_stencil_ax-direct_ax)**2)
        local_reference = sum(direct_ax**2)
        call MPI_Allreduce(local_error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(local_reference, global_reference, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        host_relative_error = sqrt(global_error/max(global_reference,tiny(1.0_c_double)))
        !$acc update device(ml_coarse_halo_x)
        call ml_coarse_halo%exchange_device(ml_coarse_halo_x)
        !$acc update self(ml_coarse_halo_x)
        local_error = sum((ml_coarse_halo_x-host_halo_reference)**2)
        local_reference = sum(host_halo_reference**2)
        call MPI_Allreduce(local_error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(local_reference, global_reference, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        halo_relative_error = sqrt(global_error/max(global_reference,tiny(1.0_c_double)))
        call ml_stencil%apply_owned_device(ml_coarse_halo_x, ml_coarse_ax)
        !$acc update self(ml_coarse_ax)
        local_error = sum((ml_coarse_ax-direct_ax)**2)
        local_reference = sum(direct_ax**2)
        call MPI_Allreduce(local_error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(local_reference, global_reference, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        relative_error = sqrt(global_error/max(global_reference,tiny(1.0_c_double)))
        deallocate(test_x, direct_ax, host_stencil_ax, host_halo_reference)
        if (solver_rank == 0) then
            write(output_unit,'(A,ES12.4,A,ES12.4,A,ES12.4)') &
                ' HICAR multilevel R A P verification: host_stencil=', host_relative_error, &
                ' device_halo=', halo_relative_error, ' device_stencil=', relative_error
            flush(output_unit)
        endif
        if (max(host_relative_error, halo_relative_error, relative_error) > 2.0e-11_c_double) then
            if (solver_rank == 0) write(output_unit,'(A)') &
                ' HICAR multilevel rejected: coarse R A P or device transport verification failed'
            call release_multilevel_preconditioner()
            multilevel_setup_attempted = .true.
            return
        endif

        call setup_recursive_multilevel_levels(deep_status)
        if (deep_status /= 0) then
            call release_multilevel_preconditioner()
            multilevel_setup_attempted = .true.
            return
        endif

        multilevel_ready = .true.
        if (solver_rank == 0) then
            write(output_unit,'(A,I0,A,I0,A,ES12.4,A,ES12.4)') &
                ' HICAR Petrov-Galerkin level ready: global coarse=', ml_transfer%nx_c_global, 'x', &
                ml_transfer%ny_c_global, ' RAP_error=', relative_error, ' min_line_pivot=', minimum_pivot
            flush(output_unit)
            write(output_unit,'(A,I0)') ' HICAR exact Galerkin hierarchy ready: total coarse levels=', &
                1+ml_deep_count
            flush(output_unit)
        endif

    contains

        subroutine apply_coarse_rap(coarse, coarse_ax)
            real(c_double), intent(in) :: coarse(:,:,:)
            real(c_double), intent(out) :: coarse_ax(:,:,:)
            integer :: ii, jj, kk

            ml_coarse_halo_x = 0.0_c_double
            ml_coarse_halo_x(1:ml_transfer%nx_c_local,:,1:ml_transfer%ny_c_local) = coarse
            !$acc update device(ml_coarse_halo_x)
            call ml_coarse_halo%exchange_device(ml_coarse_halo_x)
            call ml_transfer%prolong_owned_device(ml_coarse_halo_x, ml_fine_owned)
            ! Reuse the solver's persistent Krylov work arrays here.  Besides
            ! avoiding two full-size allocations, these arrays are the exact
            ! allocation class exercised by the production NCCL halo path.
            !$acc parallel loop gang vector collapse(3) present(x_sol)
            do jj = j_s-1, j_e+1
                do kk = k_s-1, k_e+1
                    do ii = i_s-1, i_e+1
                        x_sol(ii,kk,jj) = 0.0_c_double
                    enddo
                enddo
            enddo
            !$acc parallel loop gang vector collapse(3) present(x_sol,ml_fine_owned)
            do jj = ys, ys+ym-1
                do kk = zs, zs+zm-1
                    do ii = xs, xs+xm-1
                        x_sol(ii,kk,jj) = ml_fine_owned(ii-xs+1,kk-zs+1,jj-ys+1)
                    enddo
                enddo
            enddo
            call exchange_krylov_halos(x_sol, domain)
            call spmv(x_sol, t_vec)
            call exchange_krylov_halos(t_vec, domain)
            ! The adjoint restriction has a 3x3 horizontal footprint.  The
            ! production exchange posts all four faces concurrently, so a
            ! second pass is required to propagate newly received face data
            ! into diagonal corners (the coarse halo exchanger instead does
            ! this explicitly as x-then-y).
            call exchange_krylov_halos(t_vec, domain)
            call copy_solver_to_multilevel_halo(t_vec)
            call ml_transfer%restrict_owned_adjoint_device(ml_fine_residual, ml_fine_weight, &
                                                            ml_coarse_weight, ml_coarse_r)
            !$acc update self(ml_coarse_r)
            coarse_ax = ml_coarse_r
        end subroutine apply_coarse_rap

    end subroutine setup_multilevel_preconditioner


    subroutine setup_recursive_multilevel_levels(status)
        implicit none
        integer, intent(out) :: status
        integer :: q, ierr, local_min, global_min, line_status
        integer :: parent_nx_global, parent_ny_global, parent_x_first, parent_y_first
        integer :: parent_nx_local, parent_ny_local
        integer :: child_x_first, child_y_first, child_nx_local, child_ny_local
        integer :: west, east, south, north
        real(c_double) :: minimum_pivot

        status = 0
        ml_deep_count = 0
        ml_assembly_child = 0
        if (allocated(ml_deep)) call release_recursive_multilevel_levels()
        allocate(ml_deep(MAX_ML_DEEP_LEVELS))

        do q = 1, MAX_ML_DEEP_LEVELS
            if (q == 1) then
                parent_nx_global = ml_transfer%nx_c_global
                parent_ny_global = ml_transfer%ny_c_global
                parent_x_first = ml_transfer%x_c_first
                parent_y_first = ml_transfer%y_c_first
                parent_nx_local = ml_transfer%nx_c_local
                parent_ny_local = ml_transfer%ny_c_local
            else
                parent_nx_global = ml_deep(q-1)%stencil%nx_global
                parent_ny_global = ml_deep(q-1)%stencil%ny_global
                parent_x_first = ml_deep(q-1)%stencil%x_first
                parent_y_first = ml_deep(q-1)%stencil%y_first
                parent_nx_local = ml_deep(q-1)%stencil%nx
                parent_ny_local = ml_deep(q-1)%stencil%ny
            endif

            ! Retain every rank until agglomeration is introduced.  Stop at
            ! the last grid for which the existing decomposition is nonempty.
            if (parent_nx_global <= 3 .or. parent_ny_global <= 3) exit
            call owned_coarse_interval(parent_nx_global, parent_x_first, parent_nx_local, &
                                       child_x_first, child_nx_local)
            call owned_coarse_interval(parent_ny_global, parent_y_first, parent_ny_local, &
                                       child_y_first, child_ny_local)
            local_min = min(child_nx_local, child_ny_local)
            call MPI_Allreduce(local_min, global_min, 1, MPI_INTEGER, MPI_MIN, solver_comm, ierr)
            if (global_min < 1) exit
            call ml_deep(q)%transfer_from_parent%init(parent_nx_global, parent_ny_global, &
                parent_x_first, parent_nx_local, parent_y_first, parent_ny_local, &
                fix_lateral_boundaries=.true., fix_vertical_boundaries=.true.)

            west = merge(west_neighbor, MPI_PROC_NULL, &
                ml_deep(q)%transfer_from_parent%x_c_first > 0)
            east = merge(east_neighbor, MPI_PROC_NULL, &
                ml_deep(q)%transfer_from_parent%x_c_first + &
                ml_deep(q)%transfer_from_parent%nx_c_local < &
                ml_deep(q)%transfer_from_parent%nx_c_global)
            south = merge(south_neighbor, MPI_PROC_NULL, &
                ml_deep(q)%transfer_from_parent%y_c_first > 0)
            north = merge(north_neighbor, MPI_PROC_NULL, &
                ml_deep(q)%transfer_from_parent%y_c_first + &
                ml_deep(q)%transfer_from_parent%ny_c_local < &
                ml_deep(q)%transfer_from_parent%ny_c_global)
            call ml_deep(q)%halo%init(ml_deep(q)%transfer_from_parent%nx_c_local, &
                ml_deep(q)%transfer_from_parent%ny_c_local, mz, solver_comm, west, east, south, north)
            call allocate_recursive_level_arrays(q)
            call ml_deep(q)%transfer_from_parent%upload_device()
            call ml_deep(q)%halo%upload_device()

            if (q == 1) then
                call ml_deep(q)%transfer_from_parent%build_owned_coarse_weights_device( &
                    ml_coarse_weight_halo, ml_deep(q)%weight)
            else
                call ml_deep(q)%transfer_from_parent%build_owned_coarse_weights_device( &
                    ml_deep(q-1)%weight_halo, ml_deep(q)%weight)
            endif
            call copy_owned_to_halo_device(ml_deep(q)%weight, ml_deep(q)%weight_halo)
            call ml_deep(q)%halo%exchange_device(ml_deep(q)%weight_halo)

            ml_assembly_child = q
            if (solver_rank == 0) then
                write(output_unit,'(A,I0,A,I0,A,I0,A)') ' HICAR multilevel assembling level ', q+1, &
                    ' exact ', ml_deep(q)%transfer_from_parent%nx_c_global, 'x', &
                    ml_deep(q)%transfer_from_parent%ny_c_global, ' R A P stencil'
                flush(output_unit)
            endif
            call assemble_colored_tile_galerkin( &
                ml_deep(q)%transfer_from_parent%nx_c_global, &
                ml_deep(q)%transfer_from_parent%ny_c_global, &
                ml_deep(q)%transfer_from_parent%x_c_first, &
                ml_deep(q)%transfer_from_parent%y_c_first, &
                ml_deep(q)%transfer_from_parent%nx_c_local, &
                ml_deep(q)%transfer_from_parent%ny_c_local, mz, .true., .true., &
                apply_recursive_coarse_rap, ml_deep(q)%stencil)
            call ml_deep(q)%line_factor%factorize(ml_deep(q)%stencil, line_status, minimum_pivot)
            if (line_status /= 0) then
                if (solver_rank == 0) write(output_unit,'(A,I0,A,ES12.4)') &
                    ' HICAR multilevel rejected: singular line on level ', q+1, &
                    ', pivot=', minimum_pivot
                status = 1
                return
            endif
            call ml_deep(q)%stencil%upload_device()
            call ml_deep(q)%line_factor%upload_device()
            call verify_recursive_level(q, status)
            if (status /= 0) return
            ml_deep_count = q
            if (solver_rank == 0) then
                write(output_unit,'(A,I0,A,I0,A,I0,A,ES12.4)') &
                    ' HICAR exact Galerkin level ', q+1, ' ready: global=', &
                    ml_deep(q)%stencil%nx_global, 'x', ml_deep(q)%stencil%ny_global, &
                    ' min_line_pivot=', minimum_pivot
                flush(output_unit)
            endif
            if (ml_deep(q)%stencil%nx_global <= 4 .or. ml_deep(q)%stencil%ny_global <= 4) exit
        enddo
        ml_assembly_child = 0
    end subroutine setup_recursive_multilevel_levels


    subroutine allocate_recursive_level_arrays(q)
        implicit none
        integer, intent(in) :: q
        integer :: nx, ny

        nx = ml_deep(q)%transfer_from_parent%nx_c_local
        ny = ml_deep(q)%transfer_from_parent%ny_c_local
        allocate(ml_deep(q)%weight(nx,mz,ny), ml_deep(q)%weight_halo(0:nx+1,mz,0:ny+1), &
                 ml_deep(q)%halo_x(0:nx+1,mz,0:ny+1), ml_deep(q)%b(nx,mz,ny), &
                 ml_deep(q)%x(nx,mz,ny), ml_deep(q)%ax(nx,mz,ny), ml_deep(q)%r(nx,mz,ny), &
                 ml_deep(q)%correction(nx,mz,ny))
        ml_deep(q)%weight = 0.0_c_double
        ml_deep(q)%weight_halo = 0.0_c_double
        ml_deep(q)%halo_x = 0.0_c_double
        ml_deep(q)%b = 0.0_c_double
        ml_deep(q)%x = 0.0_c_double
        ml_deep(q)%ax = 0.0_c_double
        ml_deep(q)%r = 0.0_c_double
        ml_deep(q)%correction = 0.0_c_double
        !$acc enter data create(ml_deep(q)%weight, ml_deep(q)%weight_halo, &
        !$acc                   ml_deep(q)%halo_x, ml_deep(q)%b, ml_deep(q)%x, &
        !$acc                   ml_deep(q)%ax, ml_deep(q)%r, ml_deep(q)%correction)
        ml_deep(q)%arrays_uploaded = .true.
    end subroutine allocate_recursive_level_arrays


    subroutine apply_recursive_coarse_rap(coarse, coarse_ax)
        implicit none
        real(c_double), intent(in) :: coarse(:,:,:)
        real(c_double), intent(out) :: coarse_ax(:,:,:)
        integer :: q

        q = ml_assembly_child
        if (q < 1 .or. q > size(ml_deep)) error stop 'invalid recursive R A P level'
        ml_deep(q)%halo_x = 0.0_c_double
        ml_deep(q)%halo_x(1:size(coarse,1),:,1:size(coarse,3)) = coarse
        !$acc update device(ml_deep(q)%halo_x)
        call ml_deep(q)%halo%exchange_device(ml_deep(q)%halo_x)
        if (q == 1) then
            call ml_deep(q)%transfer_from_parent%prolong_owned_device(ml_deep(q)%halo_x, ml_coarse_b)
            call copy_owned_to_halo_device(ml_coarse_b, ml_coarse_halo_x)
            call ml_coarse_halo%exchange_device(ml_coarse_halo_x)
            call ml_stencil%apply_owned_device(ml_coarse_halo_x, ml_coarse_ax)
            call copy_owned_to_halo_device(ml_coarse_ax, ml_coarse_halo_x)
            call ml_coarse_halo%exchange_device(ml_coarse_halo_x)
            call ml_deep(q)%transfer_from_parent%restrict_owned_adjoint_device( &
                ml_coarse_halo_x, ml_coarse_weight_halo, ml_deep(q)%weight, ml_deep(q)%r)
        else
            call ml_deep(q)%transfer_from_parent%prolong_owned_device(ml_deep(q)%halo_x, ml_deep(q-1)%b)
            call copy_owned_to_halo_device(ml_deep(q-1)%b, ml_deep(q-1)%halo_x)
            call ml_deep(q-1)%halo%exchange_device(ml_deep(q-1)%halo_x)
            call ml_deep(q-1)%stencil%apply_owned_device(ml_deep(q-1)%halo_x, ml_deep(q-1)%ax)
            call copy_owned_to_halo_device(ml_deep(q-1)%ax, ml_deep(q-1)%halo_x)
            call ml_deep(q-1)%halo%exchange_device(ml_deep(q-1)%halo_x)
            call ml_deep(q)%transfer_from_parent%restrict_owned_adjoint_device( &
                ml_deep(q-1)%halo_x, ml_deep(q-1)%weight_halo, &
                ml_deep(q)%weight, ml_deep(q)%r)
        endif
        !$acc update self(ml_deep(q)%r)
        coarse_ax = ml_deep(q)%r
    end subroutine apply_recursive_coarse_rap


    subroutine verify_recursive_level(q, status)
        implicit none
        integer, intent(in) :: q
        integer, intent(out) :: status
        real(c_double), allocatable :: test_x(:,:,:), direct_ax(:,:,:), host_ax(:,:,:)
        real(c_double), allocatable :: host_halo(:,:,:)
        real(c_double) :: local_error, global_error, local_reference, global_reference
        real(c_double) :: host_error, halo_error, device_error
        integer :: i, j, k, gi, gj, ierr, nx, ny

        status = 0
        nx = ml_deep(q)%stencil%nx
        ny = ml_deep(q)%stencil%ny
        allocate(test_x(nx,mz,ny), direct_ax(nx,mz,ny), host_ax(nx,mz,ny), &
                 host_halo(0:nx+1,mz,0:ny+1))
        do j = 1, ny
            gj = ml_deep(q)%stencil%y_first+j-1
            do k = 1, mz
                do i = 1, nx
                    gi = ml_deep(q)%stencil%x_first+i-1
                    test_x(i,k,j) = sin(0.137_c_double*real(3*gi+5*(k-1)+7*gj,c_double))
                    if (gi == 0 .or. gi == ml_deep(q)%stencil%nx_global-1 .or. &
                        gj == 0 .or. gj == ml_deep(q)%stencil%ny_global-1 .or. &
                        k == 1 .or. k == mz) test_x(i,k,j) = 0.0_c_double
                enddo
            enddo
        enddo
        ml_assembly_child = q
        call apply_recursive_coarse_rap(test_x, direct_ax)
        ml_deep(q)%halo_x = 0.0_c_double
        ml_deep(q)%halo_x(1:nx,:,1:ny) = test_x
        call ml_deep(q)%halo%exchange(ml_deep(q)%halo_x)
        host_halo = ml_deep(q)%halo_x
        call ml_deep(q)%stencil%apply_owned(ml_deep(q)%halo_x, host_ax)
        local_error = sum((host_ax-direct_ax)**2)
        local_reference = sum(direct_ax**2)
        call MPI_Allreduce(local_error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(local_reference, global_reference, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        host_error = sqrt(global_error/max(global_reference,tiny(1.0_c_double)))
        !$acc update device(ml_deep(q)%halo_x)
        call ml_deep(q)%halo%exchange_device(ml_deep(q)%halo_x)
        !$acc update self(ml_deep(q)%halo_x)
        local_error = sum((ml_deep(q)%halo_x-host_halo)**2)
        local_reference = sum(host_halo**2)
        call MPI_Allreduce(local_error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(local_reference, global_reference, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        halo_error = sqrt(global_error/max(global_reference,tiny(1.0_c_double)))
        call ml_deep(q)%stencil%apply_owned_device(ml_deep(q)%halo_x, ml_deep(q)%ax)
        !$acc update self(ml_deep(q)%ax)
        local_error = sum((ml_deep(q)%ax-direct_ax)**2)
        local_reference = sum(direct_ax**2)
        call MPI_Allreduce(local_error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        call MPI_Allreduce(local_reference, global_reference, 1, MPI_DOUBLE_PRECISION, MPI_SUM, solver_comm, ierr)
        device_error = sqrt(global_error/max(global_reference,tiny(1.0_c_double)))
        if (solver_rank == 0) then
            write(output_unit,'(A,I0,A,ES12.4,A,ES12.4,A,ES12.4)') &
                ' HICAR recursive R A P verification level ', q+1, ': host_stencil=', host_error, &
                ' device_halo=', halo_error, ' device_stencil=', device_error
            flush(output_unit)
        endif
        if (max(host_error, halo_error, device_error) > 2.0e-11_c_double) status = 1
        deallocate(test_x, direct_ax, host_ax, host_halo)
    end subroutine verify_recursive_level


    subroutine copy_owned_to_halo_device(owned, halo)
        implicit none
        real(c_double), intent(in) :: owned(:,:,:)
        real(c_double), intent(inout) :: halo(0:,:,0:)
        integer :: i, j, k

        !$acc parallel loop gang vector collapse(3) present(owned,halo)
        do j = 1, size(owned,3)
            do k = 1, size(owned,2)
                do i = 1, size(owned,1)
                    halo(i,k,j) = owned(i,k,j)
                enddo
            enddo
        enddo
    end subroutine copy_owned_to_halo_device


    subroutine release_multilevel_preconditioner()
        implicit none

        call release_recursive_multilevel_levels()
        if (multilevel_arrays_uploaded) then
            !$acc exit data delete(ml_fine_owned, ml_fine_residual, ml_fine_weight, &
            !$acc                  ml_coarse_halo_x, ml_coarse_weight_halo, ml_coarse_weight, &
            !$acc                  ml_coarse_b, ml_coarse_x, ml_coarse_ax, ml_coarse_r, &
            !$acc                  ml_coarse_correction)
        endif
        multilevel_arrays_uploaded = .false.
        call ml_line_factor%release()
        call ml_stencil%release()
        call ml_transfer%release()
        call ml_fine_halo%release()
        call ml_coarse_halo%release()
        if (allocated(ml_fine_owned)) deallocate(ml_fine_owned)
        if (allocated(ml_fine_residual)) deallocate(ml_fine_residual)
        if (allocated(ml_fine_weight)) deallocate(ml_fine_weight)
        if (allocated(ml_coarse_halo_x)) deallocate(ml_coarse_halo_x)
        if (allocated(ml_coarse_weight_halo)) deallocate(ml_coarse_weight_halo)
        if (allocated(ml_coarse_weight)) deallocate(ml_coarse_weight)
        if (allocated(ml_coarse_b)) deallocate(ml_coarse_b)
        if (allocated(ml_coarse_x)) deallocate(ml_coarse_x)
        if (allocated(ml_coarse_ax)) deallocate(ml_coarse_ax)
        if (allocated(ml_coarse_r)) deallocate(ml_coarse_r)
        if (allocated(ml_coarse_correction)) deallocate(ml_coarse_correction)
        multilevel_ready = .false.
    end subroutine release_multilevel_preconditioner


    subroutine release_recursive_multilevel_levels()
        implicit none
        integer :: q

        if (.not. allocated(ml_deep)) then
            ml_deep_count = 0
            ml_assembly_child = 0
            return
        endif
        do q = 1, size(ml_deep)
            if (ml_deep(q)%arrays_uploaded) then
                !$acc exit data delete(ml_deep(q)%weight, ml_deep(q)%weight_halo, &
                !$acc                  ml_deep(q)%halo_x, ml_deep(q)%b, ml_deep(q)%x, &
                !$acc                  ml_deep(q)%ax, ml_deep(q)%r, ml_deep(q)%correction)
            endif
            ml_deep(q)%arrays_uploaded = .false.
            call ml_deep(q)%line_factor%release()
            call ml_deep(q)%stencil%release()
            call ml_deep(q)%transfer_from_parent%release()
            call ml_deep(q)%halo%release()
            if (allocated(ml_deep(q)%weight)) deallocate(ml_deep(q)%weight)
            if (allocated(ml_deep(q)%weight_halo)) deallocate(ml_deep(q)%weight_halo)
            if (allocated(ml_deep(q)%halo_x)) deallocate(ml_deep(q)%halo_x)
            if (allocated(ml_deep(q)%b)) deallocate(ml_deep(q)%b)
            if (allocated(ml_deep(q)%x)) deallocate(ml_deep(q)%x)
            if (allocated(ml_deep(q)%ax)) deallocate(ml_deep(q)%ax)
            if (allocated(ml_deep(q)%r)) deallocate(ml_deep(q)%r)
            if (allocated(ml_deep(q)%correction)) deallocate(ml_deep(q)%correction)
        enddo
        deallocate(ml_deep)
        ml_deep_count = 0
        ml_assembly_child = 0
    end subroutine release_recursive_multilevel_levels


    subroutine copy_solver_to_multilevel_halo(vec)
        implicit none
        real(c_double), dimension(i_s-1:i_e+1,k_s-1:k_e+1,j_s-1:j_e+1), intent(in) :: vec
        integer :: i, j, k, gi, gj, gk

        !$acc parallel loop gang vector collapse(3) present(vec,ml_fine_residual) private(gi,gj,gk)
        do j = 0, ym+1
            do k = 1, mz
                do i = 0, xm+1
                    gi = xs+i-1
                    gk = zs+k-1
                    gj = ys+j-1
                    if (gi >= i_s-1 .and. gi <= i_e+1 .and. &
                        gk >= k_s-1 .and. gk <= k_e+1 .and. &
                        gj >= j_s-1 .and. gj <= j_e+1) then
                        ml_fine_residual(i,k,j) = vec(gi,gk,gj)
                    else
                        ml_fine_residual(i,k,j) = 0.0_c_double
                    endif
                enddo
            enddo
        enddo
    end subroutine copy_solver_to_multilevel_halo


    subroutine apply_multilevel_preconditioner(in_vec, out_vec, domain)
        implicit none
        real(c_double), dimension(i_s-1:i_e+1,k_s-1:k_e+1,j_s-1:j_e+1), intent(in) :: in_vec
        real(c_double), dimension(i_s-1:i_e+1,k_s-1:k_e+1,j_s-1:j_e+1), intent(inout) :: out_vec
        type(domain_t), intent(in) :: domain
        integer :: i, j, k, nxc, nyc
        real(c_double), parameter :: coarse_omega = 0.8_c_double
        real(c_double), parameter :: post_omega = 0.8_c_double

        nxc = ml_transfer%nx_c_local
        nyc = ml_transfer%ny_c_local

        !$acc parallel loop gang vector collapse(3) present(out_vec)
        do j = j_s-1, j_e+1
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    out_vec(i,k,j) = 0.0_c_double
                enddo
            enddo
        enddo
        call apply_vertical_lines(in_vec, out_vec)
        call exchange_krylov_halos(out_vec, domain)
        call spmv(out_vec, prec_res)
        !$acc parallel loop gang vector collapse(3) present(in_vec,prec_res,ml_fine_residual)
        do j = ys, ys+ym-1
            do k = zs, zs+zm-1
                do i = xs, xs+xm-1
                    prec_res(i,k,j) = in_vec(i,k,j)-prec_res(i,k,j)
                enddo
            enddo
        enddo
        call exchange_krylov_halos(prec_res, domain)
        call exchange_krylov_halos(prec_res, domain)
        call copy_solver_to_multilevel_halo(prec_res)
        call ml_transfer%restrict_owned_adjoint_device(ml_fine_residual, ml_fine_weight, &
                                                        ml_coarse_weight, ml_coarse_b)
        call apply_level_one_vcycle(coarse_omega)
        !$acc parallel loop gang vector collapse(3) present(ml_coarse_halo_x,ml_coarse_x)
        do j = 1, nyc
            do k = 1, mz
                do i = 1, nxc
                    ml_coarse_halo_x(i,k,j) = ml_coarse_x(i,k,j)
                enddo
            enddo
        enddo
        call ml_coarse_halo%exchange_device(ml_coarse_halo_x)
        call ml_transfer%prolong_owned_device(ml_coarse_halo_x, ml_fine_owned)
        !$acc parallel loop gang vector collapse(3) present(out_vec,ml_fine_owned)
        do j = ys, ys+ym-1
            do k = zs, zs+zm-1
                do i = xs, xs+xm-1
                    out_vec(i,k,j) = out_vec(i,k,j)+ml_fine_owned(i-xs+1,k-zs+1,j-ys+1)
                enddo
            enddo
        enddo
        call exchange_krylov_halos(out_vec, domain)
        call spmv(out_vec, prec_res)
        !$acc parallel loop gang vector collapse(3) present(in_vec,prec_res)
        do j = ys, ys+ym-1
            do k = zs, zs+zm-1
                do i = xs, xs+xm-1
                    prec_res(i,k,j) = in_vec(i,k,j)-prec_res(i,k,j)
                enddo
            enddo
        enddo
        call apply_vertical_lines(prec_res, t_vec)
        !$acc parallel loop gang vector collapse(3) present(out_vec,t_vec)
        do j = ys, ys+ym-1
            do k = zs, zs+zm-1
                do i = xs, xs+xm-1
                    out_vec(i,k,j) = out_vec(i,k,j)+post_omega*t_vec(i,k,j)
                enddo
            enddo
        enddo
    end subroutine apply_multilevel_preconditioner


    subroutine apply_level_one_vcycle(omega)
        implicit none
        real(c_double), intent(in) :: omega
        integer :: sweep

        call ml_line_factor%apply_device(ml_coarse_b, ml_coarse_x)
        if (ml_deep_count > 0) then
            call compute_level_one_residual()
            call copy_owned_to_halo_device(ml_coarse_r, ml_coarse_halo_x)
            call ml_coarse_halo%exchange_device(ml_coarse_halo_x)
            call ml_deep(1)%transfer_from_parent%restrict_owned_adjoint_device( &
                ml_coarse_halo_x, ml_coarse_weight_halo, ml_deep(1)%weight, ml_deep(1)%b)
            call apply_recursive_vcycle(1, omega)
            call copy_owned_to_halo_device(ml_deep(1)%x, ml_deep(1)%halo_x)
            call ml_deep(1)%halo%exchange_device(ml_deep(1)%halo_x)
            call ml_deep(1)%transfer_from_parent%prolong_owned_device( &
                ml_deep(1)%halo_x, ml_coarse_correction)
            call owned_axpy_device(ml_coarse_x, ml_coarse_correction, 1.0_c_double)
            call compute_level_one_residual()
            call ml_line_factor%apply_device(ml_coarse_r, ml_coarse_correction)
            call owned_axpy_device(ml_coarse_x, ml_coarse_correction, omega)
        else
            do sweep = 2, 4
                call compute_level_one_residual()
                call ml_line_factor%apply_device(ml_coarse_r, ml_coarse_correction)
                call owned_axpy_device(ml_coarse_x, ml_coarse_correction, omega)
            enddo
        endif
    end subroutine apply_level_one_vcycle


    recursive subroutine apply_recursive_vcycle(q, omega)
        implicit none
        integer, intent(in) :: q
        real(c_double), intent(in) :: omega
        integer :: sweep

        call ml_deep(q)%line_factor%apply_device(ml_deep(q)%b, ml_deep(q)%x)
        if (q < ml_deep_count) then
            call compute_recursive_level_residual(q)
            call copy_owned_to_halo_device(ml_deep(q)%r, ml_deep(q)%halo_x)
            call ml_deep(q)%halo%exchange_device(ml_deep(q)%halo_x)
            call ml_deep(q+1)%transfer_from_parent%restrict_owned_adjoint_device( &
                ml_deep(q)%halo_x, ml_deep(q)%weight_halo, &
                ml_deep(q+1)%weight, ml_deep(q+1)%b)
            call apply_recursive_vcycle(q+1, omega)
            call copy_owned_to_halo_device(ml_deep(q+1)%x, ml_deep(q+1)%halo_x)
            call ml_deep(q+1)%halo%exchange_device(ml_deep(q+1)%halo_x)
            call ml_deep(q+1)%transfer_from_parent%prolong_owned_device( &
                ml_deep(q+1)%halo_x, ml_deep(q)%correction)
            call owned_axpy_device(ml_deep(q)%x, ml_deep(q)%correction, 1.0_c_double)
            call compute_recursive_level_residual(q)
            call ml_deep(q)%line_factor%apply_device(ml_deep(q)%r, ml_deep(q)%correction)
            call owned_axpy_device(ml_deep(q)%x, ml_deep(q)%correction, omega)
        else
            do sweep = 2, 4
                call compute_recursive_level_residual(q)
                call ml_deep(q)%line_factor%apply_device(ml_deep(q)%r, ml_deep(q)%correction)
                call owned_axpy_device(ml_deep(q)%x, ml_deep(q)%correction, omega)
            enddo
        endif
    end subroutine apply_recursive_vcycle


    subroutine compute_level_one_residual()
        implicit none

        call copy_owned_to_halo_device(ml_coarse_x, ml_coarse_halo_x)
        call ml_coarse_halo%exchange_device(ml_coarse_halo_x)
        call ml_stencil%apply_owned_device(ml_coarse_halo_x, ml_coarse_ax)
        call owned_residual_device(ml_coarse_b, ml_coarse_ax, ml_coarse_r)
    end subroutine compute_level_one_residual


    subroutine compute_recursive_level_residual(q)
        implicit none
        integer, intent(in) :: q

        call copy_owned_to_halo_device(ml_deep(q)%x, ml_deep(q)%halo_x)
        call ml_deep(q)%halo%exchange_device(ml_deep(q)%halo_x)
        call ml_deep(q)%stencil%apply_owned_device(ml_deep(q)%halo_x, ml_deep(q)%ax)
        call owned_residual_device(ml_deep(q)%b, ml_deep(q)%ax, ml_deep(q)%r)
    end subroutine compute_recursive_level_residual


    subroutine owned_residual_device(b, ax, r)
        implicit none
        real(c_double), intent(in) :: b(:,:,:), ax(:,:,:)
        real(c_double), intent(out) :: r(:,:,:)
        integer :: i, j, k

        !$acc parallel loop gang vector collapse(3) present(b,ax,r)
        do j = 1, size(b,3)
            do k = 1, size(b,2)
                do i = 1, size(b,1)
                    r(i,k,j) = b(i,k,j)-ax(i,k,j)
                enddo
            enddo
        enddo
    end subroutine owned_residual_device


    subroutine owned_axpy_device(x, correction, omega)
        implicit none
        real(c_double), intent(inout) :: x(:,:,:)
        real(c_double), intent(in) :: correction(:,:,:), omega
        integer :: i, j, k

        !$acc parallel loop gang vector collapse(3) present(x,correction)
        do j = 1, size(x,3)
            do k = 1, size(x,2)
                do i = 1, size(x,1)
                    x(i,k,j) = x(i,k,j)+omega*correction(i,k,j)
                enddo
            enddo
        enddo
    end subroutine owned_axpy_device


    !>------------------------------------------------------------
    !! Vertical-line block-Jacobi preconditioner apply with multi-sweep Richardson.
    !!
    !! Applies M^{-1} to in_vec, with M^{-1} approximated by `precond_n_sweeps`
    !! Richardson sweeps using one exact vertical tridiagonal block per (i,j):
    !!   y_0 = 0
    !!   y_{k+1} = y_k + M_line^{-1} (in - A y_k)
    !!
    !! Sweep 1 is an exact column solve. Subsequent sweeps add 1 SpMV and
    !! another column solve. Horizontal rank boundaries remain homogeneous
    !! during the local smoothing operation, as required by block-Jacobi.
    !!
    !! Block-Jacobi convention: no halo exchange between inner sweeps.
    !! Halo cells of out_vec stay at zero (initialised at first solve, never
    !! written by any solver routine) → inner SpMV uses zero halo, equivalent
    !! to homogeneous Dirichlet at the rank boundary. This is the standard
    !! approximation used by AMGX BLOCK_JACOBI and Hypre's PCJacobi.
    !!------------------------------------------------------------
    subroutine apply_precond(in_vec, out_vec, domain)
        implicit none
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: in_vec
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: out_vec
        type(domain_t), intent(in) :: domain
        integer :: i, j, k, sweep

        if (multilevel_ready) then
            call apply_multilevel_preconditioner(in_vec, out_vec, domain)
            return
        endif

        ! Zero entire out_vec INCLUDING halo cells before sweeps. The inner SpMV
        ! in sweep 2+ reads out_vec at halo positions; for block-Jacobi convention
        ! those must be zero (homogeneous Dirichlet at rank boundary). Without
        ! this, halo cells leak previous outer-iter's exchange_krylov_halos data
        ! into the smoother and convergence drifts on multi-rank runs (bug
        ! manifested as slow residual growth at ~0.25%/iter, 2-GPU 2026-05-10).
        ! Halo cells of out_vec stay at zero throughout the precond apply since
        ! sweeps only write to owned cells.
        !$acc parallel loop gang vector collapse(3) present(out_vec)
        do j = j_s-1, j_e+1
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    out_vec(i,k,j) = 0.0_c_double
                enddo
            enddo
        enddo

        ! Sweep 1: y_1 = M_line^{-1} in (starting from y_0 = 0).
        call apply_vertical_lines(in_vec, out_vec)

        ! Sweeps 2..n: y += M_line^{-1} (in - A y). Inner SpMV uses out_vec's
        ! zero halo (no MPI exchange — block-Jacobi convention). t_vec is safe
        ! scratch at both preconditioner call sites and avoids a further large
        ! per-rank allocation.
        do sweep = 2, precond_n_sweeps
            call spmv(out_vec, prec_res)
            !$acc parallel loop gang vector collapse(3) present(in_vec, prec_res)
            do j = ys, ys + ym - 1
                do k = zs, zs + zm - 1
                    do i = xs, xs + xm - 1
                        prec_res(i,k,j) = in_vec(i,k,j) - prec_res(i,k,j)
                    enddo
                enddo
            enddo
            call apply_vertical_lines(prec_res, t_vec)
            !$acc parallel loop gang vector collapse(3) &
            !$acc present(in_vec, out_vec, t_vec)
            do j = ys, ys + ym - 1
                do k = zs, zs + zm - 1
                    do i = xs, xs + xm - 1
                        out_vec(i,k,j) = out_vec(i,k,j) + t_vec(i,k,j)
                    enddo
                enddo
            enddo
        enddo
    end subroutine apply_precond


    !> Apply the factored vertical-line blocks to `in_vec`.
    !! The forward result is stored in out_vec and overwritten in place by
    !! backward substitution.  All halos are cleared so subsequent local SpMV
    !! operations retain block-Jacobi boundary semantics.
    subroutine apply_vertical_lines(in_vec, out_vec)
        implicit none
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: in_vec
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: out_vec
        integer :: i, j, k
        real(c_double) :: lower

        !$acc parallel loop gang vector collapse(3) present(out_vec)
        do j = j_s-1, j_e+1
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    out_vec(i,k,j) = 0.0_c_double
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(2) &
        !$acc present(in_vec, out_vec, D_inv, line_cprime, C_coef, dz_if) private(lower)
        do j = ys, ys + ym - 1
            do i = xs, xs + xm - 1
                out_vec(i,0,j) = D_inv(i,0,j) * in_vec(i,0,j)
                do k = 1, mz - 1
                    if (i <= 0 .or. j <= 0 .or. i >= mx-1 .or. j >= my-1) then
                        lower = 0.0_c_double
                    else if (k == mz-1) then
                        if (operator_probed) then
                            lower = 0.0_c_double
                        else
                            lower = -1.0_c_double / real(dz_if(i,k,j), c_double)
                        endif
                    else
                        lower = real(C_coef(i,k,j), c_double)
                    endif
                    out_vec(i,k,j) = D_inv(i,k,j) * (in_vec(i,k,j) - lower * out_vec(i,k-1,j))
                enddo
                do k = mz - 2, 0, -1
                    out_vec(i,k,j) = out_vec(i,k,j) - line_cprime(i,k,j) * out_vec(i,k+1,j)
                enddo
            enddo
        enddo
    end subroutine apply_vertical_lines


    !>------------------------------------------------------------
    !! Vector ops over owned cell range (xs..xs+xm-1 etc.).
    !! Local reductions only — global allreduce is the caller's job.
    !!------------------------------------------------------------

    subroutine vec_zero(v)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: v
        integer :: i, j, k
        !$acc parallel loop gang vector collapse(3) present(v)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    v(i,k,j) = 0.0_c_double
                enddo
            enddo
        enddo
    end subroutine vec_zero

    subroutine vec_copy(dst, src)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: src
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: dst
        integer :: i, j, k
        !$acc parallel loop gang vector collapse(3) present(dst, src)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    dst(i,k,j) = src(i,k,j)
                enddo
            enddo
        enddo
    end subroutine vec_copy

    !> y = a*x + y
    subroutine vec_axpy(y, a, x)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: y
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: x
        real(c_double), intent(in) :: a
        integer :: i, j, k
        !$acc parallel loop gang vector collapse(3) present(y, x)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    y(i,k,j) = a * x(i,k,j) + y(i,k,j)
                enddo
            enddo
        enddo
    end subroutine vec_axpy

    !> y = y + a*x + b*z  (fused 3-term update — used for the x_sol += alpha*p_hat + omega*s_hat
    !> step in BiCGStab, replaces two sequential vec_axpy calls and saves one full pass through y).
    subroutine vec_axpy2(y, a, x, b, z)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: y
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: x, z
        real(c_double), intent(in) :: a, b
        integer :: i, j, k
        !$acc parallel loop gang vector collapse(3) present(y, x, z)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    y(i,k,j) = y(i,k,j) + a * x(i,k,j) + b * z(i,k,j)
                enddo
            enddo
        enddo
    end subroutine vec_axpy2

    !> p = r + beta * (p - omega * v)  (classic BiCGStab search-direction update; kept for reference)
    subroutine vec_p_update(p, r, beta_s, omega_s, v)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: p
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: r, v
        real(c_double), intent(in) :: beta_s, omega_s
        integer :: i, j, k
        !$acc parallel loop gang vector collapse(3) present(p, r, v)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    p(i,k,j) = r(i,k,j) + beta_s * (p(i,k,j) - omega_s * v(i,k,j))
                enddo
            enddo
        enddo
    end subroutine vec_p_update

    !> p_hat = u + beta * (p_hat - omega * q)  (pipelined-BiCGStab fused recurrence — replaces
    !> two AXPYs and avoids fresh M^{-1} apply on p)
    subroutine vec_p_hat_update(p_hat_v, u, beta_s, omega_s, q)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: p_hat_v
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: u, q
        real(c_double), intent(in) :: beta_s, omega_s
        integer :: i, j, k
        !$acc parallel loop gang vector collapse(3) present(p_hat_v, u, q)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    p_hat_v(i,k,j) = u(i,k,j) + beta_s * (p_hat_v(i,k,j) - omega_s * q(i,k,j))
                enddo
            enddo
        enddo
    end subroutine vec_p_hat_update

    !> y = a*x + b*y  (axpby — needed for the s_hat = u - alpha*q recurrence)
    subroutine vec_axpby_into(out_v, a, x, b, y)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: out_v
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in)    :: x, y
        real(c_double), intent(in) :: a, b
        integer :: i, j, k
        !$acc parallel loop gang vector collapse(3) present(out_v, x, y)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    out_v(i,k,j) = a * x(i,k,j) + b * y(i,k,j)
                enddo
            enddo
        enddo
    end subroutine vec_axpby_into

    !> result = sum_i x_i * y_i  (LOCAL — caller MPI_Allreduces)
    subroutine vec_dot_local(x, y, result_local)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in) :: x, y
        real(c_double), intent(out) :: result_local
        integer :: i, j, k
        real(c_double) :: s
        s = 0.0_c_double
        !$acc parallel loop gang vector collapse(3) reduction(+:s) present(x, y)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    s = s + x(i,k,j) * y(i,k,j)
                enddo
            enddo
        enddo
        result_local = s
    end subroutine vec_dot_local

    !> Single-pass fused inner products used in BiCGStab's omega allreduce:
    !>   ts = <t, s>, tt = <t, t>, rs = <r_hat, s>, rt = <r_hat, t>, ss = <s, s>
    !> Reads each of t, s, r_hat once instead of 5× sequentially — saves both kernel
    !> launch overhead and ~3× the bandwidth of separate dot products.
    subroutine vec_dots_fused(t_v, s_v, r_h, ts_l, tt_l, rs_l, rt_l, ss_l)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in) :: t_v, s_v, r_h
        real(c_double), intent(out) :: ts_l, tt_l, rs_l, rt_l, ss_l
        integer :: i, j, k
        real(c_double) :: ts_acc, tt_acc, rs_acc, rt_acc, ss_acc
        real(c_double) :: tv, sv, rv

        ts_acc = 0.0_c_double; tt_acc = 0.0_c_double
        rs_acc = 0.0_c_double; rt_acc = 0.0_c_double
        ss_acc = 0.0_c_double

        !$acc parallel loop gang vector collapse(3) &
        !$acc reduction(+:ts_acc,tt_acc,rs_acc,rt_acc,ss_acc) &
        !$acc private(tv, sv, rv) present(t_v, s_v, r_h)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    tv = t_v(i,k,j)
                    sv = s_v(i,k,j)
                    rv = r_h(i,k,j)
                    ts_acc = ts_acc + tv * sv
                    tt_acc = tt_acc + tv * tv
                    rs_acc = rs_acc + rv * sv
                    rt_acc = rt_acc + rv * tv
                    ss_acc = ss_acc + sv * sv
                enddo
            enddo
        enddo

        ts_l = ts_acc
        tt_l = tt_acc
        rs_l = rs_acc
        rt_l = rt_acc
        ss_l = ss_acc
    end subroutine vec_dots_fused

    !> result = sum_i x_i^2  (LOCAL — caller MPI_Allreduces, then sqrt)
    subroutine vec_norm2_local(x, result_local)
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(in) :: x
        real(c_double), intent(out) :: result_local
        integer :: i, j, k
        real(c_double) :: s
        s = 0.0_c_double
        !$acc parallel loop gang vector collapse(3) reduction(+:s) present(x)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    s = s + x(i,k,j) * x(i,k,j)
                enddo
            enddo
        enddo
        result_local = s
    end subroutine vec_norm2_local


    !>------------------------------------------------------------
    !! RHS vector: rhs(i,k,j) = -2*div(i,k,j) for interior; 0 at BCs.
    !! Halo cells are left at their previous value (won't be read).
    !!------------------------------------------------------------
    subroutine compute_rhs_3d()
        implicit none
        integer :: i, j, k

        !$acc parallel loop gang vector collapse(3) present(rhs, div)
        do j = ys, ys + ym - 1
            do k = zs, zs + zm - 1
                do i = xs, xs + xm - 1
                    if (i <= 0 .or. j <= 0 .or. i >= mx-1 .or. j >= my-1 .or. &
                        k <= 0 .or. k >= mz-1) then
                        rhs(i,k,j) = 0.0_c_double
                    else
                        rhs(i,k,j) = real(-2.0 * div(i,k,j), c_double)
                    endif
                enddo
            enddo
        enddo
    end subroutine compute_rhs_3d


    !>------------------------------------------------------------
    !! Halo exchange for a Krylov vector.
    !!
    !! Pack all 4 faces on GPU into persistent buffers, post Irecv
    !! on all non-boundary directions, post Isend on the same, wait
    !! for completion of all 8 (or fewer) requests, unpack into
    !! halo cells on GPU.
    !!
    !! The Isend/Irecv/Waitall pattern lets all 4 directions of comm
    !! happen in parallel — total halo time ~ max(4 directions) instead
    !! of sum. Combined with face-only PCIe and persistent buffers,
    !! halo cost is now latency-bound rather than allocator-bound.
    !!
    !! Tag convention:
    !!   tag_ew (200): data going east (R sends east face, R+1 recvs as west halo)
    !!   tag_we (201): data going west (R sends west face, R-1 recvs as east halo)
    !!   tag_ns (202): data going north (j+ direction)
    !!   tag_sn (203): data going south (j- direction)
    !!------------------------------------------------------------
    subroutine exchange_krylov_halos(v, domain)
        implicit none
        real(c_double), dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: v
        type(domain_t), intent(in) :: domain
        integer :: ierr, i, j, k, nreq
        integer :: reqs(8)
        integer :: nz_w, nx_w, ny_w
        integer, parameter :: tag_ew = 200, tag_we = 201, tag_ns = 202, tag_sn = 203

        nz_w = (k_e+1) - (k_s-1) + 1
        nx_w = (i_e+1) - (i_s-1) + 1
        ny_w = (j_e+1) - (j_s-1) + 1

        ! ===== Pack all four faces on GPU =====
        if (.not. domain%east_boundary) then
            !$acc parallel loop collapse(2) present(v, east_send)
            do j = j_s-1, j_e+1
                do k = k_s-1, k_e+1
                    east_send(k, j) = v(i_e, k, j)
                enddo
            enddo
        endif
        if (.not. domain%west_boundary) then
            !$acc parallel loop collapse(2) present(v, west_send)
            do j = j_s-1, j_e+1
                do k = k_s-1, k_e+1
                    west_send(k, j) = v(i_s, k, j)
                enddo
            enddo
        endif
        if (.not. domain%north_boundary) then
            !$acc parallel loop collapse(2) present(v, north_send)
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    north_send(i, k) = v(i, k, j_e)
                enddo
            enddo
        endif
        if (.not. domain%south_boundary) then
            !$acc parallel loop collapse(2) present(v, south_send)
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    south_send(i, k) = v(i, k, j_s)
                enddo
            enddo
        endif

#ifdef USE_NCCL
        ! ===== NCCL path: device-to-device send/recv (no PCIe) =====
        ! Pack already happened on GPU into east_send/west_send/etc. above.
        ! Wrap all 4 directions of send+recv in nccl_group_start/end so NCCL
        ! fuses them into a single stream-ordered op.
        ! host_data use_device gives us device pointers for the buffers.
        !$acc host_data use_device(east_send, east_recv, west_send, west_recv, &
        !$acc                      north_send, north_recv, south_send, south_recv)
        call nccl_group_start()
        if (.not. domain%east_boundary) then
            ierr = nccl_recv_double(c_loc(east_recv),  nz_w*ny_w, east_neighbor,  nccl_comm, nccl_stream)
            ierr = nccl_send_double(c_loc(east_send),  nz_w*ny_w, east_neighbor,  nccl_comm, nccl_stream)
        endif
        if (.not. domain%west_boundary) then
            ierr = nccl_recv_double(c_loc(west_recv),  nz_w*ny_w, west_neighbor,  nccl_comm, nccl_stream)
            ierr = nccl_send_double(c_loc(west_send),  nz_w*ny_w, west_neighbor,  nccl_comm, nccl_stream)
        endif
        if (.not. domain%north_boundary) then
            ierr = nccl_recv_double(c_loc(north_recv), nx_w*nz_w, north_neighbor, nccl_comm, nccl_stream)
            ierr = nccl_send_double(c_loc(north_send), nx_w*nz_w, north_neighbor, nccl_comm, nccl_stream)
        endif
        if (.not. domain%south_boundary) then
            ierr = nccl_recv_double(c_loc(south_recv), nx_w*nz_w, south_neighbor, nccl_comm, nccl_stream)
            ierr = nccl_send_double(c_loc(south_send), nx_w*nz_w, south_neighbor, nccl_comm, nccl_stream)
        endif
        call nccl_group_end()
        !$acc end host_data
        ! NCCL ops are stream-ordered — subsequent OpenACC kernels on the same
        ! sync queue will see the recv buffers populated.

#else
        ! ===== MPI path: face data via host (face-only PCIe) =====
        ! H2D for send buffers (face-only — small)
        if (.not. domain%east_boundary)  then; !$acc update host(east_send)
        endif
        if (.not. domain%west_boundary)  then; !$acc update host(west_send)
        endif
        if (.not. domain%north_boundary) then; !$acc update host(north_send)
        endif
        if (.not. domain%south_boundary) then; !$acc update host(south_send)
        endif

        ! Post all Irecvs first (to maximise overlap with Isend)
        nreq = 0
        if (.not. domain%east_boundary) then
            nreq = nreq + 1
            call MPI_Irecv(east_recv(k_s-1, j_s-1),  nz_w*ny_w, MPI_DOUBLE_PRECISION, east_neighbor,  tag_we, &
                           solver_comm, reqs(nreq), ierr)
        endif
        if (.not. domain%west_boundary) then
            nreq = nreq + 1
            call MPI_Irecv(west_recv(k_s-1, j_s-1),  nz_w*ny_w, MPI_DOUBLE_PRECISION, west_neighbor,  tag_ew, &
                           solver_comm, reqs(nreq), ierr)
        endif
        if (.not. domain%north_boundary) then
            nreq = nreq + 1
            call MPI_Irecv(north_recv(i_s-1, k_s-1), nx_w*nz_w, MPI_DOUBLE_PRECISION, north_neighbor, tag_sn, &
                           solver_comm, reqs(nreq), ierr)
        endif
        if (.not. domain%south_boundary) then
            nreq = nreq + 1
            call MPI_Irecv(south_recv(i_s-1, k_s-1), nx_w*nz_w, MPI_DOUBLE_PRECISION, south_neighbor, tag_ns, &
                           solver_comm, reqs(nreq), ierr)
        endif

        ! Post all Isends
        if (.not. domain%east_boundary) then
            nreq = nreq + 1
            call MPI_Isend(east_send(k_s-1, j_s-1),  nz_w*ny_w, MPI_DOUBLE_PRECISION, east_neighbor,  tag_ew, &
                           solver_comm, reqs(nreq), ierr)
        endif
        if (.not. domain%west_boundary) then
            nreq = nreq + 1
            call MPI_Isend(west_send(k_s-1, j_s-1),  nz_w*ny_w, MPI_DOUBLE_PRECISION, west_neighbor,  tag_we, &
                           solver_comm, reqs(nreq), ierr)
        endif
        if (.not. domain%north_boundary) then
            nreq = nreq + 1
            call MPI_Isend(north_send(i_s-1, k_s-1), nx_w*nz_w, MPI_DOUBLE_PRECISION, north_neighbor, tag_ns, &
                           solver_comm, reqs(nreq), ierr)
        endif
        if (.not. domain%south_boundary) then
            nreq = nreq + 1
            call MPI_Isend(south_send(i_s-1, k_s-1), nx_w*nz_w, MPI_DOUBLE_PRECISION, south_neighbor, tag_sn, &
                           solver_comm, reqs(nreq), ierr)
        endif

        ! Wait for all
        if (nreq > 0) call MPI_Waitall(nreq, reqs(1:nreq), MPI_STATUSES_IGNORE, ierr)

        ! D2H for recv buffers (face-only)
        if (.not. domain%east_boundary)  then; !$acc update device(east_recv)
        endif
        if (.not. domain%west_boundary)  then; !$acc update device(west_recv)
        endif
        if (.not. domain%north_boundary) then; !$acc update device(north_recv)
        endif
        if (.not. domain%south_boundary) then; !$acc update device(south_recv)
        endif
#endif

        ! ===== Unpack on GPU =====
        if (.not. domain%east_boundary) then
            !$acc parallel loop collapse(2) present(v, east_recv)
            do j = j_s-1, j_e+1
                do k = k_s-1, k_e+1
                    v(i_e+1, k, j) = east_recv(k, j)
                enddo
            enddo
        endif
        if (.not. domain%west_boundary) then
            !$acc parallel loop collapse(2) present(v, west_recv)
            do j = j_s-1, j_e+1
                do k = k_s-1, k_e+1
                    v(i_s-1, k, j) = west_recv(k, j)
                enddo
            enddo
        endif
        if (.not. domain%north_boundary) then
            !$acc parallel loop collapse(2) present(v, north_recv)
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    v(i, k, j_e+1) = north_recv(i, k)
                enddo
            enddo
        endif
        if (.not. domain%south_boundary) then
            !$acc parallel loop collapse(2) present(v, south_recv)
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    v(i, k, j_s-1) = south_recv(i, k)
                enddo
            enddo
        endif
    end subroutine exchange_krylov_halos


    !>------------------------------------------------------------
    !! Apply the BiCGStab solution `x_sol` to the wind field.
    !!
    !! Direct port of calc_updated_winds from wind_iterative_amgx.F90:1117-1423,
    !! but consumes x_sol in 3D form (no flat→3D reshape needed).
    !!------------------------------------------------------------
    subroutine calc_updated_winds(domain, adv_den)
        implicit none
        type(domain_t), intent(inout) :: domain
        logical,        intent(in)    :: adv_den
        real, allocatable, dimension(:,:,:) :: u_dlambdz, v_dlambdz, u_temp, v_temp
        real, allocatable, dimension(:,:,:) :: lambda_3d
        real, allocatable, dimension(:,:,:) :: rho, rho_u, rho_v, rho_w
        integer :: i, j, k, i_start, i_end, j_start, j_end

        i_start = i_s
        i_end   = i_e + 1
        j_start = j_s
        j_end   = j_e + 1

        allocate(u_temp   (i_start:i_end,    k_s-1:k_e+1, j_s:j_e))
        allocate(v_temp   (i_s:i_e,          k_s-1:k_e+1, j_start:j_end))
        allocate(lambda_3d(i_s-1:i_e+1,      k_s-1:k_e+1, j_s-1:j_e+1))

        allocate(u_dlambdz(i_start:i_end,    k_s:k_e, j_s:j_e))
        allocate(v_dlambdz(i_s:i_e,          k_s:k_e, j_start:j_end))

        allocate(rho   (domain%ims:domain%ime, k_s:k_e, domain%jms:domain%jme))
        allocate(rho_u (i_start:i_end, k_s:k_e, j_s:j_e))
        allocate(rho_v (i_s:i_e,       k_s:k_e, j_start:j_end))
        allocate(rho_w (i_s:i_e,       k_s:k_e, j_s:j_e))

        !$acc enter data create(lambda_3d)

        ! Cast x_sol (double) into the single-precision lambda_3d on GPU
        !$acc parallel loop gang vector collapse(3) present(x_sol, lambda_3d)
        do j = j_s-1, j_e+1
            do k = k_s-1, k_e+1
                do i = i_s-1, i_e+1
                    lambda_3d(i,k,j) = real(x_sol(i,k,j))
                enddo
            enddo
        enddo

        ! Halo exchange for lambda_3d (single precision — reuse AMGX path's pattern)
        !$acc update host(lambda_3d)
        call exchange_lambda_halos_real(lambda_3d, domain)
        !$acc update device(lambda_3d)

        associate(density => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                  u       => domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d, &
                  v       => domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d, &
                  w       => domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d, &
                  dz      => domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                  alpha_d => domain%vars_3d(domain%var_indx(kVARS%wind_alpha)%v)%data_3d, &
                  jaco_u_domain => domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, &
                  jaco_v_domain => domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d, &
                  jaco_domain   => domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d, &
                  dzdx_u => domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d, &
                  dzdy_v => domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d, &
                  mf_mx_u => domain%mapfac_mx_u, &
                  mf_my_v => domain%mapfac_my_v)

        !$acc enter data create(u_temp, v_temp, u_dlambdz, v_dlambdz, rho, rho_u, rho_v, rho_w)
        !$acc data present(density, u, v, jaco_u_domain, jaco_v_domain, jaco_domain, dzdx_u, dzdy_v, &
        !$acc              u_temp, v_temp, u_dlambdz, v_dlambdz, rho, rho_u, rho_v, rho_w, &
        !$acc              alpha, dz_if, lambda_3d, mf_mx_u, mf_my_v)

        !$acc kernels
        rho   = 1.0
        rho_w = 1.0
        rho_u = 1.0
        rho_v = 1.0
        !$acc end kernels

        if (adv_den) then
            !$acc parallel loop gang vector collapse(3)
            do j = jms, jme
                do k = k_s, k_e
                    do i = ims, ime
                        rho(i,k,j) = density(i,k,j)
                    enddo
                enddo
            enddo
        endif

        !$acc parallel
        if (i_s == ids .and. i_e == ide) then
            !$acc loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_start+1, i_end-1
                        rho_u(i,k,j) = 0.5 * (rho(i,k,j) + rho(i-1,k,j))
                    enddo
                enddo
            enddo
            !$acc loop gang vector collapse(2)
            do j = j_s, j_e
                do k = k_s, k_e
                    rho_u(i_end,k,j)   = rho(i_end-1,k,j)
                    rho_u(i_start,k,j) = rho(i_start,k,j)
                enddo
            enddo
        else if (i_s == ids) then
            !$acc loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_start+1, i_end
                        rho_u(i,k,j) = 0.5 * (rho(i,k,j) + rho(i-1,k,j))
                    enddo
                enddo
            enddo
            !$acc loop gang vector collapse(2)
            do j = j_s, j_e
                do k = k_s, k_e
                    rho_u(i_start,k,j) = rho(i_start,k,j)
                enddo
            enddo
        else if (i_e == ide) then
            !$acc loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_start, i_end-1
                        rho_u(i,k,j) = 0.5 * (rho(i,k,j) + rho(i-1,k,j))
                    enddo
                enddo
            enddo
            !$acc loop gang vector collapse(2)
            do j = j_s, j_e
                do k = k_s, k_e
                    rho_u(i_end,k,j) = rho(i_end-1,k,j)
                enddo
            enddo
        else
            !$acc loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_start, i_end
                        rho_u(i,k,j) = 0.5 * (rho(i,k,j) + rho(i-1,k,j))
                    enddo
                enddo
            enddo
        endif
        !$acc end parallel

        !$acc parallel
        if (j_s == jds .and. j_e == jde) then
            !$acc loop gang vector collapse(3)
            do j = j_start+1, j_end-1
                do k = k_s, k_e
                    do i = i_s, i_e
                        rho_v(i,k,j) = 0.5 * (rho(i,k,j) + rho(i,k,j-1))
                    enddo
                enddo
            enddo
            !$acc loop gang vector collapse(2)
            do k = k_s, k_e
                do i = i_s, i_e
                    rho_v(i,k,j_start) = rho(i,k,j_start)
                    rho_v(i,k,j_end)   = rho(i,k,j_end-1)
                enddo
            enddo
        else if (j_s == jds) then
            !$acc loop gang vector collapse(3)
            do j = j_start+1, j_end
                do k = k_s, k_e
                    do i = i_s, i_e
                        rho_v(i,k,j) = 0.5 * (rho(i,k,j) + rho(i,k,j-1))
                    enddo
                enddo
            enddo
            !$acc loop gang vector collapse(2)
            do k = k_s, k_e
                do i = i_s, i_e
                    rho_v(i,k,j_start) = rho(i,k,j_start)
                enddo
            enddo
        else if (j_e == jde) then
            !$acc loop gang vector collapse(3)
            do j = j_start, j_end-1
                do k = k_s, k_e
                    do i = i_s, i_e
                        rho_v(i,k,j) = 0.5 * (rho(i,k,j) + rho(i,k,j-1))
                    enddo
                enddo
            enddo
            !$acc loop gang vector collapse(2)
            do k = k_s, k_e
                do i = i_s, i_e
                    rho_v(i,k,j_end) = rho(i,k,j_end-1)
                enddo
            enddo
        else
            !$acc loop gang vector collapse(3)
            do j = j_start, j_end
                do k = k_s, k_e
                    do i = i_s, i_e
                        rho_v(i,k,j) = 0.5 * (rho(i,k,j) + rho(i,k,j-1))
                    enddo
                enddo
            enddo
        endif
        !$acc end parallel

        !$acc parallel loop gang vector tile(32,2,1)
        do j = j_s, j_e
            do k = k_s, k_e-1
                do i = i_s, i_e
                    rho_w(i,k,j) = ( rho(i,k,j)*dz(i,k+1,j) + rho(i,k+1,j)*dz(i,k,j) ) / &
                                   (dz(i,k,j) + dz(i,k+1,j))
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(2)
        do j = j_s, j_e
            do i = i_s, i_e
                rho_w(i,k_e,j) = rho(i,k_e,j)
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = j_s, j_e
            do k = k_s-1, k_e+1
                do i = i_start, i_end
                    u_temp(i,k,j) = (lambda_3d(i,k,j) + lambda_3d(i-1,k,j)) / 2
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = j_start, j_end
            do k = k_s-1, k_e+1
                do i = i_s, i_e
                    v_temp(i,k,j) = (lambda_3d(i,k,j) + lambda_3d(i,k,j-1)) / 2
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_start, i_end
                    u_dlambdz(i,k,j) = u_temp(i,k+1,j) - u_temp(i,k-1,j)
                    u_dlambdz(i,k,j) = u_dlambdz(i,k,j) / (dz_if(i_s,k+1,j_s) + dz_if(i_s,k,j_s))
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = j_start, j_end
            do k = k_s, k_e
                do i = i_s, i_e
                    v_dlambdz(i,k,j) = v_temp(i,k+1,j) - v_temp(i,k-1,j)
                    v_dlambdz(i,k,j) = v_dlambdz(i,k,j) / (dz_if(i_s,k+1,j_s) + dz_if(i_s,k,j_s))
                enddo
            enddo
        enddo

        ! One-sided vertical gradients at the bottom/top levels so the
        ! correction operator never reads ghost-plane lambda (the ghost
        ! rows are identity in the probed operator, lambda = 0).
        !$acc parallel loop gang vector collapse(2)
        do j = j_s, j_e
            do i = i_start, i_end
                u_dlambdz(i,k_s,j) = (u_temp(i,k_s+1,j) - u_temp(i,k_s,j)) / dz_if(i_s,k_s+1,j_s)
                u_dlambdz(i,k_e,j) = (u_temp(i,k_e,j) - u_temp(i,k_e-1,j)) / dz_if(i_s,k_e,j_s)
            enddo
        enddo
        !$acc parallel loop gang vector collapse(2)
        do j = j_start, j_end
            do i = i_s, i_e
                v_dlambdz(i,k_s,j) = (v_temp(i,k_s+1,j) - v_temp(i,k_s,j)) / dz_if(i_s,k_s+1,j_s)
                v_dlambdz(i,k_e,j) = (v_temp(i,k_e,j) - v_temp(i,k_e-1,j)) / dz_if(i_s,k_e,j_s)
            enddo
        enddo

        ! Map factors: the horizontal lambda gradient in TRUE distance is
        ! m_x * (terrain-following bracket) — the whole bracket scales,
        ! including the dzdx cross term (true slope = m_x * grid slope).
        ! Factors are exactly 1.0 when use_map_factors is off. The probed
        ! operator A = D∘G absorbs these on the next (re-)probe.
        !$acc parallel loop gang vector collapse(3)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_start, i_end
                    u(i,k,j) = u(i,k,j) + 0.5 * mf_mx_u(i,j) * &
                                          ( (lambda_3d(i,k,j) - lambda_3d(i-1,k,j)) / dx - &
                                            dzdx_u(i,k,j) * (u_dlambdz(i,k,j)) / jaco_u_domain(i,k,j) ) &
                                          / (rho_u(i,k,j))
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = j_start, j_end
            do k = k_s, k_e
                do i = i_s, i_e
                    v(i,k,j) = v(i,k,j) + 0.5 * mf_my_v(i,j) * &
                                          ( (lambda_3d(i,k,j) - lambda_3d(i,k,j-1)) / dx - &
                                            dzdy_v(i,k,j) * (v_dlambdz(i,k,j)) / jaco_v_domain(i,k,j) ) &
                                          / (rho_v(i,k,j))
                enddo
            enddo
        enddo

        ! Vertical velocity correction: applied to the grid-relative w
        ! predictor (w%dqdt) so the corrected (u, v, w_grid) triplet
        ! jointly satisfies continuity — the unified form for both the
        ! RANS and diagnostic paths (the diagnostic path then derives its
        ! final w via balance_uvw as before, from a better-converged
        ! horizontal field). The gradient is the COMPACT one centred on
        ! interface k ((lambda(k+1)-lambda(k))/dz_if(k+1)): accurate for a
        ! face-located quantity, no ghost-plane reads, and a composition
        ! that fits the probed 15-point stencil. The lid interface (k_e)
        ! gets no correction (rigid lid).
        associate(w_g => domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d)
        !$acc parallel loop gang vector collapse(3) present(w_g, lambda_3d, dz_if)
        do j = j_s, j_e
            do k = k_s, k_e-1
                do i = i_s, i_e
                    w_g(i,k,j) = w_g(i,k,j) + 0.5 * (alpha(i,k,j)**2) * &
                                 (lambda_3d(i,k+1,j) - lambda_3d(i,k,j)) / dz_if(i,k+1,j) &
                                 / jaco_domain(i,k,j) / rho_w(i,k,j)
                enddo
            enddo
        enddo
        end associate

        !$acc end data
        !$acc exit data delete(u_temp, v_temp, u_dlambdz, v_dlambdz, rho, rho_u, rho_v, rho_w, lambda_3d)
        end associate

        deallocate(u_temp, v_temp, lambda_3d, u_dlambdz, v_dlambdz)
        deallocate(rho, rho_u, rho_v, rho_w)
    end subroutine calc_updated_winds


    !>------------------------------------------------------------
    !! Single-precision halo exchange for lambda_3d (used by
    !! calc_updated_winds). Direct mirror of
    !! exchange_lambda_halos in wind_iterative_amgx.F90:2214-2310.
    !!------------------------------------------------------------
    subroutine exchange_lambda_halos_real(lambda_3d, domain)
        implicit none
        real, dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: lambda_3d
        type(domain_t), intent(in) :: domain
        integer :: ierr, my_rank
        integer :: mpi_stat(MPI_STATUS_SIZE)
        integer :: nz, nx, ny
        integer :: east_rank, west_rank, north_rank, south_rank
        integer, parameter :: tag_ew = 100, tag_we = 101, tag_ns = 102, tag_sn = 103
        real, allocatable :: send_buff(:,:), rec_buff(:,:)

        nz = (k_e+1) - (k_s-1) + 1
        nx = (i_e+1) - (i_s-1) + 1
        ny = (j_e+1) - (j_s-1) + 1

        call MPI_Comm_rank(domain%compute_comms, my_rank, ierr)
        east_rank  = my_rank + 1
        west_rank  = my_rank - 1
        north_rank = my_rank + domain%grid%ximages
        south_rank = my_rank - domain%grid%ximages

        if (.not. domain%east_boundary) then
            allocate(send_buff(nz, ny), rec_buff(nz, ny))
            send_buff = lambda_3d(i_e, k_s-1:k_e+1, j_s-1:j_e+1)
            call MPI_Sendrecv(send_buff, nz*ny, MPI_REAL, east_rank, tag_ew, &
                              rec_buff,  nz*ny, MPI_REAL, east_rank, tag_we, &
                              domain%compute_comms, mpi_stat, ierr)
            lambda_3d(i_e+1, k_s-1:k_e+1, j_s-1:j_e+1) = rec_buff
            deallocate(send_buff, rec_buff)
        endif
        if (.not. domain%west_boundary) then
            allocate(send_buff(nz, ny), rec_buff(nz, ny))
            send_buff = lambda_3d(i_s, k_s-1:k_e+1, j_s-1:j_e+1)
            call MPI_Sendrecv(send_buff, nz*ny, MPI_REAL, west_rank, tag_we, &
                              rec_buff,  nz*ny, MPI_REAL, west_rank, tag_ew, &
                              domain%compute_comms, mpi_stat, ierr)
            lambda_3d(i_s-1, k_s-1:k_e+1, j_s-1:j_e+1) = rec_buff
            deallocate(send_buff, rec_buff)
        endif
        if (.not. domain%north_boundary) then
            allocate(send_buff(nx, nz), rec_buff(nx, nz))
            send_buff = lambda_3d(i_s-1:i_e+1, k_s-1:k_e+1, j_e)
            call MPI_Sendrecv(send_buff, nz*nx, MPI_REAL, north_rank, tag_ns, &
                              rec_buff,  nz*nx, MPI_REAL, north_rank, tag_sn, &
                              domain%compute_comms, mpi_stat, ierr)
            lambda_3d(i_s-1:i_e+1, k_s-1:k_e+1, j_e+1) = rec_buff
            deallocate(send_buff, rec_buff)
        endif
        if (.not. domain%south_boundary) then
            allocate(send_buff(nx, nz), rec_buff(nx, nz))
            send_buff = lambda_3d(i_s-1:i_e+1, k_s-1:k_e+1, j_s)
            call MPI_Sendrecv(send_buff, nz*nx, MPI_REAL, south_rank, tag_sn, &
                              rec_buff,  nz*nx, MPI_REAL, south_rank, tag_ns, &
                              domain%compute_comms, mpi_stat, ierr)
            lambda_3d(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1) = rec_buff
            deallocate(send_buff, rec_buff)
        endif
    end subroutine exchange_lambda_halos_real


    !>------------------------------------------------------------
    !! Cleanup: device delete + host deallocate.
    !!------------------------------------------------------------
    subroutine finalize_iter_winds()
        implicit none
        integer :: n

        if (.not. initialized_iter_winds) return

        ! Finalize all cached nests by restoring each one and cleaning it up.
        ! Mirrors finalize_amgx pattern.
        do n = 1, MAX_NESTS
            if (n /= active_nest_indx .and. domain_cache(n)%valid) then
                call restore_from_cache(n)
                call finalize_active_state()
            endif
        enddo

        ! Finalize the currently-active state
        call finalize_active_state()

#ifdef USE_NCCL
        if (nccl_initialized) then
            call nccl_comm_destroy(nccl_comm)
            nccl_comm = c_null_ptr
            ! nccl_stream is aliased to OpenACC's sync queue — do NOT destroy
            nccl_stream = c_null_ptr
            nccl_initialized = .false.
        endif
#endif

        initialized_iter_winds  = .false.
        structure_uploaded = .false.
        active_nest_indx   = -1
    end subroutine finalize_iter_winds


    !>------------------------------------------------------------
    !! Clean up the currently-active state (device data + host alloc).
    !! Used by finalize_iter_winds both for the live state and for each
    !! cached nest after restoring it.
    !!------------------------------------------------------------
    subroutine finalize_active_state()
        implicit none

        call release_multilevel_preconditioner()

        if (structure_uploaded) then
            !$acc exit data delete(x_sol, r_vec, r_hat, p_vec, p_hat, v_vec, &
            !$acc                  s_vec, s_hat, t_vec, rhs, D_inv, line_cprime, prec_res, &
            !$acc                  east_send, east_recv, west_send, west_recv, &
            !$acc                  north_send, north_recv, south_send, south_recv)
            !$acc exit data delete(A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, &
            !$acc                  H_coef, I_coef, J_coef, K_coef, L_coef, M_coef, N_coef, O_coef)
        endif

        if (allocated(x_sol))       deallocate(x_sol)
        if (allocated(r_vec))       deallocate(r_vec)
        if (allocated(r_hat))       deallocate(r_hat)
        if (allocated(p_vec))       deallocate(p_vec)
        if (allocated(p_hat))       deallocate(p_hat)
        if (allocated(v_vec))       deallocate(v_vec)
        if (allocated(s_vec))       deallocate(s_vec)
        if (allocated(s_hat))       deallocate(s_hat)
        if (allocated(t_vec))       deallocate(t_vec)
        if (allocated(rhs))         deallocate(rhs)
        if (allocated(D_inv))       deallocate(D_inv)
        if (allocated(line_cprime)) deallocate(line_cprime)
        if (allocated(prec_res))    deallocate(prec_res)
#ifdef USE_NCCL
        if (allocated(sigma_dev) .or. allocated(red5_dev) .or. allocated(rho0_dev)) then
            !$acc exit data delete(sigma_dev, red5_dev, rho0_dev)
        endif
        if (allocated(sigma_dev)) deallocate(sigma_dev)
        if (allocated(red5_dev))  deallocate(red5_dev)
        if (allocated(rho0_dev))  deallocate(rho0_dev)

        if (nccl_initialized) then
            call nccl_comm_destroy(nccl_comm)
            nccl_comm = c_null_ptr
            ! nccl_stream is aliased to OpenACC's sync queue — do NOT destroy
            nccl_stream = c_null_ptr
            nccl_initialized = .false.
        endif
#endif
        if (allocated(east_send))   deallocate(east_send)
        if (allocated(east_recv))   deallocate(east_recv)
        if (allocated(west_send))   deallocate(west_send)
        if (allocated(west_recv))   deallocate(west_recv)
        if (allocated(north_send))  deallocate(north_send)
        if (allocated(north_recv))  deallocate(north_recv)
        if (allocated(south_send))  deallocate(south_send)
        if (allocated(south_recv))  deallocate(south_recv)

        if (allocated(A_coef)) deallocate(A_coef)
        if (allocated(B_coef)) deallocate(B_coef)
        if (allocated(C_coef)) deallocate(C_coef)
        if (allocated(D_coef)) deallocate(D_coef)
        if (allocated(E_coef)) deallocate(E_coef)
        if (allocated(F_coef)) deallocate(F_coef)
        if (allocated(G_coef)) deallocate(G_coef)
        if (allocated(H_coef)) deallocate(H_coef)
        if (allocated(I_coef)) deallocate(I_coef)
        if (allocated(J_coef)) deallocate(J_coef)
        if (allocated(K_coef)) deallocate(K_coef)
        if (allocated(L_coef)) deallocate(L_coef)
        if (allocated(M_coef)) deallocate(M_coef)
        if (allocated(N_coef)) deallocate(N_coef)
        if (allocated(O_coef)) deallocate(O_coef)

        if (allocated(dz_if)) then
            !$acc exit data delete(dz_if, div, alpha, dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, sigma)
            !$acc exit data delete(jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag, adv_dz_col)
            deallocate(dz_if)
        endif
        if (allocated(jaco))         deallocate(jaco)
        if (allocated(dzdx))         deallocate(dzdx)
        if (allocated(dzdy))         deallocate(dzdy)
        if (allocated(dzdx_surf))    deallocate(dzdx_surf)
        if (allocated(dzdy_surf))    deallocate(dzdy_surf)
        if (allocated(sigma))        deallocate(sigma)
        if (allocated(alpha))        deallocate(alpha)
        if (allocated(div))          deallocate(div)
        if (allocated(jaco_w))       deallocate(jaco_w)
        if (allocated(dzdx_u_stag))  deallocate(dzdx_u_stag)
        if (allocated(jaco_u_stag))  deallocate(jaco_u_stag)
        if (allocated(dzdy_v_stag))  deallocate(dzdy_v_stag)
        if (allocated(jaco_v_stag))  deallocate(jaco_v_stag)
        if (allocated(adv_dz_col))   deallocate(adv_dz_col)

        structure_uploaded = .false.
    end subroutine finalize_active_state


    !>------------------------------------------------------------
    !! Physics setup — copied verbatim from wind_iterative_amgx.F90
    !! (same source-of-truth for matrix coefficients).
    !!------------------------------------------------------------
    subroutine init_module_vars(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        integer :: i, j, k, i_s_bnd, i_e_bnd, j_s_bnd, j_e_bnd

        i_s = domain%its;  i_e = domain%ite
        k_s = domain%kts;  k_e = domain%kte
        j_s = domain%jts;  j_e = domain%jte
        ims = domain%ims;  ime = domain%ime
        jms = domain%jms;  jme = domain%jme
        ids = domain%grid%ids;  ide = domain%grid%ide
        jds = domain%grid%jds;  jde = domain%grid%jde

        if (ims == ids) i_s = ids
        if (ime == ide) i_e = ide
        if (jms == jds) j_s = jds
        if (jme == jde) j_e = jde

        xs = i_s - merge(1, 0, i_s == ids)
        ys = j_s - merge(1, 0, j_s == jds)
        zs = 0

        mx = ide + 2
        my = jde + 2
        mz = domain%kde + 2

        xm = (i_e + merge(1, 0, i_e == ide)) - xs + 1
        ym = (j_e + merge(1, 0, j_e == jde)) - ys + 1
        zm = mz - zs

        i_s_bnd = i_s; i_e_bnd = i_e; j_s_bnd = j_s; j_e_bnd = j_e
        if (i_s == ids) i_s_bnd = i_s + 1
        if (i_e == ide) i_e_bnd = i_e - 1
        if (j_s == jds) j_s_bnd = j_s + 1
        if (j_e == jde) j_e_bnd = j_e - 1

        hs = domain%grid%halo_size

        if (.not. allocated(dzdx)) then
            allocate(dzdx     (i_s:i_e, k_s:k_e,   j_s:j_e))
            allocate(dzdy     (i_s:i_e, k_s:k_e,   j_s:j_e))
            allocate(jaco     (i_s:i_e, k_s:k_e,   j_s:j_e))
            allocate(dzdx_surf(i_s:i_e,            j_s:j_e))
            allocate(dzdy_surf(i_s:i_e,            j_s:j_e))
            allocate(sigma    (i_s:i_e, k_s:k_e,   j_s:j_e))
            allocate(dz_if    (i_s:i_e, k_s:k_e+1, j_s:j_e))
            allocate(alpha    (i_s:i_e, k_s:k_e,   j_s:j_e))
            allocate(div      (i_s:i_e, k_s:k_e,   j_s:j_e))
            allocate(jaco_w     (i_s:i_e,   k_s:k_e, j_s:j_e))
            allocate(dzdx_u_stag(i_s:i_e+1, k_s:k_e, j_s:j_e))
            allocate(jaco_u_stag(i_s:i_e+1, k_s:k_e, j_s:j_e))
            allocate(dzdy_v_stag(i_s:i_e,   k_s:k_e, j_s:j_e+1))
            allocate(jaco_v_stag(i_s:i_e,   k_s:k_e, j_s:j_e+1))
            allocate(adv_dz_col (k_s:k_e))

            !$acc enter data create(dz_if, div, alpha, dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, sigma)
            !$acc enter data create(jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag, adv_dz_col)

            dx = domain%dx

            associate(advection_dz_var    => domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                      neighbor_terrain_var=> domain%vars_2d(domain%var_indx(kVARS%neighbor_terrain)%v)%data_2d, &
                      adv_dz_dom         => domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                      dzdx_domain         => domain%vars_3d(domain%var_indx(kVARS%dzdx)%v)%data_3d, &
                      dzdy_domain         => domain%vars_3d(domain%var_indx(kVARS%dzdy)%v)%data_3d, &
                      dzdx_u_domain      => domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d, &
                      dzdy_v_domain      => domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d, &
                      jaco_u_domain      => domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, &
                      jaco_v_domain      => domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d, &
                      jaco_w_domain      => domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d, &
                      jaco_domain         => domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d)

            !$acc data present(advection_dz_var, neighbor_terrain_var, adv_dz_dom, dzdx_domain, dzdy_domain, dzdx_u_domain, dzdy_v_domain, &
            !$acc                  jaco_u_domain, jaco_v_domain, jaco_w_domain, jaco_domain, dz_if, div, alpha, dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, sigma, &
            !$acc                  jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag, adv_dz_col)

            !$acc parallel loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s+1, k_e
                    do i = i_s, i_e
                        dz_if(i,k,j) = (advection_dz_var(i,k,j) + advection_dz_var(i,k-1,j)) / 2
                    enddo
                enddo
            enddo
            !$acc parallel loop gang vector collapse(2)
            do j = j_s, j_e
                do i = i_s, i_e
                    dz_if(i,k_s,j)   = advection_dz_var(i,k_s,j)
                    dz_if(i,k_e+1,j) = advection_dz_var(i,k_e,j)
                    dzdx_surf(i,j) = dzdx_domain(i,k_s,j)
                    dzdy_surf(i,j) = dzdy_domain(i,k_s,j)
                enddo
            enddo
            !$acc parallel loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_s, i_e
                        dzdx(i,k,j)  = dzdx_domain(i,k,j)
                        dzdy(i,k,j)  = dzdy_domain(i,k,j)
                        jaco(i,k,j)  = jaco_domain(i,k,j)
                        jaco_w(i,k,j)  = jaco_w_domain(i,k,j)
                        sigma(i,k,j) = dz_if(i,k,j) / dz_if(i,k+1,j)
                    enddo
                enddo
            enddo
            !$acc parallel loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_s, i_e+1
                        dzdx_u_stag(i,k,j)  = dzdx_u_domain(i,k,j)
                        jaco_u_stag(i,k,j)  = jaco_u_domain(i,k,j)
                    enddo
                enddo
            enddo
            !$acc parallel loop gang vector collapse(3)
            do j = j_s, j_e+1
                do k = k_s, k_e
                    do i = i_s, i_e
                        dzdy_v_stag(i,k,j)  = dzdy_v_domain(i,k,j)
                        jaco_v_stag(i,k,j)  = jaco_v_domain(i,k,j)
                    enddo
                enddo
            enddo

            !$acc parallel loop gang vector collapse(2)
            do j = j_s, j_e
                do i = i_s_bnd, i_e_bnd
                    dzdx_surf(i,j) = (neighbor_terrain_var(i+1,j) - neighbor_terrain_var(i-1,j)) / (2*dx)
                enddo
            enddo
            !$acc parallel loop gang vector collapse(2)
            do j = j_s_bnd, j_e_bnd
                do i = i_s, i_e
                    dzdy_surf(i,j) = (neighbor_terrain_var(i,j+1) - neighbor_terrain_var(i,j-1)) / (2*dx)
                enddo
            enddo
            !$acc parallel loop
            do k = k_s, k_e
                adv_dz_col(k) = adv_dz_dom(i_s, k, j_s)
            enddo
            !$acc end data
            end associate

        endif
        !$acc update host(dz_if, div, alpha, dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, sigma)
        !$acc update host(jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag, adv_dz_col)

    end subroutine init_module_vars


    subroutine initialize_coefs(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        integer :: k
        real :: mixed_denom

        allocate(A_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(B_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(C_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(D_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(E_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(F_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(G_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(H_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(I_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(J_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(K_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(L_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(M_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(N_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(O_coef(i_s:i_e, k_s:k_e, j_s:j_e))

        A_coef = 1; B_coef = 0; C_coef = 0; D_coef = 0; E_coef = 0
        F_coef = 0; G_coef = 0; H_coef = 0; I_coef = 0; J_coef = 0
        K_coef = 0; L_coef = 0; M_coef = 0; N_coef = 0; O_coef = 0

        D_coef = 1.0/(domain%dx**2)
        E_coef = 1.0/(domain%dx**2)
        F_coef = 1.0/(domain%dx**2)
        G_coef = 1.0/(domain%dx**2)

        do k = k_s, k_e
            mixed_denom = 2*domain%dx*(dz_if(i_s,k+1,j_s) + dz_if(i_s,k,j_s))
            H_coef(:,k,:) = -(dzdx(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s+1:i_e+1,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s+1:i_e+1,k,j_s:j_e))/mixed_denom
            I_coef(:,k,:) =  (dzdx(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s:i_e,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s:i_e,k,j_s:j_e))/mixed_denom
            L_coef(:,k,:) = -(dzdy(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,k,j_s+1:j_e+1)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,k,j_s+1:j_e+1))/mixed_denom
            M_coef(:,k,:) =  (dzdy(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,k,j_s:j_e))/mixed_denom
            J_coef(:,k,:) =  (dzdx(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s+1:i_e+1,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s+1:i_e+1,k,j_s:j_e))/mixed_denom
            K_coef(:,k,:) = -(dzdx(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s:i_e,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s:i_e,k,j_s:j_e))/mixed_denom
            N_coef(:,k,:) =  (dzdy(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,k,j_s+1:j_e+1)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,k,j_s+1:j_e+1))/mixed_denom
            O_coef(:,k,:) = -(dzdy(:,k,:)/jaco(:,k,:) + domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,k,j_s:j_e))/mixed_denom
        enddo

        ! B / C / A — host computation (used on first call before update_coefs_gpu)
        call update_coefs_host(domain)
    end subroutine initialize_coefs


    subroutine update_coefs_host(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        real, allocatable, dimension(:,:,:) :: mixed_denom, X_coef
        real, allocatable, dimension(:,:)   :: M_up, M_dwn
        integer :: k

        allocate(mixed_denom(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(X_coef     (i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(M_up       (i_s:i_e,          j_s:j_e))
        allocate(M_dwn      (i_s:i_e,          j_s:j_e))

        mixed_denom = (dz_if(:,k_s+1:k_e+1,:) + dz_if(:,k_s:k_e,:)) * 2 * domain%dx

        do k = k_s, k_e
            if (k == k_s) then
                M_up = dzdx(:,k,:)*(dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        dzdx(:,k+1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = dzdx(:,k,:)*dzdx_surf
                M_up = M_up + dzdy(:,k,:)*(dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        dzdy(:,k+1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + dzdy(:,k,:)*dzdy_surf
                M_up = M_up + (alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        alpha(:,k+1,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + alpha(:,k,:)**2
            else if (k == k_e) then
                M_up = dzdx(:,k,:)*(dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s) + &
                        dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))
                M_dwn = dzdx(:,k,:)*(dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        dzdx(:,k-1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
                M_up = M_up + dzdy(:,k,:)*(dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s) + &
                        dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))
                M_dwn = M_dwn + dzdy(:,k,:)*(dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        dzdy(:,k-1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
                M_up = M_up + (alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s) + &
                        alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))
                M_dwn = M_dwn + (alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        alpha(:,k-1,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
            else
                M_up = dzdx(:,k,:)*(dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        dzdx(:,k+1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = dzdx(:,k,:)*(dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        dzdx(:,k-1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
                M_up = M_up + dzdy(:,k,:)*(dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        dzdy(:,k+1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + dzdy(:,k,:)*(dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        dzdy(:,k-1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
                M_up = M_up + (alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        alpha(:,k+1,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + (alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        alpha(:,k-1,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
            endif
            B_coef(:,k,:) = M_up /(jaco(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)*dz_if(i_s,k+1,j_s))
            C_coef(:,k,:) = M_dwn/(jaco(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)*dz_if(i_s,k,  j_s))
        enddo

        A_coef = -4/(domain%dx**2) - B_coef - C_coef

        X_coef = -((domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s+1:i_e+1,:,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s+1:i_e+1,:,j_s:j_e) - domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s:i_e,:,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s:i_e,:,j_s:j_e)) + &
                   (domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,:,j_s+1:j_e+1)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,:,j_s+1:j_e+1) - domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,:,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,:,j_s:j_e)))/mixed_denom

        B_coef = B_coef + X_coef
        C_coef = C_coef - X_coef

        deallocate(mixed_denom, X_coef, M_up, M_dwn)
    end subroutine update_coefs_host


    subroutine update_coefs_gpu()
        implicit none
        integer :: i, j, k
        real :: M_up_val, M_dwn_val, X_coef_val, mixed_denom_val
        real :: adv_dz_k, adv_dz_kp1, adv_dz_km1, sum_dz_up, sum_dz_dwn

        !$acc parallel loop gang vector collapse(3) &
        !$acc present(B_coef, C_coef, A_coef, dzdx, dzdy, jaco, alpha, dz_if, &
        !$acc         dzdx_surf, dzdy_surf, jaco_w, adv_dz_col, &
        !$acc         dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_s, i_e
                    adv_dz_k = adv_dz_col(k)

                    if (k == k_e) then
                        M_up_val = dzdx(i,k,j)**2 + dzdy(i,k,j)**2 + alpha(i,k,j)**2
                    else
                        adv_dz_kp1 = adv_dz_col(k+1)
                        sum_dz_up  = adv_dz_k + adv_dz_kp1
                        M_up_val = dzdx(i,k,j) * (dzdx(i,k,j)*adv_dz_kp1 + dzdx(i,k+1,j)*adv_dz_k) / sum_dz_up / jaco_w(i,k,j) &
                                 + dzdy(i,k,j) * (dzdy(i,k,j)*adv_dz_kp1 + dzdy(i,k+1,j)*adv_dz_k) / sum_dz_up / jaco_w(i,k,j) &
                                 + (alpha(i,k,j)**2*adv_dz_kp1 + alpha(i,k+1,j)**2*adv_dz_k) / sum_dz_up / jaco_w(i,k,j)
                    endif

                    if (k == k_s) then
                        M_dwn_val = dzdx(i,k,j)*dzdx_surf(i,j) + dzdy(i,k,j)*dzdy_surf(i,j) + alpha(i,k,j)**2
                    else
                        adv_dz_km1 = adv_dz_col(k-1)
                        sum_dz_dwn = adv_dz_k + adv_dz_km1
                        M_dwn_val = dzdx(i,k,j) * (dzdx(i,k,j)*adv_dz_km1 + dzdx(i,k-1,j)*adv_dz_k) / sum_dz_dwn / jaco_w(i,k-1,j) &
                                  + dzdy(i,k,j) * (dzdy(i,k,j)*adv_dz_km1 + dzdy(i,k-1,j)*adv_dz_k) / sum_dz_dwn / jaco_w(i,k-1,j) &
                                  + (alpha(i,k,j)**2*adv_dz_km1 + alpha(i,k-1,j)**2*adv_dz_k) / sum_dz_dwn / jaco_w(i,k-1,j)
                    endif

                    B_coef(i,k,j) = M_up_val  / (jaco(i,k,j) * adv_dz_k * dz_if(i_s,k+1,j_s))
                    C_coef(i,k,j) = M_dwn_val / (jaco(i,k,j) * adv_dz_k * dz_if(i_s,k,  j_s))

                    mixed_denom_val = (dz_if(i,k+1,j) + dz_if(i,k,j)) * 2.0 * dx
                    X_coef_val = -((dzdx_u_stag(i+1,k,j)/jaco_u_stag(i+1,k,j) - dzdx_u_stag(i,k,j)/jaco_u_stag(i,k,j)) + &
                                   (dzdy_v_stag(i,k,j+1)/jaco_v_stag(i,k,j+1) - dzdy_v_stag(i,k,j)/jaco_v_stag(i,k,j))) / mixed_denom_val

                    B_coef(i,k,j) = B_coef(i,k,j) + X_coef_val
                    C_coef(i,k,j) = C_coef(i,k,j) - X_coef_val
                    A_coef(i,k,j) = -4.0/(dx**2) - B_coef(i,k,j) - C_coef(i,k,j)
                enddo
            enddo
        enddo
    end subroutine update_coefs_gpu


    !>------------------------------------------------------------
    !! Save current module state into a cache slot via MOVE_ALLOC.
    !! Host arrays are transferred (no copy); device data must be
    !! removed first because OpenACC tracks attachments by host
    !! address. After this routine, the cache slot owns all the
    !! state; the module's allocatables are unallocated.
    !!------------------------------------------------------------
    subroutine save_to_cache(slot)
        implicit none
        integer, intent(in) :: slot

        ! The hierarchy is derived state.  Rebuild it from the cached exact
        ! operator when this nest becomes active again instead of moving
        ! device-attached derived types through host cache slots.
        call release_multilevel_preconditioner()

        ! Grid scalars
        domain_cache(slot)%i_s = i_s; domain_cache(slot)%i_e = i_e
        domain_cache(slot)%k_s = k_s; domain_cache(slot)%k_e = k_e
        domain_cache(slot)%j_s = j_s; domain_cache(slot)%j_e = j_e
        domain_cache(slot)%ims = ims; domain_cache(slot)%ime = ime
        domain_cache(slot)%jms = jms; domain_cache(slot)%jme = jme
        domain_cache(slot)%ids = ids; domain_cache(slot)%ide = ide
        domain_cache(slot)%jds = jds; domain_cache(slot)%jde = jde
        domain_cache(slot)%xs  = xs;  domain_cache(slot)%ys  = ys;  domain_cache(slot)%zs  = zs
        domain_cache(slot)%xm  = xm;  domain_cache(slot)%ym  = ym;  domain_cache(slot)%zm  = zm
        domain_cache(slot)%mx  = mx;  domain_cache(slot)%my  = my;  domain_cache(slot)%mz  = mz
        domain_cache(slot)%hs  = hs;  domain_cache(slot)%dx  = dx
        domain_cache(slot)%n_rows        = n_rows
        domain_cache(slot)%n_rows_global = n_rows_global
        domain_cache(slot)%solver_rank   = solver_rank
        domain_cache(slot)%east_neighbor  = east_neighbor
        domain_cache(slot)%west_neighbor  = west_neighbor
        domain_cache(slot)%north_neighbor = north_neighbor
        domain_cache(slot)%south_neighbor = south_neighbor
        domain_cache(slot)%solver_comm    = solver_comm
        domain_cache(slot)%structure_uploaded   = structure_uploaded
        domain_cache(slot)%operator_probed           = operator_probed
        domain_cache(slot)%operator_audit_enabled    = operator_audit_enabled
        domain_cache(slot)%operator_audit_done       = operator_audit_done
        domain_cache(slot)%krylov_audit_written      = krylov_audit_written
        domain_cache(slot)%bootstrap_rhs_saved       = bootstrap_rhs_saved
        domain_cache(slot)%multilevel_requested      = multilevel_requested
        domain_cache(slot)%wind_solver_max_iters = wind_solver_max_iters
        domain_cache(slot)%precond_n_sweeps      = precond_n_sweeps

        ! Sync device → host before removing device attachments
        !$acc wait
        if (structure_uploaded) then
            !$acc update host(A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, &
            !$acc             H_coef, I_coef, J_coef, K_coef, L_coef, M_coef, N_coef, O_coef)
            !$acc update host(x_sol, r_vec, r_hat, p_vec, p_hat, v_vec, &
            !$acc             s_vec, s_hat, t_vec, rhs, D_inv, line_cprime, prec_res)
            !$acc update host(east_send, east_recv, west_send, west_recv, &
            !$acc             north_send, north_recv, south_send, south_recv)
            !$acc exit data delete(A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, &
            !$acc                  H_coef, I_coef, J_coef, K_coef, L_coef, M_coef, N_coef, O_coef)
            !$acc exit data delete(x_sol, r_vec, r_hat, p_vec, p_hat, v_vec, &
            !$acc                  s_vec, s_hat, t_vec, rhs, D_inv, line_cprime, prec_res)
            !$acc exit data delete(east_send, east_recv, west_send, west_recv, &
            !$acc                  north_send, north_recv, south_send, south_recv)
#ifdef USE_NCCL
            if (allocated(sigma_dev)) then
                !$acc update host(sigma_dev, red5_dev, rho0_dev)
                !$acc exit data delete(sigma_dev, red5_dev, rho0_dev)
            endif
#endif
        endif
        if (allocated(dz_if)) then
            !$acc exit data delete(dz_if, div, alpha, dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, sigma)
            !$acc exit data delete(jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag, adv_dz_col)
        endif

        ! Geometry / coefficient arrays
        if (allocated(A_coef))       call move_alloc(A_coef,       domain_cache(slot)%A_coef)
        if (allocated(B_coef))       call move_alloc(B_coef,       domain_cache(slot)%B_coef)
        if (allocated(C_coef))       call move_alloc(C_coef,       domain_cache(slot)%C_coef)
        if (allocated(D_coef))       call move_alloc(D_coef,       domain_cache(slot)%D_coef)
        if (allocated(E_coef))       call move_alloc(E_coef,       domain_cache(slot)%E_coef)
        if (allocated(F_coef))       call move_alloc(F_coef,       domain_cache(slot)%F_coef)
        if (allocated(G_coef))       call move_alloc(G_coef,       domain_cache(slot)%G_coef)
        if (allocated(H_coef))       call move_alloc(H_coef,       domain_cache(slot)%H_coef)
        if (allocated(I_coef))       call move_alloc(I_coef,       domain_cache(slot)%I_coef)
        if (allocated(J_coef))       call move_alloc(J_coef,       domain_cache(slot)%J_coef)
        if (allocated(K_coef))       call move_alloc(K_coef,       domain_cache(slot)%K_coef)
        if (allocated(L_coef))       call move_alloc(L_coef,       domain_cache(slot)%L_coef)
        if (allocated(M_coef))       call move_alloc(M_coef,       domain_cache(slot)%M_coef)
        if (allocated(N_coef))       call move_alloc(N_coef,       domain_cache(slot)%N_coef)
        if (allocated(O_coef))       call move_alloc(O_coef,       domain_cache(slot)%O_coef)
        if (allocated(div))          call move_alloc(div,          domain_cache(slot)%div)
        if (allocated(dz_if))        call move_alloc(dz_if,        domain_cache(slot)%dz_if)
        if (allocated(jaco))         call move_alloc(jaco,         domain_cache(slot)%jaco)
        if (allocated(dzdx))         call move_alloc(dzdx,         domain_cache(slot)%dzdx)
        if (allocated(dzdy))         call move_alloc(dzdy,         domain_cache(slot)%dzdy)
        if (allocated(sigma))        call move_alloc(sigma,        domain_cache(slot)%sigma)
        if (allocated(alpha))        call move_alloc(alpha,        domain_cache(slot)%alpha)
        if (allocated(dzdx_surf))    call move_alloc(dzdx_surf,    domain_cache(slot)%dzdx_surf)
        if (allocated(dzdy_surf))    call move_alloc(dzdy_surf,    domain_cache(slot)%dzdy_surf)
        if (allocated(jaco_w))       call move_alloc(jaco_w,       domain_cache(slot)%jaco_w)
        if (allocated(dzdx_u_stag))  call move_alloc(dzdx_u_stag,  domain_cache(slot)%dzdx_u_stag)
        if (allocated(jaco_u_stag))  call move_alloc(jaco_u_stag,  domain_cache(slot)%jaco_u_stag)
        if (allocated(dzdy_v_stag))  call move_alloc(dzdy_v_stag,  domain_cache(slot)%dzdy_v_stag)
        if (allocated(jaco_v_stag))  call move_alloc(jaco_v_stag,  domain_cache(slot)%jaco_v_stag)
        if (allocated(adv_dz_col))   call move_alloc(adv_dz_col,   domain_cache(slot)%adv_dz_col)

        ! Krylov / RHS / precond
        if (allocated(x_sol))    call move_alloc(x_sol,    domain_cache(slot)%x_sol)
        if (allocated(r_vec))    call move_alloc(r_vec,    domain_cache(slot)%r_vec)
        if (allocated(r_hat))    call move_alloc(r_hat,    domain_cache(slot)%r_hat)
        if (allocated(p_vec))    call move_alloc(p_vec,    domain_cache(slot)%p_vec)
        if (allocated(p_hat))    call move_alloc(p_hat,    domain_cache(slot)%p_hat)
        if (allocated(v_vec))    call move_alloc(v_vec,    domain_cache(slot)%v_vec)
        if (allocated(s_vec))    call move_alloc(s_vec,    domain_cache(slot)%s_vec)
        if (allocated(s_hat))    call move_alloc(s_hat,    domain_cache(slot)%s_hat)
        if (allocated(t_vec))    call move_alloc(t_vec,    domain_cache(slot)%t_vec)
        if (allocated(rhs))      call move_alloc(rhs,      domain_cache(slot)%rhs)
        if (allocated(D_inv))    call move_alloc(D_inv,    domain_cache(slot)%D_inv)
        if (allocated(line_cprime)) call move_alloc(line_cprime, domain_cache(slot)%line_cprime)
        if (allocated(prec_res)) call move_alloc(prec_res, domain_cache(slot)%prec_res)

        ! Halo buffers
        if (allocated(east_send))   call move_alloc(east_send,   domain_cache(slot)%east_send)
        if (allocated(east_recv))   call move_alloc(east_recv,   domain_cache(slot)%east_recv)
        if (allocated(west_send))   call move_alloc(west_send,   domain_cache(slot)%west_send)
        if (allocated(west_recv))   call move_alloc(west_recv,   domain_cache(slot)%west_recv)
        if (allocated(north_send))  call move_alloc(north_send,  domain_cache(slot)%north_send)
        if (allocated(north_recv))  call move_alloc(north_recv,  domain_cache(slot)%north_recv)
        if (allocated(south_send))  call move_alloc(south_send,  domain_cache(slot)%south_send)
        if (allocated(south_recv))  call move_alloc(south_recv,  domain_cache(slot)%south_recv)

#ifdef USE_NCCL
        if (allocated(sigma_dev))   call move_alloc(sigma_dev,   domain_cache(slot)%sigma_dev)
        if (allocated(red5_dev))    call move_alloc(red5_dev,    domain_cache(slot)%red5_dev)
        if (allocated(rho0_dev))    call move_alloc(rho0_dev,    domain_cache(slot)%rho0_dev)
#endif

        domain_cache(slot)%valid = .true.
    end subroutine save_to_cache


    !>------------------------------------------------------------
    !! Restore module state from a cache slot via MOVE_ALLOC.
    !! After this routine, the cache slot is empty and the module
    !! owns all the state. Device data is re-registered on GPU.
    !!------------------------------------------------------------
    subroutine restore_from_cache(slot)
        implicit none
        integer, intent(in) :: slot

        ! Grid scalars
        i_s = domain_cache(slot)%i_s; i_e = domain_cache(slot)%i_e
        k_s = domain_cache(slot)%k_s; k_e = domain_cache(slot)%k_e
        j_s = domain_cache(slot)%j_s; j_e = domain_cache(slot)%j_e
        ims = domain_cache(slot)%ims; ime = domain_cache(slot)%ime
        jms = domain_cache(slot)%jms; jme = domain_cache(slot)%jme
        ids = domain_cache(slot)%ids; ide = domain_cache(slot)%ide
        jds = domain_cache(slot)%jds; jde = domain_cache(slot)%jde
        xs  = domain_cache(slot)%xs;  ys  = domain_cache(slot)%ys;  zs  = domain_cache(slot)%zs
        xm  = domain_cache(slot)%xm;  ym  = domain_cache(slot)%ym;  zm  = domain_cache(slot)%zm
        mx  = domain_cache(slot)%mx;  my  = domain_cache(slot)%my;  mz  = domain_cache(slot)%mz
        hs  = domain_cache(slot)%hs;  dx  = domain_cache(slot)%dx
        n_rows        = domain_cache(slot)%n_rows
        n_rows_global = domain_cache(slot)%n_rows_global
        solver_rank   = domain_cache(slot)%solver_rank
        east_neighbor  = domain_cache(slot)%east_neighbor
        west_neighbor  = domain_cache(slot)%west_neighbor
        north_neighbor = domain_cache(slot)%north_neighbor
        south_neighbor = domain_cache(slot)%south_neighbor
        solver_comm    = domain_cache(slot)%solver_comm
        structure_uploaded   = domain_cache(slot)%structure_uploaded
        operator_probed          = domain_cache(slot)%operator_probed
        operator_audit_enabled   = domain_cache(slot)%operator_audit_enabled
        operator_audit_done      = domain_cache(slot)%operator_audit_done
        krylov_audit_written     = domain_cache(slot)%krylov_audit_written
        bootstrap_rhs_saved      = domain_cache(slot)%bootstrap_rhs_saved
        multilevel_requested     = domain_cache(slot)%multilevel_requested
        multilevel_setup_attempted = .false.
        multilevel_ready = .false.
        wind_solver_max_iters = domain_cache(slot)%wind_solver_max_iters
        precond_n_sweeps      = domain_cache(slot)%precond_n_sweeps

        ! Geometry / coefficient arrays
        if (allocated(domain_cache(slot)%A_coef))       call move_alloc(domain_cache(slot)%A_coef,       A_coef)
        if (allocated(domain_cache(slot)%B_coef))       call move_alloc(domain_cache(slot)%B_coef,       B_coef)
        if (allocated(domain_cache(slot)%C_coef))       call move_alloc(domain_cache(slot)%C_coef,       C_coef)
        if (allocated(domain_cache(slot)%D_coef))       call move_alloc(domain_cache(slot)%D_coef,       D_coef)
        if (allocated(domain_cache(slot)%E_coef))       call move_alloc(domain_cache(slot)%E_coef,       E_coef)
        if (allocated(domain_cache(slot)%F_coef))       call move_alloc(domain_cache(slot)%F_coef,       F_coef)
        if (allocated(domain_cache(slot)%G_coef))       call move_alloc(domain_cache(slot)%G_coef,       G_coef)
        if (allocated(domain_cache(slot)%H_coef))       call move_alloc(domain_cache(slot)%H_coef,       H_coef)
        if (allocated(domain_cache(slot)%I_coef))       call move_alloc(domain_cache(slot)%I_coef,       I_coef)
        if (allocated(domain_cache(slot)%J_coef))       call move_alloc(domain_cache(slot)%J_coef,       J_coef)
        if (allocated(domain_cache(slot)%K_coef))       call move_alloc(domain_cache(slot)%K_coef,       K_coef)
        if (allocated(domain_cache(slot)%L_coef))       call move_alloc(domain_cache(slot)%L_coef,       L_coef)
        if (allocated(domain_cache(slot)%M_coef))       call move_alloc(domain_cache(slot)%M_coef,       M_coef)
        if (allocated(domain_cache(slot)%N_coef))       call move_alloc(domain_cache(slot)%N_coef,       N_coef)
        if (allocated(domain_cache(slot)%O_coef))       call move_alloc(domain_cache(slot)%O_coef,       O_coef)
        if (allocated(domain_cache(slot)%div))          call move_alloc(domain_cache(slot)%div,          div)
        if (allocated(domain_cache(slot)%dz_if))        call move_alloc(domain_cache(slot)%dz_if,        dz_if)
        if (allocated(domain_cache(slot)%jaco))         call move_alloc(domain_cache(slot)%jaco,         jaco)
        if (allocated(domain_cache(slot)%dzdx))         call move_alloc(domain_cache(slot)%dzdx,         dzdx)
        if (allocated(domain_cache(slot)%dzdy))         call move_alloc(domain_cache(slot)%dzdy,         dzdy)
        if (allocated(domain_cache(slot)%sigma))        call move_alloc(domain_cache(slot)%sigma,        sigma)
        if (allocated(domain_cache(slot)%alpha))        call move_alloc(domain_cache(slot)%alpha,        alpha)
        if (allocated(domain_cache(slot)%dzdx_surf))    call move_alloc(domain_cache(slot)%dzdx_surf,    dzdx_surf)
        if (allocated(domain_cache(slot)%dzdy_surf))    call move_alloc(domain_cache(slot)%dzdy_surf,    dzdy_surf)
        if (allocated(domain_cache(slot)%jaco_w))       call move_alloc(domain_cache(slot)%jaco_w,       jaco_w)
        if (allocated(domain_cache(slot)%dzdx_u_stag))  call move_alloc(domain_cache(slot)%dzdx_u_stag,  dzdx_u_stag)
        if (allocated(domain_cache(slot)%jaco_u_stag))  call move_alloc(domain_cache(slot)%jaco_u_stag,  jaco_u_stag)
        if (allocated(domain_cache(slot)%dzdy_v_stag))  call move_alloc(domain_cache(slot)%dzdy_v_stag,  dzdy_v_stag)
        if (allocated(domain_cache(slot)%jaco_v_stag))  call move_alloc(domain_cache(slot)%jaco_v_stag,  jaco_v_stag)
        if (allocated(domain_cache(slot)%adv_dz_col))   call move_alloc(domain_cache(slot)%adv_dz_col,   adv_dz_col)

        ! Krylov / RHS / precond
        if (allocated(domain_cache(slot)%x_sol))    call move_alloc(domain_cache(slot)%x_sol,    x_sol)
        if (allocated(domain_cache(slot)%r_vec))    call move_alloc(domain_cache(slot)%r_vec,    r_vec)
        if (allocated(domain_cache(slot)%r_hat))    call move_alloc(domain_cache(slot)%r_hat,    r_hat)
        if (allocated(domain_cache(slot)%p_vec))    call move_alloc(domain_cache(slot)%p_vec,    p_vec)
        if (allocated(domain_cache(slot)%p_hat))    call move_alloc(domain_cache(slot)%p_hat,    p_hat)
        if (allocated(domain_cache(slot)%v_vec))    call move_alloc(domain_cache(slot)%v_vec,    v_vec)
        if (allocated(domain_cache(slot)%s_vec))    call move_alloc(domain_cache(slot)%s_vec,    s_vec)
        if (allocated(domain_cache(slot)%s_hat))    call move_alloc(domain_cache(slot)%s_hat,    s_hat)
        if (allocated(domain_cache(slot)%t_vec))    call move_alloc(domain_cache(slot)%t_vec,    t_vec)
        if (allocated(domain_cache(slot)%rhs))      call move_alloc(domain_cache(slot)%rhs,      rhs)
        if (allocated(domain_cache(slot)%D_inv))    call move_alloc(domain_cache(slot)%D_inv,    D_inv)
        if (allocated(domain_cache(slot)%line_cprime)) call move_alloc(domain_cache(slot)%line_cprime, line_cprime)
        if (allocated(domain_cache(slot)%prec_res)) call move_alloc(domain_cache(slot)%prec_res, prec_res)

        ! Halo buffers
        if (allocated(domain_cache(slot)%east_send))   call move_alloc(domain_cache(slot)%east_send,   east_send)
        if (allocated(domain_cache(slot)%east_recv))   call move_alloc(domain_cache(slot)%east_recv,   east_recv)
        if (allocated(domain_cache(slot)%west_send))   call move_alloc(domain_cache(slot)%west_send,   west_send)
        if (allocated(domain_cache(slot)%west_recv))   call move_alloc(domain_cache(slot)%west_recv,   west_recv)
        if (allocated(domain_cache(slot)%north_send))  call move_alloc(domain_cache(slot)%north_send,  north_send)
        if (allocated(domain_cache(slot)%north_recv))  call move_alloc(domain_cache(slot)%north_recv,  north_recv)
        if (allocated(domain_cache(slot)%south_send))  call move_alloc(domain_cache(slot)%south_send,  south_send)
        if (allocated(domain_cache(slot)%south_recv))  call move_alloc(domain_cache(slot)%south_recv,  south_recv)

#ifdef USE_NCCL
        if (allocated(domain_cache(slot)%sigma_dev))   call move_alloc(domain_cache(slot)%sigma_dev,   sigma_dev)
        if (allocated(domain_cache(slot)%red5_dev))    call move_alloc(domain_cache(slot)%red5_dev,    red5_dev)
        if (allocated(domain_cache(slot)%rho0_dev))    call move_alloc(domain_cache(slot)%rho0_dev,    rho0_dev)
#endif

        ! Re-register on device
        if (allocated(dz_if)) then
            !$acc enter data copyin(dz_if, div, alpha, dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, sigma)
            !$acc enter data copyin(jaco_w, dzdx_u_stag, jaco_u_stag, dzdy_v_stag, jaco_v_stag, adv_dz_col)
        endif
        if (structure_uploaded) then
            !$acc enter data copyin(A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, &
            !$acc                   H_coef, I_coef, J_coef, K_coef, L_coef, M_coef, N_coef, O_coef)
            !$acc enter data copyin(x_sol, r_vec, r_hat, p_vec, p_hat, v_vec, &
            !$acc                   s_vec, s_hat, t_vec, rhs, D_inv, line_cprime, prec_res)
            !$acc enter data copyin(east_send, east_recv, west_send, west_recv, &
            !$acc                   north_send, north_recv, south_send, south_recv)
#ifdef USE_NCCL
            if (allocated(sigma_dev)) then
                !$acc enter data copyin(sigma_dev, red5_dev, rho0_dev)
            endif
#endif
        endif

        domain_cache(slot)%valid = .false.
    end subroutine restore_from_cache


end module wind_iterative
