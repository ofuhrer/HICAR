!> ----------------------------------------------------------------------------
!!  Driver to call different advection schemes
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module advection
    use icar_constants
    use adv_std,                    only : adv_std_init, adv_std_var_request, adv_std_advect3d, adv_std_compute_wind, &
                                           adv_std_clean_wind_arrays, adv_std_clean_wind_arrays_fm, &
                                           adv_theta_ref
    use adv_fluxcorr,               only : init_fluxcorr, set_sign_arrays, compute_upwind_fluxes_async
    ! use debug_module,               only : domain_fix
    use options_interface,          only: options_t
    use domain_interface,           only: domain_t
    use variable_dict_interface,  only : var_dict_t
    use variable_interface,       only : variable_t
    use timer_interface,          only : timer_t
    use data_structures,          only : index_type
    use mpi_f08,                  only : MPI_Allreduce, MPI_IN_PLACE, MPI_REAL, MPI_MIN, MPI_MAX

    implicit none
    private
    integer :: ims, ime, kms, kme, jms, jme, its, ite, jts, jte

    ! Persistent module-level wind/denom buffers, populated by adv_std_compute_wind
    ! on first call of each nest and reused across timesteps. Promoted from `save`
    ! locals inside `advect` so the nest-switch path can free them — otherwise their
    ! shape (and device mapping) lags behind the new nest's bounds and triggers OOB.
    real, allocatable :: U_m(:,:,:), V_m(:,:,:), W_m(:,:,:), denom(:,:,:)

    public :: advect, adv_init, adv_var_request, adv_bound_theta_perturbation
contains

    subroutine adv_init(domain,options, context_chng)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in) :: options
        logical, optional, intent(in) :: context_chng

        logical :: context_change

        if (present(context_chng)) then
            context_change = context_chng
        else
            context_change = .false.
        endif

        if (STD_OUT_PE .and. .not.context_change) write(*,*) ""
        if (STD_OUT_PE .and. .not.context_change) write(*,*) "Initializing Advection"
        if (options%physics%advection==kADV_STD) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Standard"
            ! Drop persistent buffers before re-bounding for the new nest. They
            ! will be re-allocated (host + device) on next adv_std_compute_wind.
            call clean_wind_arrays()
            call adv_std_clean_wind_arrays()
            ! Drop the fine-mesh wind/flux buffers too so they are re-bounded for
            ! the new nest on the next adv_std_compute_wind_2d_fm (snow_drift).
            call adv_std_clean_wind_arrays_fm()
            call adv_std_init(domain,options)
            if (options%adv%flux_corr > 0) call init_fluxcorr(domain)
        endif
        
        ims = domain%ims; ime = domain%ime
        kms = domain%kms; kme = domain%kme
        jms = domain%jms; jme = domain%jme
        its = domain%its; ite = domain%ite
        jts = domain%jts; jte = domain%jte
        
    end subroutine adv_init

    subroutine clean_wind_arrays()
        implicit none

        ! `finalize` zeros the dynamic refcount in one shot. adv_std_compute_wind
        ! calls `enter data create` on these every timestep, so the refcount has
        ! grown to ~Nsteps; without finalize the device mapping survives
        ! deallocate() and later collides with other arrays in the freed host
        ! range (partial-present fatal).
        if (allocated(U_m)) then
            !$acc exit data delete(U_m) finalize
            deallocate(U_m)
        endif
        if (allocated(V_m)) then
            !$acc exit data delete(V_m) finalize
            deallocate(V_m)
        endif
        if (allocated(W_m)) then
            !$acc exit data delete(W_m) finalize
            deallocate(W_m)
        endif
        if (allocated(denom)) then
            !$acc exit data delete(denom) finalize
            deallocate(denom)
        endif
    end subroutine clean_wind_arrays

    subroutine adv_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        !if (options%physics%advection==kADV_UPWIND) then
        !    call upwind_var_request(options)
        if (options%physics%advection==kADV_STD) then
            call adv_std_var_request(options)
        endif
    end subroutine
    
    subroutine advect(domain, options, dt,flux_time, flux_corr_time, sum_time, adv_wind_time)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,intent(in) :: dt
        type(timer_t), intent(inout) :: flux_time, flux_corr_time, sum_time, adv_wind_time
        type(variable_t) :: var_to_advect
        integer :: n

        if (options%physics%advection==kADV_STD) then

            call adv_wind_time%start()
            call adv_std_compute_wind(domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d, &
                                      domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                                      domain%mapfac_my_u, domain%mapfac_mx_v, domain%mapfac_mxy, &
                                      domain%dx,options,dt, U_m, V_m, W_m, denom)
            if (options%adv%flux_corr==kFLUXCOR_MONO) call set_sign_arrays(U_m,V_m,W_m)
            call adv_wind_time%stop()
        else
            return
        endif

        if (options%time%RK3 .and. options%physics%advection==kADV_STD) then
            call RK3_adv(domain, options, U_m, V_m, W_m, denom, dt, flux_time, flux_corr_time, sum_time, adv_wind_time)
        elseif (options%physics%advection==kADV_STD) then
            do n = 1, size(domain%adv_vars)
            if (options%physics%advection==kADV_STD) then
                call adv_std_advect3d(domain%vars_3d(domain%adv_vars(n)%v)%data_3d,domain%vars_3d(domain%adv_vars(n)%v)%data_3d, &
                    U_m, V_m, W_m, denom, domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d,flux_time, flux_corr_time, sum_time)
            endif
            enddo
        endif

    end subroutine advect

    subroutine RK3_adv(domain, options, U_m, V_m, W_m, denom, dt, flux_time, flux_corr_time, sum_time, adv_wind_time)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(in) :: options
        real, allocatable, intent(in) :: U_m(:,:,:), V_m(:,:,:), W_m(:,:,:), denom(:,:,:)
        real, intent(in) :: dt
        type(timer_t), intent(inout) :: flux_time, flux_corr_time, sum_time, adv_wind_time

        integer :: RK3_step, n, n_adv, flux_corr, flux_corr_n, q_id, theta_bound_radius, kd, ku
        logical :: apply_cz_diff_n
        real :: t_fac, theta_min, theta_max
        real, allocatable :: temp_all(:,:,:,:), theta_ref_stage(:,:,:)
        real, allocatable :: theta_min_k(:), theta_max_k(:), theta_bound_min_k(:), theta_bound_max_k(:)

        integer :: i, j, k

        n_adv = size(domain%adv_vars)

        ! A passive monotone advection step should not create new extrema in
        ! full theta. The split-theta correction below keeps a same-level lower
        ! floor so stable upper-level theta cannot ratchet downward over many
        ! RK3 steps, while the upper cap uses a vertical-neighborhood envelope
        ! to allow local import of high-theta air.
        allocate(theta_min_k(kms:kme), theta_max_k(kms:kme))
        allocate(theta_bound_min_k(kms:kme), theta_bound_max_k(kms:kme))
        associate(th => domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d)
        do k = kms, kme
            theta_min =  huge(1.0)
            theta_max = -huge(1.0)
            !$acc parallel loop gang vector collapse(2) default(present) &
            !$acc          reduction(min:theta_min) reduction(max:theta_max)
            do j = jts, jte
                do i = its, ite
                    theta_min = min(theta_min, th(i,k,j))
                    theta_max = max(theta_max, th(i,k,j))
                enddo
            enddo
            theta_min_k(k) = theta_min
            theta_max_k(k) = theta_max
        enddo
        end associate
        call MPI_Allreduce(MPI_IN_PLACE, theta_min_k, kme-kms+1, MPI_REAL, MPI_MIN, domain%compute_comms)
        call MPI_Allreduce(MPI_IN_PLACE, theta_max_k, kme-kms+1, MPI_REAL, MPI_MAX, domain%compute_comms)
        theta_bound_radius = max(1, (options%adv%v_order + 1) / 2)
        do k = kms, kme
            kd = max(kms, k - theta_bound_radius)
            ku = min(kme, k + theta_bound_radius)
            theta_bound_min_k(k) = theta_min_k(k)
            theta_bound_max_k(k) = maxval(theta_max_k(kd:ku))
        enddo

        ! ----  convert theta -> theta' = theta - theta_bar(z) before the
        ! RK3 save, so the whole RK3 (and its inter-stage halo exchanges) acts
        ! on the smooth, weakly-stratified perturbation. theta is rebuilt after
        ! the RK3 loop, so no other code ever sees the perturbation.
        associate(th => domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d)
        !$acc parallel loop gang vector collapse(3) default(present)
        do j = jms, jme
            do k = kms, kme
                do i = ims, ime
                    th(i,k,j) = th(i,k,j) - adv_theta_ref(i,k,j)
                enddo
            enddo
        enddo
        end associate

        ! Allocate 4D temporary array: one 3D slab per advected variable
        allocate(temp_all(ims:ime, kms:kme, jms:jme, n_adv))
        allocate(theta_ref_stage(ims:ime, kms:kme, jms:jme))
        !$acc data create(temp_all, theta_ref_stage) copyin(theta_bound_min_k, theta_bound_max_k)

        ! Save initial states for all advected variables
        do n = 1, n_adv
            associate( curr_var_data3d => domain%vars_3d(domain%adv_vars(n)%v)%data_3d)
            !$acc kernels present(curr_var_data3d, temp_all) async(n)
            temp_all(:,:,:,n) = curr_var_data3d(:,:,:)
            !$acc end kernels
            end associate
        enddo
        !$acc wait

        ! RK3 loop (outer) over variables (inner), with one batch exchange per step
        do RK3_step = 1, 3

            select case(RK3_step)
            case (1)
                t_fac = 1.0/3.0
                flux_corr = 0
            case (2)
                t_fac = 0.5
                flux_corr = 0
            case (3)
                t_fac = 1.0
                flux_corr = options%adv%flux_corr
            end select

            do n = 1, n_adv
                q_id = n

                ! Phase E: under zb_diff, theta BYPASSES the FCT limiter (it
                ! relies on centered advection + the explicit constant-z
                ! diffusion + its minmod limiter; FCT would otherwise re-inject
                ! along-eta upwind diffusion on clipped slope cells). All other
                ! advected vars keep FCT. Skipping FCT for theta also skips its
                ! upwind-reference precompute.
                flux_corr_n = flux_corr
                if (domain%adv_vars(n)%id == domain%var_indx(kVARS%potential_temperature)%id) flux_corr_n = 0

                ! Per-variable cz_diff dispatch: theta always gets it (when
                ! cz_diff_order > 0); qv optionally (cz_diff_qv namelist); other
                ! species never — they aren't strongly stratified and would
                ! waste the constant-z reconstruction cost on every flux3 call.
                apply_cz_diff_n = (domain%adv_vars(n)%id == domain%var_indx(kVARS%potential_temperature)%id .or. &
                                   domain%adv_vars(n)%id == domain%var_indx(kVARS%water_vapor)%id)

                if (flux_corr_n==kFLUXCOR_MONO) then
                    call flux_corr_time%start()
                    call compute_upwind_fluxes_async(temp_all(:,:,:,n), U_m, V_m, W_m, &
                        domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, denom, (q_id+100))
                    call flux_corr_time%stop()
                endif


                call adv_std_advect3d(domain%vars_3d(domain%adv_vars(n)%v)%data_3d, &
                    temp_all(:,:,:,n), U_m, V_m, W_m, denom, &
                    domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, flux_time, flux_corr_time, sum_time, &
                    t_factor_in=t_fac, q_id_in=q_id, flux_corr_in=flux_corr_n, apply_cz_diff_in=apply_cz_diff_n)

                if (domain%adv_vars(n)%id == domain%var_indx(kVARS%potential_temperature)%id) then
                    ! theta is advected as theta' = theta - theta_bar(z).
                    ! Account for transport of theta_bar with the same discrete
                    ! advection operator and RK3 stage factor, then add only
                    ! advected(theta_bar)-theta_bar to theta'. This replaces the
                    ! unstable post-RK3 pointwise -w*dtheta_bar/dz correction.
                    !$acc kernels present(theta_ref_stage, adv_theta_ref)
                    theta_ref_stage(:,:,:) = adv_theta_ref(:,:,:)
                    !$acc end kernels
                    call adv_std_advect3d(theta_ref_stage, adv_theta_ref, U_m, V_m, W_m, denom, &
                        domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, flux_time, flux_corr_time, sum_time, &
                        t_factor_in=t_fac, q_id_in=(q_id+n_adv), flux_corr_in=0, apply_cz_diff_in=.false.)
                    associate(th => domain%vars_3d(domain%adv_vars(n)%v)%data_3d)
                    !$acc parallel loop gang vector collapse(3) default(present)
                    do j = jms, jme
                        do k = kms, kme
                            do i = ims, ime
                                th(i,k,j) = th(i,k,j) + (theta_ref_stage(i,k,j) - adv_theta_ref(i,k,j))
                            enddo
                        enddo
                    enddo
                    end associate
                    call adv_bound_theta_perturbation( &
                        domain%vars_3d(domain%adv_vars(n)%v)%data_3d, &
                        adv_theta_ref, theta_bound_min_k, theta_bound_max_k)
                endif
            enddo

            ! if (RK3_step < 3) then
            ! Single batch exchange for all advected variables after each RK3 step
            call domain%halo_3d_send()
            call domain%halo_3d_retrieve()
            ! endif
        enddo

        ! ---- theta currently holds the RK3-advected perturbation theta'.
        ! Rebuild full theta = theta_bar + theta' over the full memory range
        ! (including boundaries, where theta' = theta_forcing - theta_bar, so
        ! the reconstruction restores theta_forcing exactly).
        !
        ! theta_bar transport has already been included in flux form inside
        ! each RK3 stage above.

        associate(th => domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d)
        !$acc parallel loop gang vector collapse(3) default(present)
        do j = jms, jme
            do k = kms, kme
                do i = ims, ime
                    th(i,k,j) = th(i,k,j) + adv_theta_ref(i,k,j)
                enddo
            enddo
        enddo
        call bound_theta_full(th, theta_bound_min_k, theta_bound_max_k)
        end associate


        !$acc end data
        deallocate(temp_all, theta_ref_stage, theta_min_k, theta_max_k, theta_bound_min_k, theta_bound_max_k)

    end subroutine RK3_adv

    subroutine adv_bound_theta_perturbation(theta_prime, theta_ref, theta_min_k, theta_max_k)
        implicit none
        real, dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: theta_prime
        real, dimension(ims:ime,kms:kme,jms:jme), intent(in) :: theta_ref
        real, dimension(kms:kme), intent(in) :: theta_min_k, theta_max_k

        integer :: i, j, k
        real :: theta_full

        !$acc parallel loop gang vector collapse(3) default(present) private(theta_full)
        do j = jts, jte
            do k = kms, kme
                do i = its, ite
                    if (theta_min_k(k) > theta_max_k(k)) cycle
                    theta_full = theta_prime(i,k,j) + theta_ref(i,k,j)
                    if (theta_full < theta_min_k(k)) then
                        theta_prime(i,k,j) = theta_min_k(k) - theta_ref(i,k,j)
                    else if (theta_full > theta_max_k(k)) then
                        theta_prime(i,k,j) = theta_max_k(k) - theta_ref(i,k,j)
                    endif
                enddo
            enddo
        enddo

    end subroutine adv_bound_theta_perturbation

    subroutine bound_theta_full(theta, theta_min_k, theta_max_k)
        implicit none
        real, dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: theta
        real, dimension(kms:kme), intent(in) :: theta_min_k, theta_max_k

        integer :: i, j, k

        !$acc parallel loop gang vector collapse(3) default(present)
        do j = jts, jte
            do k = kms, kme
                do i = its, ite
                    if (theta_min_k(k) <= theta_max_k(k)) then
                        theta(i,k,j) = min(theta_max_k(k), max(theta_min_k(k), theta(i,k,j)))
                    endif
                enddo
            enddo
        enddo

    end subroutine bound_theta_full

end module advection
