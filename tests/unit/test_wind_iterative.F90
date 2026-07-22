! test_wind_iterative.F90
! Unit tests for the hand-rolled BiCGStab iterative wind solver (wind_iterative.F90).
!
! The motivating regression: the solver's Krylov halo-exchange buffers
! (east_send/east_recv/...) were allocated with non-1-based bounds (k_s-1:k_e+1,
! ...) and passed whole to the mpi_f08 MPI_Isend/Irecv assumed-rank async buffer
! dummy. gfortran rebinds such an actual's descriptor to 1-based across the call,
! so the pack/unpack then ran off the ends of the buffers. The bug only fires when
! the domain is decomposed (>1 compute rank) so a halo exchange actually runs, and
! it is silent except under -fcheck=bounds (debug). This suite therefore EXISTS to
! be run multi-rank (the CI full-test runs `mpiexec -np 4 HICAR-tester`): with >=2
! ranks at least one of the four face exchanges executes and the bug crashes here.
!
! The test also checks the solve is functional (finite winds, reduced divergence),
! so it doubles as a guard against silent halo-correctness regressions on GPU/release
! where there is no bounds check.

module test_wind_iterative

    use, intrinsic :: iso_c_binding, only : c_double
    use mpi
    use icar_constants
    use testdrive,          only : new_unittest, unittest_type, error_type, test_failed
    use domain_interface,   only : domain_t
    use options_interface,  only : options_t
    use wind,               only : wind_var_request, init_winds, calc_divergence
    use wind_iterative,     only : calc_iter_winds, finalize_iter_winds
    use wind_multilevel,    only : horizontal_transfer_t, horizontal_tile_transfer_t, horizontal_coarse_extent, &
                                    horizontal_coarse_coordinate, owned_coarse_interval, &
                                    galerkin_stencil_t, galerkin_tile_stencil_t, vertical_line_factor_t, &
                                    assemble_colored_galerkin, assemble_colored_tile_galerkin, &
                                    relax_with_vertical_lines
    use wind_multilevel_mpi, only : horizontal_halo_exchange_t
    use advection,          only : adv_var_request
    use io_routines,        only : check_file_exists
    implicit none
    private

    real, parameter :: PI = 3.14159265358979

    public :: collect_wind_iterative_suite

contains

    !> Collect all exported unit tests
    subroutine collect_wind_iterative_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("multilevel_transfer", test_multilevel_transfer), &
            new_unittest("multilevel_device", test_multilevel_device), &
            new_unittest("multilevel_halo", test_multilevel_halo), &
            new_unittest("iter_wind_solve_decomp", test_iter_wind_solve) &
            ]
    end subroutine collect_wind_iterative_suite


    !> Build the iterative wind solver on the flat test domain, drive one solve of a
    !> seeded divergent wind field, and verify it runs across the MPI decomposition
    !> (exercising exchange_krylov_halos) and reduces the divergence to a finite field.
    subroutine test_iter_wind_solve(error)
        type(error_type), allocatable, intent(out) :: error

        type(domain_t)  :: domain
        type(options_t) :: options
        real, allocatable :: div(:,:,:)
        real    :: div0_max, div1_max, div0_max_g, div1_max_g
        integer :: ierr
        integer :: ims, ime, jms, jme, kms, kme, its, ite, jts, jte
        logical :: ok
        character(len=256) :: msg

        STD_OUT_PE = .False.

        ! --- options: iterative wind solver with a CONSTANT alpha so the solve does
        !     not depend on the Froude / calc_alpha machinery, and with the optional
        !     wind add-ons (Sx, thermal, linear) disabled. ----------------------------
        call options%init()
        options%domain%init_conditions_file = '../tests/Test_Cases/domains/flat_plane_250m.nc'
        options%domain%hgt_hi = 'topo'
        options%domain%lat_hi = 'lat'
        options%domain%lon_hi = 'lon'
        if (trim(options%domain%init_conditions_file) /= '') then
            call check_file_exists(trim(options%domain%init_conditions_file), &
                message='The test domain file does not exist. Ensure the HICAR Test-Data repo was '// &
                        'installed and ../tests/Test_Cases/domains/flat_plane_250m.nc exists.')
        endif
        options%domain%dx              = 250.0
        options%domain%nz              = 20
        options%domain%sleve           = .True.
        options%domain%use_map_factors = .False.
        options%physics%advection      = kADV_STD
        options%physics%windtype       = kITERATIVE_WINDS
        options%wind%alpha_const       = 1.0       ! constant alpha -> skip Froude / calc_alpha
        options%wind%Sx                = .False.
        options%wind%thermal           = .False.
        options%wind%linear_theory     = .False.

        ! --- build the domain over MPI_COMM_WORLD (this is what splits it into tiles,
        !     so the solver's halo exchange actually runs) -----------------------------
        domain%compute_comms = MPI_COMM_WORLD
        call adv_var_request(options)              ! core dynamics vars (u/v/w/density/jaco/dz)
        call wind_var_request(options)             ! wind_alpha / w_real
        call domain%init(options, 1)
        call init_winds(domain, options)           ! -> init_iter_winds: solver state + neighbours

        ims = domain%ims; ime = domain%ime; jms = domain%jms; jme = domain%jme
        kms = domain%kms; kme = domain%kme
        its = domain%its; ite = domain%ite; jts = domain%jts; jte = domain%jte

        allocate(div(ims:ime, kms:kme, jms:jme)); div = 0.0

        ! The iterative solver operates on the forcing-tendency winds (dqdt_3d), which
        ! are auto-allocated only for variables flagged as forcing inputs (i.e. when
        ! the forcing %uvar/%vvar/%wvar names are set). This test drives the solver
        ! directly with no forcing read, so allocate u/v/w dqdt_3d to mirror data_3d.
        call ensure_dqdt(domain%var_indx(kVARS%u)%v)
        call ensure_dqdt(domain%var_indx(kVARS%v)%v)
        call ensure_dqdt(domain%var_indx(kVARS%w)%v)

        ! --- seed a smooth, tile-continuous, DIVERGENT forcing wind in dqdt_3d. HICAR
        !     memory indices are global, so a function of (i,j) is automatically
        !     continuous across tile boundaries and the halo cells are well-defined. --
        call seed_divergent_winds()
        domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d    = 1.0
        domain%vars_3d(domain%var_indx(kVARS%wind_alpha)%v)%data_3d = 1.0

        ! These tendency arrays are allocated locally by this fixture rather
        ! than by the normal forcing reader.  The production reader enters
        ! them into OpenACC data; mirror that contract before the correction
        ! kernel's PRESENT clauses are evaluated.
        !$acc enter data copyin(domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d, &
        !$acc                       domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d)

        !$acc data copy(kVARS) copy(div)
        ! divergence of the seeded field
        call calc_divergence(div, domain, advect_density=.False., horz_only=.False., use_dqdt=.True.)
        !$acc update host(div)
        div0_max = maxval(abs(div(its:ite, kms:kme, jts:jte)))

        ! THE solve: drives bicgstab_solve -> exchange_krylov_halos on every iteration
        call calc_iter_winds(domain, &
            domain%vars_3d(domain%var_indx(kVARS%wind_alpha)%v)%data_3d, div, .False.)

        ! divergence of the corrected field
        call calc_divergence(div, domain, advect_density=.False., horz_only=.False., use_dqdt=.True.)
        !$acc update host(div)
        div1_max = maxval(abs(div(its:ite, kms:kme, jts:jte)))
        !$acc end data
        !$acc exit data delete(domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d, &
        !$acc                      domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d)

        ! reduce to a global max so the assertion is decomposition-independent
        call MPI_Allreduce(div0_max, div0_max_g, 1, MPI_REAL, MPI_MAX, domain%compute_comms, ierr)
        call MPI_Allreduce(div1_max, div1_max_g, 1, MPI_REAL, MPI_MAX, domain%compute_comms, ierr)

        ! --- evaluate BEFORE tearing down (so cleanup always runs) -------------------
        ok = .True.; msg = ''
        if (div1_max_g /= div1_max_g) then                       ! NaN
            ok = .False.
            msg = 'corrected wind divergence is NaN'
        else if (.not. (div1_max_g < div0_max_g)) then
            ! The single uncalibrated analytic-operator solve should still reduce the
            ! peak divergence. If this proves marginal in CI, relax to a factor
            ! (e.g. div1 < 0.9*div0) rather than removing it -- a broken halo crashes
            ! under -fcheck before reaching here anyway.
            ok = .False.
            write(msg,'(A,ES12.4,A,ES12.4)') &
                'iterative solver did not reduce divergence: initial=', div0_max_g, &
                ' final=', div1_max_g
        endif

        ! reset the solver's module-level state so it does not leak into other suites
        call finalize_iter_winds()
        call domain%release()

        if (.not. ok) call test_failed(error, "test_iter_wind_solve", trim(msg))

    contains

        !> Allocate a variable's dqdt_3d to exactly mirror its data_3d bounds (incl.
        !> the staggered +1 face) and zero it, if not already allocated.
        subroutine ensure_dqdt(vidx)
            integer, intent(in) :: vidx
            integer :: l1,u1,l2,u2,l3,u3
            if (.not. allocated(domain%vars_3d(vidx)%dqdt_3d)) then
                l1 = lbound(domain%vars_3d(vidx)%data_3d,1); u1 = ubound(domain%vars_3d(vidx)%data_3d,1)
                l2 = lbound(domain%vars_3d(vidx)%data_3d,2); u2 = ubound(domain%vars_3d(vidx)%data_3d,2)
                l3 = lbound(domain%vars_3d(vidx)%data_3d,3); u3 = ubound(domain%vars_3d(vidx)%data_3d,3)
                allocate(domain%vars_3d(vidx)%dqdt_3d(l1:u1, l2:u2, l3:u3))
            endif
            domain%vars_3d(vidx)%dqdt_3d = 0.0
        end subroutine ensure_dqdt

        !> Fill u/v/w dqdt with a smooth, globally-continuous field whose horizontal
        !> divergence is non-zero. Loop over each array's OWN bounds (the staggered
        !> u/v carry an extra face) so the seeding itself is bounds-clean under -fcheck.
        subroutine seed_divergent_winds()
            integer :: ii, jj, kk
            real    :: lx, ly

            lx = real(max(domain%grid%ide - domain%grid%ids, 1))
            ly = real(max(domain%grid%jde - domain%grid%jds, 1))

            associate(ud => domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d)
                do jj = lbound(ud,3), ubound(ud,3)
                    do kk = lbound(ud,2), ubound(ud,2)
                        do ii = lbound(ud,1), ubound(ud,1)
                            ud(ii,kk,jj) = 8.0 * cos(2.0*PI*real(ii - domain%grid%ids)/lx)
                        enddo
                    enddo
                enddo
            end associate

            associate(vd => domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d)
                do jj = lbound(vd,3), ubound(vd,3)
                    do kk = lbound(vd,2), ubound(vd,2)
                        do ii = lbound(vd,1), ubound(vd,1)
                            vd(ii,kk,jj) = 6.0 * sin(2.0*PI*real(jj - domain%grid%jds)/ly)
                        enddo
                    enddo
                enddo
            end associate

            domain%vars_3d(domain%var_indx(kVARS%w)%v)%dqdt_3d = 0.0
        end subroutine seed_divergent_winds

    end subroutine test_iter_wind_solve


    subroutine test_multilevel_transfer(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: nx = 8, ny = 7, nz = 5
        type(horizontal_transfer_t) :: transfer, free_transfer
        type(horizontal_tile_transfer_t) :: tile_transfer
        type(galerkin_stencil_t) :: stencil
        type(galerkin_tile_stencil_t) :: tile_stencil
        type(vertical_line_factor_t) :: line_factor
        real(c_double), allocatable :: coarse_x(:,:,:), coarse_y(:,:,:), coarse_r(:,:,:), coarse_stencil(:,:,:)
        real(c_double), allocatable :: line_rhs(:,:,:), line_x(:,:,:), line_ax(:,:,:)
        real(c_double), allocatable :: coarse_weight(:,:,:), free_weight(:,:,:), free_coarse(:,:,:)
        real(c_double), allocatable :: fine_x(:,:,:), fine_y(:,:,:), fine_v(:,:,:), fine_a_x(:,:,:)
        real(c_double), allocatable :: fine_weight(:,:,:), free_fine(:,:,:)
        real(c_double), allocatable :: tile_coarse_halo(:,:,:), tile_fine_halo(:,:,:), tile_weight_halo(:,:,:)
        real(c_double), allocatable :: tile_fine(:,:,:), tile_coarse_weight(:,:,:), tile_coarse_r(:,:,:)
        real(c_double), allocatable :: stencil_tile_halo(:,:,:), stencil_tile_result(:,:,:)
        real(c_double) :: lhs, rhs, scale, boundary_max, constant_error, minimum_pivot
        integer :: i, j, k, dk, line_status, rank_case, coarse_first, coarse_count
        integer :: gi, gj, gc_i, gc_j
        integer :: owner_count(0:4)
        integer, parameter :: partition_first(3) = [0, 2, 5]
        integer, parameter :: partition_count(3) = [2, 3, 3]

        if (horizontal_coarse_extent(nx) /= 5 .or. horizontal_coarse_extent(ny) /= 4) then
            call test_failed(error, 'test_multilevel_transfer', 'coarse extent does not retain both boundaries')
            return
        endif

        owner_count = 0
        do rank_case = 1, size(partition_first)
            call owned_coarse_interval(nx, partition_first(rank_case), partition_count(rank_case), &
                                       coarse_first, coarse_count)
            do i = coarse_first, coarse_first+coarse_count-1
                if (horizontal_coarse_coordinate(i,nx) < partition_first(rank_case) .or. &
                    horizontal_coarse_coordinate(i,nx) >= partition_first(rank_case)+partition_count(rank_case)) then
                    call test_failed(error, 'test_multilevel_transfer', 'coarse point assigned outside owning fine tile')
                    return
                endif
                owner_count(i) = owner_count(i) + 1
            enddo
        enddo
        if (any(owner_count /= 1)) then
            call test_failed(error, 'test_multilevel_transfer', 'distributed coarse ownership is not unique and complete')
            return
        endif

        call transfer%init(nx, ny, fix_lateral_boundaries=.true., fix_vertical_boundaries=.true.)
        allocate(coarse_x(transfer%nx_c,nz,transfer%ny_c), coarse_y(transfer%nx_c,nz,transfer%ny_c), &
                 coarse_r(transfer%nx_c,nz,transfer%ny_c), coarse_stencil(transfer%nx_c,nz,transfer%ny_c), &
                 coarse_weight(transfer%nx_c,nz,transfer%ny_c), &
                 line_rhs(transfer%nx_c,nz,transfer%ny_c), line_x(transfer%nx_c,nz,transfer%ny_c), &
                 line_ax(transfer%nx_c,nz,transfer%ny_c), &
                 stencil_tile_halo(0:transfer%nx_c+1,nz,0:transfer%ny_c+1), &
                 stencil_tile_result(transfer%nx_c,nz,transfer%ny_c))
        allocate(fine_x(nx,nz,ny), fine_y(nx,nz,ny), fine_v(nx,nz,ny), fine_a_x(nx,nz,ny), &
                 fine_weight(nx,nz,ny))

        do j = 1, transfer%ny_c
            do k = 1, nz
                do i = 1, transfer%nx_c
                    coarse_x(i,k,j) = sin(0.31_c_double*real(i,c_double)) + &
                                      0.07_c_double*real(k*j,c_double)
                    coarse_y(i,k,j) = cos(0.23_c_double*real(i*j,c_double)) - &
                                      0.04_c_double*real(k,c_double)
                enddo
            enddo
        enddo
        call zero_fixed_coarse(coarse_x)
        call zero_fixed_coarse(coarse_y)

        do j = 1, ny
            do k = 1, nz
                do i = 1, nx
                    fine_weight(i,k,j) = 1.0_c_double + 0.03_c_double*real(i,c_double) + &
                                         0.02_c_double*real(j,c_double) + 0.01_c_double*real(k,c_double)
                    fine_v(i,k,j) = sin(0.17_c_double*real(i+2*j+3*k,c_double))
                enddo
            enddo
        enddo

        call transfer%build_coarse_weights(fine_weight, coarse_weight)
        call transfer%prolong(coarse_x, fine_x)
        call transfer%prolong(coarse_y, fine_y)
        call transfer%restrict_adjoint(fine_v, fine_weight, coarse_weight, coarse_r)

        lhs = sum(fine_weight * fine_x * fine_v)
        rhs = sum(coarse_weight * coarse_x * coarse_r)
        scale = max(1.0_c_double, abs(lhs), abs(rhs))
        if (abs(lhs-rhs) > 5.0e-13_c_double*scale) then
            call test_failed(error, 'test_multilevel_transfer', 'metric-weighted transfer adjointness failed')
            call transfer%release()
            return
        endif

        ! Decompose the same transfer into three uneven x tiles.  Halo-filled
        ! local P and R must reproduce the global result exactly and every
        ! coarse row must be formed only by its unique owner.
        do rank_case = 1, size(partition_first)
            call tile_transfer%init(nx, ny, partition_first(rank_case), partition_count(rank_case), 0, ny, &
                                    fix_lateral_boundaries=.true., fix_vertical_boundaries=.true.)
            allocate(tile_coarse_halo(0:tile_transfer%nx_c_local+1,nz,0:tile_transfer%ny_c_local+1), &
                     tile_fine_halo(0:tile_transfer%nx_f_local+1,nz,0:tile_transfer%ny_f_local+1), &
                     tile_weight_halo(0:tile_transfer%nx_f_local+1,nz,0:tile_transfer%ny_f_local+1), &
                     tile_fine(tile_transfer%nx_f_local,nz,tile_transfer%ny_f_local), &
                     tile_coarse_weight(tile_transfer%nx_c_local,nz,tile_transfer%ny_c_local), &
                     tile_coarse_r(tile_transfer%nx_c_local,nz,tile_transfer%ny_c_local))
            tile_coarse_halo = 0.0_c_double
            do j = 0, tile_transfer%ny_c_local+1
                gc_j = tile_transfer%y_c_first+j-1
                if (gc_j < 0 .or. gc_j >= transfer%ny_c) cycle
                do i = 0, tile_transfer%nx_c_local+1
                    gc_i = tile_transfer%x_c_first+i-1
                    if (gc_i < 0 .or. gc_i >= transfer%nx_c) cycle
                    tile_coarse_halo(i,:,j) = coarse_x(gc_i+1,:,gc_j+1)
                enddo
            enddo
            call tile_transfer%prolong_owned(tile_coarse_halo, tile_fine)
            scale = max(1.0_c_double, maxval(abs(fine_x(partition_first(rank_case)+1: &
                        partition_first(rank_case)+partition_count(rank_case),:,:))))
            if (maxval(abs(tile_fine-fine_x(partition_first(rank_case)+1: &
                    partition_first(rank_case)+partition_count(rank_case),:,:))) > 5.0e-14_c_double*scale) then
                call test_failed(error, 'test_multilevel_transfer', 'tile prolongation differs from global prolongation')
                call tile_transfer%release()
                return
            endif

            tile_fine_halo = 0.0_c_double
            tile_weight_halo = 1.0_c_double
            do j = 0, tile_transfer%ny_f_local+1
                gj = tile_transfer%y_f_first+j-1
                if (gj < 0 .or. gj >= ny) cycle
                do i = 0, tile_transfer%nx_f_local+1
                    gi = tile_transfer%x_f_first+i-1
                    if (gi < 0 .or. gi >= nx) cycle
                    tile_fine_halo(i,:,j) = fine_v(gi+1,:,gj+1)
                    tile_weight_halo(i,:,j) = fine_weight(gi+1,:,gj+1)
                enddo
            enddo
            call tile_transfer%build_owned_coarse_weights(tile_weight_halo, tile_coarse_weight)
            call tile_transfer%restrict_owned_adjoint(tile_fine_halo, tile_weight_halo, tile_coarse_weight, tile_coarse_r)
            do j = 1, tile_transfer%ny_c_local
                gc_j = tile_transfer%y_c_first+j
                do i = 1, tile_transfer%nx_c_local
                    gc_i = tile_transfer%x_c_first+i
                    if (maxval(abs(tile_coarse_weight(i,:,j)-coarse_weight(gc_i,:,gc_j))) > 5.0e-14_c_double .or. &
                        maxval(abs(tile_coarse_r(i,:,j)-coarse_r(gc_i,:,gc_j))) > 5.0e-14_c_double) then
                        call test_failed(error, 'test_multilevel_transfer', 'tile restriction differs from global adjoint')
                        call tile_transfer%release()
                        return
                    endif
                enddo
            enddo
            deallocate(tile_coarse_halo, tile_fine_halo, tile_weight_halo, tile_fine, &
                       tile_coarse_weight, tile_coarse_r)
            call tile_transfer%release()
        enddo

        boundary_max = max(maxval(abs(fine_x(1,:,:))), maxval(abs(fine_x(nx,:,:))), &
                           maxval(abs(fine_x(:,:,1))), maxval(abs(fine_x(:,:,ny))), &
                           maxval(abs(fine_x(:,1,:))), maxval(abs(fine_x(:,nz,:))))
        if (boundary_max /= 0.0_c_double) then
            call test_failed(error, 'test_multilevel_transfer', 'prolongation changed a fixed identity boundary')
            call transfer%release()
            return
        endif

        ! Verify the actual Petrov-Galerkin composition, including a deliberately
        ! nonsymmetric fine-grid operator: <Py,A Px>_f = <y,R A P x>_c.
        call apply_test_operator(fine_x, fine_a_x)
        call transfer%restrict_adjoint(fine_a_x, fine_weight, coarse_weight, coarse_r)
        lhs = sum(fine_weight * fine_y * fine_a_x)
        rhs = sum(coarse_weight * coarse_y * coarse_r)
        scale = max(1.0_c_double, abs(lhs), abs(rhs))
        if (abs(lhs-rhs) > 5.0e-13_c_double*scale) then
            call test_failed(error, 'test_multilevel_transfer', 'Petrov-Galerkin composition identity failed')
            call transfer%release()
            return
        endif

        ! Assemble R A P without coefficient rediscretization.  The colored
        ! probe must reproduce a direct matrix-free application even though
        ! the test operator is nonsymmetric and the grid has one even extent.
        call assemble_colored_galerkin(transfer, fine_weight, coarse_weight, apply_test_operator, stencil)
        call stencil%apply(coarse_x, coarse_stencil)
        scale = max(1.0_c_double, maxval(abs(coarse_r)), maxval(abs(coarse_stencil)))
        if (maxval(abs(coarse_r-coarse_stencil)) > 2.0e-12_c_double*scale) then
            call test_failed(error, 'test_multilevel_transfer', 'colored Galerkin stencil does not equal R A P')
            call stencil%release()
            call transfer%release()
            return
        endif
        call assemble_colored_tile_galerkin(transfer%nx_c, transfer%ny_c, 0, 0, transfer%nx_c, &
            transfer%ny_c, nz, .true., .true., apply_test_galerkin, tile_stencil)
        if (maxval(abs(tile_stencil%value-stencil%value)) > 2.0e-12_c_double*scale) then
            call test_failed(error, 'test_multilevel_transfer', 'tile colored assembly differs from global assembly')
            call tile_stencil%release()
            call stencil%release()
            call transfer%release()
            return
        endif
        stencil_tile_halo = 0.0_c_double
        stencil_tile_halo(1:transfer%nx_c,:,1:transfer%ny_c) = coarse_x
        call tile_stencil%apply_owned(stencil_tile_halo, stencil_tile_result)
        if (maxval(abs(stencil_tile_result-coarse_r)) > 2.0e-12_c_double*scale) then
            call test_failed(error, 'test_multilevel_transfer', 'tile Galerkin application differs from R A P')
            call tile_stencil%release()
            call stencil%release()
            call transfer%release()
            return
        endif

        call line_factor%factorize(stencil, line_status, minimum_pivot)
        if (line_status /= 0 .or. minimum_pivot <= 0.0_c_double) then
            call test_failed(error, 'test_multilevel_transfer', 'coarse vertical block factorization failed')
            call stencil%release()
            call transfer%release()
            return
        endif
        do j = 1, transfer%ny_c
            do k = 1, nz
                do i = 1, transfer%nx_c
                    line_rhs(i,k,j) = cos(0.13_c_double*real(3*i+2*k+j,c_double))
                enddo
            enddo
        enddo
        call line_factor%apply(line_rhs, line_x)
        line_ax = 0.0_c_double
        do j = 1, transfer%ny_c
            do k = 1, nz
                do i = 1, transfer%nx_c
                    do dk = -1, 1
                        if (k+dk < 1 .or. k+dk > nz) cycle
                        line_ax(i,k,j) = line_ax(i,k,j) + &
                            stencil%value(0,dk,0,i,k,j)*line_x(i,k+dk,j)
                    enddo
                enddo
            enddo
        enddo
        scale = max(1.0_c_double, maxval(abs(line_rhs)))
        if (maxval(abs(line_ax-line_rhs)) > 2.0e-12_c_double*scale) then
            call test_failed(error, 'test_multilevel_transfer', 'exact coarse vertical line solve failed')
            call line_factor%release()
            call stencil%release()
            call transfer%release()
            return
        endif
        line_x = 0.0_c_double
        call relax_with_vertical_lines(stencil, line_factor, line_rhs, line_x, line_ax, coarse_stencil, &
                                       n_sweeps=4, omega=0.7_c_double)
        call stencil%apply(line_x, line_ax)
        lhs = sqrt(sum((line_rhs-line_ax)**2))
        rhs = sqrt(sum(line_rhs**2))
        if (lhs >= 0.5_c_double*rhs) then
            call test_failed(error, 'test_multilevel_transfer', 'vertical-line relaxation did not reduce residual')
            call line_factor%release()
            call stencil%release()
            call transfer%release()
            return
        endif

        ! Without fixed identity rows the normalized adjoint restriction must
        ! preserve constants on both odd and even horizontal extents.
        call free_transfer%init(nx, ny, fix_lateral_boundaries=.false., fix_vertical_boundaries=.false.)
        allocate(free_coarse(free_transfer%nx_c,nz,free_transfer%ny_c), &
                 free_weight(free_transfer%nx_c,nz,free_transfer%ny_c), free_fine(nx,nz,ny))
        free_coarse = 1.0_c_double
        call free_transfer%prolong(free_coarse, free_fine)
        constant_error = maxval(abs(free_fine - 1.0_c_double))
        call free_transfer%build_coarse_weights(fine_weight, free_weight)
        call free_transfer%restrict_adjoint(free_fine, fine_weight, free_weight, free_coarse)
        constant_error = max(constant_error, maxval(abs(free_coarse - 1.0_c_double)))
        if (constant_error > 5.0e-14_c_double) then
            call test_failed(error, 'test_multilevel_transfer', 'unconstrained transfer does not preserve constants')
        endif

        call free_transfer%release()
        call line_factor%release()
        call tile_stencil%release()
        call stencil%release()
        call transfer%release()

    contains

        subroutine zero_fixed_coarse(field)
            real(c_double), intent(inout) :: field(:,:,:)
            field(1,:,:) = 0.0_c_double
            field(size(field,1),:,:) = 0.0_c_double
            field(:,:,1) = 0.0_c_double
            field(:,:,size(field,3)) = 0.0_c_double
            field(:,1,:) = 0.0_c_double
            field(:,size(field,2),:) = 0.0_c_double
        end subroutine zero_fixed_coarse

        subroutine apply_test_operator(x, ax)
            real(c_double), intent(in) :: x(:,:,:)
            real(c_double), intent(out) :: ax(:,:,:)
            integer :: ii, jj, kk

            ax = 0.0_c_double
            do jj = 2, size(x,3)-1
                do kk = 2, size(x,2)-1
                    do ii = 2, size(x,1)-1
                        ax(ii,kk,jj) = 3.7_c_double*x(ii,kk,jj) &
                            - 0.8_c_double*x(ii-1,kk,jj) + 0.35_c_double*x(ii+1,kk,jj) &
                            - 0.6_c_double*x(ii,kk,jj-1) + 0.15_c_double*x(ii,kk,jj+1) &
                            - 1.1_c_double*x(ii,kk-1,jj) + 0.45_c_double*x(ii,kk+1,jj)
                    enddo
                enddo
            enddo
        end subroutine apply_test_operator

        subroutine apply_test_galerkin(x, ax)
            real(c_double), intent(in) :: x(:,:,:)
            real(c_double), intent(out) :: ax(:,:,:)

            call transfer%prolong(x, fine_x)
            call apply_test_operator(fine_x, fine_a_x)
            call transfer%restrict_adjoint(fine_a_x, fine_weight, coarse_weight, ax)
        end subroutine apply_test_galerkin

    end subroutine test_multilevel_transfer


    subroutine test_multilevel_device(error)
        type(error_type), allocatable, intent(out) :: error
#ifdef _OPENACC
        integer, parameter :: nx = 8, ny = 7, nz = 5
        type(horizontal_tile_transfer_t) :: transfer
        type(galerkin_tile_stencil_t) :: stencil
        type(vertical_line_factor_t) :: line_factor
        real(c_double), allocatable :: coarse_halo(:,:,:), fine(:,:,:), fine_reference(:,:,:)
        real(c_double), allocatable :: fine_halo(:,:,:), fine_weight(:,:,:)
        real(c_double), allocatable :: coarse_weight(:,:,:), coarse_weight_reference(:,:,:)
        real(c_double), allocatable :: coarse_r(:,:,:), coarse_r_reference(:,:,:)
        real(c_double), allocatable :: stencil_x(:,:,:), stencil_ax(:,:,:), stencil_reference(:,:,:)
        real(c_double), allocatable :: line_rhs(:,:,:), line_x(:,:,:), line_reference(:,:,:)
        real(c_double) :: scale, minimum_pivot
        integer :: i, j, k, line_status

        call transfer%init(nx, ny, 0, nx, 0, ny, &
                           fix_lateral_boundaries=.true., fix_vertical_boundaries=.true.)
        allocate(coarse_halo(0:transfer%nx_c_local+1,nz,0:transfer%ny_c_local+1), &
                 fine(nx,nz,ny), fine_reference(nx,nz,ny), &
                 fine_halo(0:nx+1,nz,0:ny+1), fine_weight(0:nx+1,nz,0:ny+1), &
                 coarse_weight(transfer%nx_c_local,nz,transfer%ny_c_local), &
                 coarse_weight_reference(transfer%nx_c_local,nz,transfer%ny_c_local), &
                 coarse_r(transfer%nx_c_local,nz,transfer%ny_c_local), &
                 coarse_r_reference(transfer%nx_c_local,nz,transfer%ny_c_local))
        do j = 0, transfer%ny_c_local+1
            do k = 1, nz
                do i = 0, transfer%nx_c_local+1
                    coarse_halo(i,k,j) = sin(0.11_c_double*real(3*i+2*k+j,c_double))
                enddo
            enddo
        enddo
        do j = 0, ny+1
            do k = 1, nz
                do i = 0, nx+1
                    fine_halo(i,k,j) = cos(0.09_c_double*real(i+3*k+2*j,c_double))
                    fine_weight(i,k,j) = 1.0_c_double + 0.01_c_double*real(i+k+j,c_double)
                enddo
            enddo
        enddo
        call transfer%prolong_owned(coarse_halo, fine_reference)
        call transfer%build_owned_coarse_weights(fine_weight, coarse_weight_reference)
        call transfer%restrict_owned_adjoint(fine_halo, fine_weight, coarse_weight_reference, coarse_r_reference)

        call transfer%upload_device()
        !$acc enter data copyin(coarse_halo, fine_halo, fine_weight) &
        !$acc            create(fine, coarse_weight, coarse_r)
        call transfer%prolong_owned_device(coarse_halo, fine)
        call transfer%build_owned_coarse_weights_device(fine_weight, coarse_weight)
        call transfer%restrict_owned_adjoint_device(fine_halo, fine_weight, coarse_weight, coarse_r)
        !$acc update self(fine, coarse_weight, coarse_r)
        !$acc exit data delete(coarse_halo, fine_halo, fine_weight, fine, coarse_weight, coarse_r)
        scale = max(1.0_c_double, maxval(abs(fine_reference)), maxval(abs(coarse_r_reference)))
        if (maxval(abs(fine-fine_reference)) > 5.0e-13_c_double*scale .or. &
            maxval(abs(coarse_weight-coarse_weight_reference)) > 5.0e-13_c_double*scale .or. &
            maxval(abs(coarse_r-coarse_r_reference)) > 5.0e-13_c_double*scale) then
            call test_failed(error, 'test_multilevel_device', 'device transfer differs from host reference')
            call transfer%release()
            return
        endif

        stencil%nx_global = transfer%nx_c_global
        stencil%ny_global = transfer%ny_c_global
        stencil%x_first = 0
        stencil%y_first = 0
        stencil%nx = transfer%nx_c_local
        stencil%ny = transfer%ny_c_local
        stencil%nz = nz
        stencil%fix_lateral_boundaries = .true.
        stencil%fix_vertical_boundaries = .true.
        allocate(stencil%value(-1:1,-1:1,-1:1,stencil%nx,nz,stencil%ny), &
                 stencil_x(0:stencil%nx+1,nz,0:stencil%ny+1), &
                 stencil_ax(stencil%nx,nz,stencil%ny), stencil_reference(stencil%nx,nz,stencil%ny), &
                 line_rhs(stencil%nx,nz,stencil%ny), line_x(stencil%nx,nz,stencil%ny), &
                 line_reference(stencil%nx,nz,stencil%ny))
        stencil%value = 0.0_c_double
        stencil%value(0,0,0,:,:,:) = 4.0_c_double
        stencil%value(0,-1,0,:,:,:) = -0.45_c_double
        stencil%value(0,1,0,:,:,:) = -0.35_c_double
        stencil%value(-1,0,0,:,:,:) = -0.20_c_double
        stencil%value(1,0,0,:,:,:) = -0.10_c_double
        stencil%value(0,0,-1,:,:,:) = -0.15_c_double
        stencil%value(0,0,1,:,:,:) = -0.05_c_double
        do j = 1, stencil%ny
            do k = 1, nz
                do i = 1, stencil%nx
                    if (i == 1 .or. i == stencil%nx .or. j == 1 .or. j == stencil%ny .or. &
                        k == 1 .or. k == nz) then
                        stencil%value(:,:,:,i,k,j) = 0.0_c_double
                        stencil%value(0,0,0,i,k,j) = 1.0_c_double
                    endif
                enddo
            enddo
        enddo
        do j = 0, stencil%ny+1
            do k = 1, nz
                do i = 0, stencil%nx+1
                    stencil_x(i,k,j) = sin(0.07_c_double*real(2*i+5*k+3*j,c_double))
                enddo
            enddo
        enddo
        call stencil%apply_owned(stencil_x, stencil_reference)
        call line_factor%factorize(stencil, line_status, minimum_pivot)
        if (line_status /= 0 .or. minimum_pivot <= 0.0_c_double) then
            call test_failed(error, 'test_multilevel_device', 'device-test line factorization failed')
            call stencil%release()
            call transfer%release()
            return
        endif
        do j = 1, stencil%ny
            do k = 1, nz
                do i = 1, stencil%nx
                    line_rhs(i,k,j) = cos(0.13_c_double*real(i+2*k+4*j,c_double))
                enddo
            enddo
        enddo
        call line_factor%apply(line_rhs, line_reference)

        call stencil%upload_device()
        call line_factor%upload_device()
        !$acc enter data copyin(stencil_x, line_rhs) create(stencil_ax, line_x)
        call stencil%apply_owned_device(stencil_x, stencil_ax)
        call line_factor%apply_device(line_rhs, line_x)
        !$acc update self(stencil_ax, line_x)
        !$acc exit data delete(stencil_x, line_rhs, stencil_ax, line_x)
        scale = max(1.0_c_double, maxval(abs(stencil_reference)), maxval(abs(line_reference)))
        if (maxval(abs(stencil_ax-stencil_reference)) > 5.0e-13_c_double*scale .or. &
            maxval(abs(line_x-line_reference)) > 5.0e-13_c_double*scale) then
            call test_failed(error, 'test_multilevel_device', 'device coarse kernel differs from host reference')
        endif
        call line_factor%release()
        call stencil%release()
        call transfer%release()
#endif
    end subroutine test_multilevel_device


    subroutine test_multilevel_halo(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: nx_global = 8, ny_global = 8, nz = 3
        integer, parameter :: nx_local = 4, ny_local = 4
        type(horizontal_halo_exchange_t) :: halo
        real(c_double), allocatable :: field(:,:,:)
        real(c_double) :: expected
        integer :: rank, nprocs, ierr, rx, ry, west, east, south, north
        integer :: i, j, k, gi, gj

        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
        if (nprocs /= 4) return

        rx = modulo(rank,2)
        ry = rank/2
        west = merge(rank-1, MPI_PROC_NULL, rx > 0)
        east = merge(rank+1, MPI_PROC_NULL, rx < 1)
        south = merge(rank-2, MPI_PROC_NULL, ry > 0)
        north = merge(rank+2, MPI_PROC_NULL, ry < 1)
        call halo%init(nx_local, ny_local, nz, MPI_COMM_WORLD, west, east, south, north)
        allocate(field(0:nx_local+1,nz,0:ny_local+1))
        field = 0.0_c_double
        do j = 1, ny_local
            gj = ry*ny_local+j-1
            do k = 1, nz
                do i = 1, nx_local
                    gi = rx*nx_local+i-1
                    field(i,k,j) = real(1000*k+100*gj+gi,c_double)
                enddo
            enddo
        enddo
        call halo%exchange(field)

        do j = 0, ny_local+1
            do k = 1, nz
                do i = 0, nx_local+1
                    if (i >= 1 .and. i <= nx_local .and. j >= 1 .and. j <= ny_local) cycle
                    gi = rx*nx_local+i-1
                    gj = ry*ny_local+j-1
                    if (gi < 0 .or. gi >= nx_global .or. gj < 0 .or. gj >= ny_global) then
                        expected = 0.0_c_double
                    else
                        expected = real(1000*k+100*gj+gi,c_double)
                    endif
                    if (field(i,k,j) /= expected) then
                        call test_failed(error, 'test_multilevel_halo', 'two-stage halo or corner value is incorrect')
                        call halo%release()
                        return
                    endif
                enddo
            enddo
        enddo
#ifdef _OPENACC
        field = 0.0_c_double
        do j = 1, ny_local
            gj = ry*ny_local+j-1
            do k = 1, nz
                do i = 1, nx_local
                    gi = rx*nx_local+i-1
                    field(i,k,j) = real(1000*k+100*gj+gi,c_double)
                enddo
            enddo
        enddo
        call halo%upload_device()
        !$acc enter data copyin(field)
        call halo%exchange_device(field)
        !$acc update self(field)
        !$acc exit data delete(field)
        do j = 0, ny_local+1
            do k = 1, nz
                do i = 0, nx_local+1
                    if (i >= 1 .and. i <= nx_local .and. j >= 1 .and. j <= ny_local) cycle
                    gi = rx*nx_local+i-1
                    gj = ry*ny_local+j-1
                    if (gi < 0 .or. gi >= nx_global .or. gj < 0 .or. gj >= ny_global) then
                        expected = 0.0_c_double
                    else
                        expected = real(1000*k+100*gj+gi,c_double)
                    endif
                    if (field(i,k,j) /= expected) then
                        call test_failed(error, 'test_multilevel_halo', 'device halo or corner value is incorrect')
                        call halo%release()
                        return
                    endif
                enddo
            enddo
        enddo
#endif
        call halo%release()
    end subroutine test_multilevel_halo

end module test_wind_iterative
