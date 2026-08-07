!>------------------------------------------------------------
!! Test suite for small utility modules, salvaged from the ICAR-era
!! standalone tests:
!!   - string round trips        (test_string.F90)
!!   - linear_space              (test_array_utilities.F90)
!!   - var_dict_t add/get        (test_variable_dict.F90)
!!   - fftshift / ifftshift      (test_fftshift.F90 — which only
!!     printed; the assertions here are new: fftshift must equal
!!     cshift(x, n/2) and ifftshift must invert it. fftshifter is
!!     still live via linear_theory_winds.)
!!------------------------------------------------------------
module test_utilities

    use iso_fortran_env,         only : real32, real64
    use icar_constants,          only : kVARS
    use string,                  only : get_double, get_real, get_integer, str
    use array_utilities,         only : linear_space
    use variable_dict_interface, only : var_dict_t
    use variable_interface,      only : variable_t
    use grid_interface,          only : grid_t
    use fftshifter,              only : fftshift, ifftshift
    use io_routines,             only : io_read, io_write
    use mod_atm_utilities,       only : cal_cldfra3, cal_cldfra3_level, horizon_azimuth_index, &
                                         terrain_direct_shortwave, terrain_diffuse_shortwave, &
                                         terrain_reflected_shortwave
    use domain_interface,        only : auto_dz
    use options_interface,       only : options_t
    use testdrive,               only : new_unittest, unittest_type, error_type, check
    use mpi,                     only : MPI_COMM_WORLD, MPI_Comm_rank
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

    implicit none
    private

    public :: collect_utilities_suite

contains

    subroutine collect_utilities_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("string_roundtrip", test_string_roundtrip), &
            new_unittest("linear_space",     test_linear_space), &
            new_unittest("variable_dict",    test_variable_dict), &
            new_unittest("io_read_extra_dimension", test_io_read_extra_dimension), &
            new_unittest("horizon_azimuth_index", test_horizon_azimuth_index), &
            new_unittest("terrain_shortwave_components", test_terrain_shortwave_components), &
            new_unittest("manual_vertical_grid", test_manual_vertical_grid), &
            new_unittest("fftshift_1d",      test_fftshift_1d), &
            new_unittest("fftshift_2d",      test_fftshift_2d), &
            new_unittest("cloud_fraction_level", test_cloud_fraction_level) &
          ]

    end subroutine collect_utilities_suite


    subroutine test_string_roundtrip(error)
        type(error_type), allocatable, intent(out) :: error

        real(real64), parameter :: real64_datum = 0.3_real64
        real(real32), parameter :: real32_datum = 0.3_real32
        integer,      parameter :: integer_datum = 1234567
        character(len=32) :: real64_string, real32_string, integer_string

        write(real64_string,*) real64_datum
        write(real32_string,*) real32_datum
        write(integer_string,*) integer_datum

        call check(error, get_double(real64_string) == real64_datum,   "get_double round trip")
        if (allocated(error)) return
        call check(error, get_real(real32_string) == real32_datum,     "get_real round trip")
        if (allocated(error)) return
        call check(error, get_integer(integer_string) == integer_datum,"get_integer round trip")
        if (allocated(error)) return
        call check(error, adjustl(real64_string) == str(real64_datum), "str(real64) matches list-directed write")
        if (allocated(error)) return
        call check(error, adjustl(real32_string) == str(real32_datum), "str(real32) matches list-directed write")
        if (allocated(error)) return
        call check(error, adjustl(integer_string) == str(integer_datum),"str(integer) matches list-directed write")
    end subroutine test_string_roundtrip


    subroutine test_linear_space(error)
        type(error_type), allocatable, intent(out) :: error

        real, allocatable :: a(:)
        real, parameter :: vmin = 1e-4, vmax = 3.0, tol = 1e-6
        integer, parameter :: n = 100
        real :: dv
        integer :: i

        call linear_space(a, vmin, vmax, n)

        call check(error, allocated(a), "linear_space must allocate the array")
        if (allocated(error)) return
        call check(error, size(a) == n, "linear_space array has the requested size")
        if (allocated(error)) return
        call check(error, abs(a(1) - vmin) < tol, "linear_space starts at vmin")
        if (allocated(error)) return
        call check(error, abs(a(n) - vmax) < tol, "linear_space ends at vmax")
        if (allocated(error)) return

        dv = (vmax - vmin) / (n - 1)
        do i = 2, n
            if (abs((a(i) - a(i-1)) - dv) > tol) then
                call check(error, .false., "linear_space step size is not uniform")
                return
            endif
        enddo
        call check(error, .true., "unreachable")
    end subroutine test_linear_space


    subroutine test_variable_dict(error)
        type(error_type), allocatable, intent(out) :: error

        type(var_dict_t) :: var_collection
        type(variable_t) :: var1, var2, var3, output_var
        type(grid_t)     :: grid

        ! initialize requires a kVARS metadata index (test_driver calls
        ! initialize_var_constants at startup)
        call grid%set_grid_dimensions(nx=10, ny=10, nz=5)
        call var1%initialize(kVARS%potential_temperature, grid)

        call grid%set_grid_dimensions(nx=2, ny=15, nz=5)
        call var2%initialize(kVARS%pressure, grid)

        call var3%initialize(kVARS%skin_temperature, [2,5])

        ! the dictionary is keyed by kVARS id (integer) in current HICAR
        call var_collection%add_var(kVARS%potential_temperature, var1)
        call var_collection%add_var(kVARS%pressure,              var2)
        call var_collection%add_var(kVARS%skin_temperature,      var3)

        ! retrieving the middle entry must return that variable's metadata
        output_var = var_collection%get_var(kVARS%pressure)

        call check(error, output_var%three_d .eqv. var2%three_d, &
                   "var_dict returned wrong dimensionality flag")
        if (allocated(error)) return
        call check(error, size(output_var%dim_len) == size(var2%dim_len), &
                   "var_dict returned wrong number of dimensions")
        if (allocated(error)) return
        call check(error, all(output_var%dim_len == var2%dim_len), &
                   "var_dict returned wrong dim_len")
    end subroutine test_variable_dict


    subroutine test_io_read_extra_dimension(error)
        type(error_type), allocatable, intent(out) :: error

        real :: source(2,3,4,2)
        real, allocatable :: selected(:,:,:)
        character(len=128) :: filename
        integer :: i, j, k, n, rank, ierr, unit

        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        write(filename, '(".tmp_io_extra_dimension_",I0,".nc")') rank

        do n = 1, size(source, 4)
            do k = 1, size(source, 3)
                do j = 1, size(source, 2)
                    do i = 1, size(source, 1)
                        source(i,j,k,n) = real(i + 10*j + 100*k + 1000*n)
                    enddo
                enddo
            enddo
        enddo

        call io_write(filename, "sample", source)
        call io_read(filename, "sample", selected, extradim_start=2)

        call check(error, allocated(selected), "io_read must allocate the selected record")
        if (allocated(error)) return
        call check(error, all(shape(selected) == shape(source(:,:,:,2))), &
                   "io_read selected record has the wrong shape")
        if (allocated(error)) return
        call check(error, all(selected == source(:,:,:,2)), &
                   "io_read must select one trailing record without explicit spatial counts")

        open(newunit=unit, file=filename, status="old")
        close(unit, status="delete")
    end subroutine test_io_read_extra_dimension


    subroutine test_horizon_azimuth_index(error)
        type(error_type), allocatable, intent(out) :: error

        real, parameter :: deg = acos(-1.0)/180.0

        call check(error, horizon_azimuth_index(0.0, 90) == 1, &
                   "zero azimuth starts the first horizon sector")
        if (allocated(error)) return
        call check(error, horizon_azimuth_index(3.9*deg, 90) == 1, &
                   "the first four-degree sector has no off-by-one shift")
        if (allocated(error)) return
        call check(error, horizon_azimuth_index(4.1*deg, 90) == 2, &
                   "azimuth past four degrees advances to sector two")
        if (allocated(error)) return
        call check(error, horizon_azimuth_index(359.9*deg, 90) == 90, &
                   "azimuth near 360 degrees uses the last sector")
        if (allocated(error)) return
        call check(error, horizon_azimuth_index(361.0*deg, 90) == 90 .and. &
                          horizon_azimuth_index(-1.0*deg, 90) == 1, &
                   "out-of-range azimuths clamp to the horizon axis")
    end subroutine test_horizon_azimuth_index


    subroutine test_terrain_shortwave_components(error)
        type(error_type), allocatable, intent(out) :: error

        real, parameter :: tol = 1.0e-5
        real :: horizontal_direct, sin_elevation

        horizontal_direct = 400.0
        sin_elevation = 0.5

        call check(error, abs(terrain_direct_shortwave(horizontal_direct, sin_elevation, &
                   sin_elevation, 1361.0, .true.) - horizontal_direct) < tol, &
                   "flat visible terrain preserves horizontal direct shortwave")
        horizontal_direct = 1000.0
        sin_elevation = 0.8
        call check(error, abs(terrain_direct_shortwave(horizontal_direct, sin_elevation, &
                   sin_elevation, 1361.0, .true.) - horizontal_direct) < tol, &
                   "flat terrain preserves high valid RRTMGP direct shortwave")
        if (allocated(error)) return
        call check(error, terrain_direct_shortwave(horizontal_direct, sin_elevation, &
                   0.8, 1361.0, .false.) == 0.0, &
                   "horizon obstruction removes direct shortwave")
        if (allocated(error)) return
        call check(error, terrain_direct_shortwave(horizontal_direct, sin_elevation, &
                   -0.2, 1361.0, .true.) == 0.0, &
                   "back-facing slopes receive no direct shortwave")
        if (allocated(error)) return

        call check(error, abs(terrain_diffuse_shortwave(200.0, 0.75) - 150.0) < tol, &
                   "diffuse shortwave scales only with sky-view factor")
        if (allocated(error)) return
        call check(error, abs(terrain_diffuse_shortwave(200.0, 1.0) - 200.0) < tol, &
                   "open sky preserves diffuse shortwave")
        if (allocated(error)) return

        call check(error, terrain_reflected_shortwave(0.0, 100.0, 0.2) == 0.0, &
                   "open sky has no terrain-reflected shortwave")
        if (allocated(error)) return
        call check(error, abs(terrain_reflected_shortwave(0.5, 100.0, 0.2) - &
                   (50.0 / 0.9)) < tol, &
                   "terrain-reflected shortwave preserves irradiance units and correction")
    end subroutine test_terrain_shortwave_components


    subroutine test_manual_vertical_grid(error)
        type(error_type), allocatable, intent(out) :: error

        type(options_t), allocatable :: options
        real, parameter :: manual(5) = [100.0, 200.0, 300.0, 400.0, 500.0]

        allocate(options)
        call options%init()
        options%domain%auto_level = 0
        options%domain%nz = size(manual)
        if (allocated(options%domain%dz_levels)) deallocate(options%domain%dz_levels)
        allocate(options%domain%dz_levels(size(manual)), source=manual)

        call auto_dz(options)

        call check(error, allocated(options%domain%dz_levels), &
                   "manual dz_levels remain allocated")
        if (allocated(error)) return
        call check(error, size(options%domain%dz_levels) == size(manual) .and. &
                          all(options%domain%dz_levels == manual), &
                   "auto_level=0 preserves supplied dz_levels")
    end subroutine test_manual_vertical_grid


    subroutine test_fftshift_1d(error)
        type(error_type), allocatable, intent(out) :: error

        real    :: r_odd(5), r_even(6), r_orig_odd(5), r_orig_even(6)
        complex :: c_odd(5), c_orig_odd(5)
        integer :: i

        r_orig_odd  = [(real(i), i = 1, 5)]
        r_orig_even = [(real(i), i = 1, 6)]
        c_orig_odd  = cmplx(r_orig_odd, -r_orig_odd)

        ! fftshift is a cyclic shift by n/2 (matlab/numpy convention)
        r_odd = r_orig_odd
        call fftshift(r_odd)
        call check(error, all(r_odd == cshift(r_orig_odd, size(r_odd)/2)), &
                   "1d real fftshift (odd n) is not cshift(x, n/2)")
        if (allocated(error)) return
        call ifftshift(r_odd)
        call check(error, all(r_odd == r_orig_odd), &
                   "1d real ifftshift does not invert fftshift (odd n)")
        if (allocated(error)) return

        r_even = r_orig_even
        call fftshift(r_even)
        call check(error, all(r_even == cshift(r_orig_even, size(r_even)/2)), &
                   "1d real fftshift (even n) is not cshift(x, n/2)")
        if (allocated(error)) return
        call ifftshift(r_even)
        call check(error, all(r_even == r_orig_even), &
                   "1d real ifftshift does not invert fftshift (even n)")
        if (allocated(error)) return

        c_odd = c_orig_odd
        call fftshift(c_odd)
        call check(error, all(c_odd == cshift(c_orig_odd, size(c_odd)/2)), &
                   "1d complex fftshift is not cshift(x, n/2)")
        if (allocated(error)) return
        call ifftshift(c_odd)
        call check(error, all(c_odd == c_orig_odd), &
                   "1d complex ifftshift does not invert fftshift")
    end subroutine test_fftshift_1d


    subroutine test_fftshift_2d(error)
        type(error_type), allocatable, intent(out) :: error

        real    :: r2(5,6), r2_orig(5,6), r2_expect(5,6)
        complex :: c2(5,6), c2_orig(5,6)
        integer :: i, j

        do j = 1, 6
            do i = 1, 5
                r2_orig(i,j) = i + j*50.0
            enddo
        enddo
        c2_orig = cmplx(r2_orig, -r2_orig)

        ! 2d fftshift = cyclic shift by n/2 along each dimension
        r2_expect = cshift(cshift(r2_orig, size(r2_orig,1)/2, dim=1), &
                           size(r2_orig,2)/2, dim=2)

        r2 = r2_orig
        call fftshift(r2)
        call check(error, all(r2 == r2_expect), &
                   "2d real fftshift is not cshift by n/2 in both dims")
        if (allocated(error)) return
        call ifftshift(r2)
        call check(error, all(r2 == r2_orig), &
                   "2d real ifftshift does not invert fftshift")
        if (allocated(error)) return

        c2 = c2_orig
        call fftshift(c2)
        call ifftshift(c2)
        call check(error, all(c2 == c2_orig), &
                   "2d complex ifftshift does not invert fftshift")
    end subroutine test_fftshift_2d


    subroutine test_cloud_fraction_level(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nlev = 5
        real, parameter :: tol = 2.0e-7
        real :: cldfra_column(nlev), cldfra_level(nlev)
        real :: qv(nlev), qc(nlev), qi(nlev), qs(nlev)
        real :: dz(nlev), pressure(nlev), temperature(nlev)
        real :: qvs, rh, rhoa
        integer :: k

        qv = [0.010, 0.006, 0.002, 0.0005, 0.004]
        qc = [0.0, 2.0e-6, 1.0e-8, 0.0, 0.0]
        qi = [0.0, 0.0, 0.0, 2.0e-7, 0.0]
        qs = 0.0
        dz = [60.0, 100.0, 180.0, 300.0, 120.0]
        pressure = [90000.0, 80000.0, 65000.0, 45000.0, 75000.0]
        temperature = [300.0, 280.0, 260.0, 240.0, 270.0]
        cldfra_column = -1.0

        call cal_cldfra3(cldfra_column, qv, qc, qi, qs, dz, pressure, temperature, &
                         1.0, 0.2, 1.5, 1, nlev, .false., .false.)

        do k = 1, nlev
            call cal_cldfra3_level(cldfra_level(k), qvs, rh, rhoa, &
                                   qv(k), qc(k), qi(k), qs(k), dz(k), &
                                   pressure(k), temperature(k), 1.0, 0.2, 1.5, .false.)
        enddo

        call check(error, all(ieee_is_finite(cldfra_level)), &
                   "level-local cloud fractions must be finite")
        if (allocated(error)) return
        call check(error, maxval(abs(cldfra_level - cldfra_column)) < tol, &
                   "level-local cloud fraction must match the column routine")
        if (allocated(error)) return
        call check(error, cldfra_level(1) == 0.0, &
                   "warm clear air must have zero cloud fraction")
        if (allocated(error)) return
        call check(error, cldfra_level(2) == 1.0 .and. cldfra_level(4) == 1.0, &
                   "resolved condensate must produce full cloud fraction")
        if (allocated(error)) return
        call check(error, cldfra_level(3) > 0.0 .and. cldfra_level(3) < 1.0, &
                   "trace condensate must produce a bounded partial cloud fraction")
    end subroutine test_cloud_fraction_level

end module test_utilities
