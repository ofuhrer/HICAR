module wind_multilevel
    use, intrinsic :: iso_c_binding, only : c_double
    implicit none
    private

    type, public :: horizontal_transfer_t
        integer :: nx_f = 0
        integer :: ny_f = 0
        integer :: nx_c = 0
        integer :: ny_c = 0
        logical :: fix_lateral_boundaries = .true.
        logical :: fix_vertical_boundaries = .true.
        integer, allocatable :: i_lo(:), i_hi(:), j_lo(:), j_hi(:)
        real(c_double), allocatable :: i_hi_weight(:), j_hi_weight(:)
    contains
        procedure :: init => init_horizontal_transfer
        procedure :: release => release_horizontal_transfer
        procedure :: prolong => prolong_horizontal
        procedure :: build_coarse_weights
        procedure :: restrict_adjoint
    end type horizontal_transfer_t

    type, public :: galerkin_stencil_t
        integer :: nx = 0
        integer :: ny = 0
        integer :: nz = 0
        logical :: fix_lateral_boundaries = .true.
        logical :: fix_vertical_boundaries = .true.
        real(c_double), allocatable :: value(:,:,:,:,:,:)
    contains
        procedure :: release => release_galerkin_stencil
        procedure :: apply => apply_galerkin_stencil
    end type galerkin_stencil_t

    type, public :: vertical_line_factor_t
        integer :: nx = 0
        integer :: ny = 0
        integer :: nz = 0
        logical :: ready = .false.
        real(c_double), allocatable :: diagonal_inverse(:,:,:)
        real(c_double), allocatable :: upper_prime(:,:,:)
        real(c_double), allocatable :: lower_coefficient(:,:,:)
    contains
        procedure :: factorize => factorize_vertical_lines
        procedure :: release => release_vertical_line_factor
        procedure :: apply => apply_vertical_line_factor
    end type vertical_line_factor_t

    abstract interface
        subroutine fine_operator_apply(x, ax)
            import c_double
            real(c_double), intent(in) :: x(:,:,:)
            real(c_double), intent(out) :: ax(:,:,:)
        end subroutine fine_operator_apply
    end interface

    public :: horizontal_coarse_extent, horizontal_coarse_coordinate
    public :: horizontal_coarse_bracket, owned_coarse_interval
    public :: assemble_colored_galerkin

contains

    pure integer function horizontal_coarse_extent(n_f) result(n_c)
        integer, intent(in) :: n_f

        if (n_f <= 1) then
            n_c = n_f
        else
            ! Retain both physical boundary points.  For an even fine-grid
            ! extent the final coarse interval has length one; this avoids
            ! dropping the east/north identity row.
            n_c = n_f / 2 + 1
        endif
    end function horizontal_coarse_extent


    pure integer function horizontal_coarse_coordinate(c, n_f) result(f)
        integer, intent(in) :: c, n_f

        f = min(2*c, n_f-1)
    end function horizontal_coarse_coordinate


    pure subroutine horizontal_coarse_bracket(f, n_f, lo, hi, hi_weight)
        integer, intent(in) :: f, n_f
        integer, intent(out) :: lo, hi
        real(c_double), intent(out) :: hi_weight
        integer :: n_c, left_coordinate, right_coordinate

        n_c = horizontal_coarse_extent(n_f)
        if (f == n_f-1) then
            lo = n_c-1
            hi = n_c-1
            hi_weight = 0.0_c_double
            return
        endif

        lo = f / 2
        hi = min(lo + 1, n_c-1)
        left_coordinate = horizontal_coarse_coordinate(lo, n_f)
        right_coordinate = horizontal_coarse_coordinate(hi, n_f)
        if (right_coordinate == left_coordinate) then
            hi_weight = 0.0_c_double
        else
            hi_weight = real(f-left_coordinate,c_double) / &
                        real(right_coordinate-left_coordinate,c_double)
        endif
    end subroutine horizontal_coarse_bracket


    pure subroutine owned_coarse_interval(n_f, fine_first, fine_count, coarse_first, coarse_count)
        integer, intent(in) :: n_f, fine_first, fine_count
        integer, intent(out) :: coarse_first, coarse_count
        integer :: c, coordinate, n_c

        n_c = horizontal_coarse_extent(n_f)
        coarse_first = 0
        coarse_count = 0
        do c = 0, n_c-1
            coordinate = horizontal_coarse_coordinate(c, n_f)
            if (coordinate >= fine_first .and. coordinate < fine_first+fine_count) then
                if (coarse_count == 0) coarse_first = c
                coarse_count = coarse_count + 1
            endif
        enddo
    end subroutine owned_coarse_interval


    subroutine init_horizontal_transfer(this, nx_f, ny_f, fix_lateral_boundaries, fix_vertical_boundaries)
        class(horizontal_transfer_t), intent(inout) :: this
        integer, intent(in) :: nx_f, ny_f
        logical, intent(in), optional :: fix_lateral_boundaries, fix_vertical_boundaries

        if (nx_f < 2 .or. ny_f < 2) error stop 'horizontal transfer requires at least 2x2 points'

        call this%release()
        this%nx_f = nx_f
        this%ny_f = ny_f
        this%nx_c = horizontal_coarse_extent(nx_f)
        this%ny_c = horizontal_coarse_extent(ny_f)
        if (present(fix_lateral_boundaries)) this%fix_lateral_boundaries = fix_lateral_boundaries
        if (present(fix_vertical_boundaries)) this%fix_vertical_boundaries = fix_vertical_boundaries

        allocate(this%i_lo(nx_f), this%i_hi(nx_f), this%i_hi_weight(nx_f))
        allocate(this%j_lo(ny_f), this%j_hi(ny_f), this%j_hi_weight(ny_f))
        call build_axis_map(nx_f, this%i_lo, this%i_hi, this%i_hi_weight)
        call build_axis_map(ny_f, this%j_lo, this%j_hi, this%j_hi_weight)
    end subroutine init_horizontal_transfer


    subroutine release_horizontal_transfer(this)
        class(horizontal_transfer_t), intent(inout) :: this

        if (allocated(this%i_lo)) deallocate(this%i_lo, this%i_hi, this%i_hi_weight)
        if (allocated(this%j_lo)) deallocate(this%j_lo, this%j_hi, this%j_hi_weight)
        this%nx_f = 0
        this%ny_f = 0
        this%nx_c = 0
        this%ny_c = 0
        this%fix_lateral_boundaries = .true.
        this%fix_vertical_boundaries = .true.
    end subroutine release_horizontal_transfer


    pure subroutine build_axis_map(n_f, lo, hi, hi_weight)
        integer, intent(in) :: n_f
        integer, intent(out) :: lo(n_f), hi(n_f)
        real(c_double), intent(out) :: hi_weight(n_f)
        integer :: f, lo_zero, hi_zero

        do f = 1, n_f
            call horizontal_coarse_bracket(f-1, n_f, lo_zero, hi_zero, hi_weight(f))
            lo(f) = lo_zero + 1
            hi(f) = hi_zero + 1
        enddo
    end subroutine build_axis_map


    subroutine prolong_horizontal(this, coarse, fine)
        class(horizontal_transfer_t), intent(in) :: this
        real(c_double), intent(in) :: coarse(:,:,:)
        real(c_double), intent(out) :: fine(:,:,:)
        integer :: i, j, k, il, ih, jl, jh
        real(c_double) :: tx, ty

        call require_transfer_shapes(this, coarse, fine)

        do j = 1, this%ny_f
            jl = this%j_lo(j)
            jh = this%j_hi(j)
            ty = this%j_hi_weight(j)
            do k = 1, size(fine, 2)
                do i = 1, this%nx_f
                    if (is_fixed_fine_point(this, i, k, j, size(fine, 2))) then
                        fine(i,k,j) = 0.0_c_double
                    else
                        il = this%i_lo(i)
                        ih = this%i_hi(i)
                        tx = this%i_hi_weight(i)
                        fine(i,k,j) = &
                            (1.0_c_double-tx) * (1.0_c_double-ty) * coarse(il,k,jl) + &
                            tx                * (1.0_c_double-ty) * coarse(ih,k,jl) + &
                            (1.0_c_double-tx) * ty                * coarse(il,k,jh) + &
                            tx                * ty                * coarse(ih,k,jh)
                    endif
                enddo
            enddo
        enddo
    end subroutine prolong_horizontal


    subroutine build_coarse_weights(this, fine_weight, coarse_weight)
        class(horizontal_transfer_t), intent(in) :: this
        real(c_double), intent(in) :: fine_weight(:,:,:)
        real(c_double), intent(out) :: coarse_weight(:,:,:)
        integer :: i, j, k, il, ih, jl, jh
        real(c_double) :: tx, ty, w

        call require_transfer_shapes(this, coarse_weight, fine_weight)
        if (any(fine_weight <= 0.0_c_double)) error stop 'fine transfer weights must be positive'

        coarse_weight = 0.0_c_double
        do j = 1, this%ny_f
            jl = this%j_lo(j)
            jh = this%j_hi(j)
            ty = this%j_hi_weight(j)
            do k = 1, size(fine_weight, 2)
                do i = 1, this%nx_f
                    if (is_fixed_fine_point(this, i, k, j, size(fine_weight, 2))) cycle
                    il = this%i_lo(i)
                    ih = this%i_hi(i)
                    tx = this%i_hi_weight(i)
                    w = fine_weight(i,k,j)
                    coarse_weight(il,k,jl) = coarse_weight(il,k,jl) + &
                        w * (1.0_c_double-tx) * (1.0_c_double-ty)
                    coarse_weight(ih,k,jl) = coarse_weight(ih,k,jl) + &
                        w * tx * (1.0_c_double-ty)
                    coarse_weight(il,k,jh) = coarse_weight(il,k,jh) + &
                        w * (1.0_c_double-tx) * ty
                    coarse_weight(ih,k,jh) = coarse_weight(ih,k,jh) + w * tx * ty
                enddo
            enddo
        enddo

        call set_fixed_coarse_points(this, coarse_weight, 1.0_c_double)
        if (any(coarse_weight <= 0.0_c_double)) error stop 'coarse transfer weight is not positive'
    end subroutine build_coarse_weights


    subroutine restrict_adjoint(this, fine, fine_weight, coarse_weight, coarse)
        class(horizontal_transfer_t), intent(in) :: this
        real(c_double), intent(in) :: fine(:,:,:), fine_weight(:,:,:), coarse_weight(:,:,:)
        real(c_double), intent(out) :: coarse(:,:,:)
        integer :: i, j, k, il, ih, jl, jh
        real(c_double) :: tx, ty, weighted_value

        call require_transfer_shapes(this, coarse, fine)
        call require_transfer_shapes(this, coarse_weight, fine_weight)
        if (any(coarse_weight <= 0.0_c_double)) error stop 'coarse transfer weights must be positive'

        coarse = 0.0_c_double
        do j = 1, this%ny_f
            jl = this%j_lo(j)
            jh = this%j_hi(j)
            ty = this%j_hi_weight(j)
            do k = 1, size(fine, 2)
                do i = 1, this%nx_f
                    if (is_fixed_fine_point(this, i, k, j, size(fine, 2))) cycle
                    il = this%i_lo(i)
                    ih = this%i_hi(i)
                    tx = this%i_hi_weight(i)
                    weighted_value = fine_weight(i,k,j) * fine(i,k,j)
                    coarse(il,k,jl) = coarse(il,k,jl) + &
                        weighted_value * (1.0_c_double-tx) * (1.0_c_double-ty)
                    coarse(ih,k,jl) = coarse(ih,k,jl) + &
                        weighted_value * tx * (1.0_c_double-ty)
                    coarse(il,k,jh) = coarse(il,k,jh) + &
                        weighted_value * (1.0_c_double-tx) * ty
                    coarse(ih,k,jh) = coarse(ih,k,jh) + weighted_value * tx * ty
                enddo
            enddo
        enddo

        coarse = coarse / coarse_weight
        call set_fixed_coarse_points(this, coarse, 0.0_c_double)
    end subroutine restrict_adjoint


    pure logical function is_fixed_fine_point(this, i, k, j, nz) result(fixed)
        class(horizontal_transfer_t), intent(in) :: this
        integer, intent(in) :: i, k, j, nz

        fixed = (this%fix_lateral_boundaries .and. &
                 (i == 1 .or. i == this%nx_f .or. j == 1 .or. j == this%ny_f)) .or. &
                (this%fix_vertical_boundaries .and. (k == 1 .or. k == nz))
    end function is_fixed_fine_point


    subroutine set_fixed_coarse_points(this, field, value)
        class(horizontal_transfer_t), intent(in) :: this
        real(c_double), intent(inout) :: field(:,:,:)
        real(c_double), intent(in) :: value

        if (this%fix_lateral_boundaries) then
            field(1,:,:) = value
            field(this%nx_c,:,:) = value
            field(:,:,1) = value
            field(:,:,this%ny_c) = value
        endif
        if (this%fix_vertical_boundaries) then
            field(:,1,:) = value
            field(:,size(field,2),:) = value
        endif
    end subroutine set_fixed_coarse_points


    subroutine require_transfer_shapes(this, coarse, fine)
        class(horizontal_transfer_t), intent(in) :: this
        real(c_double), intent(in) :: coarse(:,:,:), fine(:,:,:)

        if (size(fine,1) /= this%nx_f .or. size(fine,3) /= this%ny_f) &
            error stop 'fine array does not match horizontal transfer'
        if (size(coarse,1) /= this%nx_c .or. size(coarse,3) /= this%ny_c) &
            error stop 'coarse array does not match horizontal transfer'
        if (size(coarse,2) /= size(fine,2)) error stop 'horizontal transfer cannot coarsen vertically'
    end subroutine require_transfer_shapes


    subroutine assemble_colored_galerkin(transfer, fine_weight, coarse_weight, apply_fine, stencil)
        type(horizontal_transfer_t), intent(in) :: transfer
        real(c_double), intent(in) :: fine_weight(:,:,:), coarse_weight(:,:,:)
        procedure(fine_operator_apply) :: apply_fine
        type(galerkin_stencil_t), intent(inout) :: stencil
        real(c_double), allocatable :: coarse_probe(:,:,:), coarse_response(:,:,:)
        real(c_double), allocatable :: fine_probe(:,:,:), fine_response(:,:,:)
        integer :: ci, cj, ck, i, j, k, di, dj, dk

        call require_transfer_shapes(transfer, coarse_weight, fine_weight)
        call stencil%release()
        stencil%nx = transfer%nx_c
        stencil%ny = transfer%ny_c
        stencil%nz = size(fine_weight, 2)
        stencil%fix_lateral_boundaries = transfer%fix_lateral_boundaries
        stencil%fix_vertical_boundaries = transfer%fix_vertical_boundaries
        allocate(stencil%value(-1:1,-1:1,-1:1,stencil%nx,stencil%nz,stencil%ny))
        allocate(coarse_probe(stencil%nx,stencil%nz,stencil%ny), &
                 coarse_response(stencil%nx,stencil%nz,stencil%ny))
        allocate(fine_probe(transfer%nx_f,stencil%nz,transfer%ny_f), &
                 fine_response(transfer%nx_f,stencil%nz,transfer%ny_f))
        stencil%value = 0.0_c_double

        ! A 3x3x3 coloring recovers an exact nearest-neighbour coarse stencil
        ! in 27 operator applications.  Bilinear horizontal P, its weighted
        ! adjoint R, and a 3x3x3 fine stencil cannot connect coarse points
        ! farther than one index apart.  Each row therefore sees at most one
        ! active probe point of a given color, including for nonsymmetric A.
        do cj = 0, 2
            do ck = 0, 2
                do ci = 0, 2
                    coarse_probe = 0.0_c_double
                    do j = 1, stencil%ny
                        do k = 1, stencil%nz
                            do i = 1, stencil%nx
                                if (is_fixed_coarse_point(transfer, i, k, j, stencil%nz)) cycle
                                if (modulo(i-1,3) == ci .and. modulo(k-1,3) == ck .and. &
                                    modulo(j-1,3) == cj) coarse_probe(i,k,j) = 1.0_c_double
                            enddo
                        enddo
                    enddo
                    call transfer%prolong(coarse_probe, fine_probe)
                    call apply_fine(fine_probe, fine_response)
                    call transfer%restrict_adjoint(fine_response, fine_weight, coarse_weight, coarse_response)

                    do j = 1, stencil%ny
                        dj = modulo(cj - modulo(j-1,3) + 1, 3) - 1
                        do k = 1, stencil%nz
                            dk = modulo(ck - modulo(k-1,3) + 1, 3) - 1
                            do i = 1, stencil%nx
                                if (is_fixed_coarse_point(transfer, i, k, j, stencil%nz)) cycle
                                di = modulo(ci - modulo(i-1,3) + 1, 3) - 1
                                if (i+di < 1 .or. i+di > stencil%nx .or. &
                                    k+dk < 1 .or. k+dk > stencil%nz .or. &
                                    j+dj < 1 .or. j+dj > stencil%ny) cycle
                                if (is_fixed_coarse_point(transfer, i+di, k+dk, j+dj, stencil%nz)) cycle
                                stencil%value(di,dk,dj,i,k,j) = coarse_response(i,k,j)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        do j = 1, stencil%ny
            do k = 1, stencil%nz
                do i = 1, stencil%nx
                    if (is_fixed_coarse_point(transfer, i, k, j, stencil%nz)) &
                        stencil%value(0,0,0,i,k,j) = 1.0_c_double
                enddo
            enddo
        enddo
    end subroutine assemble_colored_galerkin


    subroutine release_galerkin_stencil(this)
        class(galerkin_stencil_t), intent(inout) :: this

        if (allocated(this%value)) deallocate(this%value)
        this%nx = 0
        this%ny = 0
        this%nz = 0
        this%fix_lateral_boundaries = .true.
        this%fix_vertical_boundaries = .true.
    end subroutine release_galerkin_stencil


    subroutine apply_galerkin_stencil(this, x, ax)
        class(galerkin_stencil_t), intent(in) :: this
        real(c_double), intent(in) :: x(:,:,:)
        real(c_double), intent(out) :: ax(:,:,:)
        integer :: i, j, k, di, dj, dk

        if (size(x,1) /= this%nx .or. size(x,2) /= this%nz .or. size(x,3) /= this%ny) &
            error stop 'input does not match Galerkin stencil'
        if (any(shape(ax) /= shape(x))) error stop 'Galerkin output shape mismatch'

        ax = 0.0_c_double
        do j = 1, this%ny
            do k = 1, this%nz
                do i = 1, this%nx
                    if (is_fixed_stencil_point(this, i, k, j)) then
                        ax(i,k,j) = x(i,k,j)
                    else
                        do dj = -1, 1
                            if (j+dj < 1 .or. j+dj > this%ny) cycle
                            do dk = -1, 1
                                if (k+dk < 1 .or. k+dk > this%nz) cycle
                                do di = -1, 1
                                    if (i+di < 1 .or. i+di > this%nx) cycle
                                    ax(i,k,j) = ax(i,k,j) + this%value(di,dk,dj,i,k,j) * &
                                                               x(i+di,k+dk,j+dj)
                                enddo
                            enddo
                        enddo
                    endif
                enddo
            enddo
        enddo
    end subroutine apply_galerkin_stencil


    subroutine factorize_vertical_lines(this, stencil, status, minimum_pivot)
        class(vertical_line_factor_t), intent(inout) :: this
        type(galerkin_stencil_t), intent(in) :: stencil
        integer, intent(out) :: status
        real(c_double), intent(out), optional :: minimum_pivot
        real(c_double), parameter :: pivot_floor = 1.0e-28_c_double
        real(c_double) :: diagonal, lower, upper, pivot, min_pivot
        integer :: i, j, k

        call this%release()
        this%nx = stencil%nx
        this%ny = stencil%ny
        this%nz = stencil%nz
        allocate(this%diagonal_inverse(this%nx,this%nz,this%ny), &
                 this%upper_prime(this%nx,this%nz,this%ny), &
                 this%lower_coefficient(this%nx,this%nz,this%ny))
        this%diagonal_inverse = 0.0_c_double
        this%upper_prime = 0.0_c_double
        this%lower_coefficient = 0.0_c_double
        min_pivot = huge(1.0_c_double)
        status = 0

        do j = 1, this%ny
            do i = 1, this%nx
                diagonal = stencil%value(0,0,0,i,1,j)
                upper = stencil%value(0,1,0,i,1,j)
                min_pivot = min(min_pivot, abs(diagonal))
                if (abs(diagonal) <= pivot_floor) then
                    status = 1
                    cycle
                endif
                this%diagonal_inverse(i,1,j) = 1.0_c_double / diagonal
                this%upper_prime(i,1,j) = upper * this%diagonal_inverse(i,1,j)

                do k = 2, this%nz
                    lower = stencil%value(0,-1,0,i,k,j)
                    diagonal = stencil%value(0,0,0,i,k,j)
                    upper = stencil%value(0,1,0,i,k,j)
                    pivot = diagonal - lower * this%upper_prime(i,k-1,j)
                    min_pivot = min(min_pivot, abs(pivot))
                    if (abs(pivot) <= pivot_floor) then
                        status = 1
                        exit
                    endif
                    this%diagonal_inverse(i,k,j) = 1.0_c_double / pivot
                    this%upper_prime(i,k,j) = upper * this%diagonal_inverse(i,k,j)
                    this%lower_coefficient(i,k,j) = lower
                enddo
            enddo
        enddo

        this%ready = (status == 0)
        if (present(minimum_pivot)) minimum_pivot = min_pivot
    end subroutine factorize_vertical_lines


    subroutine release_vertical_line_factor(this)
        class(vertical_line_factor_t), intent(inout) :: this

        if (allocated(this%diagonal_inverse)) deallocate(this%diagonal_inverse)
        if (allocated(this%upper_prime)) deallocate(this%upper_prime)
        if (allocated(this%lower_coefficient)) deallocate(this%lower_coefficient)
        this%nx = 0
        this%ny = 0
        this%nz = 0
        this%ready = .false.
    end subroutine release_vertical_line_factor


    subroutine apply_vertical_line_factor(this, rhs, x)
        class(vertical_line_factor_t), intent(in) :: this
        real(c_double), intent(in) :: rhs(:,:,:)
        real(c_double), intent(out) :: x(:,:,:)
        integer :: i, j, k

        if (.not. this%ready) error stop 'vertical line factor is not ready'
        if (size(rhs,1) /= this%nx .or. size(rhs,2) /= this%nz .or. size(rhs,3) /= this%ny) &
            error stop 'vertical line right-hand side shape mismatch'
        if (any(shape(x) /= shape(rhs))) error stop 'vertical line output shape mismatch'

        do j = 1, this%ny
            do i = 1, this%nx
                x(i,1,j) = this%diagonal_inverse(i,1,j) * rhs(i,1,j)
                do k = 2, this%nz
                    x(i,k,j) = this%diagonal_inverse(i,k,j) * &
                        (rhs(i,k,j) - this%lower_coefficient(i,k,j) * x(i,k-1,j))
                enddo
                do k = this%nz-1, 1, -1
                    x(i,k,j) = x(i,k,j) - this%upper_prime(i,k,j) * x(i,k+1,j)
                enddo
            enddo
        enddo
    end subroutine apply_vertical_line_factor


    pure logical function is_fixed_coarse_point(transfer, i, k, j, nz) result(fixed)
        type(horizontal_transfer_t), intent(in) :: transfer
        integer, intent(in) :: i, k, j, nz

        fixed = (transfer%fix_lateral_boundaries .and. &
                 (i == 1 .or. i == transfer%nx_c .or. j == 1 .or. j == transfer%ny_c)) .or. &
                (transfer%fix_vertical_boundaries .and. (k == 1 .or. k == nz))
    end function is_fixed_coarse_point


    pure logical function is_fixed_stencil_point(stencil, i, k, j) result(fixed)
        class(galerkin_stencil_t), intent(in) :: stencil
        integer, intent(in) :: i, k, j

        fixed = (stencil%fix_lateral_boundaries .and. &
                 (i == 1 .or. i == stencil%nx .or. j == 1 .or. j == stencil%ny)) .or. &
                (stencil%fix_vertical_boundaries .and. (k == 1 .or. k == stencil%nz))
    end function is_fixed_stencil_point

end module wind_multilevel
