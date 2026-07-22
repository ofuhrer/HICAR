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

    public :: horizontal_coarse_extent

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
        call build_axis_map(nx_f, this%nx_c, this%i_lo, this%i_hi, this%i_hi_weight)
        call build_axis_map(ny_f, this%ny_c, this%j_lo, this%j_hi, this%j_hi_weight)
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


    pure subroutine build_axis_map(n_f, n_c, lo, hi, hi_weight)
        integer, intent(in) :: n_f, n_c
        integer, intent(out) :: lo(n_f), hi(n_f)
        real(c_double), intent(out) :: hi_weight(n_f)
        integer :: f, left_coordinate, right_coordinate

        do f = 1, n_f
            if (f == n_f) then
                lo(f) = n_c
                hi(f) = n_c
                hi_weight(f) = 0.0_c_double
            else
                lo(f) = (f - 1) / 2 + 1
                hi(f) = min(lo(f) + 1, n_c)
                left_coordinate = 2 * (lo(f) - 1)
                if (hi(f) == n_c) then
                    right_coordinate = n_f - 1
                else
                    right_coordinate = 2 * (hi(f) - 1)
                endif
                if (right_coordinate == left_coordinate) then
                    hi_weight(f) = 0.0_c_double
                else
                    hi_weight(f) = real((f - 1) - left_coordinate, c_double) / &
                                   real(right_coordinate - left_coordinate, c_double)
                endif
            endif
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

end module wind_multilevel
