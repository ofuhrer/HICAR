module sparse_lbc_reader
    use iso_fortran_env, only : real64
    use netcdf
    use time_object, only : Time_type, canonical_time_seconds
    use icar_constants, only : kMAX_FILE_LENGTH
    implicit none
    private

    integer, parameter, public :: LBC_MASS = 1, LBC_U = 2, LBC_V = 3
    integer, parameter :: NAME_LENGTH = 16
    real(real64), parameter :: TIME_TOLERANCE = 1.0e-3_real64

    type, public :: sparse_lbc_points_t
        integer :: global_count = 0
        integer, allocatable :: i(:), j(:), file_position(:)
        real, allocatable :: weight(:)
    contains
        procedure :: release => release_points
    end type sparse_lbc_points_t

    type, public :: sparse_lbc_field_t
        character(len=NAME_LENGTH) :: name = ""
        integer :: grid_kind = LBC_MASS
        real, allocatable :: left(:,:), right(:,:)
    contains
        procedure :: release => release_field
    end type sparse_lbc_field_t

    type, public :: sparse_lbc_reader_t
        logical :: active = .False.
        integer :: left_index = 0
        integer :: right_index = 0
        integer :: nz = 0
        real :: relaxation_timescale_seconds = 3600.0
        character(len=256) :: static_sha256 = ""
        character(len=256) :: target_grid_fingerprint = ""
        character(len=256) :: relaxation_profile = ""
        character(len=256) :: lateral_w_policy = ""
        character(len=kMAX_FILE_LENGTH), allocatable :: files(:)
        real(real64), allocatable :: valid_seconds(:)
        type(sparse_lbc_points_t) :: mass, u, v
        type(sparse_lbc_field_t), allocatable :: fields(:)
    contains
        procedure :: init
        procedure :: ensure_right_time
        procedure :: field_index
        procedure :: release => release_reader
    end type sparse_lbc_reader_t

contains

    subroutine nc_check(status, context)
        integer, intent(in) :: status
        character(len=*), intent(in) :: context
        if (status /= nf90_noerr) then
            write(*,*) "Sparse LBC NetCDF error: ", trim(context)
            write(*,*) trim(nf90_strerror(status))
            error stop "Invalid sparse LBC input"
        endif
    end subroutine nc_check

    function global_attribute(ncid, name) result(value)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        character(len=256) :: value
        value = ""
        call nc_check(nf90_get_att(ncid, nf90_global, trim(name), value), &
                      "reading global attribute "//trim(name))
        value = trim(value)
    end function global_attribute

    function iso_time_seconds(value) result(seconds)
        character(len=*), intent(in) :: value
        real(real64) :: seconds
        type(Time_type) :: parsed
        integer :: year, month, day, hour, minute, second, ios
        character(len=64) :: normalized

        normalized = adjustl(value)
        if (len_trim(normalized) < 19) error stop "Sparse LBC valid_time is not ISO-8601"
        read(normalized(1:4), *, iostat=ios) year
        if (ios /= 0) error stop "Sparse LBC valid_time has an invalid year"
        read(normalized(6:7), *, iostat=ios) month
        if (ios /= 0) error stop "Sparse LBC valid_time has an invalid month"
        read(normalized(9:10), *, iostat=ios) day
        if (ios /= 0) error stop "Sparse LBC valid_time has an invalid day"
        read(normalized(12:13), *, iostat=ios) hour
        if (ios /= 0) error stop "Sparse LBC valid_time has an invalid hour"
        read(normalized(15:16), *, iostat=ios) minute
        if (ios /= 0) error stop "Sparse LBC valid_time has an invalid minute"
        read(normalized(18:19), *, iostat=ios) second
        if (ios /= 0) error stop "Sparse LBC valid_time has an invalid second"
        if (normalized(5:5) /= "-" .or. normalized(8:8) /= "-" .or. &
            (normalized(11:11) /= "T" .and. normalized(11:11) /= " ") .or. &
            normalized(14:14) /= ":" .or. normalized(17:17) /= ":") then
            error stop "Sparse LBC valid_time is not ISO-8601"
        endif
        if (len_trim(normalized) > 19) then
            if (trim(normalized(20:)) /= "Z") then
                error stop "Sparse LBC valid_time must be UTC (Z)"
            endif
        endif
        call parsed%set(year, month, day, hour, minute, second)
        seconds = canonical_time_seconds(parsed%seconds())
    end function iso_time_seconds

    subroutine validate_contract(ncid, filename)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: filename
        character(len=256) :: value

        value = global_attribute(ncid, "product_type")
        if (trim(value) /= "hicar_lateral_boundary_state") then
            write(*,*) trim(filename), ": product_type is ", trim(value)
            error stop "Not a HICAR sparse LBC product"
        endif
        value = global_attribute(ncid, "hicar_water_conversion")
        if (trim(value) /= "APPLIED_JOINT_ALL_WATER_SPECIES") then
            error stop "Sparse LBC water species are not HICAR dry-air mixing ratios"
        endif
    end subroutine validate_contract

    subroutine read_relaxation_timescale(ncid, value)
        integer, intent(in) :: ncid
        real, intent(out) :: value
        call nc_check(nf90_get_att(ncid, nf90_global, &
                      "relaxation_timescale_seconds", value), &
                      "reading relaxation_timescale_seconds")
        if (value <= 0.0) error stop "Sparse LBC relaxation timescale must be positive"
    end subroutine read_relaxation_timescale

    subroutine read_index_variable(ncid, name, values)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, allocatable, intent(out) :: values(:)
        integer :: varid, ndims, dimids(nf90_max_var_dims), n

        call nc_check(nf90_inq_varid(ncid, trim(name), varid), "finding "//trim(name))
        call nc_check(nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids), &
                      "inquiring "//trim(name))
        if (ndims /= 1) error stop "Sparse LBC index variable is not one-dimensional"
        call nc_check(nf90_inquire_dimension(ncid, dimids(1), len=n), &
                      "inquiring sparse point dimension")
        allocate(values(n))
        call nc_check(nf90_get_var(ncid, varid, values), "reading "//trim(name))
    end subroutine read_index_variable

    subroutine initialize_points(ncid, prefix, ids, ide, jds, jde, &
                                 global_nx, global_ny, points)
        integer, intent(in) :: ncid, ids, ide, jds, jde, global_nx, global_ny
        character(len=*), intent(in) :: prefix
        type(sparse_lbc_points_t), intent(inout) :: points
        integer, allocatable :: rows(:), columns(:)
        real, allocatable :: weights(:)
        integer :: varid, n, p
        character(len=64) :: row_name, column_name, weight_name

        row_name = trim(prefix)//"row"
        column_name = trim(prefix)//"column"
        weight_name = trim(prefix)//"relaxation_weight"
        call read_index_variable(ncid, trim(row_name), rows)
        call read_index_variable(ncid, trim(column_name), columns)
        if (size(rows) /= size(columns)) error stop "Sparse LBC row/column lengths differ"
        call nc_check(nf90_inq_varid(ncid, trim(weight_name), varid), &
                      "finding "//trim(weight_name))
        allocate(weights(size(rows)))
        call nc_check(nf90_get_var(ncid, varid, weights), "reading "//trim(weight_name))
        if (any(columns < 0) .or. any(columns >= global_nx) .or. &
            any(rows < 0) .or. any(rows >= global_ny)) then
            error stop "Sparse LBC support index lies outside its global grid"
        endif
        if (any(weights < 0.0) .or. any(weights > 1.0)) then
            error stop "Sparse LBC relaxation weights must lie in [0,1]"
        endif
        points%global_count = size(rows)
        n = count(columns + 1 >= ids .and. columns + 1 <= ide .and. &
                  rows + 1 >= jds .and. rows + 1 <= jde .and. weights > 0.0)
        allocate(points%i(n), points%j(n), points%file_position(n), points%weight(n))
        n = 0
        do p = 1, size(rows)
            if (columns(p) + 1 < ids .or. columns(p) + 1 > ide .or. &
                rows(p) + 1 < jds .or. rows(p) + 1 > jde .or. weights(p) <= 0.0) cycle
            n = n + 1
            points%i(n) = columns(p) + 1
            points%j(n) = rows(p) + 1
            points%file_position(n) = p
            points%weight(n) = weights(p)
        enddo
        deallocate(rows, columns, weights)
    end subroutine initialize_points

    logical function same_points(left, right) result(same)
        type(sparse_lbc_points_t), intent(in) :: left, right
        same = left%global_count == right%global_count .and. &
               size(left%i) == size(right%i)
        if (.not. same) return
        same = all(left%i == right%i) .and. all(left%j == right%j) .and. &
               all(left%file_position == right%file_position) .and. &
               all(abs(left%weight - right%weight) <= 1.0e-7)
    end function same_points

    logical function variable_present(ncid, name)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer :: varid
        variable_present = (nf90_inq_varid(ncid, trim(name), varid) == nf90_noerr)
    end function variable_present

    subroutine configure_fields(this, ncid)
        class(sparse_lbc_reader_t), intent(inout) :: this
        integer, intent(in) :: ncid
        character(len=NAME_LENGTH), parameter :: mandatory_names(7) = &
            [character(len=NAME_LENGTH) :: "T", "P", "QV", "QC", "QI", "U", "V"]
        character(len=NAME_LENGTH), parameter :: optional_names(4) = &
            [character(len=NAME_LENGTH) :: "QR", "QS", "QG", "W"]
        integer :: n, p

        do p = 1, size(mandatory_names)
            if (.not. variable_present(ncid, mandatory_names(p))) then
                write(*,*) "Missing mandatory sparse LBC field ", trim(mandatory_names(p))
                error stop "Incomplete sparse LBC state"
            endif
        enddo
        n = size(mandatory_names)
        do p = 1, size(optional_names)
            if (variable_present(ncid, optional_names(p))) n = n + 1
        enddo
        allocate(this%fields(n))
        n = 0
        do p = 1, size(mandatory_names)
            n = n + 1
            this%fields(n)%name = mandatory_names(p)
        enddo
        do p = 1, size(optional_names)
            if (.not. variable_present(ncid, optional_names(p))) cycle
            n = n + 1
            this%fields(n)%name = optional_names(p)
        enddo
        do p = 1, size(this%fields)
            select case(trim(this%fields(p)%name))
            case("U")
                this%fields(p)%grid_kind = LBC_U
            case("V")
                this%fields(p)%grid_kind = LBC_V
            case default
                this%fields(p)%grid_kind = LBC_MASS
            end select
        enddo
    end subroutine configure_fields

    subroutine read_selected_field(ncid, name, points, nz, values)
        integer, intent(in) :: ncid, nz
        character(len=*), intent(in) :: name
        type(sparse_lbc_points_t), intent(in) :: points
        real, allocatable, intent(out) :: values(:,:)
        integer :: varid, ndims, dimids(nf90_max_var_dims), dim1, dim2
        integer :: first, last, run, local_first
        real, allocatable :: buffer(:,:)

        call nc_check(nf90_inq_varid(ncid, trim(name), varid), "finding field "//trim(name))
        call nc_check(nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids), &
                      "inquiring field "//trim(name))
        if (ndims /= 2) error stop "Sparse LBC atmospheric field must be two-dimensional"
        call nc_check(nf90_inquire_dimension(ncid, dimids(1), len=dim1), &
                      "inquiring first field dimension")
        call nc_check(nf90_inquire_dimension(ncid, dimids(2), len=dim2), &
                      "inquiring second field dimension")
        if (dim2 /= nz) then
            write(*,*) trim(name), ": expected ", nz, " levels but found ", dim2
            error stop "Sparse LBC vertical dimension does not match HICAR"
        endif
        if (dim1 /= points%global_count) then
            error stop "Sparse LBC field point dimension does not match its support"
        endif
        allocate(values(size(points%file_position), nz))
        if (size(points%file_position) == 0) return
        first = 1
        do while (first <= size(points%file_position))
            last = first
            do while (last < size(points%file_position))
                if (points%file_position(last + 1) /= points%file_position(last) + 1) exit
                last = last + 1
            enddo
            run = last - first + 1
            local_first = points%file_position(first)
            allocate(buffer(run, nz))
            call nc_check(nf90_get_var(ncid, varid, buffer, &
                          start=[local_first, 1], count=[run, nz]), &
                          "reading local runs for "//trim(name))
            values(first:last, :) = buffer
            deallocate(buffer)
            first = last + 1
        enddo
        if (any(.not. (values == values))) error stop "Sparse LBC contains NaN values"
    end subroutine read_selected_field

    subroutine read_frame(this, file_index, right_side)
        class(sparse_lbc_reader_t), intent(inout) :: this
        integer, intent(in) :: file_index
        logical, intent(in) :: right_side
        integer :: ncid, p, expected
        real, allocatable :: values(:,:)

        call nc_check(nf90_open(trim(this%files(file_index)), nf90_nowrite, ncid), &
                      "opening "//trim(this%files(file_index)))
        call validate_contract(ncid, this%files(file_index))
        do p = 1, size(this%fields)
            if (.not. variable_present(ncid, this%fields(p)%name)) then
                error stop "Sparse LBC field schema changed between frames"
            endif
            select case(this%fields(p)%grid_kind)
            case(LBC_U)
                call read_selected_field(ncid, this%fields(p)%name, this%u, &
                                         this%nz, values)
            case(LBC_V)
                call read_selected_field(ncid, this%fields(p)%name, this%v, &
                                         this%nz, values)
            case default
            call read_selected_field(ncid, this%fields(p)%name, this%mass, &
                                         this%nz, values)
            end select
            if (right_side) then
                if (allocated(this%fields(p)%right)) then
                    if (any(shape(this%fields(p)%right) /= shape(values))) then
                        error stop "Sparse LBC local field shape changed between frames"
                    endif
                    this%fields(p)%right = values
                    deallocate(values)
                else
                    call move_alloc(values, this%fields(p)%right)
                endif
            else
                if (allocated(this%fields(p)%left)) then
                    if (any(shape(this%fields(p)%left) /= shape(values))) then
                        error stop "Sparse LBC local field shape changed between frames"
                    endif
                    this%fields(p)%left = values
                    deallocate(values)
                else
                    call move_alloc(values, this%fields(p)%left)
                endif
            endif
        enddo
        expected = size(this%fields)
        if (variable_present(ncid, "QR")) expected = expected - 1
        if (variable_present(ncid, "QS")) expected = expected - 1
        if (variable_present(ncid, "QG")) expected = expected - 1
        if (variable_present(ncid, "W")) expected = expected - 1
        if (expected /= 7) error stop "Sparse LBC optional-field schema changed between frames"
        call nc_check(nf90_close(ncid), "closing sparse LBC frame")
    end subroutine read_frame

    subroutine init(this, files, nz, domain_nx, domain_ny, &
                    mass_bounds, u_bounds, v_bounds, &
                    start_seconds, end_seconds, interval_seconds)
        class(sparse_lbc_reader_t), intent(inout) :: this
        character(len=*), intent(in) :: files(:)
        integer, intent(in) :: nz, domain_nx, domain_ny
        integer, intent(in) :: mass_bounds(4), u_bounds(4), v_bounds(4)
        real(real64), intent(in) :: start_seconds, end_seconds, interval_seconds
        integer :: ncid, n, bracket, frame_nx, frame_ny
        character(len=256) :: valid_time
        real(real64) :: delta
        real :: frame_timescale
        type(sparse_lbc_points_t) :: frame_mass, frame_u, frame_v

        call this%release()
        if (size(files) < 2) error stop "Sparse LBC requires at least two frames"
        allocate(this%files(size(files)), this%valid_seconds(size(files)))
        this%files = files
        this%nz = nz
        do n = 1, size(files)
            call nc_check(nf90_open(trim(files(n)), nf90_nowrite, ncid), &
                          "opening sparse LBC metadata")
            call validate_contract(ncid, files(n))
            call nc_check(nf90_get_att(ncid, nf90_global, "domain_nx", frame_nx), &
                          "reading domain_nx")
            call nc_check(nf90_get_att(ncid, nf90_global, "domain_ny", frame_ny), &
                          "reading domain_ny")
            if (frame_nx /= domain_nx .or. frame_ny /= domain_ny) then
                error stop "Sparse LBC target dimensions do not match HICAR"
            endif
            call read_relaxation_timescale(ncid, frame_timescale)
            if (n == 1) then
                this%relaxation_timescale_seconds = frame_timescale
                this%static_sha256 = global_attribute(ncid, "static_sha256")
                this%target_grid_fingerprint = &
                    global_attribute(ncid, "target_grid_fingerprint")
                this%relaxation_profile = global_attribute(ncid, "relaxation_profile")
                this%lateral_w_policy = global_attribute(ncid, "lateral_w_policy")
            else if (abs(frame_timescale - this%relaxation_timescale_seconds) > 1.0e-6) then
                error stop "Sparse LBC relaxation timescale changed between frames"
            endif
            if (n > 1) then
                if (trim(global_attribute(ncid, "static_sha256")) /= &
                    trim(this%static_sha256) .or. &
                    trim(global_attribute(ncid, "target_grid_fingerprint")) /= &
                    trim(this%target_grid_fingerprint) .or. &
                    trim(global_attribute(ncid, "relaxation_profile")) /= &
                    trim(this%relaxation_profile) .or. &
                    trim(global_attribute(ncid, "lateral_w_policy")) /= &
                    trim(this%lateral_w_policy)) then
                    error stop "Sparse LBC product contract changed between frames"
                endif
            endif
            valid_time = global_attribute(ncid, "valid_time")
            this%valid_seconds(n) = iso_time_seconds(valid_time)
            if (n == 1) then
                call initialize_points(ncid, "", mass_bounds(1), mass_bounds(2), &
                                       mass_bounds(3), mass_bounds(4), &
                                       domain_nx, domain_ny, this%mass)
                call initialize_points(ncid, "u_", u_bounds(1), u_bounds(2), &
                                       u_bounds(3), u_bounds(4), &
                                       domain_nx + 1, domain_ny, this%u)
                call initialize_points(ncid, "v_", v_bounds(1), v_bounds(2), &
                                       v_bounds(3), v_bounds(4), &
                                       domain_nx, domain_ny + 1, this%v)
                call configure_fields(this, ncid)
            else
                call initialize_points(ncid, "", mass_bounds(1), mass_bounds(2), &
                                       mass_bounds(3), mass_bounds(4), &
                                       domain_nx, domain_ny, frame_mass)
                call initialize_points(ncid, "u_", u_bounds(1), u_bounds(2), &
                                       u_bounds(3), u_bounds(4), &
                                       domain_nx + 1, domain_ny, frame_u)
                call initialize_points(ncid, "v_", v_bounds(1), v_bounds(2), &
                                       v_bounds(3), v_bounds(4), &
                                       domain_nx, domain_ny + 1, frame_v)
                if (.not. same_points(this%mass, frame_mass) .or. &
                    .not. same_points(this%u, frame_u) .or. &
                    .not. same_points(this%v, frame_v)) then
                    error stop "Sparse LBC support or weights changed between frames"
                endif
                call frame_mass%release()
                call frame_u%release()
                call frame_v%release()
            endif
            call nc_check(nf90_close(ncid), "closing sparse LBC metadata")
            if (n > 1) then
                delta = this%valid_seconds(n) - this%valid_seconds(n - 1)
                if (delta <= 0.0_real64) error stop "Sparse LBC valid times are not increasing"
                if (abs(delta - interval_seconds) > TIME_TOLERANCE) then
                    error stop "Sparse LBC cadence does not match forcing inputinterval"
                endif
            endif
        enddo
        if (start_seconds < this%valid_seconds(1) - TIME_TOLERANCE .or. &
            end_seconds > this%valid_seconds(size(files)) + TIME_TOLERANCE) then
            error stop "Sparse LBC sequence does not cover the simulation"
        endif
        bracket = 0
        do n = 1, size(files) - 1
            if (start_seconds >= this%valid_seconds(n) - TIME_TOLERANCE .and. &
                start_seconds < this%valid_seconds(n + 1) - TIME_TOLERANCE) then
                bracket = n
                exit
            endif
        enddo
        if (bracket == 0) error stop "Could not bracket HICAR start time with sparse LBCs"
        this%left_index = bracket
        this%right_index = bracket + 1
        call read_frame(this, this%left_index, .False.)
        call read_frame(this, this%right_index, .True.)
        this%active = .True.
        !$acc enter data copyin(this%mass%i, this%mass%j, this%mass%weight)
        !$acc enter data copyin(this%u%i, this%u%j, this%u%weight)
        !$acc enter data copyin(this%v%i, this%v%j, this%v%weight)
        do n = 1, size(this%fields)
            !$acc enter data copyin(this%fields(n)%left, this%fields(n)%right)
        enddo
    end subroutine init

    subroutine ensure_right_time(this, target_seconds)
        class(sparse_lbc_reader_t), intent(inout) :: this
        real(real64), intent(in) :: target_seconds
        integer :: target, n

        if (.not. this%active) return
        target = 0
        do n = 2, size(this%valid_seconds)
            if (abs(this%valid_seconds(n) - canonical_time_seconds(target_seconds)) <= &
                TIME_TOLERANCE) then
                target = n
                exit
            endif
        enddo
        if (target == 0) error stop "Forcing event has no exact sparse LBC timestamp"
        if (target == this%right_index) return
        if (target /= this%right_index + 1) then
            error stop "Sparse LBC runtime attempted to skip or reverse a bracket"
        endif
        do n = 1, size(this%fields)
            !$acc update self(this%fields(n)%right)
            this%fields(n)%left = this%fields(n)%right
            !$acc update device(this%fields(n)%left)
        enddo
        this%left_index = this%right_index
        this%right_index = target
        call read_frame(this, target, .True.)
        do n = 1, size(this%fields)
            !$acc update device(this%fields(n)%right)
        enddo
    end subroutine ensure_right_time

    integer function field_index(this, name) result(index)
        class(sparse_lbc_reader_t), intent(in) :: this
        character(len=*), intent(in) :: name
        integer :: n
        index = 0
        if (.not. allocated(this%fields)) return
        do n = 1, size(this%fields)
            if (trim(this%fields(n)%name) == trim(name)) then
                index = n
                return
            endif
        enddo
    end function field_index

    subroutine release_points(this)
        class(sparse_lbc_points_t), intent(inout) :: this
        if (allocated(this%i)) deallocate(this%i)
        if (allocated(this%j)) deallocate(this%j)
        if (allocated(this%file_position)) deallocate(this%file_position)
        if (allocated(this%weight)) deallocate(this%weight)
        this%global_count = 0
    end subroutine release_points

    subroutine release_field(this)
        class(sparse_lbc_field_t), intent(inout) :: this
        if (allocated(this%left)) deallocate(this%left)
        if (allocated(this%right)) deallocate(this%right)
        this%name = ""
    end subroutine release_field

    subroutine release_reader(this)
        class(sparse_lbc_reader_t), intent(inout) :: this
        integer :: n
        if (this%active) then
            !$acc exit data delete(this%mass%i, this%mass%j, this%mass%weight)
            !$acc exit data delete(this%u%i, this%u%j, this%u%weight)
            !$acc exit data delete(this%v%i, this%v%j, this%v%weight)
            if (allocated(this%fields)) then
                do n = 1, size(this%fields)
                    !$acc exit data delete(this%fields(n)%left, this%fields(n)%right)
                enddo
            endif
        endif
        if (allocated(this%fields)) then
            do n = 1, size(this%fields)
                call this%fields(n)%release()
            enddo
            deallocate(this%fields)
        endif
        call this%mass%release()
        call this%u%release()
        call this%v%release()
        if (allocated(this%files)) deallocate(this%files)
        if (allocated(this%valid_seconds)) deallocate(this%valid_seconds)
        this%active = .False.
        this%left_index = 0
        this%right_index = 0
        this%relaxation_timescale_seconds = 3600.0
        this%static_sha256 = ""
        this%target_grid_fingerprint = ""
        this%relaxation_profile = ""
        this%lateral_w_policy = ""
    end subroutine release_reader

end module sparse_lbc_reader
