module wind_multilevel_mpi
    use, intrinsic :: iso_c_binding, only : c_double
    use mpi
    implicit none
    private

    type, public :: horizontal_halo_exchange_t
        integer :: nx = 0
        integer :: ny = 0
        integer :: nz = 0
        integer :: communicator = MPI_COMM_NULL
        integer :: west = MPI_PROC_NULL
        integer :: east = MPI_PROC_NULL
        integer :: south = MPI_PROC_NULL
        integer :: north = MPI_PROC_NULL
        logical :: device_uploaded = .false.
        real(c_double), allocatable :: west_send(:,:), west_recv(:,:)
        real(c_double), allocatable :: east_send(:,:), east_recv(:,:)
        real(c_double), allocatable :: south_send(:,:), south_recv(:,:)
        real(c_double), allocatable :: north_send(:,:), north_recv(:,:)
    contains
        procedure :: init => init_horizontal_halo_exchange
        procedure :: release => release_horizontal_halo_exchange
        procedure :: upload_device => upload_horizontal_halo_exchange
        procedure :: delete_device => delete_horizontal_halo_exchange
        procedure :: exchange => exchange_horizontal_halos
        procedure :: exchange_device => exchange_horizontal_halos_device
    end type horizontal_halo_exchange_t

contains

    subroutine init_horizontal_halo_exchange(this, nx, ny, nz, communicator, west, east, south, north)
        class(horizontal_halo_exchange_t), intent(inout) :: this
        integer, intent(in) :: nx, ny, nz, communicator, west, east, south, north

        if (nx < 1 .or. ny < 1 .or. nz < 1) error stop 'halo exchange requires a nonempty tile'
        call this%release()
        this%nx = nx; this%ny = ny; this%nz = nz
        this%communicator = communicator
        this%west = west; this%east = east; this%south = south; this%north = north
        allocate(this%west_send(nz,ny), this%west_recv(nz,ny), &
                 this%east_send(nz,ny), this%east_recv(nz,ny))
        allocate(this%south_send(0:nx+1,nz), this%south_recv(0:nx+1,nz), &
                 this%north_send(0:nx+1,nz), this%north_recv(0:nx+1,nz))
        this%west_send = 0.0_c_double; this%west_recv = 0.0_c_double
        this%east_send = 0.0_c_double; this%east_recv = 0.0_c_double
        this%south_send = 0.0_c_double; this%south_recv = 0.0_c_double
        this%north_send = 0.0_c_double; this%north_recv = 0.0_c_double
    end subroutine init_horizontal_halo_exchange


    subroutine release_horizontal_halo_exchange(this)
        class(horizontal_halo_exchange_t), intent(inout) :: this

        call this%delete_device()
        if (allocated(this%west_send)) deallocate(this%west_send, this%west_recv, this%east_send, this%east_recv)
        if (allocated(this%south_send)) deallocate(this%south_send, this%south_recv, this%north_send, this%north_recv)
        this%nx = 0; this%ny = 0; this%nz = 0
        this%communicator = MPI_COMM_NULL
        this%west = MPI_PROC_NULL; this%east = MPI_PROC_NULL
        this%south = MPI_PROC_NULL; this%north = MPI_PROC_NULL
    end subroutine release_horizontal_halo_exchange


    subroutine upload_horizontal_halo_exchange(this)
        class(horizontal_halo_exchange_t), intent(inout) :: this

        if (this%device_uploaded) return
        !$acc enter data copyin(this%west_send, this%west_recv, this%east_send, this%east_recv, &
        !$acc                   this%south_send, this%south_recv, this%north_send, this%north_recv)
        this%device_uploaded = .true.
    end subroutine upload_horizontal_halo_exchange


    subroutine delete_horizontal_halo_exchange(this)
        class(horizontal_halo_exchange_t), intent(inout) :: this

        if (.not. this%device_uploaded) return
        !$acc exit data delete(this%west_send, this%west_recv, this%east_send, this%east_recv, &
        !$acc                  this%south_send, this%south_recv, this%north_send, this%north_recv)
        this%device_uploaded = .false.
    end subroutine delete_horizontal_halo_exchange


    subroutine exchange_horizontal_halos(this, field)
        class(horizontal_halo_exchange_t), intent(inout) :: this
        real(c_double), intent(inout) :: field(0:,:,0:)
        integer, parameter :: tag_to_west = 731, tag_to_east = 732
        integer, parameter :: tag_to_south = 733, tag_to_north = 734
        integer :: i, j, k, ierr, n_requests
        integer :: requests(4), statuses(MPI_STATUS_SIZE,4)

        if (size(field,1) /= this%nx+2 .or. size(field,2) /= this%nz .or. &
            size(field,3) /= this%ny+2) error stop 'field shape does not match halo exchange'

        ! East/west first.  The subsequent north/south messages include these
        ! freshly received x halos, which propagates all four corner values
        ! without requiring diagonal-neighbour ranks.  Galerkin levels have
        ! a full 3x3 horizontal stencil, unlike the fine 15-point operator.
        do j = 1, this%ny
            do k = 1, this%nz
                this%west_send(k,j) = field(1,k,j)
                this%east_send(k,j) = field(this%nx,k,j)
            enddo
        enddo
        this%west_recv = 0.0_c_double
        this%east_recv = 0.0_c_double
        n_requests = 0
        call post_receive(this%west_recv, size(this%west_recv), this%west, tag_to_east)
        call post_receive(this%east_recv, size(this%east_recv), this%east, tag_to_west)
        call post_send(this%west_send, size(this%west_send), this%west, tag_to_west)
        call post_send(this%east_send, size(this%east_send), this%east, tag_to_east)
        if (n_requests > 0) call MPI_Waitall(n_requests, requests, statuses, ierr)
        if (this%west /= MPI_PROC_NULL) field(0,:,1:this%ny) = this%west_recv
        if (this%east /= MPI_PROC_NULL) field(this%nx+1,:,1:this%ny) = this%east_recv
        if (this%west == MPI_PROC_NULL) field(0,:,1:this%ny) = 0.0_c_double
        if (this%east == MPI_PROC_NULL) field(this%nx+1,:,1:this%ny) = 0.0_c_double

        do k = 1, this%nz
            do i = 0, this%nx+1
                this%south_send(i,k) = field(i,k,1)
                this%north_send(i,k) = field(i,k,this%ny)
            enddo
        enddo
        this%south_recv = 0.0_c_double
        this%north_recv = 0.0_c_double
        n_requests = 0
        call post_receive(this%south_recv, size(this%south_recv), this%south, tag_to_north)
        call post_receive(this%north_recv, size(this%north_recv), this%north, tag_to_south)
        call post_send(this%south_send, size(this%south_send), this%south, tag_to_south)
        call post_send(this%north_send, size(this%north_send), this%north, tag_to_north)
        if (n_requests > 0) call MPI_Waitall(n_requests, requests, statuses, ierr)
        if (this%south /= MPI_PROC_NULL) field(:,:,0) = this%south_recv
        if (this%north /= MPI_PROC_NULL) field(:,:,this%ny+1) = this%north_recv
        if (this%south == MPI_PROC_NULL) field(:,:,0) = 0.0_c_double
        if (this%north == MPI_PROC_NULL) field(:,:,this%ny+1) = 0.0_c_double

    contains

        subroutine post_receive(buffer, count, source, tag)
            real(c_double), intent(inout) :: buffer(*)
            integer, intent(in) :: count, source, tag
            if (source == MPI_PROC_NULL) return
            n_requests = n_requests + 1
            call MPI_Irecv(buffer, count, MPI_DOUBLE_PRECISION, source, tag, this%communicator, &
                           requests(n_requests), ierr)
        end subroutine post_receive

        subroutine post_send(buffer, count, destination, tag)
            real(c_double), intent(in) :: buffer(*)
            integer, intent(in) :: count, destination, tag
            if (destination == MPI_PROC_NULL) return
            n_requests = n_requests + 1
            call MPI_Isend(buffer, count, MPI_DOUBLE_PRECISION, destination, tag, this%communicator, &
                           requests(n_requests), ierr)
        end subroutine post_send

    end subroutine exchange_horizontal_halos


    subroutine exchange_horizontal_halos_device(this, field)
        class(horizontal_halo_exchange_t), intent(inout) :: this
        real(c_double), intent(inout) :: field(0:,:,0:)
        integer, parameter :: tag_to_west = 741, tag_to_east = 742
        integer, parameter :: tag_to_south = 743, tag_to_north = 744
        integer :: i, j, k, ierr, n_requests
        integer :: requests(4), statuses(MPI_STATUS_SIZE,4)

        if (size(field,1) /= this%nx+2 .or. size(field,2) /= this%nz .or. &
            size(field,3) /= this%ny+2) error stop 'field shape does not match halo exchange'
        if (.not. this%device_uploaded) error stop 'halo buffers are not on the device'

        !$acc parallel loop gang vector collapse(2) &
        !$acc present(field,this%west_send,this%east_send)
        do j = 1, this%ny
            do k = 1, this%nz
                this%west_send(k,j) = field(1,k,j)
                this%east_send(k,j) = field(this%nx,k,j)
            enddo
        enddo
        !$acc update self(this%west_send,this%east_send)
        this%west_recv = 0.0_c_double
        this%east_recv = 0.0_c_double
        n_requests = 0
        call post_receive(this%west_recv, size(this%west_recv), this%west, tag_to_east)
        call post_receive(this%east_recv, size(this%east_recv), this%east, tag_to_west)
        call post_send(this%west_send, size(this%west_send), this%west, tag_to_west)
        call post_send(this%east_send, size(this%east_send), this%east, tag_to_east)
        if (n_requests > 0) call MPI_Waitall(n_requests, requests, statuses, ierr)
        !$acc update device(this%west_recv,this%east_recv)
        !$acc parallel loop gang vector collapse(2) &
        !$acc present(field,this%west_recv,this%east_recv)
        do j = 1, this%ny
            do k = 1, this%nz
                if (this%west == MPI_PROC_NULL) then
                    field(0,k,j) = 0.0_c_double
                else
                    field(0,k,j) = this%west_recv(k,j)
                endif
                if (this%east == MPI_PROC_NULL) then
                    field(this%nx+1,k,j) = 0.0_c_double
                else
                    field(this%nx+1,k,j) = this%east_recv(k,j)
                endif
            enddo
        enddo

        !$acc parallel loop gang vector collapse(2) &
        !$acc present(field,this%south_send,this%north_send)
        do k = 1, this%nz
            do i = 0, this%nx+1
                this%south_send(i,k) = field(i,k,1)
                this%north_send(i,k) = field(i,k,this%ny)
            enddo
        enddo
        !$acc update self(this%south_send,this%north_send)
        this%south_recv = 0.0_c_double
        this%north_recv = 0.0_c_double
        n_requests = 0
        call post_receive(this%south_recv, size(this%south_recv), this%south, tag_to_north)
        call post_receive(this%north_recv, size(this%north_recv), this%north, tag_to_south)
        call post_send(this%south_send, size(this%south_send), this%south, tag_to_south)
        call post_send(this%north_send, size(this%north_send), this%north, tag_to_north)
        if (n_requests > 0) call MPI_Waitall(n_requests, requests, statuses, ierr)
        !$acc update device(this%south_recv,this%north_recv)
        !$acc parallel loop gang vector collapse(2) &
        !$acc present(field,this%south_recv,this%north_recv)
        do k = 1, this%nz
            do i = 0, this%nx+1
                if (this%south == MPI_PROC_NULL) then
                    field(i,k,0) = 0.0_c_double
                else
                    field(i,k,0) = this%south_recv(i,k)
                endif
                if (this%north == MPI_PROC_NULL) then
                    field(i,k,this%ny+1) = 0.0_c_double
                else
                    field(i,k,this%ny+1) = this%north_recv(i,k)
                endif
            enddo
        enddo

    contains

        subroutine post_receive(buffer, count, source, tag)
            real(c_double), intent(inout) :: buffer(*)
            integer, intent(in) :: count, source, tag
            if (source == MPI_PROC_NULL) return
            n_requests = n_requests+1
            call MPI_Irecv(buffer,count,MPI_DOUBLE_PRECISION,source,tag,this%communicator,requests(n_requests),ierr)
        end subroutine post_receive

        subroutine post_send(buffer, count, destination, tag)
            real(c_double), intent(in) :: buffer(*)
            integer, intent(in) :: count, destination, tag
            if (destination == MPI_PROC_NULL) return
            n_requests = n_requests+1
            call MPI_Isend(buffer,count,MPI_DOUBLE_PRECISION,destination,tag,this%communicator,requests(n_requests),ierr)
        end subroutine post_send

    end subroutine exchange_horizontal_halos_device

end module wind_multilevel_mpi
