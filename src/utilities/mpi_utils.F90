module mpi_utils_module
    use mpi, only: MPI_ADDRESS_KIND, MPI_Comm_Rank, MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, &
                   MPI_IN_PLACE, MPI_INT, MPI_MIN, MPI_SUM
    implicit none

contains

    subroutine allreduce_double_sum_in_place(values, count, comm)
        implicit none
        double precision, intent(inout) :: values(:)
        integer, intent(in) :: count, comm
        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, values, count, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    end subroutine allreduce_double_sum_in_place

    subroutine allreduce_integer_min_in_place(value, comm)
        implicit none
        integer, intent(inout) :: value
        integer, intent(in) :: comm
        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, value, 1, MPI_INT, MPI_MIN, comm, ierr)
    end subroutine allreduce_integer_min_in_place

    subroutine put_integer(origin, origin_count, origin_type, target_rank, target_disp, target_count, target_type, win)
        implicit none
        integer, intent(in) :: origin
        integer, intent(in) :: origin_count, origin_type, target_rank, target_count, target_type, win
        integer(kind=MPI_ADDRESS_KIND), intent(in) :: target_disp
        integer :: ierr

        call MPI_Put(origin, origin_count, origin_type, target_rank, target_disp, target_count, target_type, win, ierr)
    end subroutine put_integer

    function get_mpi_global_rank()
        implicit none
        integer :: get_mpi_global_rank, ierr

        ! Get the rank of this MPI process on the global communicator
        call MPI_Comm_rank(MPI_COMM_WORLD, get_mpi_global_rank, ierr)

    end function get_mpi_global_rank

end module mpi_utils_module
