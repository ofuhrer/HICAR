module wind_coarse_solve
    ! Terminal solver for the distributed exact Galerkin hierarchy.  Once the
    ! horizontal grid is small, agglomerate its interior rows on one rank and
    ! solve the assembled nonsymmetric operator with unrestarted, two-pass
    ! MGS GMRES.  Vertical-line right preconditioning preserves the strong
    ! column coupling without imposing symmetry on the coarse operator.
    use, intrinsic :: iso_c_binding, only : c_double
    use, intrinsic :: iso_fortran_env, only : output_unit
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use mpi
    use wind_multilevel, only : galerkin_tile_stencil_t
    implicit none
    private

    integer, parameter :: STENCIL_WIDTH = 27
    real(c_double), parameter :: COARSE_RELATIVE_TOLERANCE = 1.0e-14_c_double
    real(c_double), parameter :: COARSE_ABSOLUTE_TOLERANCE = 1.0e-15_c_double
    ! Runtime coarse right-hand sides can amplify roundoff in this conditioned
    ! operator.  This remains three orders tighter than the outer relative
    ! wind tolerance; the deterministic setup proof separately requires 1e-10.
    real(c_double), parameter :: COARSE_ACCEPTANCE_TOLERANCE = 1.0e-8_c_double
    real(c_double), parameter :: COARSE_BREAKDOWN_TOLERANCE = 1.0e-30_c_double

    type, public :: collective_coarse_solver_t
        logical :: ready = .false.
        integer :: communicator = MPI_COMM_NULL
        integer :: rank = -1
        integer :: n_ranks = 0
        integer :: root = 0
        integer :: nx_global = 0
        integer :: ny_global = 0
        integer :: nz = 0
        integer :: nx_local = 0
        integer :: ny_local = 0
        integer :: n_global = 0
        integer :: n_local = 0
        integer :: solve_count = 0
        integer, allocatable :: counts(:), displacements(:)
        integer, allocatable :: local_i(:), local_k(:), local_j(:)
        integer, allocatable :: packed_row_ids(:)
        integer, allocatable :: matrix_nnz(:), matrix_columns(:,:)
        real(c_double), allocatable :: matrix_values(:,:)
        real(c_double), allocatable :: diagonal_inverse(:)
        real(c_double), allocatable :: upper_prime(:), lower_coefficient(:)
        real(c_double), allocatable :: packed_rhs(:), packed_solution(:)
        real(c_double), allocatable :: global_rhs(:), global_solution(:)
        real(c_double), allocatable :: arnoldi_v(:,:), arnoldi_z(:,:)
        real(c_double), allocatable :: hessenberg(:,:), givens_c(:), givens_s(:)
        real(c_double), allocatable :: least_squares_rhs(:), coefficients(:)
    contains
        procedure :: setup => setup_collective_coarse_solver
        procedure :: solve => solve_collective_coarse_system
        procedure :: release => release_collective_coarse_solver
    end type collective_coarse_solver_t

contains

    subroutine setup_collective_coarse_solver(this, stencil, communicator, status)
        class(collective_coarse_solver_t), intent(inout) :: this
        type(galerkin_tile_stencil_t), intent(in) :: stencil
        integer, intent(in) :: communicator
        integer, intent(out) :: status
        integer, allocatable :: coefficient_counts(:), coefficient_displacements(:)
        integer, allocatable :: local_ids(:), packed_ids(:), row_seen(:)
        real(c_double), allocatable :: local_coefficients(:,:), packed_coefficients(:,:)
        integer :: ierr, p, q, i, j, k, di, dk, dj, gx, gy, neighbor
        integer :: local_count, total_count, coefficient_count, matrix_status

        call this%release()
        status = 0
        this%communicator = communicator
        call MPI_Comm_rank(communicator, this%rank, ierr)
        call MPI_Comm_size(communicator, this%n_ranks, ierr)
        this%nx_global = stencil%nx_global
        this%ny_global = stencil%ny_global
        this%nz = stencil%nz
        this%nx_local = stencil%nx
        this%ny_local = stencil%ny

        if (this%nx_global <= 2 .or. this%ny_global <= 2 .or. this%nz <= 2) then
            status = 1
            return
        endif
        this%n_global = (this%nx_global-2)*(this%ny_global-2)*(this%nz-2)

        local_count = 0
        do j = 1, stencil%ny
            gy = stencil%y_first+j-1
            if (gy <= 0 .or. gy >= this%ny_global-1) cycle
            do k = 2, this%nz-1
                do i = 1, stencil%nx
                    gx = stencil%x_first+i-1
                    if (gx <= 0 .or. gx >= this%nx_global-1) cycle
                    local_count = local_count+1
                enddo
            enddo
        enddo
        this%n_local = local_count
        allocate(this%counts(this%n_ranks), this%displacements(this%n_ranks))
        call MPI_Allgather(local_count, 1, MPI_INTEGER, this%counts, 1, MPI_INTEGER, communicator, ierr)
        this%displacements(1) = 0
        do p = 2, this%n_ranks
            this%displacements(p) = this%displacements(p-1)+this%counts(p-1)
        enddo
        total_count = sum(this%counts)
        if (total_count /= this%n_global) status = 1
        call MPI_Allreduce(MPI_IN_PLACE, status, 1, MPI_INTEGER, MPI_MAX, communicator, ierr)
        if (status /= 0) then
            call this%release()
            return
        endif

        allocate(this%local_i(max(1,local_count)), this%local_k(max(1,local_count)), &
                 this%local_j(max(1,local_count)))
        allocate(local_ids(max(1,local_count)), local_coefficients(STENCIL_WIDTH,max(1,local_count)))
        local_ids = 0
        local_coefficients = 0.0_c_double
        p = 0
        do j = 1, stencil%ny
            gy = stencil%y_first+j-1
            if (gy <= 0 .or. gy >= this%ny_global-1) cycle
            do k = 2, this%nz-1
                do i = 1, stencil%nx
                    gx = stencil%x_first+i-1
                    if (gx <= 0 .or. gx >= this%nx_global-1) cycle
                    p = p+1
                    this%local_i(p) = i
                    this%local_k(p) = k
                    this%local_j(p) = j
                    local_ids(p) = interior_row(this, gx, k, gy)
                    q = 0
                    do dj = -1, 1
                        do dk = -1, 1
                            do di = -1, 1
                                q = q+1
                                local_coefficients(q,p) = stencil%value(di,dk,dj,i,k,j)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        allocate(this%packed_row_ids(max(1,total_count)), packed_ids(max(1,total_count)))
        this%packed_row_ids = 0
        packed_ids = 0
        call MPI_Gatherv(local_ids, local_count, MPI_INTEGER, packed_ids, this%counts, &
                         this%displacements, MPI_INTEGER, this%root, communicator, ierr)

        allocate(coefficient_counts(this%n_ranks), coefficient_displacements(this%n_ranks))
        coefficient_counts = STENCIL_WIDTH*this%counts
        coefficient_displacements = STENCIL_WIDTH*this%displacements
        coefficient_count = STENCIL_WIDTH*local_count
        allocate(packed_coefficients(STENCIL_WIDTH,max(1,total_count)))
        packed_coefficients = 0.0_c_double
        call MPI_Gatherv(local_coefficients, coefficient_count, MPI_DOUBLE_PRECISION, &
                         packed_coefficients, coefficient_counts, coefficient_displacements, &
                         MPI_DOUBLE_PRECISION, this%root, communicator, ierr)

        matrix_status = 0
        if (this%rank == this%root) then
            this%packed_row_ids(1:total_count) = packed_ids(1:total_count)
            allocate(this%matrix_nnz(this%n_global), &
                     this%matrix_columns(STENCIL_WIDTH,this%n_global), &
                     this%matrix_values(STENCIL_WIDTH,this%n_global), row_seen(this%n_global))
            this%matrix_nnz = 0
            this%matrix_columns = 0
            this%matrix_values = 0.0_c_double
            row_seen = 0
            do p = 1, total_count
                i = packed_ids(p)
                if (i < 1 .or. i > this%n_global) then
                    matrix_status = 1
                    cycle
                endif
                row_seen(i) = row_seen(i)+1
                call decode_interior_row(this, i, gx, k, gy)
                q = 0
                do dj = -1, 1
                    do dk = -1, 1
                        do di = -1, 1
                            q = q+1
                            if (gx+di <= 0 .or. gx+di >= this%nx_global-1 .or. &
                                gy+dj <= 0 .or. gy+dj >= this%ny_global-1 .or. &
                                k+dk <= 1 .or. k+dk >= this%nz) cycle
                            if (packed_coefficients(q,p) == 0.0_c_double) cycle
                            neighbor = interior_row(this, gx+di, k+dk, gy+dj)
                            this%matrix_nnz(i) = this%matrix_nnz(i)+1
                            this%matrix_columns(this%matrix_nnz(i),i) = neighbor
                            this%matrix_values(this%matrix_nnz(i),i) = packed_coefficients(q,p)
                        enddo
                    enddo
                enddo
            enddo
            if (any(row_seen /= 1) .or. any(this%matrix_nnz < 1)) matrix_status = 1
            if (matrix_status == 0) call factorize_root_vertical_lines(this, matrix_status)
            deallocate(row_seen)
        endif
        call MPI_Bcast(matrix_status, 1, MPI_INTEGER, this%root, communicator, ierr)
        if (matrix_status /= 0) then
            status = 1
            deallocate(local_ids, packed_ids, local_coefficients, packed_coefficients, &
                       coefficient_counts, coefficient_displacements)
            call this%release()
            return
        endif

        allocate(this%packed_rhs(max(1,total_count)), this%packed_solution(max(1,total_count)))
        this%packed_rhs = 0.0_c_double
        this%packed_solution = 0.0_c_double
        if (this%rank == this%root) then
            allocate(this%global_rhs(this%n_global), this%global_solution(this%n_global), &
                     this%arnoldi_v(this%n_global,this%n_global+1), &
                     this%arnoldi_z(this%n_global,this%n_global), &
                     this%hessenberg(this%n_global+1,this%n_global), &
                     this%givens_c(this%n_global), this%givens_s(this%n_global), &
                     this%least_squares_rhs(this%n_global+1), this%coefficients(this%n_global))
            this%global_rhs = 0.0_c_double
            this%global_solution = 0.0_c_double
        endif
        this%ready = .true.

        deallocate(local_ids, packed_ids, local_coefficients, packed_coefficients, &
                   coefficient_counts, coefficient_displacements)
    end subroutine setup_collective_coarse_solver


    subroutine solve_collective_coarse_system(this, local_rhs, local_solution, status, iterations, relative_residual)
        class(collective_coarse_solver_t), intent(inout) :: this
        real(c_double), intent(in) :: local_rhs(:,:,:)
        real(c_double), intent(out) :: local_solution(:,:,:)
        integer, intent(out) :: status, iterations
        real(c_double), intent(out) :: relative_residual
        real(c_double), allocatable :: local_values(:)
        integer :: ierr, p

        status = 1
        iterations = 0
        relative_residual = huge(1.0_c_double)
        local_solution = 0.0_c_double
        if (.not. this%ready) return
        if (size(local_rhs,1) /= this%nx_local .or. size(local_rhs,2) /= this%nz .or. &
            size(local_rhs,3) /= this%ny_local .or. any(shape(local_solution) /= shape(local_rhs))) return

        allocate(local_values(max(1,this%n_local)))
        local_values = 0.0_c_double
        do p = 1, this%n_local
            local_values(p) = local_rhs(this%local_i(p),this%local_k(p),this%local_j(p))
        enddo
        this%packed_rhs = 0.0_c_double
        call MPI_Gatherv(local_values, this%n_local, MPI_DOUBLE_PRECISION, this%packed_rhs, &
                         this%counts, this%displacements, MPI_DOUBLE_PRECISION, &
                         this%root, this%communicator, ierr)

        if (this%rank == this%root) then
            this%global_rhs = 0.0_c_double
            do p = 1, this%n_global
                this%global_rhs(this%packed_row_ids(p)) = this%packed_rhs(p)
            enddo
            call root_refined_gmres(this, status, iterations, relative_residual)
            if (this%solve_count == 1 .or. status /= 0) then
                write(output_unit,'(A,I0,A,ES12.4,A,I0)') &
                    ' HICAR terminal physical solve: iterations=', iterations, &
                    ' relative_residual=', relative_residual, ' status=', status
                flush(output_unit)
            endif
            do p = 1, this%n_global
                this%packed_solution(p) = this%global_solution(this%packed_row_ids(p))
            enddo
        endif
        call MPI_Bcast(status, 1, MPI_INTEGER, this%root, this%communicator, ierr)
        call MPI_Bcast(iterations, 1, MPI_INTEGER, this%root, this%communicator, ierr)
        call MPI_Bcast(relative_residual, 1, MPI_DOUBLE_PRECISION, this%root, this%communicator, ierr)
        call MPI_Scatterv(this%packed_solution, this%counts, this%displacements, MPI_DOUBLE_PRECISION, &
                          local_values, this%n_local, MPI_DOUBLE_PRECISION, this%root, this%communicator, ierr)
        do p = 1, this%n_local
            local_solution(this%local_i(p),this%local_k(p),this%local_j(p)) = local_values(p)
        enddo
        this%solve_count = this%solve_count+1
        deallocate(local_values)
    end subroutine solve_collective_coarse_system


    subroutine root_refined_gmres(this, status, iterations, relative_residual)
        class(collective_coarse_solver_t), intent(inout) :: this
        integer, intent(out) :: status, iterations
        real(c_double), intent(out) :: relative_residual
        real(c_double), allocatable :: original_rhs(:), accumulated_solution(:)
        real(c_double), allocatable :: residual(:), matvec_input(:), matvec_output(:)
        real(c_double) :: original_norm, previous_residual
        integer :: refinement, stage_status, stage_iterations
        real(c_double) :: stage_residual

        allocate(original_rhs(this%n_global), accumulated_solution(this%n_global), &
                 residual(this%n_global), matvec_input(this%n_global), matvec_output(this%n_global))
        original_rhs = this%global_rhs
        accumulated_solution = 0.0_c_double
        residual = original_rhs
        original_norm = sqrt(max(dot_product(original_rhs,original_rhs),0.0_c_double))
        iterations = 0
        relative_residual = 0.0_c_double
        status = 0
        if (original_norm == 0.0_c_double) then
            this%global_solution = 0.0_c_double
            return
        endif

        status = 1
        previous_residual = huge(1.0_c_double)
        do refinement = 1, 3
            this%global_rhs = residual
            call root_full_gmres(this, stage_status, stage_iterations, stage_residual)
            iterations = iterations+stage_iterations
            if (.not. ieee_is_finite(stage_residual)) exit
            if (stage_status /= 0 .and. stage_iterations == 0) exit
            accumulated_solution = accumulated_solution+this%global_solution
            matvec_input = accumulated_solution
            call root_sparse_matvec(this, matvec_input, matvec_output)
            residual = original_rhs-matvec_output
            relative_residual = sqrt(max(dot_product(residual,residual),0.0_c_double))/original_norm
            if (.not. ieee_is_finite(relative_residual)) exit
            if (relative_residual <= COARSE_ACCEPTANCE_TOLERANCE) then
                status = 0
                exit
            endif
            if (relative_residual >= 0.95_c_double*previous_residual) exit
            previous_residual = relative_residual
        enddo
        this%global_rhs = original_rhs
        this%global_solution = accumulated_solution
    end subroutine root_refined_gmres


    subroutine root_full_gmres(this, status, iterations, relative_residual)
        class(collective_coarse_solver_t), intent(inout) :: this
        integer, intent(out) :: status, iterations
        real(c_double), intent(out) :: relative_residual
        real(c_double), allocatable :: preconditioned(:), work_vector(:)
        real(c_double) :: beta, rhs_norm, target, projection, denominator, temporary
        integer :: i, j, pass, used

        status = 1
        iterations = 0
        relative_residual = huge(1.0_c_double)
        this%global_solution = 0.0_c_double
        rhs_norm = sqrt(max(dot_product(this%global_rhs,this%global_rhs),0.0_c_double))
        if (rhs_norm == 0.0_c_double) then
            status = 0
            relative_residual = 0.0_c_double
            return
        endif
        target = max(COARSE_ABSOLUTE_TOLERANCE, COARSE_RELATIVE_TOLERANCE*rhs_norm)
        beta = rhs_norm
        this%arnoldi_v = 0.0_c_double
        this%arnoldi_z = 0.0_c_double
        this%hessenberg = 0.0_c_double
        this%givens_c = 0.0_c_double
        this%givens_s = 0.0_c_double
        this%least_squares_rhs = 0.0_c_double
        this%coefficients = 0.0_c_double
        this%arnoldi_v(:,1) = this%global_rhs/beta
        this%least_squares_rhs(1) = beta
        used = this%n_global
        allocate(preconditioned(this%n_global), work_vector(this%n_global))

        do j = 1, this%n_global
            preconditioned = this%arnoldi_v(:,j)
            call apply_root_vertical_lines(this, preconditioned, work_vector)
            this%arnoldi_z(:,j) = work_vector
            preconditioned = work_vector
            call root_sparse_matvec(this, preconditioned, work_vector)
            do pass = 1, 2
                do i = 1, j
                    projection = dot_product(this%arnoldi_v(:,i),work_vector)
                    this%hessenberg(i,j) = this%hessenberg(i,j)+projection
                    work_vector = work_vector-projection*this%arnoldi_v(:,i)
                enddo
            enddo
            this%hessenberg(j+1,j) = sqrt(max(dot_product(work_vector,work_vector),0.0_c_double))
            if (this%hessenberg(j+1,j) > COARSE_BREAKDOWN_TOLERANCE .and. j < this%n_global) &
                this%arnoldi_v(:,j+1) = work_vector/this%hessenberg(j+1,j)

            do i = 1, j-1
                temporary = this%givens_c(i)*this%hessenberg(i,j) + &
                            this%givens_s(i)*this%hessenberg(i+1,j)
                this%hessenberg(i+1,j) = -this%givens_s(i)*this%hessenberg(i,j) + &
                                         this%givens_c(i)*this%hessenberg(i+1,j)
                this%hessenberg(i,j) = temporary
            enddo
            denominator = sqrt(this%hessenberg(j,j)**2+this%hessenberg(j+1,j)**2)
            if (denominator <= COARSE_BREAKDOWN_TOLERANCE) then
                used = j
                exit
            endif
            this%givens_c(j) = this%hessenberg(j,j)/denominator
            this%givens_s(j) = this%hessenberg(j+1,j)/denominator
            this%hessenberg(j,j) = denominator
            this%hessenberg(j+1,j) = 0.0_c_double
            temporary = this%givens_c(j)*this%least_squares_rhs(j) + &
                        this%givens_s(j)*this%least_squares_rhs(j+1)
            this%least_squares_rhs(j+1) = -this%givens_s(j)*this%least_squares_rhs(j) + &
                                          this%givens_c(j)*this%least_squares_rhs(j+1)
            this%least_squares_rhs(j) = temporary
            iterations = j
            if (abs(this%least_squares_rhs(j+1)) <= target) then
                used = j
                exit
            endif
        enddo

        do i = used, 1, -1
            this%coefficients(i) = this%least_squares_rhs(i)
            do j = i+1, used
                this%coefficients(i) = this%coefficients(i)- &
                                       this%hessenberg(i,j)*this%coefficients(j)
            enddo
            if (abs(this%hessenberg(i,i)) <= COARSE_BREAKDOWN_TOLERANCE) return
            this%coefficients(i) = this%coefficients(i)/this%hessenberg(i,i)
        enddo
        do i = 1, used
            this%global_solution = this%global_solution + &
                                   this%coefficients(i)*this%arnoldi_z(:,i)
        enddo
        preconditioned = this%global_solution
        call root_sparse_matvec(this, preconditioned, work_vector)
        work_vector = this%global_rhs-work_vector
        relative_residual = sqrt(max(dot_product(work_vector,work_vector),0.0_c_double))/rhs_norm
        iterations = used
        if (relative_residual <= COARSE_ACCEPTANCE_TOLERANCE .or. &
            sqrt(max(dot_product(work_vector,work_vector),0.0_c_double)) <= COARSE_ABSOLUTE_TOLERANCE) status = 0
    end subroutine root_full_gmres


    subroutine root_sparse_matvec(this, x, ax)
        class(collective_coarse_solver_t), intent(in) :: this
        real(c_double), intent(in) :: x(:)
        real(c_double), intent(out) :: ax(:)
        integer :: row, q

        ax = 0.0_c_double
        do row = 1, this%n_global
            do q = 1, this%matrix_nnz(row)
                ax(row) = ax(row)+this%matrix_values(q,row)*x(this%matrix_columns(q,row))
            enddo
        enddo
    end subroutine root_sparse_matvec


    subroutine factorize_root_vertical_lines(this, status)
        class(collective_coarse_solver_t), intent(inout) :: this
        integer, intent(out) :: status
        real(c_double) :: diagonal, lower, upper, pivot
        integer :: gx, gy, k, row

        status = 0
        allocate(this%diagonal_inverse(this%n_global), this%upper_prime(this%n_global), &
                 this%lower_coefficient(this%n_global))
        this%diagonal_inverse = 0.0_c_double
        this%upper_prime = 0.0_c_double
        this%lower_coefficient = 0.0_c_double
        do gy = 1, this%ny_global-2
            do gx = 1, this%nx_global-2
                do k = 2, this%nz-1
                    row = interior_row(this,gx,k,gy)
                    diagonal = matrix_entry(this,row,row)
                    lower = 0.0_c_double
                    upper = 0.0_c_double
                    if (k > 2) lower = matrix_entry(this,row,row-1)
                    if (k < this%nz-1) upper = matrix_entry(this,row,row+1)
                    if (k == 2) then
                        pivot = diagonal
                    else
                        pivot = diagonal-lower*this%upper_prime(row-1)
                    endif
                    if (abs(pivot) <= 1.0e-28_c_double) then
                        status = 1
                        return
                    endif
                    this%diagonal_inverse(row) = 1.0_c_double/pivot
                    this%upper_prime(row) = upper*this%diagonal_inverse(row)
                    this%lower_coefficient(row) = lower
                enddo
            enddo
        enddo
    end subroutine factorize_root_vertical_lines


    subroutine apply_root_vertical_lines(this, rhs, solution)
        class(collective_coarse_solver_t), intent(in) :: this
        real(c_double), intent(in) :: rhs(:)
        real(c_double), intent(out) :: solution(:)
        integer :: gx, gy, k, row

        solution = 0.0_c_double
        do gy = 1, this%ny_global-2
            do gx = 1, this%nx_global-2
                do k = 2, this%nz-1
                    row = interior_row(this,gx,k,gy)
                    if (k == 2) then
                        solution(row) = rhs(row)*this%diagonal_inverse(row)
                    else
                        solution(row) = (rhs(row)-this%lower_coefficient(row)*solution(row-1))* &
                                        this%diagonal_inverse(row)
                    endif
                enddo
                do k = this%nz-2, 2, -1
                    row = interior_row(this,gx,k,gy)
                    solution(row) = solution(row)-this%upper_prime(row)*solution(row+1)
                enddo
            enddo
        enddo
    end subroutine apply_root_vertical_lines


    real(c_double) function matrix_entry(this, row, column) result(value)
        class(collective_coarse_solver_t), intent(in) :: this
        integer, intent(in) :: row, column
        integer :: q

        value = 0.0_c_double
        do q = 1, this%matrix_nnz(row)
            if (this%matrix_columns(q,row) == column) then
                value = this%matrix_values(q,row)
                return
            endif
        enddo
    end function matrix_entry


    pure integer function interior_row(this, gx, k, gy) result(row)
        class(collective_coarse_solver_t), intent(in) :: this
        integer, intent(in) :: gx, k, gy

        row = ((gy-1)*(this%nx_global-2)+(gx-1))*(this%nz-2)+(k-2)+1
    end function interior_row


    pure subroutine decode_interior_row(this, row, gx, k, gy)
        class(collective_coarse_solver_t), intent(in) :: this
        integer, intent(in) :: row
        integer, intent(out) :: gx, k, gy
        integer :: horizontal_index

        k = modulo(row-1,this%nz-2)+2
        horizontal_index = (row-1)/(this%nz-2)
        gx = modulo(horizontal_index,this%nx_global-2)+1
        gy = horizontal_index/(this%nx_global-2)+1
    end subroutine decode_interior_row


    subroutine release_collective_coarse_solver(this)
        class(collective_coarse_solver_t), intent(inout) :: this

        if (allocated(this%counts)) deallocate(this%counts, this%displacements)
        if (allocated(this%local_i)) deallocate(this%local_i, this%local_k, this%local_j)
        if (allocated(this%packed_row_ids)) deallocate(this%packed_row_ids)
        if (allocated(this%matrix_nnz)) deallocate(this%matrix_nnz, this%matrix_columns, this%matrix_values)
        if (allocated(this%diagonal_inverse)) deallocate(this%diagonal_inverse, this%upper_prime, &
                                                         this%lower_coefficient)
        if (allocated(this%packed_rhs)) deallocate(this%packed_rhs, this%packed_solution)
        if (allocated(this%global_rhs)) deallocate(this%global_rhs, this%global_solution)
        if (allocated(this%arnoldi_v)) deallocate(this%arnoldi_v, this%arnoldi_z, this%hessenberg, &
                                                  this%givens_c, this%givens_s, &
                                                  this%least_squares_rhs, this%coefficients)
        this%ready = .false.
        this%communicator = MPI_COMM_NULL
        this%rank = -1
        this%n_ranks = 0
        this%nx_global = 0
        this%ny_global = 0
        this%nz = 0
        this%nx_local = 0
        this%ny_local = 0
        this%n_global = 0
        this%n_local = 0
        this%solve_count = 0
    end subroutine release_collective_coarse_solver

end module wind_coarse_solve
