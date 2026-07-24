!> Discretely adjoint mass-consistent wind projection primitives.
!!
!! The constraint is the finite-volume mass-flux balance
!!
!!     B q = 0,
!!
!! where B is the cell-volume-weighted divergence and q contains the
!! staggered velocity corrections.  Given a positive face energy M, the
!! correction and Schur complement are defined algebraically by
!!
!!     delta q = -M^{-1} B^T lambda,
!!     K lambda = B M^{-1} B^T lambda.
!!
!! This makes K symmetric positive definite (subject to the selected boundary
!! conditions) by construction.  It also guarantees that a converged solve
!! removes the *same* discrete mass imbalance measured by B.
!!
!! The initial diagonal metric implemented here is the exact flat-grid limit
!! of HICAR's existing correction.  Its face energies penalize squared mass
!! flux.  The factors include physical dual volumes, density, map factors,
!! the vertical Jacobian, and alpha.  Terrain cross-coupling belongs in an
!! explicit positive metric extension; it must not be introduced through a
!! separately discretized, non-adjoint gradient.
module wind_adjoint_projection
    use, intrinsic :: iso_c_binding, only : c_double
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    public :: adjoint_projection_t

    type :: adjoint_projection_t
        integer :: nx = 0
        integer :: ny = 0
        integer :: nz = 0
        real(c_double), allocatable :: cell_volume(:,:,:)
        real(c_double), allocatable :: b_u(:,:,:), inv_m_u(:,:,:)
        real(c_double), allocatable :: b_v(:,:,:), inv_m_v(:,:,:)
        real(c_double), allocatable :: b_w(:,:,:), inv_m_w(:,:,:)
    contains
        procedure :: initialize_diagonal_metric
        procedure :: apply_constraint
        procedure :: apply_correction
        procedure :: apply_schur
        procedure :: release
    end type adjoint_projection_t

contains

    !> Construct B and the diagonal inverse face energy from HICAR geometry.
    !!
    !! Array ordering follows HICAR: x, vertical, y.  The w arrays contain
    !! only the nz-1 correctable interior interfaces; the ground flux is
    !! identically zero and the top interface is the rigid lid.
    subroutine initialize_diagonal_metric(this, dx, dz, dz_w, jaco_c, &
                                          jaco_u, rho_u, map_mx_u, map_my_u, &
                                          jaco_v, rho_v, map_mx_v, map_my_v, &
                                          jaco_w, rho_w, alpha_w, map_mxy, status)
        class(adjoint_projection_t), intent(inout) :: this
        real(c_double), intent(in) :: dx
        real(c_double), intent(in) :: dz(:), dz_w(:)
        real(c_double), intent(in) :: jaco_c(:,:,:)
        real(c_double), intent(in) :: jaco_u(:,:,:), rho_u(:,:,:)
        real(c_double), intent(in) :: map_mx_u(:,:), map_my_u(:,:)
        real(c_double), intent(in) :: jaco_v(:,:,:), rho_v(:,:,:)
        real(c_double), intent(in) :: map_mx_v(:,:), map_my_v(:,:)
        real(c_double), intent(in) :: jaco_w(:,:,:), rho_w(:,:,:), alpha_w(:,:,:)
        real(c_double), intent(in) :: map_mxy(:,:)
        integer, intent(out) :: status
        integer :: i, j, k, nx, ny, nz

        call this%release()
        status = 1
        nx = size(jaco_c,1)
        nz = size(jaco_c,2)
        ny = size(jaco_c,3)
        if (nx < 1 .or. ny < 1 .or. nz < 2) return
        if (dx <= 0.0_c_double .or. .not. ieee_is_finite(dx)) return
        if (size(dz) /= nz .or. size(dz_w) /= nz-1) return
        if (any(shape(jaco_u) /= [nx+1,nz,ny]) .or. any(shape(rho_u) /= [nx+1,nz,ny])) return
        if (any(shape(map_mx_u) /= [nx+1,ny]) .or. any(shape(map_my_u) /= [nx+1,ny])) return
        if (any(shape(jaco_v) /= [nx,nz,ny+1]) .or. any(shape(rho_v) /= [nx,nz,ny+1])) return
        if (any(shape(map_mx_v) /= [nx,ny+1]) .or. any(shape(map_my_v) /= [nx,ny+1])) return
        if (any(shape(jaco_w) /= [nx,nz-1,ny]) .or. any(shape(rho_w) /= [nx,nz-1,ny])) return
        if (any(shape(alpha_w) /= [nx,nz-1,ny]) .or. any(shape(map_mxy) /= [nx,ny])) return
        if (any(dz <= 0.0_c_double) .or. any(dz_w <= 0.0_c_double)) return
        if (any(.not. ieee_is_finite(dz)) .or. any(.not. ieee_is_finite(dz_w))) return
        if (any(.not. ieee_is_finite(jaco_c)) .or. any(.not. ieee_is_finite(jaco_u)) .or. &
            any(.not. ieee_is_finite(jaco_v)) .or. any(.not. ieee_is_finite(jaco_w))) return
        if (any(.not. ieee_is_finite(rho_u)) .or. any(.not. ieee_is_finite(rho_v)) .or. &
            any(.not. ieee_is_finite(rho_w)) .or. any(.not. ieee_is_finite(alpha_w))) return
        if (any(.not. ieee_is_finite(map_mx_u)) .or. any(.not. ieee_is_finite(map_my_u)) .or. &
            any(.not. ieee_is_finite(map_mx_v)) .or. any(.not. ieee_is_finite(map_my_v)) .or. &
            any(.not. ieee_is_finite(map_mxy))) return
        if (any(jaco_c <= 0.0_c_double) .or. any(jaco_u <= 0.0_c_double) .or. &
            any(jaco_v <= 0.0_c_double) .or. any(jaco_w <= 0.0_c_double)) return
        if (any(rho_u <= 0.0_c_double) .or. any(rho_v <= 0.0_c_double) .or. &
            any(rho_w <= 0.0_c_double) .or. any(alpha_w <= 0.0_c_double)) return
        if (any(map_mx_u <= 0.0_c_double) .or. any(map_my_u <= 0.0_c_double) .or. &
            any(map_mx_v <= 0.0_c_double) .or. any(map_my_v <= 0.0_c_double) .or. &
            any(map_mxy <= 0.0_c_double)) return

        this%nx = nx
        this%ny = ny
        this%nz = nz
        allocate(this%cell_volume(nx,nz,ny))
        allocate(this%b_u(nx+1,nz,ny), this%inv_m_u(nx+1,nz,ny))
        allocate(this%b_v(nx,nz,ny+1), this%inv_m_v(nx,nz,ny+1))
        allocate(this%b_w(nx,nz-1,ny), this%inv_m_w(nx,nz-1,ny))

        do j = 1, ny
            do k = 1, nz
                do i = 1, nx
                    this%cell_volume(i,k,j) = dx*dx*dz(k)*jaco_c(i,k,j)/map_mxy(i,j)
                enddo
            enddo
        enddo

        do j = 1, ny
            do k = 1, nz
                do i = 1, nx+1
                    ! B face coefficient = physical face area times mass density.
                    this%b_u(i,k,j) = dx*dz(k)*jaco_u(i,k,j)*rho_u(i,k,j)/map_my_u(i,j)
                    ! M_u = 2 rho^2 V_u.  B/M reduces exactly to
                    ! map_mx_u/(2 rho dx) on the diagonal metric.
                    this%inv_m_u(i,k,j) = map_mx_u(i,j)*map_my_u(i,j) / &
                        (2.0_c_double*rho_u(i,k,j)**2*dx*dx*dz(k)*jaco_u(i,k,j))
                enddo
            enddo
        enddo

        do j = 1, ny+1
            do k = 1, nz
                do i = 1, nx
                    this%b_v(i,k,j) = dx*dz(k)*jaco_v(i,k,j)*rho_v(i,k,j)/map_mx_v(i,j)
                    this%inv_m_v(i,k,j) = map_mx_v(i,j)*map_my_v(i,j) / &
                        (2.0_c_double*rho_v(i,k,j)**2*dx*dx*dz(k)*jaco_v(i,k,j))
                enddo
            enddo
        enddo

        do j = 1, ny
            do k = 1, nz-1
                do i = 1, nx
                    this%b_w(i,k,j) = dx*dx*jaco_w(i,k,j)*rho_w(i,k,j)/map_mxy(i,j)
                    ! Grid-relative w converts to physical vertical velocity
                    ! with a factor jaco_w, hence jaco_w^2 in its energy.
                    this%inv_m_w(i,k,j) = alpha_w(i,k,j)**2*map_mxy(i,j) / &
                        (2.0_c_double*rho_w(i,k,j)**2*dx*dx*dz_w(k)*jaco_w(i,k,j)**2)
                enddo
            enddo
        enddo
        status = 0
    end subroutine initialize_diagonal_metric


    !> Apply the volume-integrated discrete mass constraint B q.
    subroutine apply_constraint(this, u, v, w, constraint)
        class(adjoint_projection_t), intent(in) :: this
        real(c_double), intent(in) :: u(this%nx+1,this%nz,this%ny)
        real(c_double), intent(in) :: v(this%nx,this%nz,this%ny+1)
        real(c_double), intent(in) :: w(this%nx,this%nz-1,this%ny)
        real(c_double), intent(out) :: constraint(this%nx,this%nz,this%ny)
        integer :: i, j, k

        constraint = 0.0_c_double
        do j = 1, this%ny
            do k = 1, this%nz
                do i = 1, this%nx
                    constraint(i,k,j) = &
                        this%b_u(i+1,k,j)*u(i+1,k,j) - this%b_u(i,k,j)*u(i,k,j) + &
                        this%b_v(i,k,j+1)*v(i,k,j+1) - this%b_v(i,k,j)*v(i,k,j)
                    if (k < this%nz) constraint(i,k,j) = constraint(i,k,j) + &
                        this%b_w(i,k,j)*w(i,k,j)
                    if (k > 1) constraint(i,k,j) = constraint(i,k,j) - &
                        this%b_w(i,k-1,j)*w(i,k-1,j)
                enddo
            enddo
        enddo
    end subroutine apply_constraint


    !> Apply -M^{-1} B^T lambda to the staggered correction fields.
    subroutine apply_correction(this, lambda, u, v, w)
        class(adjoint_projection_t), intent(in) :: this
        real(c_double), intent(in) :: lambda(this%nx,this%nz,this%ny)
        real(c_double), intent(out) :: u(this%nx+1,this%nz,this%ny)
        real(c_double), intent(out) :: v(this%nx,this%nz,this%ny+1)
        real(c_double), intent(out) :: w(this%nx,this%nz-1,this%ny)
        integer :: i, j, k

        do j = 1, this%ny
            do k = 1, this%nz
                u(1,k,j) = this%inv_m_u(1,k,j)*this%b_u(1,k,j)*lambda(1,k,j)
                do i = 2, this%nx
                    u(i,k,j) = this%inv_m_u(i,k,j)*this%b_u(i,k,j)* &
                        (lambda(i,k,j)-lambda(i-1,k,j))
                enddo
                u(this%nx+1,k,j) = -this%inv_m_u(this%nx+1,k,j)* &
                    this%b_u(this%nx+1,k,j)*lambda(this%nx,k,j)
            enddo
        enddo

        do j = 1, this%ny+1
            do k = 1, this%nz
                do i = 1, this%nx
                    if (j == 1) then
                        v(i,k,j) = this%inv_m_v(i,k,j)*this%b_v(i,k,j)*lambda(i,k,1)
                    else if (j == this%ny+1) then
                        v(i,k,j) = -this%inv_m_v(i,k,j)*this%b_v(i,k,j)*lambda(i,k,this%ny)
                    else
                        v(i,k,j) = this%inv_m_v(i,k,j)*this%b_v(i,k,j)* &
                            (lambda(i,k,j)-lambda(i,k,j-1))
                    endif
                enddo
            enddo
        enddo

        do j = 1, this%ny
            do k = 1, this%nz-1
                do i = 1, this%nx
                    w(i,k,j) = this%inv_m_w(i,k,j)*this%b_w(i,k,j)* &
                        (lambda(i,k+1,j)-lambda(i,k,j))
                enddo
            enddo
        enddo
    end subroutine apply_correction


    !> Apply K = B M^{-1} B^T without assembling a sparse matrix.
    subroutine apply_schur(this, lambda, result)
        class(adjoint_projection_t), intent(in) :: this
        real(c_double), intent(in) :: lambda(this%nx,this%nz,this%ny)
        real(c_double), intent(out) :: result(this%nx,this%nz,this%ny)
        real(c_double), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)

        allocate(u(this%nx+1,this%nz,this%ny))
        allocate(v(this%nx,this%nz,this%ny+1))
        allocate(w(this%nx,this%nz-1,this%ny))
        call this%apply_correction(lambda,u,v,w)
        call this%apply_constraint(u,v,w,result)
        result = -result
        deallocate(u,v,w)
    end subroutine apply_schur


    subroutine release(this)
        class(adjoint_projection_t), intent(inout) :: this

        if (allocated(this%cell_volume)) deallocate(this%cell_volume)
        if (allocated(this%b_u)) deallocate(this%b_u)
        if (allocated(this%inv_m_u)) deallocate(this%inv_m_u)
        if (allocated(this%b_v)) deallocate(this%b_v)
        if (allocated(this%inv_m_v)) deallocate(this%inv_m_v)
        if (allocated(this%b_w)) deallocate(this%b_w)
        if (allocated(this%inv_m_w)) deallocate(this%inv_m_w)
        this%nx = 0
        this%ny = 0
        this%nz = 0
    end subroutine release

end module wind_adjoint_projection
