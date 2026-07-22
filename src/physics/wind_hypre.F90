!> Exact calibrated wind operator backend using HYPRE ParCSR/BoomerAMG.
!! The adapter deliberately stages HICAR's OpenACC fields through host buffers
!! for correctness first; wind_iterative retains ownership of the mandatory
!! independent true-residual check before a correction can be applied.
module wind_hypre
    use iso_c_binding
    implicit none
    private
    public :: wind_hypre_available, wind_hypre_invalidate, wind_hypre_solve, wind_hypre_apply

    logical :: matrix_ready = .false.

#ifdef USE_HYPRE
    interface
        integer(c_int) function hicar_hypre_build(comm, xs, ys, zs, xm, ym, zm, mx, my, mz, &
                ci_s, ci_e, ck_s, ck_e, cj_s, cj_e, &
                a, b, c, d, e, f, g, h, ii, jj, kk, l, m, n, o) bind(C)
            import :: c_int, c_float
            integer(c_int), value :: comm, xs, ys, zs, xm, ym, zm, mx, my, mz
            integer(c_int), value :: ci_s, ci_e, ck_s, ck_e, cj_s, cj_e
            real(c_float), intent(in) :: a(*), b(*), c(*), d(*), e(*), f(*), g(*), h(*)
            real(c_float), intent(in) :: ii(*), jj(*), kk(*), l(*), m(*), n(*), o(*)
        end function
        integer(c_int) function hicar_hypre_initialize() bind(C)
            import :: c_int
        end function
        integer(c_int) function hicar_hypre_solve(rhs, x, max_iter, tol, iterations, residual) bind(C)
            import :: c_int, c_double
            real(c_double), intent(in) :: rhs(*)
            real(c_double), intent(inout) :: x(*)
            integer(c_int), value :: max_iter
            real(c_double), value :: tol
            integer(c_int), intent(out) :: iterations
            real(c_double), intent(out) :: residual
        end function
        subroutine hicar_hypre_destroy() bind(C)
        end subroutine
        integer(c_int) function hicar_hypre_apply(x, y) bind(C)
            import :: c_int, c_double
            real(c_double), intent(in) :: x(*)
            real(c_double), intent(out) :: y(*)
        end function
    end interface
#endif

contains

    logical function wind_hypre_available()
#ifdef USE_HYPRE
        wind_hypre_available = .true.
#else
        wind_hypre_available = .false.
#endif
    end function wind_hypre_available

    subroutine wind_hypre_invalidate()
        matrix_ready = .false.
#ifdef USE_HYPRE
        call hicar_hypre_destroy()
#endif
    end subroutine wind_hypre_invalidate

    subroutine wind_hypre_solve(comm, xs, ys, zs, xm, ym, zm, mx, my, mz, i_s, i_e, k_s, k_e, j_s, j_e, &
                                a, b, c, d, e, f, g, h, ii, jj, kk, l, m, n, o, rhs, x, max_iter, tol, &
                                status, iterations, residual)
        integer, intent(in) :: comm, xs, ys, zs, xm, ym, zm, mx, my, mz, i_s, i_e, k_s, k_e, j_s, j_e, max_iter
        real(c_float), intent(in) :: a(i_s:i_e,k_s:k_e,j_s:j_e), b(i_s:i_e,k_s:k_e,j_s:j_e)
        real(c_float), intent(in) :: c(i_s:i_e,k_s:k_e,j_s:j_e), d(i_s:i_e,k_s:k_e,j_s:j_e)
        real(c_float), intent(in) :: e(i_s:i_e,k_s:k_e,j_s:j_e), f(i_s:i_e,k_s:k_e,j_s:j_e)
        real(c_float), intent(in) :: g(i_s:i_e,k_s:k_e,j_s:j_e), h(i_s:i_e,k_s:k_e,j_s:j_e)
        real(c_float), intent(in) :: ii(i_s:i_e,k_s:k_e,j_s:j_e), jj(i_s:i_e,k_s:k_e,j_s:j_e)
        real(c_float), intent(in) :: kk(i_s:i_e,k_s:k_e,j_s:j_e), l(i_s:i_e,k_s:k_e,j_s:j_e)
        real(c_float), intent(in) :: m(i_s:i_e,k_s:k_e,j_s:j_e), n(i_s:i_e,k_s:k_e,j_s:j_e), o(i_s:i_e,k_s:k_e,j_s:j_e)
        real(c_double), intent(in) :: rhs(i_s-1:i_e+1,k_s-1:k_e+1,j_s-1:j_e+1)
        real(c_double), intent(inout) :: x(i_s-1:i_e+1,k_s-1:k_e+1,j_s-1:j_e+1)
        real(c_double), intent(in) :: tol
        integer, intent(out) :: status, iterations
        real(c_double), intent(out) :: residual
        real(c_double), allocatable :: rhs_pack(:), x_pack(:)
        integer :: i, j, k, q
#ifdef USE_HYPRE
        integer(c_int) :: rc, its

        status = 1; iterations = 0; residual = huge(1.0_c_double)
        rc = hicar_hypre_initialize()
        if (rc /= 0_c_int) then
            status = int(rc)
            return
        endif
        if (.not. matrix_ready) then
            rc = hicar_hypre_build(int(comm,c_int), int(xs,c_int), int(ys,c_int), int(zs,c_int), &
                 int(xm,c_int), int(ym,c_int), int(zm,c_int), int(mx,c_int), int(my,c_int), int(mz,c_int), &
                 int(i_s,c_int), int(i_e,c_int), int(k_s,c_int), int(k_e,c_int), int(j_s,c_int), int(j_e,c_int), &
                 a, b, c, d, e, f, g, h, ii, jj, kk, l, m, n, o)
            if (rc /= 0_c_int) then
                status = int(rc)
                return
            endif
            matrix_ready = .true.
        endif

        allocate(rhs_pack(xm*ym*zm), x_pack(xm*ym*zm))
        q = 0
        do j = ys, ys+ym-1
            do k = zs, zs+zm-1
                do i = xs, xs+xm-1
                    q = q + 1
                    rhs_pack(q) = rhs(i,k,j)
                    x_pack(q) = x(i,k,j)
                enddo
            enddo
        enddo
        rc = hicar_hypre_solve(rhs_pack, x_pack, int(max_iter,c_int), tol, its, residual)
        iterations = int(its)
        if (rc == 0_c_int) then
            q = 0
            do j = ys, ys+ym-1
                do k = zs, zs+zm-1
                    do i = xs, xs+xm-1
                        q = q + 1
                        x(i,k,j) = x_pack(q)
                    enddo
                enddo
            enddo
            status = 0
        else
            status = int(rc)
        endif
        deallocate(rhs_pack, x_pack)
#else
        status = -999; iterations = 0; residual = huge(1.0_c_double)
#endif
    end subroutine wind_hypre_solve

    subroutine wind_hypre_apply(xs, ys, zs, xm, ym, zm, i_s, i_e, k_s, k_e, j_s, j_e, x, y, status)
        integer, intent(in) :: xs, ys, zs, xm, ym, zm, i_s, i_e, k_s, k_e, j_s, j_e
        real(c_double), intent(in) :: x(i_s-1:i_e+1,k_s-1:k_e+1,j_s-1:j_e+1)
        real(c_double), intent(out) :: y(i_s-1:i_e+1,k_s-1:k_e+1,j_s-1:j_e+1)
        integer, intent(out) :: status
        real(c_double), allocatable :: x_pack(:), y_pack(:)
        integer :: i, j, k, q
#ifdef USE_HYPRE
        integer(c_int) :: rc
        status = -10
        if (.not. matrix_ready) return
        allocate(x_pack(xm*ym*zm), y_pack(xm*ym*zm))
        q = 0
        do j = ys, ys+ym-1
            do k = zs, zs+zm-1
                do i = xs, xs+xm-1
                    q = q + 1
                    x_pack(q) = x(i,k,j)
                enddo
            enddo
        enddo
        rc = hicar_hypre_apply(x_pack, y_pack)
        if (rc == 0_c_int) then
            q = 0
            do j = ys, ys+ym-1
                do k = zs, zs+zm-1
                    do i = xs, xs+xm-1
                        q = q + 1
                        y(i,k,j) = y_pack(q)
                    enddo
                enddo
            enddo
            status = 0
        else
            status = int(rc)
        endif
        deallocate(x_pack, y_pack)
#else
        status = -999
#endif
    end subroutine wind_hypre_apply

end module wind_hypre
