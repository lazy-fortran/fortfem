module meridian_test_fields
    !! Analytic physical fields u(R, Z) exp(i n phi), their per-n operators
    !! (right-handed (R, phi, Z) basis) and sparse helpers for the meridian
    !! FEEC tests.
    use fortfem_kinds, only: dp
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_load, &
        meridian_real_to_complex
    use fortfem_sparse_direct, only: sparse_direct_solve_once
    use fortsparse, only: csc_from_triplet, csc_z_t, csc_t, fortsparse_status_t, &
        status_ok
    implicit none
    private

    integer, public :: mode_n = 0
    real(dp), save :: ann_c = 3.0_dp, ann_r1 = 0.4_dp, ann_r2 = 1.0_dp
    real(dp), parameter :: fd_h = 2.0e-3_dp
    real(dp), public :: wscale = 1.0_dp

    public :: scalar_f, grad_f, vec_a, curl_a, vec_b, div_b
    public :: project, solve_z, zero_rows, zadd
    logical, public :: rect_chi = .false.
    public :: ann_c, ann_r1, ann_r2, fun_astar, fun_bstar, fun_f1, fun_bs, &
        fun_dbs, fun_f2, fun_fs, fun_gfs, fun_f0, fun_a0star, fun_b0star, &
        fun_j0star

contains

    subroutine scalar_f(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        v = (0.0_dp, 0.0_dp)
        v(1) = cmplx(sin(1.3_dp*x(1))*cos(0.9_dp*x(2) + 0.2_dp), &
            0.3_dp*x(1)*sin(x(2)), dp)
    end subroutine scalar_f

    subroutine grad_f(x, v)
        !! (d_R f, i n f / R, d_Z f)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        complex(dp) :: f
        real(dp) :: r, z
        r = x(1)
        z = x(2)
        f = cmplx(sin(1.3_dp*r)*cos(0.9_dp*z + 0.2_dp), 0.3_dp*r*sin(z), dp)
        v(1) = cmplx(1.3_dp*cos(1.3_dp*r)*cos(0.9_dp*z + 0.2_dp), &
            0.3_dp*sin(z), dp)
        v(2) = cmplx(0.0_dp, real(mode_n, dp), dp)*f/r
        v(3) = cmplx(-0.9_dp*sin(1.3_dp*r)*sin(0.9_dp*z + 0.2_dp), &
            0.3_dp*r*cos(z), dp)
    end subroutine grad_f

    subroutine vec_a(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        real(dp) :: r, z
        r = x(1)
        z = x(2)
        v(1) = cmplx(sin(r + 0.5_dp*z), 0.2_dp*cos(z), dp)
        v(2) = cmplx(cos(0.7_dp*r - z), 0.1_dp*r, dp)
        v(3) = cmplx(exp(0.3_dp*r)*cos(z), sin(0.4_dp*r*z), dp)
    end subroutine vec_a

    subroutine curl_a(x, v)
        !! curl_n of vec_a in physical components.
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        complex(dp) :: ar, ap, az, dzap, drrap, draz, dzar
        complex(dp), parameter :: ii = (0.0_dp, 1.0_dp)
        real(dp) :: r, z
        r = x(1)
        z = x(2)
        ar = cmplx(sin(r + 0.5_dp*z), 0.2_dp*cos(z), dp)
        ap = cmplx(cos(0.7_dp*r - z), 0.1_dp*r, dp)
        az = cmplx(exp(0.3_dp*r)*cos(z), sin(0.4_dp*r*z), dp)
        dzar = cmplx(0.5_dp*cos(r + 0.5_dp*z), -0.2_dp*sin(z), dp)
        dzap = cmplx(sin(0.7_dp*r - z), 0.0_dp, dp)
        drrap = ap + r*cmplx(-0.7_dp*sin(0.7_dp*r - z), 0.1_dp, dp)
        draz = cmplx(0.3_dp*exp(0.3_dp*r)*cos(z), 0.4_dp*z*cos(0.4_dp*r*z), dp)
        v(1) = ii*real(mode_n, dp)*az/r - dzap
        v(2) = dzar - draz
        v(3) = (drrap - ii*real(mode_n, dp)*ar)/r
    end subroutine curl_a


    subroutine vec_b(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        real(dp) :: r, z
        r = x(1)
        z = x(2)
        v(1) = cmplx(cos(0.8_dp*r + z), 0.1_dp*sin(r), dp)
        v(2) = cmplx(sin(r*z*0.5_dp), 0.3_dp*cos(z), dp)
        v(3) = cmplx(exp(-0.2_dp*z)*sin(r), 0.0_dp, dp)
    end subroutine vec_b

    subroutine div_b(x, v)
        !! div_n of vec_b: (1/R) d_R(R B_R) + i n B_phi / R + d_Z B_Z.
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        complex(dp) :: br, bp, drbr, dzbz
        complex(dp), parameter :: ii = (0.0_dp, 1.0_dp)
        real(dp) :: r, z
        r = x(1)
        z = x(2)
        br = cmplx(cos(0.8_dp*r + z), 0.1_dp*sin(r), dp)
        bp = cmplx(sin(r*z*0.5_dp), 0.3_dp*cos(z), dp)
        drbr = br + r*cmplx(-0.8_dp*sin(0.8_dp*r + z), 0.1_dp*cos(r), dp)
        dzbz = cmplx(-0.2_dp*exp(-0.2_dp*z)*sin(r), 0.0_dp, dp)
        v = (0.0_dp, 0.0_dp)
        v(1) = drbr/r + ii*real(mode_n, dp)*bp/r + dzbz
    end subroutine div_b

    ! ---- manufactured fields with zero trace on the annulus boundary ----

    pure function chi(x) result(c)
        real(dp), intent(in) :: x(2)
        real(dp) :: c, rho2
        if (rect_chi) then
            ! vanishes on the boundary of [2, 3] x [-0.5, 0.5]
            c = (x(1) - 2.0_dp)*(3.0_dp - x(1))*(x(2) + 0.5_dp)*(0.5_dp - x(2))
            return
        end if
        rho2 = (x(1) - ann_c)**2 + x(2)**2
        c = (rho2 - ann_r1**2)*(ann_r2**2 - rho2)
    end function chi

    recursive subroutine eval_id(id, x, v)
        !! 1 A*, 2 curl A*, 3 curl curl A* + A*, 4 B*, 5 div B*,
        !! 6 -grad div B* + B*, 7 f*, 8 grad f*, 9 -div grad f* + f*.
        integer, intent(in) :: id
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        complex(dp) :: w(3)
        select case (id)
        case (1)
            call vec_a(x, w)
            v = chi(x)*w
        case (2)
            call fd_curl(1, x, v)
        case (3)
            call fd_curl(2, x, v)
            call eval_id(1, x, w)
            v = v + w
        case (4)
            call vec_b(x, w)
            v = chi(x)*w
        case (5)
            call fd_div(4, x, v)
        case (6)
            call fd_grad(5, x, v)
            call eval_id(4, x, w)
            v = -v + w
        case (7)
            call scalar_f(x, w)
            v = chi(x)*w
        case (8)
            call fd_grad(7, x, v)
        case (9)
            call fd_div(8, x, v)
            call eval_id(7, x, w)
            v = -v + w
        case (10)
            ! gauge-fixed (A_phi = 0) potential, its curl and curl curl
            call vec_a(x, w)
            v = chi(x)*[w(1), (0.0_dp, 0.0_dp), w(3)]
        case (11)
            call fd_curl(10, x, v)
        case (12)
            call fd_curl(11, x, v)
        case default
            v = (0.0_dp, 0.0_dp)
        end select
    end subroutine eval_id

    recursive subroutine dpart(id, x, dir, d)
        !! Fourth-order central derivative of all three components.
        integer, intent(in) :: id, dir
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: d(3)
        complex(dp) :: f1(3), f2(3), f3(3), f4(3)
        real(dp) :: e(2)
        e = 0.0_dp
        e(dir) = fd_h
        call eval_id(id, x + e, f1)
        call eval_id(id, x - e, f2)
        call eval_id(id, x + 2.0_dp*e, f3)
        call eval_id(id, x - 2.0_dp*e, f4)
        d = (8.0_dp*(f1 - f2) - (f3 - f4))/(12.0_dp*fd_h)
    end subroutine dpart

    recursive subroutine fd_grad(id, x, v)
        integer, intent(in) :: id
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        complex(dp) :: dr(3), dz(3), f(3)
        call dpart(id, x, 1, dr)
        call dpart(id, x, 2, dz)
        call eval_id(id, x, f)
        v(1) = dr(1)
        v(2) = cmplx(0.0_dp, real(mode_n, dp), dp)*f(1)/x(1)
        v(3) = dz(1)
    end subroutine fd_grad

    recursive subroutine fd_curl(id, x, v)
        integer, intent(in) :: id
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        complex(dp) :: dr(3), dz(3), u(3)
        complex(dp), parameter :: ii = (0.0_dp, 1.0_dp)
        call dpart(id, x, 1, dr)
        call dpart(id, x, 2, dz)
        call eval_id(id, x, u)
        v(1) = ii*real(mode_n, dp)*u(3)/x(1) - dz(2)
        v(2) = dz(1) - dr(3)
        v(3) = (u(2) + x(1)*dr(2) - ii*real(mode_n, dp)*u(1))/x(1)
    end subroutine fd_curl

    recursive subroutine fd_div(id, x, v)
        integer, intent(in) :: id
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        complex(dp) :: dr(3), dz(3), u(3)
        call dpart(id, x, 1, dr)
        call dpart(id, x, 2, dz)
        call eval_id(id, x, u)
        v = (0.0_dp, 0.0_dp)
        v(1) = (u(1) + x(1)*dr(1))/x(1) &
            + cmplx(0.0_dp, real(mode_n, dp), dp)*u(2)/x(1) + dz(3)
    end subroutine fd_div

    subroutine fun_astar(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        call eval_id(1, x, v)
    end subroutine fun_astar
    subroutine fun_bstar(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        call eval_id(2, x, v)
    end subroutine fun_bstar
    subroutine fun_f1(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        call eval_id(3, x, v)
    end subroutine fun_f1
    subroutine fun_bs(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        call eval_id(4, x, v)
    end subroutine fun_bs
    subroutine fun_dbs(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        call eval_id(5, x, v)
    end subroutine fun_dbs
    subroutine fun_f2(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        call eval_id(6, x, v)
    end subroutine fun_f2
    subroutine fun_fs(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        call eval_id(7, x, v)
    end subroutine fun_fs
    subroutine fun_gfs(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        call eval_id(8, x, v)
    end subroutine fun_gfs
    subroutine fun_a0star(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        call eval_id(10, x, v)
    end subroutine fun_a0star
    subroutine fun_b0star(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        call eval_id(11, x, v)
    end subroutine fun_b0star
    subroutine fun_j0star(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        call eval_id(12, x, v)
    end subroutine fun_j0star
    subroutine fun_f0(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        call eval_id(9, x, v)
    end subroutine fun_f0

    subroutine solve_z(a, b, x, ok)
        type(csc_z_t), intent(in) :: a
        complex(dp), intent(in) :: b(:)
        complex(dp), intent(out) :: x(:)
        logical, intent(out) :: ok
        type(fortsparse_status_t) :: st

        call sparse_direct_solve_once(a, b, x, st)
        ok = status_ok(st)
    end subroutine solve_z

    subroutine project(fs, space, fun, coeff, ok)
        !! L2 projection onto the given space of an analytic field.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: space
        interface
            subroutine fun(x, v)
                import :: dp
                real(dp), intent(in) :: x(2)
                complex(dp), intent(out) :: v(3)
            end subroutine fun
        end interface
        complex(dp), allocatable, intent(out) :: coeff(:)
        logical, intent(out) :: ok
        complex(dp), allocatable :: b(:)
        type(csc_z_t) :: mz
        integer :: st

        call meridian_load(fs, space, fun, b, st)
        ok = st == 0
        if (.not. ok) return
        call meridian_real_to_complex(fs%mass(space), mz)
        allocate (coeff(size(b)))
        call solve_z(mz, b, coeff, ok)
    end subroutine project

    subroutine zadd(a, b, c)
        type(csc_z_t), intent(in) :: a, b
        type(csc_z_t), intent(out) :: c
        integer, allocatable :: r(:), cc(:)
        complex(dp), allocatable :: v(:)
        type(fortsparse_status_t) :: st
        integer :: j, p, n

        allocate (r(a%nnz + b%nnz), cc(a%nnz + b%nnz), v(a%nnz + b%nnz))
        n = 0
        do j = 1, a%ncol
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                n = n + 1
                r(n) = a%row_idx(p)
                cc(n) = j
                v(n) = a%val(p)
            end do
            do p = b%col_ptr(j), b%col_ptr(j + 1) - 1
                n = n + 1
                r(n) = b%row_idx(p)
                cc(n) = j
                v(n) = b%val(p)
            end do
        end do
        call csc_from_triplet(a%nrow, a%ncol, r, cc, v, c, st)
    end subroutine zadd

    subroutine zero_rows(a, mask)
        !! Replace rows and columns of masked dofs by the identity.
        type(csc_z_t), intent(inout) :: a
        logical, intent(in) :: mask(:)
        integer :: j, p

        do j = 1, a%ncol
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                if (mask(j) .or. mask(a%row_idx(p))) then
                    if (j == a%row_idx(p)) then
                        a%val(p) = (1.0_dp, 0.0_dp)
                    else
                        a%val(p) = (0.0_dp, 0.0_dp)
                    end if
                end if
            end do
        end do
    end subroutine zero_rows

end module meridian_test_fields
