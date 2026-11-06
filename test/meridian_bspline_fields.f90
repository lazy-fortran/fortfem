module meridian_bspline_fields
    !! Fields and caller quadratures for the per-mode B-spline tests.
    !!
    !! The manufactured potential A, its field B = curl_n A and the current
    !! curl_n B are FortSym leaves (test/generated); all vanish outside the
    !! disc |x - xc| < width.  The n = 0 toroidal field T = (0, b, 0) is
    !! divergence free with nonzero toroidal flux, so its potential is not
    !! compactly supported.  grad_n chi is a smooth curl-free field.
    use fortfem_kinds, only: dp
    use fortfem_generated_meridian_bspline_potential, only: &
        generated_meridian_bspline_potential
    use fortfem_generated_meridian_bspline_gauge, only: &
        generated_meridian_bspline_gauge
    use fortfem_generated_meridian_bspline_toroidal, only: &
        generated_meridian_bspline_toroidal
    use fortfem_meridian_bspline_n, only: meridian_bspline_t
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_special_elliptic, only: ellip_k, ellip_e
    implicit none
    private

    real(dp), parameter, public :: pi = acos(-1.0_dp)
    real(dp), parameter, public :: centre(2) = [2.5_dp, 0.05_dp]
    real(dp), parameter, public :: width = 0.3_dp

    public :: bump_fields, gauge_gradient, toroidal_fields
    public :: disc_quadrature, circle_quadrature, patch_quadrature
    public :: loop_potential, cross
    public :: manufactured_solve, loop_solve

contains

    subroutine manufactured_solve(nr, p, n, with_gauge, err_b, err_out, &
            ext, pmax, status)
        !! Recover B = curl_n A (plus T for n = 0) on [2, 3] x [-0.5, 0.5]
        !! from the source f = (B_s, curl_n A')_V + <n x B_s, A'>_dV on the
        !! disc V = {|x - centre| < 0.4}.  with_gauge adds grad_n chi to B_s,
        !! which changes both terms but not f.  Returns the relative L2
        !! errors of B on the patch and the L2 norm of B_h outside the
        !! support (relative to |B|), the ratio of exterior to interior
        !! energy and max |p| / max |u|.
        use fortfem_meridian_bspline_n, only: meridian_bspline_build, &
            meridian_bspline_exterior_build, meridian_bspline_solve, &
            meridian_bspline_load_curl, meridian_bspline_load_value, &
            meridian_bspline_eval, meridian_bspline_field_energy
        use fortfem_meridian_magnetic_coupling, only: meridian_exterior_t, &
            meridian_exterior_energy
        integer, intent(in) :: nr, p, n
        logical, intent(in) :: with_gauge
        real(dp), intent(out) :: err_b, err_out, ext, pmax
        integer, intent(out) :: status

        type(meridian_bspline_t) :: sp
        type(meridian_exterior_t) :: ex
        real(dp), allocatable :: x(:, :), w(:), nrm(:, :)
        complex(dp), allocatable :: v(:, :), f(:), u(:), pm(:), ah(:, :), &
            bh(:, :)
        complex(dp) :: b(3), g(3), nv(3)
        real(dp) :: num, den, out, eint
        integer :: k

        call meridian_bspline_build(2.0_dp, 3.0_dp, -0.5_dp, 0.5_dp, nr, nr, &
            p, sp, status)
        if (status /= 0) return
        call meridian_bspline_exterior_build(sp, n, ex, status)
        if (status /= 0) return
        allocate (f(sp%n1))
        f = (0.0_dp, 0.0_dp)
        call disc_quadrature(centre, 0.4_dp, 64, 128, x, w)
        allocate (v(3, size(w)))
        do k = 1, size(w)
            call exact_field(x(:, k), n, b)
            v(:, k) = b
            if (with_gauge) then
                call gauge_gradient(x(:, k), n, g)
                v(:, k) = v(:, k) + g
            end if
        end do
        call meridian_bspline_load_curl(sp, n, x, w, v, f, status)
        if (status /= 0) return
        if (with_gauge) then
            call circle_quadrature(centre, 0.4_dp, 1024, x, w, nrm)
            deallocate (v)
            allocate (v(3, size(w)))
            do k = 1, size(w)
                call gauge_gradient(x(:, k), n, g)
                nv = [cmplx(nrm(1, k), 0.0_dp, dp), (0.0_dp, 0.0_dp), &
                    cmplx(nrm(2, k), 0.0_dp, dp)]
                v(:, k) = cross(nv, g)
            end do
            call meridian_bspline_load_value(sp, n, x, w, v, f, status)
            if (status /= 0) return
        end if
        call meridian_bspline_solve(sp, ex, f, u, pm, status)
        if (status /= 0) return

        call patch_quadrature(sp, p + 3, x, w)
        allocate (ah(3, size(w)), bh(3, size(w)))
        call meridian_bspline_eval(sp, n, u, x, ah, bh, status)
        if (status /= 0) return
        num = 0.0_dp
        den = 0.0_dp
        out = 0.0_dp
        do k = 1, size(w)
            call exact_field(x(:, k), n, b)
            num = num + w(k)*sum(abs(bh(:, k) - b)**2)
            den = den + w(k)*sum(abs(b)**2)
            if (norm2(x(:, k) - centre) >= width) &
                out = out + w(k)*sum(abs(bh(:, k))**2)
        end do
        err_b = sqrt(num/den)
        err_out = sqrt(out/den)
        call meridian_bspline_field_energy(sp, n, u, eint, status)
        if (status /= 0) return
        ext = meridian_exterior_energy(ex, u)/eint
        pmax = maxval(abs(pm))/maxval(abs(u))
    end subroutine manufactured_solve

    subroutine exact_field(x, n, b)
        !! curl_n A, plus T for n = 0.
        real(dp), intent(in) :: x(2)
        integer, intent(in) :: n
        complex(dp), intent(out) :: b(3)
        complex(dp) :: a(3), c(3), t(3), ct(3)

        call bump_fields(x, n, a, b, c)
        if (n /= 0) return
        call toroidal_fields(x, t, ct)
        b = b + t
    end subroutine exact_field

    subroutine loop_solve(nr, p, a0, err, status)
        !! n = 0 field of the filament (a0, 0), mu0 I = 1, on
        !! [2, 3] x [-0.5, 0.5]; err = max |A_phi,h - A_phi| / max |A_phi|
        !! over points at distance >= 0.2 from the filament.
        use fortfem_meridian_bspline_n, only: meridian_bspline_build, &
            meridian_bspline_exterior_build, meridian_bspline_solve, &
            meridian_bspline_load_value, meridian_bspline_eval
        use fortfem_meridian_magnetic_coupling, only: meridian_exterior_t
        integer, intent(in) :: nr, p
        real(dp), intent(in) :: a0
        real(dp), intent(out) :: err
        integer, intent(out) :: status

        type(meridian_bspline_t) :: sp
        type(meridian_exterior_t) :: ex
        real(dp), allocatable :: x(:, :), w(:)
        complex(dp), allocatable :: f(:), u(:), pm(:), ah(:, :), bh(:, :)
        complex(dp) :: g(3, 1)
        real(dp) :: xs(2, 1), ws(1), ex_val, big
        integer :: i, j, k

        call meridian_bspline_build(2.0_dp, 3.0_dp, -0.5_dp, 0.5_dp, nr, nr, &
            p, sp, status)
        if (status /= 0) return
        call meridian_bspline_exterior_build(sp, 0, ex, status)
        if (status /= 0) return
        allocate (f(sp%n1))
        f = (0.0_dp, 0.0_dp)
        xs(:, 1) = [a0, 0.0_dp]
        ws(1) = a0
        g(:, 1) = [(0.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)]
        call meridian_bspline_load_value(sp, 0, xs, ws, g, f, status)
        if (status /= 0) return
        call meridian_bspline_solve(sp, ex, f, u, pm, status)
        if (status /= 0) return
        allocate (x(2, 21*21), w(21*21))
        k = 0
        do j = 0, 20
            do i = 0, 20
                if (norm2([2.0_dp + 0.05_dp*i, -0.5_dp + 0.05_dp*j] &
                    - [a0, 0.0_dp]) < 0.2_dp) cycle
                k = k + 1
                x(:, k) = [2.0_dp + 0.05_dp*i, -0.5_dp + 0.05_dp*j]
            end do
        end do
        allocate (ah(3, k), bh(3, k))
        call meridian_bspline_eval(sp, 0, u, x(:, :k), ah, bh, status)
        if (status /= 0) return
        err = 0.0_dp
        big = 0.0_dp
        do i = 1, k
            ex_val = loop_potential(a0, x(:, i))
            err = max(err, abs(ah(2, i) - ex_val))
            big = max(big, abs(ex_val))
        end do
        err = err/big
    end subroutine loop_solve

    subroutine bump_fields(x, n, a, b, c)
        !! A, B = curl_n A and curl_n B at x for mode n.
        real(dp), intent(in) :: x(2)
        integer, intent(in) :: n
        complex(dp), intent(out) :: a(3), b(3), c(3)
        real(dp) :: v(15)

        a = (0.0_dp, 0.0_dp)
        b = a
        c = a
        if (norm2(x - centre) >= width) return
        call generated_meridian_bspline_potential(x(1), x(2), real(n, dp), &
            centre(1), centre(2), width, v(1), v(2), v(3), v(4), v(5), v(6), &
            v(7), v(8), v(9), v(10), v(11), v(12), v(13), v(14), v(15))
        a = [cmplx(v(1), 0.0_dp, dp), cmplx(0.0_dp, v(2), dp), &
            cmplx(v(3), 0.0_dp, dp)]
        b = [cmplx(v(4), v(5), dp), cmplx(v(6), v(7), dp), &
            cmplx(v(8), v(9), dp)]
        c = [cmplx(v(10), v(11), dp), cmplx(v(12), v(13), dp), &
            cmplx(v(14), v(15), dp)]
    end subroutine bump_fields

    subroutine toroidal_fields(x, t, c)
        !! T = (0, b, 0) and curl_0 T.
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: t(3), c(3)
        real(dp) :: tp, cr, cz

        t = (0.0_dp, 0.0_dp)
        c = t
        if (norm2(x - centre) >= width) return
        call generated_meridian_bspline_toroidal(x(1), x(2), centre(1), &
            centre(2), width, tp, cr, cz)
        t(2) = tp
        c(1) = cr
        c(3) = cz
    end subroutine toroidal_fields

    subroutine gauge_gradient(x, n, g)
        real(dp), intent(in) :: x(2)
        integer, intent(in) :: n
        complex(dp), intent(out) :: g(3)
        real(dp) :: gr, gi, gz

        call generated_meridian_bspline_gauge(x(1), x(2), real(n, dp), gr, gi, &
            gz)
        g = [cmplx(gr, 0.0_dp, dp), cmplx(0.0_dp, gi, dp), cmplx(gz, 0.0_dp, dp)]
    end subroutine gauge_gradient

    pure function cross(u, v) result(w)
        !! Cross product in the right-handed frame (R, phi, Z).
        complex(dp), intent(in) :: u(3), v(3)
        complex(dp) :: w(3)

        w = [u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), &
            u(1)*v(2) - u(2)*v(1)]
    end function cross

    subroutine disc_quadrature(xc, radius, nrad, nang, x, w)
        !! Polar Gauss rule on the disc, weights of R dR dZ.
        real(dp), intent(in) :: xc(2), radius
        integer, intent(in) :: nrad, nang
        real(dp), allocatable, intent(out) :: x(:, :), w(:)
        real(dp) :: rho(nrad), wr(nrad), th
        integer :: i, j, k

        call gauss_legendre_ab(nrad, 0.0_dp, radius, rho, wr)
        allocate (x(2, nrad*nang), w(nrad*nang))
        k = 0
        do j = 1, nang
            th = 2.0_dp*pi*(real(j, dp) - 0.5_dp)/real(nang, dp)
            do i = 1, nrad
                k = k + 1
                x(:, k) = xc + rho(i)*[cos(th), sin(th)]
                w(k) = wr(i)*rho(i)*(2.0_dp*pi/real(nang, dp))*x(1, k)
            end do
        end do
    end subroutine disc_quadrature

    subroutine circle_quadrature(xc, radius, nang, x, w, normal)
        !! Trapezoidal rule on the circle, weights of R ds, outward normal.
        real(dp), intent(in) :: xc(2), radius
        integer, intent(in) :: nang
        real(dp), allocatable, intent(out) :: x(:, :), w(:), normal(:, :)
        real(dp) :: th
        integer :: j

        allocate (x(2, nang), w(nang), normal(2, nang))
        do j = 1, nang
            th = 2.0_dp*pi*(real(j, dp) - 0.5_dp)/real(nang, dp)
            normal(:, j) = [cos(th), sin(th)]
            x(:, j) = xc + radius*normal(:, j)
            w(j) = radius*(2.0_dp*pi/real(nang, dp))*x(1, j)
        end do
    end subroutine circle_quadrature

    subroutine patch_quadrature(sp, nq, x, w)
        !! Tensor Gauss rule on the knot cells of the patch, weights of
        !! R dR dZ.
        type(meridian_bspline_t), intent(in) :: sp
        integer, intent(in) :: nq
        real(dp), allocatable, intent(out) :: x(:, :), w(:)
        real(dp) :: xr(nq), wr(nq), xz(nq), wz(nq)
        integer :: sr, sz, i, j, k, ncell

        ncell = (sp%mr - sp%p)*(sp%mz - sp%p)
        allocate (x(2, ncell*nq*nq), w(ncell*nq*nq))
        k = 0
        do sz = sp%p + 1, sp%mz
            call gauss_legendre_ab(nq, sp%kz(sz), sp%kz(sz + 1), xz, wz)
            do sr = sp%p + 1, sp%mr
                call gauss_legendre_ab(nq, sp%kr(sr), sp%kr(sr + 1), xr, wr)
                do j = 1, nq
                    do i = 1, nq
                        k = k + 1
                        x(:, k) = [xr(i), xz(j)]
                        w(k) = wr(i)*wz(j)*xr(i)
                    end do
                end do
            end do
        end do
    end subroutine patch_quadrature

    function loop_potential(a0, x) result(aphi)
        !! A_phi of a circular filament of radius a0 in Z = 0 with
        !! mu0 I = 1: (sqrt(a0 / R) / (pi k)) ((1 - k^2/2) K(k) - E(k)),
        !! k^2 = 4 a0 R / ((a0 + R)^2 + Z^2) (Jackson, Classical
        !! Electrodynamics, Sec. 5.5).
        real(dp), intent(in) :: a0, x(2)
        real(dp) :: aphi, m

        m = 4.0_dp*a0*x(1)/((a0 + x(1))**2 + x(2)**2)
        aphi = sqrt(a0/x(1))/(pi*sqrt(m)) &
            *((1.0_dp - 0.5_dp*m)*ellip_k(m) - ellip_e(m))
    end function loop_potential

end module meridian_bspline_fields
