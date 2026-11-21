module fortfem_curve_chart
    !! Charts (s1, s2, t) -> x in R^3 for toroidal extensions of a section
    !! (doc/design/multipatch_spline_spaces.md, toroidal direction).
    !!
    !! Cylindrical: x = (s1 cos phi, s1 sin phi, s2) with phi = 2 pi t/nfp,
    !! so t in [0, 1) covers one field period.
    !! Guide curve: x = r(t) + s1 e1(t) + s2 e2(t) on the planes normal to a
    !! closed curve r. (e1, e2) is the rotation-minimising frame of the
    !! double reflection method (Wang, Juttler, Zheng, Liu, ACM TOG 27, 2008)
    !! with its holonomy angle spread linearly over t, so the frame closes.
    !! r, e1 and e2 are trigonometric interpolants of the uniformly sampled
    !! curve and frames; the curve may cross itself in projection (figure-8).
    use fortfem_kinds, only: dp
    implicit none
    private

    integer, parameter, public :: CHART_CYLINDRICAL = 0
    integer, parameter, public :: CHART_CURVE = 1

    type, public :: curve_chart_t
        integer :: kind = CHART_CYLINDRICAL
        integer :: nfp = 1
        !> Fourier coefficients of r, e1, e2: (3, 0:nmax) cosine, sine.
        real(dp), allocatable :: rc(:, :), rs(:, :), e1c(:, :), e1s(:, :)
        real(dp), allocatable :: e2c(:, :), e2s(:, :)
    end type curve_chart_t

    public :: curve_chart_cylindrical
    public :: curve_chart_from_samples
    public :: curve_chart_map

    real(dp), parameter :: two_pi = 2.0_dp*acos(-1.0_dp)

contains

    subroutine curve_chart_cylindrical(chart, nfp)
        type(curve_chart_t), intent(out) :: chart
        integer, intent(in), optional :: nfp

        chart%kind = CHART_CYLINDRICAL
        chart%nfp = 1
        if (present(nfp)) chart%nfp = max(nfp, 1)
    end subroutine curve_chart_cylindrical

    subroutine curve_chart_from_samples(points, chart, status)
        !! points(3, n): the closed guide curve at t = (k - 1)/n, n >= 8, not
        !! repeated at the end. The frame is computed on 4096 points of the
        !! interpolated curve and kept to 64 harmonics. Status 1: invalid
        !! input, 2: a degenerate tangent.
        real(dp), intent(in) :: points(:, :)
        type(curve_chart_t), intent(out) :: chart
        integer, intent(out) :: status

        integer, parameter :: fine = 4096, frame_modes = 64
        real(dp), allocatable :: x(:, :), tangent(:, :), e1(:, :), e2(:, :)
        real(dp) :: v1(3), v2(3), rl(3), tl(3), c1, c2, angle, d(3), cs, sn
        integer :: n, k

        status = 1
        n = size(points, 2)
        if (size(points, 1) /= 3 .or. n < 8) return
        chart%kind = CHART_CURVE
        call fit(points, (n - 1)/2, chart%rc, chart%rs)
        allocate(x(3, fine), tangent(3, fine), e1(3, fine + 1), e2(3, fine))
        status = 2
        do k = 1, fine
            call evaluate(chart%rc, chart%rs, real(k - 1, dp)/fine, x(:, k), d)
            if (norm2(d) == 0.0_dp) return
            tangent(:, k) = d/norm2(d)
        end do
        e1(:, 1) = any_normal(tangent(:, 1))
        do k = 1, fine
            v1 = x(:, modulo(k, fine) + 1) - x(:, k)
            c1 = dot_product(v1, v1)
            rl = e1(:, k) - (2.0_dp/c1)*dot_product(v1, e1(:, k))*v1
            tl = tangent(:, k) - (2.0_dp/c1)*dot_product(v1, tangent(:, k))*v1
            v2 = tangent(:, modulo(k, fine) + 1) - tl
            c2 = dot_product(v2, v2)
            if (c2 > 0.0_dp) then
                e1(:, k + 1) = rl - (2.0_dp/c2)*dot_product(v2, rl)*v2
            else
                e1(:, k + 1) = rl
            end if
        end do
        angle = atan2(dot_product(cross(e1(:, 1), e1(:, fine + 1)), tangent(:, 1)), &
            dot_product(e1(:, 1), e1(:, fine + 1)))
        do k = 1, fine
            cs = cos(-angle*(k - 1)/fine)
            sn = sin(-angle*(k - 1)/fine)
            e1(:, k) = cs*e1(:, k) + sn*cross(tangent(:, k), e1(:, k))
            e1(:, k) = e1(:, k) - dot_product(e1(:, k), tangent(:, k))*tangent(:, k)
            e1(:, k) = e1(:, k)/norm2(e1(:, k))
            e2(:, k) = cross(tangent(:, k), e1(:, k))
        end do
        call fit(e1(:, :fine), frame_modes, chart%e1c, chart%e1s)
        call fit(e2, frame_modes, chart%e2c, chart%e2s)
        status = 0
    end subroutine curve_chart_from_samples

    pure subroutine curve_chart_map(chart, s, t, x, jacobian)
        !! x and jacobian(:, j) = dx/d(s1, s2, t)_j.
        type(curve_chart_t), intent(in) :: chart
        real(dp), intent(in) :: s(2), t
        real(dp), intent(out) :: x(3), jacobian(3, 3)

        real(dp) :: phi, r(3), dr(3), a(3), da(3), b(3), db(3)

        if (chart%kind == CHART_CYLINDRICAL) then
            phi = two_pi*t/chart%nfp
            x = [s(1)*cos(phi), s(1)*sin(phi), s(2)]
            jacobian(:, 1) = [cos(phi), sin(phi), 0.0_dp]
            jacobian(:, 2) = [0.0_dp, 0.0_dp, 1.0_dp]
            jacobian(:, 3) = two_pi/chart%nfp*s(1)*[-sin(phi), cos(phi), 0.0_dp]
            return
        end if
        call evaluate(chart%rc, chart%rs, t, r, dr)
        call evaluate(chart%e1c, chart%e1s, t, a, da)
        call evaluate(chart%e2c, chart%e2s, t, b, db)
        x = r + s(1)*a + s(2)*b
        jacobian(:, 1) = a
        jacobian(:, 2) = b
        jacobian(:, 3) = dr + s(1)*da + s(2)*db
    end subroutine curve_chart_map

    subroutine fit(values, nmax, c, s)
        !! Trigonometric interpolant of values(3, n) at t = (k - 1)/n for odd
        !! n (for even n the Nyquist mode is dropped).
        real(dp), intent(in) :: values(:, :)
        integer, intent(in) :: nmax
        real(dp), allocatable, intent(out) :: c(:, :), s(:, :)

        integer :: n, k, m
        real(dp) :: t

        n = size(values, 2)
        allocate(c(3, 0:nmax), s(3, 0:nmax))
        c = 0.0_dp
        s = 0.0_dp
        do k = 1, n
            t = real(k - 1, dp)/n
            do m = 0, nmax
                c(:, m) = c(:, m) + values(:, k)*cos(two_pi*m*t)
                s(:, m) = s(:, m) + values(:, k)*sin(two_pi*m*t)
            end do
        end do
        c = 2.0_dp*c/n
        s = 2.0_dp*s/n
        c(:, 0) = 0.5_dp*c(:, 0)
    end subroutine fit

    pure subroutine evaluate(c, s, t, x, dx)
        real(dp), intent(in) :: c(:, 0:), s(:, 0:), t
        real(dp), intent(out) :: x(3), dx(3)
        integer :: m

        x = c(:, 0)
        dx = 0.0_dp
        do m = 1, ubound(c, 2)
            x = x + c(:, m)*cos(two_pi*m*t) + s(:, m)*sin(two_pi*m*t)
            dx = dx + two_pi*m*(-c(:, m)*sin(two_pi*m*t) + s(:, m)*cos(two_pi*m*t))
        end do
    end subroutine evaluate

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)

        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross

    pure function any_normal(t) result(n)
        real(dp), intent(in) :: t(3)
        real(dp) :: n(3), e(3)

        e = [1.0_dp, 0.0_dp, 0.0_dp]
        if (abs(t(1)) > 0.9_dp) e = [0.0_dp, 1.0_dp, 0.0_dp]
        n = e - dot_product(e, t)*t
        n = n/norm2(n)
    end function any_normal

end module fortfem_curve_chart
