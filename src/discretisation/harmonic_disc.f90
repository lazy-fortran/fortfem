module fortfem_harmonic_disc
    !! Harmonic map of a plane domain onto the unit disc and the pullback of
    !! a disc O-grid (doc/HARMONIC_DISC_MAP.md).
    !!
    !! G = (xi, eta) solves Laplace's equation in the domain bounded by a
    !! smooth simple closed curve, with the curve mapped onto the unit circle
    !! by normalised arc length, or (grading > 0) by the normalised integral
    !! of the arc length weighted with (1 + (grading kappa L/2 pi)^2)^(1/3),
    !! which places more of the circle, and so more cells of a pulled-back
    !! mesh, where the curvature kappa is large. The target is convex and
    !! either correspondence is a homeomorphism of the curves, so G is a
    !! diffeomorphism for every such curve, nonconvex ones included
    !! (Rado-Kneser-Choquet). G is computed without a mesh: xi = Re v, v the
    !! Cauchy integral of a real density on the curve (the double-layer
    !! potential), discretised by the trapezoidal rule on samples that are
    !! equally spaced in a smooth periodic parameter. Inside, v and v' follow
    !! from the boundary values of v by the barycentric form of Cauchy's
    !! formula, which stays accurate up to the curve.
    !!
    !! harmonic_disc_pull_back replaces the block maps of a C0 spline space on
    !! the O-grid of the unit disc by the tensor Greville interpolants of
    !! G^-1 o (exact disc block map): a mesh of the domain that is valid by
    !! construction in the limit, with the certificate and the wall errors of
    !! spline_jacobian_certificate and harmonic_disc_wall_error for the
    !! discrete one.
    use fortfem_kinds, only: dp
    use fortfem_bspline_feec, only: evaluate_bspline_basis
    use fortfem_spline_layout, only: layout_side_point, LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_greville, SPLINE_C0
    use fortfem_spline_eval, only: spline_geometry_at
    use fortfem_generated_harmonic_disc_cauchy_terms, only: harmonic_disc_cauchy_terms
    use fortfem_generated_harmonic_disc_nystrom_offdiag, only: &
        harmonic_disc_nystrom_offdiag
    use fortfem_generated_harmonic_disc_nystrom_self, only: harmonic_disc_nystrom_self
    use fortfem_generated_harmonic_disc_interior_limit_self, only: &
        harmonic_disc_interior_limit_self
    use fortfem_generated_harmonic_disc_boundary_weight, only: &
        harmonic_disc_boundary_weight
    use fortfem_generated_spline_jacobian_determinant, only: spline_jacobian_determinant
    use fortfem_generated_spline_bernstein_interpolation, only: &
        spline_bernstein_interpolation
    implicit none
    private

    type, public :: harmonic_disc_t
        integer :: n = 0
        complex(dp), allocatable :: y(:) ! curve nodes x1 + i x2
        complex(dp), allocatable :: w(:) ! y' h
        complex(dp), allocatable :: v(:, :) ! (n, 2) boundary values of v
        complex(dp), allocatable :: coef(:) ! Fourier coefficients of y
        complex(dp), allocatable :: arc_coef(:) ! and of |y'|
        !> Fourier coefficients of the density of the boundary
        !> correspondence in the curve parameter (|y'| without grading).
        complex(dp), allocatable :: weight_coef(:)
        real(dp) :: grading = 0.0_dp
        complex(dp), allocatable :: table(:, :) ! G^-1 on a polar grid
        real(dp) :: length = 0.0_dp, area = 0.0_dp
        !> Largest error of the boundary values of G at the nodes.
        real(dp) :: residual = 0.0_dp
    end type harmonic_disc_t

    integer, parameter :: nring = 32, nangle = 128
    real(dp), parameter :: two_pi = 2.0_dp*acos(-1.0_dp)
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)

    public :: harmonic_disc_build
    public :: harmonic_disc_forward
    public :: harmonic_disc_inverse
    public :: harmonic_disc_boundary
    public :: harmonic_disc_pull_back
    public :: harmonic_disc_wall_error
    public :: spline_jacobian_certificate

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

contains

    subroutine harmonic_disc_build(map, samples, status, grading)
        !! samples(2, n): the curve, counter-clockwise, equally spaced in a
        !! smooth periodic parameter. The first sample is mapped to the disc
        !! angle 0. grading (default 0: arc length) weights the boundary
        !! correspondence with the curvature (module header); the samples
        !! must then resolve the weight too (map%residual). Status 1: invalid
        !! input (fewer than 8 samples, clockwise, negative grading), 2:
        !! singular system, 3: the inversion on the polar grid failed (the
        !! curve is not resolved by the samples or not simple).
        type(harmonic_disc_t), intent(out) :: map
        real(dp), intent(in) :: samples(:, :)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: grading

        complex(dp), allocatable :: dy(:), ddy(:), cm(:, :), work(:)
        real(dp), allocatable :: a(:, :), rhs(:, :)
        complex(dp) :: e, z, zeta
        real(dp) :: h, angle, cr, ci, weight
        integer, allocatable :: pivot(:)
        integer :: n, j, k, q, info

        status = 1
        n = size(samples, 2)
        if (size(samples, 1) /= 2 .or. n < 8) return
        if (present(grading)) map%grading = grading
        if (map%grading < 0.0_dp) return
        map%n = n
        h = two_pi/real(n, dp)
        allocate(map%y(n), map%w(n), map%v(n, 2), dy(n), ddy(n), work(n))
        map%y = cmplx(samples(1, :), samples(2, :), dp)
        call fourier_coefficients(map%y, map%coef)
        call derivative(map%coef, 1, dy)
        call derivative(map%coef, 2, ddy)
        map%w = dy*h
        map%area = 0.5_dp*sum(aimag(conjg(map%y)*map%w))
        if (map%area <= 0.0_dp) return
        work = cmplx(abs(dy), 0.0_dp, dp)
        call fourier_coefficients(work, map%arc_coef)
        map%length = two_pi*real(map%arc_coef(0), dp)
        if (map%grading > 0.0_dp) then
            do k = 1, n
                call harmonic_disc_boundary_weight(real(dy(k), dp), aimag(dy(k)), &
                    real(ddy(k), dp), aimag(ddy(k)), map%grading*map%length/two_pi, &
                    weight)
                work(k) = cmplx(weight, 0.0_dp, dp)
            end do
            call fourier_coefficients(work, map%weight_coef)
        else
            map%weight_coef = map%arc_coef
        end if
        allocate(a(n, n), rhs(n, 2), pivot(n), cm(n, n))
        do k = 1, n
            do j = 1, n
                if (j == k) then
                    call harmonic_disc_nystrom_self(real(ddy(k), dp), aimag(ddy(k)), &
                        real(dy(k), dp), aimag(dy(k)), h, cr, ci)
                    cm(k, j) = cmplx(cr, ci, dp)
                else
                    call harmonic_disc_nystrom_offdiag(real(map%w(j), dp), &
                        aimag(map%w(j)), real(map%y(j), dp), aimag(map%y(j)), &
                        real(map%y(k), dp), aimag(map%y(k)), cr, ci)
                    cm(k, j) = cmplx(cr, ci, dp)
                end if
            end do
        end do
        a = real(cm, dp)
        do k = 1, n
            a(k, k) = a(k, k) + 0.5_dp
            angle = boundary_measure(map, h*real(k - 1, dp))/ &
                real(map%weight_coef(0), dp)
            rhs(k, :) = [cos(angle), sin(angle)]
        end do
        ddy = cmplx(rhs(:, 1), rhs(:, 2), dp) ! boundary values, kept
        call dgesv(n, 2, a, n, pivot, rhs, n, info)
        status = 2
        if (info /= 0) return
        ! interior limit of the Cauchy integral of the density at the nodes
        map%residual = 0.0_dp
        do q = 1, 2
            block
                complex(dp), allocatable :: mu_coef(:)

                work = cmplx(rhs(:, q), 0.0_dp, dp)
                call fourier_coefficients(work, mu_coef)
                call derivative(mu_coef, 1, dy)
            end block
            do k = 1, n
                call harmonic_disc_interior_limit_self(real(dy(k), dp), h, cr, ci)
                e = cmplx(cr, ci, dp)
                do j = 1, n
                    if (j /= k) e = e + cm(k, j)*(rhs(j, q) - rhs(k, q))
                end do
                map%v(k, q) = rhs(k, q) + e
                if (q == 1) then
                    map%residual = max(map%residual, &
                        abs(real(map%v(k, q), dp) - real(ddy(k), dp)))
                else
                    map%residual = max(map%residual, &
                        abs(real(map%v(k, q), dp) - aimag(ddy(k))))
                end if
            end do
        end do
        ! G^-1 on the polar grid rho = 1 - j/nring, marching inward
        status = 3
        allocate(map%table(nangle, 0:nring - 1))
        do k = 1, nangle
            angle = two_pi*real(k - 1, dp)/real(nangle, dp)
            call harmonic_disc_boundary(map, angle, z)
            map%table(k, 0) = z
            do j = 1, nring - 1
                zeta = (1.0_dp - real(j, dp)/real(nring, dp))*exp(iu*angle)
                call newton(map, zeta, z, info)
                if (info /= 0) return
                map%table(k, j) = z
            end do
        end do
        status = 0
    end subroutine harmonic_disc_build

    subroutine fourier_coefficients(f, coef)
        !! coef(m), |m| <= (n - 1)/2, of the trigonometric interpolant of the
        !! nodal values f (direct sums; n is a few hundred).
        complex(dp), intent(in) :: f(:)
        complex(dp), allocatable, intent(out) :: coef(:)

        integer :: n, m, j

        n = size(f)
        allocate(coef(-(n - 1)/2:(n - 1)/2))
        do m = lbound(coef, 1), ubound(coef, 1)
            coef(m) = (0.0_dp, 0.0_dp)
            do j = 1, n
                coef(m) = coef(m) + f(j)*exp(-iu*two_pi*real(m*(j - 1), dp)/ &
                    real(n, dp))
            end do
            coef(m) = coef(m)/real(n, dp)
        end do
    end subroutine fourier_coefficients

    subroutine derivative(coef, order, df)
        !! Derivative of the given order at the n nodes.
        complex(dp), intent(out) :: df(:)
        complex(dp), intent(in) :: coef(-(size(df) - 1)/2:)
        integer, intent(in) :: order

        integer :: n, k, m

        n = size(df)
        df = (0.0_dp, 0.0_dp)
        do m = lbound(coef, 1), ubound(coef, 1)
            do k = 1, n
                df(k) = df(k) + (iu*real(m, dp))**order*coef(m)* &
                    exp(iu*two_pi*real(m*(k - 1), dp)/real(n, dp))
            end do
        end do
    end subroutine derivative

    pure subroutine curve_point(map, t, y, dy, ddy)
        !! The curve and its first two derivatives at the parameter t.
        type(harmonic_disc_t), intent(in) :: map
        real(dp), intent(in) :: t
        complex(dp), intent(out) :: y, dy, ddy

        complex(dp) :: e
        integer :: m

        y = (0.0_dp, 0.0_dp)
        dy = y
        ddy = y
        do m = lbound(map%coef, 1), ubound(map%coef, 1)
            e = map%coef(m)*exp(iu*real(m, dp)*t)
            y = y + e
            dy = dy + iu*real(m, dp)*e
            ddy = ddy - real(m*m, dp)*e
        end do
    end subroutine curve_point

    pure real(dp) function boundary_measure(map, t) result(s)
        !! Integral from 0 to t of the interpolant of the density of the
        !! boundary correspondence (the arc length without grading); the
        !! whole curve has 2 pi weight_coef(0).
        type(harmonic_disc_t), intent(in) :: map
        real(dp), intent(in) :: t

        integer :: m

        s = real(map%weight_coef(0), dp)*t
        do m = lbound(map%weight_coef, 1), ubound(map%weight_coef, 1)
            if (m == 0) cycle
            s = s + real(map%weight_coef(m)*(exp(iu*real(m, dp)*t) - 1.0_dp)/ &
                (iu*real(m, dp)), dp)
        end do
    end function boundary_measure

    pure real(dp) function boundary_density(map, t) result(s)
        type(harmonic_disc_t), intent(in) :: map
        real(dp), intent(in) :: t

        integer :: m

        s = 0.0_dp
        do m = lbound(map%weight_coef, 1), ubound(map%weight_coef, 1)
            s = s + real(map%weight_coef(m)*exp(iu*real(m, dp)*t), dp)
        end do
    end function boundary_density

    pure subroutine harmonic_disc_boundary(map, angle, point, parameter)
        !! The curve point with the disc angle `angle` (and its curve
        !! parameter): the boundary correspondence inverted by Newton's
        !! method.
        type(harmonic_disc_t), intent(in) :: map
        real(dp), intent(in) :: angle
        complex(dp), intent(out) :: point
        real(dp), intent(out), optional :: parameter

        complex(dp) :: dy, ddy
        real(dp) :: target, t, step
        integer :: it

        target = modulo(angle, two_pi)*real(map%weight_coef(0), dp)
        t = modulo(angle, two_pi)
        do it = 1, 50
            step = (boundary_measure(map, t) - target)/boundary_density(map, t)
            t = t - step
            if (abs(step) < 1.0e-14_dp) exit
        end do
        call curve_point(map, t, point, dy, ddy)
        if (present(parameter)) parameter = t
    end subroutine harmonic_disc_boundary

    pure recursive subroutine harmonic_disc_forward(map, z, g, dg)
        !! G(z) = xi + i eta at the interior point z, and the derivatives of
        !! the two holomorphic functions: grad xi = (Re, -Im) dg(1), grad eta
        !! the same of dg(2).
        type(harmonic_disc_t), intent(in) :: map
        complex(dp), intent(in) :: z
        complex(dp), intent(out) :: g, dg(2)

        complex(dp) :: q, den, f(2), d(2), r
        real(dp) :: qr, qi, pr, pim
        integer :: j

        den = (0.0_dp, 0.0_dp)
        f = den
        do j = 1, map%n
            r = map%y(j) - z
            if (abs(r) < 1.0e-15_dp*abs(map%w(j))*real(map%n, dp)) then
                ! on a node: its values, the derivative from just inside
                call harmonic_disc_forward(map, z + 1.0e-7_dp*iu*map%w(j), g, dg)
                g = cmplx(real(map%v(j, 1), dp), real(map%v(j, 2), dp), dp)
                return
            end if
            call harmonic_disc_cauchy_terms(real(map%y(j), dp), aimag(map%y(j)), &
                real(z, dp), aimag(z), real(map%w(j), dp), aimag(map%w(j)), qr, qi, pr, pim)
            q = cmplx(qr, qi, dp)
            den = den + q
            f = f + q*map%v(j, :)
        end do
        f = f/den
        d = (0.0_dp, 0.0_dp)
        do j = 1, map%n
            r = map%y(j) - z
            call harmonic_disc_cauchy_terms(real(map%y(j), dp), aimag(map%y(j)), &
                real(z, dp), aimag(z), real(map%w(j), dp), aimag(map%w(j)), qr, qi, pr, pim)
            d = d + cmplx(pr, pim, dp)*(map%v(j, :) - f)
        end do
        dg = d/den
        g = cmplx(real(f(1), dp), real(f(2), dp), dp)
    end subroutine harmonic_disc_forward

    pure subroutine newton(map, zeta, z, status)
        !! z with G(z) = zeta from the guess z; a step is accepted only when
        !! the residual falls.
        type(harmonic_disc_t), intent(in) :: map
        complex(dp), intent(in) :: zeta
        complex(dp), intent(inout) :: z
        integer, intent(out) :: status

        complex(dp) :: g, dg(2), step, zn, gn, dgn(2)
        real(dp) :: a, b, c, d, det, r(2), lambda
        integer :: it, half

        call harmonic_disc_forward(map, z, g, dg)
        status = 0
        do it = 1, 60
            r = [real(g - zeta, dp), aimag(g - zeta)]
            if (norm2(r) < 1.0e-13_dp) return
            a = real(dg(1), dp)
            b = -aimag(dg(1))
            c = real(dg(2), dp)
            d = -aimag(dg(2))
            det = a*d - b*c
            if (det == 0.0_dp) exit
            step = -cmplx(d*r(1) - b*r(2), -c*r(1) + a*r(2), dp)/det
            lambda = 1.0_dp
            do half = 1, 30
                zn = z + lambda*step
                call harmonic_disc_forward(map, zn, gn, dgn)
                if (abs(gn - zeta) < norm2(r) .and. abs(gn) < 1.5_dp) exit
                lambda = 0.5_dp*lambda
            end do
            if (half > 30) exit
            z = zn
            g = gn
            dg = dgn
        end do
        status = 1
        if (abs(g - zeta) < 1.0e-10_dp) status = 0
    end subroutine newton

    pure subroutine harmonic_disc_inverse(map, zeta, z, status)
        !! z = G^-1(zeta) for |zeta| <= 1: on the circle the curve point,
        !! inside Newton's method from the nearest point of the polar table,
        !! on failure by continuation along the segment from it. Status 1:
        !! not found.
        type(harmonic_disc_t), intent(in) :: map
        complex(dp), intent(in) :: zeta
        complex(dp), intent(out) :: z
        integer, intent(out) :: status

        complex(dp) :: z0, zeta0
        real(dp) :: rho, angle
        integer :: j, k, m

        rho = abs(zeta)
        angle = atan2(aimag(zeta), real(zeta, dp))
        status = 0
        if (rho >= 1.0_dp - 1.0e-13_dp) then
            call harmonic_disc_boundary(map, angle, z)
            return
        end if
        j = min(max(nint((1.0_dp - rho)*real(nring, dp)), 1), nring - 1)
        k = modulo(nint(angle/two_pi*real(nangle, dp)), nangle) + 1
        z0 = map%table(k, j)
        zeta0 = (1.0_dp - real(j, dp)/real(nring, dp))* &
            exp(iu*two_pi*real(k - 1, dp)/real(nangle, dp))
        z = z0
        call newton(map, zeta, z, status)
        if (status == 0) return
        z = z0
        do m = 1, 16
            call newton(map, zeta0 + (zeta - zeta0)*real(m, dp)/16.0_dp, z, status)
            if (status /= 0) return
        end do
    end subroutine harmonic_disc_inverse

    subroutine harmonic_disc_pull_back(space, map, status, preserve_area)
        !! Replace the block maps of a C0 space on layout_ogrid of the unit
        !! disc (centre 0, radius 1, no annulus) by the tensor Greville
        !! interpolants of G^-1 o (disc block map), the disc block maps being
        !! the exact ruled maps between the circular arcs and the core. The
        !! wall of the mesh is one periodic spline that interpolates the
        !! curve (wall_spline_data). With preserve_area the wall data are
        !! moved along the curve normal by one constant so that the mesh
        !! encloses the area of the curve.
        !! Status 1: not such a space, 2: an inversion failed, 3: singular
        !! interpolation.
        type(spline_space_t), intent(inout) :: space
        type(harmonic_disc_t), intent(in) :: map
        integer, intent(out) :: status
        logical, intent(in), optional :: preserve_area

        real(dp), allocatable :: data(:, :), normal(:, :), g1(:), g2(:)
        real(dp), allocatable :: wall_angle(:)
        logical, allocatable :: on_wall(:)
        real(dp) :: corner(2, 4), p(2), q(2), tangent(2), tau, t, offset
        real(dp) :: distance, area_error, previous_offset, previous_error, step
        complex(dp) :: z, dy, ddy
        integer :: k, i, j, n1, n2, side, info, pass, index
        logical :: keep_area

        status = 1
        if (space%regularity /= SPLINE_C0 .or. space%block_count /= 5) return
        if (any(space%layout%periodic) .or. any(space%layout%coons)) return
        if (space%layout%arc_side(1) /= 0 .or. &
            any(space%layout%arc_side(2:5) == 0)) return
        keep_area = .false.
        if (present(preserve_area)) keep_area = preserve_area
        allocate(data(2, space%broken_count), normal(2, space%broken_count))
        allocate(wall_angle(space%broken_count), on_wall(space%broken_count))
        normal = 0.0_dp
        on_wall = .false.
        do k = 1, space%block_count
            associate (blk => space%blocks(k), layout => space%layout)
                n1 = blk%n(1)
                n2 = blk%n(2)
                allocate(g1(n1), g2(n2))
                call spline_greville(blk%knots1, space%degree, g1)
                call spline_greville(blk%knots2, space%degree, g2)
                corner = layout%vertices(:, layout%corners(:, k))
                side = layout%arc_side(k)
                do j = 1, n2
                    do i = 1, n1
                        index = blk%offset + i + n1*(j - 1)
                        select case (side)
                        case (0)
                            p = bilinear(corner, g1(i), g2(j))
                            tau = 1.0_dp
                        case (1, 3)
                            call layout_side_point(layout, k, side, g1(i), p, tangent)
                            call layout_side_point(layout, k, 4 - side, g1(i), q, &
                                tangent)
                            tau = merge(g2(j), 1.0_dp - g2(j), side == 1)
                            p = (1.0_dp - tau)*p + tau*q
                        case default
                            call layout_side_point(layout, k, side, g2(j), p, tangent)
                            call layout_side_point(layout, k, 6 - side, g2(j), q, &
                                tangent)
                            tau = merge(g1(i), 1.0_dp - g1(i), side == 4)
                            p = (1.0_dp - tau)*p + tau*q
                        end select
                        if (tau == 0.0_dp) then
                            ! filled by wall_spline_data below
                            on_wall(index) = .true.
                            wall_angle(index) = atan2(p(2), p(1))
                            z = (0.0_dp, 0.0_dp)
                        else
                            call harmonic_disc_inverse(map, cmplx(p(1), p(2), dp), z, &
                                info)
                            if (info /= 0) then
                                status = 2
                                return
                            end if
                        end if
                        data(:, index) = [real(z, dp), aimag(z)]
                    end do
                end do
                deallocate(g1, g2)
            end associate
        end do
        call wall_spline_data(space, map, wall_angle, on_wall, data, normal, status)
        if (status /= 0) return
        ! secant iteration on the constant normal offset of the wall data
        offset = 0.0_dp
        previous_offset = 0.0_dp
        previous_error = 0.0_dp
        do pass = 1, 12
            call interpolate(space, data + offset*normal, status)
            if (status /= 0) return
            if (.not. keep_area) exit
            call harmonic_disc_wall_error(space, map, distance, area_error, info)
            if (info /= 0 .or. abs(area_error) < 1.0e-13_dp) exit
            if (pass == 1) then
                step = -area_error*map%area/map%length
            else
                if (area_error == previous_error) exit
                step = -area_error*(offset - previous_offset)/ &
                    (area_error - previous_error)
            end if
            previous_offset = offset
            previous_error = area_error
            offset = offset + step
        end do
    end subroutine harmonic_disc_pull_back

    subroutine wall_spline_data(space, map, wall_angle, on_wall, data, normal, status)
        !! Wall data of the block interpolation from one periodic spline
        !! around the whole wall: its breaks are the cell breaks of the wall
        !! sides of all blocks (as disc angles), it interpolates the curve at
        !! its Greville angles, and data = its values at the wall Greville
        !! points of the blocks. On a block the periodic spline lies in the
        !! spline space of that block side, so the block interpolation
        !! reproduces it and the mesh wall is that one spline: as smooth at
        !! the block corners as inside a block (C^(p-1)). Interpolating the
        !! curve per block instead leaves a kink of the wall at every block
        !! corner. normal: the same spline of the curve normals at the
        !! Greville angles, so that data + c normal is the spline of the
        !! curve moved by c along its normal.
        !! Status 3: singular interpolation.
        type(spline_space_t), intent(in) :: space
        type(harmonic_disc_t), intent(in) :: map
        real(dp), intent(in) :: wall_angle(:)
        logical, intent(in) :: on_wall(:)
        real(dp), intent(inout) :: data(:, :), normal(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: angles(:), breaks(:), knots(:), a(:, :), rhs(:, :)
        real(dp), allocatable :: values(:), derivatives(:)
        integer, allocatable :: pivot(:)
        real(dp) :: point(2), tangent(2), site, t, x, tolerance
        complex(dp) :: z, dy, ddy
        integer :: k, side, e, i, j, b, n, p, count

        p = space%degree
        allocate(angles(0))
        do k = 1, space%block_count
            side = space%layout%arc_side(k)
            if (side == 0) cycle
            if (side == 1 .or. side == 3) then
                breaks = space%blocks(k)%breaks1
            else
                breaks = space%blocks(k)%breaks2
            end if
            do e = 1, size(breaks)
                call layout_side_point(space%layout, k, side, breaks(e), point, tangent)
                angles = [angles, modulo(atan2(point(2), point(1)), two_pi)]
            end do
        end do
        ! sorted distinct break angles
        do i = 2, size(angles)
            x = angles(i)
            j = i - 1
            do while (j >= 1)
                if (angles(j) <= x) exit
                angles(j + 1) = angles(j)
                j = j - 1
            end do
            angles(j + 1) = x
        end do
        tolerance = 1.0e-10_dp
        count = 1
        do i = 2, size(angles)
            if (angles(i) - angles(count) > tolerance .and. &
                angles(1) + two_pi - angles(i) > tolerance) then
                count = count + 1
                angles(count) = angles(i)
            end if
        end do
        n = count
        status = 3
        if (n < p + 1) return
        ! periodic knots k_j = angle_(j mod n) + 2 pi floor(j/n), j = -p .. n + p
        allocate(knots(n + 2*p + 1))
        do j = -p, n + p
            knots(j + p + 1) = angles(modulo(j, n) + 1) + &
                two_pi*real(floor(real(j, dp)/real(n, dp)), dp)
        end do
        allocate(a(n, n), rhs(n, 4), pivot(n))
        a = 0.0_dp
        do i = 1, n
            ! Greville angle of the periodic basis function i
            site = sum(knots(i + p + 1:i + 2*p))/real(max(p, 1), dp)
            x = knots(p + 1) + modulo(site - knots(p + 1), two_pi)
            call evaluate_bspline_basis(knots, p, x, values, derivatives, status)
            if (status /= 0) then
                status = 3
                return
            end if
            do b = 1, n + p
                j = modulo(b - 1, n) + 1
                a(i, j) = a(i, j) + values(b)
            end do
            call harmonic_disc_boundary(map, site, z, t)
            call curve_point(map, t, z, dy, ddy)
            rhs(i, 1:2) = [real(z, dp), aimag(z)]
            rhs(i, 3:4) = [aimag(dy), -real(dy, dp)]/abs(dy)
        end do
        call dgesv(n, 4, a, n, pivot, rhs, n, status)
        if (status /= 0) then
            status = 3
            return
        end if
        do i = 1, size(on_wall)
            if (.not. on_wall(i)) cycle
            x = knots(p + 1) + modulo(wall_angle(i) - knots(p + 1), two_pi)
            call evaluate_bspline_basis(knots, p, x, values, derivatives, status)
            if (status /= 0) then
                status = 3
                return
            end if
            data(:, i) = 0.0_dp
            normal(:, i) = 0.0_dp
            do b = 1, n + p
                j = modulo(b - 1, n) + 1
                data(:, i) = data(:, i) + values(b)*rhs(j, 1:2)
                normal(:, i) = normal(:, i) + values(b)*rhs(j, 3:4)
            end do
        end do
        status = 0
    end subroutine wall_spline_data

    subroutine interpolate(space, data, status)
        !! Block control points from the values at the tensor Greville points.
        type(spline_space_t), intent(inout) :: space
        real(dp), intent(in) :: data(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: m1(:, :), m2(:, :), g1(:), g2(:), work(:, :)
        real(dp), allocatable :: values(:), derivatives(:), a(:, :)
        integer, allocatable :: pivot(:)
        integer :: k, i, j, d, n1, n2

        do k = 1, space%block_count
            associate (blk => space%blocks(k))
                n1 = blk%n(1)
                n2 = blk%n(2)
                allocate(g1(n1), g2(n2), m1(n1, n1), m2(n2, n2), pivot(max(n1, n2)))
                call spline_greville(blk%knots1, space%degree, g1)
                call spline_greville(blk%knots2, space%degree, g2)
                do i = 1, n1
                    call evaluate_bspline_basis(blk%knots1, space%degree, g1(i), &
                        values, derivatives, status)
                    if (status /= 0) return
                    m1(i, :) = values
                end do
                do j = 1, n2
                    call evaluate_bspline_basis(blk%knots2, space%degree, g2(j), &
                        values, derivatives, status)
                    if (status /= 0) return
                    m2(j, :) = values
                end do
                do d = 1, 2
                    allocate(work(n1, n2))
                    do j = 1, n2
                        work(:, j) = data(d, blk%offset + 1 + n1*(j - 1): &
                            blk%offset + n1*j)
                    end do
                    a = m1
                    call dgesv(n1, n2, a, n1, pivot, work, n1, status)
                    if (status == 0) then
                        work = transpose(work)
                        a = m2
                        call dgesv(n2, n1, a, n2, pivot, work, n2, status)
                    end if
                    if (status /= 0) then
                        status = 3
                        return
                    end if
                    work = transpose(work)
                    do j = 1, n2
                        space%geometry(d, blk%offset + 1 + n1*(j - 1): &
                            blk%offset + n1*j) = work(:, j)
                    end do
                    deallocate(work)
                end do
                deallocate(g1, g2, m1, m2, pivot)
            end associate
        end do
        status = 0
    end subroutine interpolate

    pure function bilinear(corner, u, v) result(x)
        real(dp), intent(in) :: corner(2, 4), u, v
        real(dp) :: x(2)

        x = (1.0_dp - u)*(1.0_dp - v)*corner(:, 1) + u*(1.0_dp - v)*corner(:, 2) + &
            u*v*corner(:, 3) + (1.0_dp - u)*v*corner(:, 4)
    end function bilinear

    pure function side_xi(side, t) result(xi)
        integer, intent(in) :: side
        real(dp), intent(in) :: t
        real(dp) :: xi(2)

        select case (side)
        case (1)
            xi = [t, 0.0_dp]
        case (2)
            xi = [1.0_dp, t]
        case (3)
            xi = [t, 1.0_dp]
        case default
            xi = [0.0_dp, t]
        end select
    end function side_xi

    subroutine harmonic_disc_wall_error(space, map, distance, area_error, status)
        !! Of the LAYOUT_WALL sides of the mesh against the curve of map:
        !! the largest distance of the mesh wall from the curve (eight points
        !! per cell side) and the relative error of the enclosed area,
        !! (mesh - curve)/curve. Status 1: no wall side or a failed evaluation.
        type(spline_space_t), intent(in) :: space
        type(harmonic_disc_t), intent(in) :: map
        real(dp), intent(out) :: distance, area_error
        integer, intent(out) :: status

        integer, parameter :: points = 8
        real(dp), parameter :: gauss(4) = [0.0694318442029737_dp, &
            0.3300094782075719_dp, 0.6699905217924281_dp, 0.9305681557970263_dp], &
            weight(4) = [0.1739274225687269_dp, 0.3260725774312731_dp, &
            0.3260725774312731_dp, 0.1739274225687269_dp]
        real(dp), allocatable :: breaks(:)
        real(dp) :: jet(2, 6), xi(2), t, area, d(2), best, s, step
        complex(dp) :: y, dy, ddy, p
        integer :: k, side, e, q, direction, it, node, info
        logical :: found

        status = 1
        distance = 0.0_dp
        area = 0.0_dp
        found = .false.
        do k = 1, space%block_count
            do side = 1, 4
                if (space%layout%side_tag(side, k) /= LAYOUT_WALL) cycle
                found = .true.
                direction = merge(1, 2, side == 1 .or. side == 3)
                if (direction == 1) then
                    breaks = space%blocks(k)%breaks1
                else
                    breaks = space%blocks(k)%breaks2
                end if
                do e = 1, size(breaks) - 1
                    do q = 1, 4
                        t = breaks(e) + gauss(q)*(breaks(e + 1) - breaks(e))
                        xi = side_xi(side, t)
                        call spline_geometry_at(space, k, xi, jet, info)
                        if (info /= 0) return
                        d = jet(:, direction + 1)
                        area = area + weight(q)*(breaks(e + 1) - breaks(e))*0.5_dp* &
                            (jet(1, 1)*d(2) - jet(2, 1)*d(1))
                    end do
                    do q = 1, points
                        t = breaks(e) + (real(q, dp) - 0.5_dp)/real(points, dp)* &
                            (breaks(e + 1) - breaks(e))
                        xi = side_xi(side, t)
                        call spline_geometry_at(space, k, xi, jet, info)
                        if (info /= 0) return
                        p = cmplx(jet(1, 1), jet(2, 1), dp)
                        best = huge(1.0_dp)
                        node = 1
                        do it = 1, map%n
                            if (abs(map%y(it) - p) < best) then
                                best = abs(map%y(it) - p)
                                node = it
                            end if
                        end do
                        s = two_pi*real(node - 1, dp)/real(map%n, dp)
                        do it = 1, 20
                            call curve_point(map, s, y, dy, ddy)
                            step = real(conjg(dy)*(y - p), dp)/ &
                                (abs(dy)**2 + real(conjg(ddy)*(y - p), dp))
                            s = s - step
                            if (abs(step) < 1.0e-14_dp) exit
                        end do
                        call curve_point(map, s, y, dy, ddy)
                        distance = max(distance, abs(y - p))
                    end do
                end do
            end do
        end do
        if (.not. found) return
        area_error = (abs(area) - map%area)/map%area
        status = 0
    end subroutine harmonic_disc_wall_error

    subroutine spline_jacobian_certificate(space, certified, ratio, status, failed)
        !! Sufficient condition for det(dX/dxi) > 0 everywhere, not only at
        !! quadrature points: on every cell the determinant is a polynomial of
        !! degree 2p - 1 per direction, and it is positive when all its
        !! Bernstein coefficients are. certified: all cells pass. ratio: the
        !! smallest coefficient over the largest (a lower bound of min/max of
        !! the determinant when certified). failed: cells that do not pass
        !! (they may still be valid; the condition is not necessary).
        !! Non-periodic blocks of degree 1 to 5 only. Status 1: unsupported
        !! space or a failed evaluation.
        type(spline_space_t), intent(in) :: space
        logical, intent(out) :: certified
        real(dp), intent(out) :: ratio
        integer, intent(out) :: status
        integer, intent(out), optional :: failed

        real(dp), allocatable :: nodes(:), inverse(:, :), f(:, :)
        real(dp) :: jet(2, 6), xi(2), low, high, cell_low
        integer :: m, i, j, k, e1, e2, info, bad

        status = 1
        certified = .false.
        ratio = 0.0_dp
        if (any(space%layout%periodic)) return
        call spline_bernstein_interpolation(space%degree, nodes, inverse, info)
        if (info /= 0) return
        m = size(nodes)
        allocate(f(m, m))
        low = huge(1.0_dp)
        high = -huge(1.0_dp)
        bad = 0
        do k = 1, space%block_count
            associate (b1 => space%blocks(k)%breaks1, b2 => space%blocks(k)%breaks2)
                do e2 = 1, size(b2) - 1
                    do e1 = 1, size(b1) - 1
                        do j = 1, m
                            do i = 1, m
                                xi = [b1(e1) + nodes(i)*(b1(e1 + 1) - b1(e1)), &
                                    b2(e2) + nodes(j)*(b2(e2 + 1) - b2(e2))]
                                call spline_geometry_at(space, k, xi, jet, info)
                                if (info /= 0) return
                                call spline_jacobian_determinant(jet(1, 2), jet(2, 2), &
                                    jet(1, 3), jet(2, 3), f(i, j))
                            end do
                        end do
                        f = matmul(inverse, matmul(f, transpose(inverse)))
                        cell_low = minval(f)
                        if (cell_low <= 0.0_dp) bad = bad + 1
                        low = min(low, cell_low)
                        high = max(high, maxval(f))
                    end do
                end do
            end associate
        end do
        certified = bad == 0
        if (high > 0.0_dp) ratio = low/high
        if (present(failed)) failed = bad
        status = 0
    end subroutine spline_jacobian_certificate

end module fortfem_harmonic_disc
