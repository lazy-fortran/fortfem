module fortfem_spline_extraction_1d
    !! Univariate spline data for multipatch spline spaces: graded
    !! breakpoints, open knot vectors with per-breakpoint multiplicity, and
    !! the Bezier extraction of a knot vector.
    !!
    !! On element e with breakpoints [x_e, x_{e+1}] and t = (x - x_e)/h_e the
    !! p + 1 nonzero B-splines are N_{f_e + a - 1} = sum_b C_e(a, b) B_b^p(t),
    !! where B_b^p are the Bernstein polynomials on [0, 1] and f_e is the
    !! first nonzero B-spline. The extraction follows Algorithm 1 of Borden,
    !! Scott, Evans, Hughes, IJNME 87 (2011) 15-47.
    use fortfem_kinds, only: dp
    implicit none
    private

    integer, parameter, public :: SPLINE_GRID_UNIFORM = 0
    integer, parameter, public :: SPLINE_GRID_SQRT_S = 1
    integer, parameter, public :: SPLINE_GRID_S2 = 2
    integer, parameter, public :: SPLINE_GRID_BUMP = 3
    integer, parameter, public :: SPLINE_GRID_BUMP_EDGE = 4

    public :: build_graded_breakpoints
    public :: build_open_knot_vector
    public :: build_bezier_extraction_1d
    public :: build_periodic_extraction_1d

contains

    subroutine build_graded_breakpoints(element_count, grading, breakpoints, status)
        !! Breakpoints on [0, 1]. The gradings are the radial grid types of
        !! GVEC (src/grid/sgrid.F90, sGrid_init, MIT licence): SQRT_S is finer
        !! at s = 1, S2 finer at s = 0, BUMP finer at both ends, BUMP_EDGE
        !! nearly uniform at s = 0 and finer towards s = 1.
        integer, intent(in) :: element_count, grading
        real(dp), allocatable, intent(out) :: breakpoints(:)
        integer, intent(out) :: status

        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: s(element_count + 1), scale
        integer :: i

        status = 1
        if (element_count < 1) return
        s = [(real(i, dp)/real(element_count, dp), i = 0, element_count)]
        select case (grading)
        case (SPLINE_GRID_UNIFORM)
        case (SPLINE_GRID_SQRT_S)
            s = sqrt(s)
        case (SPLINE_GRID_S2)
            s = s*s
        case (SPLINE_GRID_BUMP)
            s = s - 0.05_dp*sin(2.0_dp*pi*s)
        case (SPLINE_GRID_BUMP_EDGE)
            scale = 1.0_dp - 0.75_dp*((1.0_dp - 0.4_dp)**3 + 0.4_dp**3)
            s = (s - 0.75_dp*((s - 0.4_dp)**3 + 0.4_dp**3))/scale
        case default
            return
        end select
        s(1) = 0.0_dp
        s(element_count + 1) = 1.0_dp
        breakpoints = s
        status = 0
    end subroutine build_graded_breakpoints

    subroutine build_open_knot_vector(degree, breakpoints, multiplicity, knots, status)
        !! Open knot vector with end multiplicity degree + 1 and interior
        !! multiplicity m_k, so the spline is C^(degree - m_k) at breakpoint k.
        integer, intent(in) :: degree
        real(dp), intent(in) :: breakpoints(:)
        integer, intent(in) :: multiplicity(:)
        real(dp), allocatable, intent(out) :: knots(:)
        integer, intent(out) :: status

        integer :: k, n

        status = 1
        n = size(breakpoints)
        if (degree < 0 .or. n < 2 .or. size(multiplicity) /= n) return
        if (any(breakpoints(2:) <= breakpoints(:n - 1))) return
        if (n > 2) then
            if (any(multiplicity(2:n - 1) < 1) .or. &
                any(multiplicity(2:n - 1) > degree)) return
        end if
        allocate(knots(0))
        knots = [knots, spread(breakpoints(1), 1, degree + 1)]
        do k = 2, n - 1
            knots = [knots, spread(breakpoints(k), 1, multiplicity(k))]
        end do
        knots = [knots, spread(breakpoints(n), 1, degree + 1)]
        status = 0
    end subroutine build_open_knot_vector

    subroutine build_bezier_extraction_1d( &
            knots, degree, extraction, first_basis, breakpoints, status)
        !! extraction(a, b, e) = C_e(a, b); first_basis(e) = f_e. The knot
        !! vector must be open (end multiplicity degree + 1) with interior
        !! multiplicities at most degree.
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree
        real(dp), allocatable, intent(out) :: extraction(:, :, :)
        integer, allocatable, intent(out) :: first_basis(:)
        real(dp), allocatable, intent(out) :: breakpoints(:)
        integer, intent(out) :: status

        real(dp), allocatable :: c(:, :, :)
        real(dp) :: alphas(max(degree, 1)), alpha, numer
        integer :: m, p, a, b, i, j, k, r, s, save, mult, nb, element_count
        integer :: ident

        status = 1
        p = degree
        m = size(knots)
        if (.not. valid_open_knots(knots, p)) return
        element_count = count(knots(2:) > knots(:m - 1))
        allocate(c(p + 1, p + 1, element_count + 1), first_basis(element_count))
        allocate(breakpoints(element_count + 1))
        c = 0.0_dp
        do ident = 1, p + 1
            c(ident, ident, :) = 1.0_dp
        end do
        a = p + 1
        b = a + 1
        nb = 1
        breakpoints(1) = knots(a)
        do while (b < m)
            first_basis(nb) = a - p
            i = b
            do while (b < m)
                if (knots(b + 1) /= knots(b)) exit
                b = b + 1
            end do
            mult = b - i + 1
            if (mult < p) then
                numer = knots(b) - knots(a)
                do j = p, mult + 1, -1
                    alphas(j - mult) = numer/(knots(a + j) - knots(a))
                end do
                r = p - mult
                do j = 1, r
                    save = r - j + 1
                    s = mult + j
                    do k = p + 1, s + 1, -1
                        alpha = alphas(k - s)
                        c(:, k, nb) = alpha*c(:, k, nb) + &
                            (1.0_dp - alpha)*c(:, k - 1, nb)
                    end do
                    if (b < m) c(save:save + j, save, nb + 1) = &
                        c(p - j + 1:p + 1, p + 1, nb)
                end do
            end if
            breakpoints(nb + 1) = knots(b)
            nb = nb + 1
            if (b < m) then
                a = b
                b = b + 1
            end if
        end do
        if (nb - 1 /= element_count) return
        extraction = c(:, :, :element_count)
        status = 0
    end subroutine build_bezier_extraction_1d

    pure logical function valid_open_knots(knots, degree) result(valid)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree

        integer :: m, k, run

        valid = .false.
        m = size(knots)
        if (degree < 1 .or. m < 2*degree + 2) return
        if (any(knots(2:) < knots(:m - 1))) return
        if (any(knots(2:degree + 1) /= knots(1))) return
        if (any(knots(m - degree:m - 1) /= knots(m))) return
        if (knots(degree + 2) == knots(1)) return
        if (knots(m - degree - 1) == knots(m)) return
        run = 1
        do k = degree + 3, m - degree - 1
            if (knots(k) == knots(k - 1)) then
                run = run + 1
                if (run > degree) return
            else
                run = 1
            end if
        end do
        valid = .true.
    end function valid_open_knots

    subroutine build_periodic_extraction_1d(degree, breakpoints, knots, extraction, &
            first_basis, status)
        !! Periodic splines of maximal smoothness on the breakpoints
        !! 0 = x_0 < ... < x_n = 1 as the n + p unclamped B-splines on the
        !! periodically extended knots x_(i - n) - 1, x_i, x_(i + n) + 1; the
        !! functions i and i + n (i <= p) coincide on [0, 1]. Element e
        !! carries the functions e .. e + p, and extraction(:, :, e) maps the
        !! Bernstein polynomials of degree p on the element to them.
        !! Status 1: fewer than p elements or invalid breakpoints.
        integer, intent(in) :: degree
        real(dp), intent(in) :: breakpoints(:)
        real(dp), allocatable, intent(out) :: knots(:), extraction(:, :, :)
        integer, allocatable, intent(out) :: first_basis(:)
        integer, intent(out) :: status

        real(dp), allocatable :: bernstein(:, :), splines(:, :), value(:), work(:)
        real(dp) :: t, x
        integer :: n, p, e, q, b

        status = 1
        p = degree
        n = size(breakpoints) - 1
        if (p < 1 .or. n < p) return
        if (abs(breakpoints(1)) > 0.0_dp .or. abs(breakpoints(n + 1) - 1.0_dp) > &
            0.0_dp) return
        if (any(breakpoints(2:) <= breakpoints(:n))) return
        knots = [breakpoints(n - p + 1:n) - 1.0_dp, breakpoints, &
            breakpoints(2:p + 1) + 1.0_dp]
        allocate(extraction(p + 1, p + 1, n), first_basis(n), &
            bernstein(p + 1, p + 1), splines(p + 1, p + 1), value(p + 1), &
            work(p + 1))
        do e = 1, n
            first_basis(e) = e
            do q = 1, p + 1
                t = 0.5_dp*(1.0_dp - cos(acos(-1.0_dp)*(q - 0.5_dp)/(p + 1)))
                x = breakpoints(e) + t*(breakpoints(e + 1) - breakpoints(e))
                do b = 0, p
                    bernstein(q, b + 1) = binomial(p, b)*t**b*(1.0_dp - t)**(p - b)
                end do
                call local_bsplines(knots, p, e + p, x, value)
                splines(q, :) = value
            end do
            ! bernstein C^T = splines, column by column
            call solve_dense(bernstein, splines, work, status)
            if (status /= 0) return
            extraction(:, :, e) = transpose(splines)
        end do
        status = 0
    end subroutine build_periodic_extraction_1d

    pure subroutine local_bsplines(knots, p, span, x, value)
        !! The p + 1 B-splines nonzero on [knots(span), knots(span + 1)),
        !! functions span - p .. span (Cox-de Boor).
        real(dp), intent(in) :: knots(:), x
        integer, intent(in) :: p, span
        real(dp), intent(out) :: value(p + 1)

        real(dp) :: left(p), right(p), saved, temp
        integer :: j, r

        value(1) = 1.0_dp
        do j = 1, p
            left(j) = x - knots(span + 1 - j)
            right(j) = knots(span + j) - x
            saved = 0.0_dp
            do r = 0, j - 1
                temp = value(r + 1)/(right(r + 1) + left(j - r))
                value(r + 1) = saved + right(r + 1)*temp
                saved = left(j - r)*temp
            end do
            value(j + 1) = saved
        end do
    end subroutine local_bsplines

    pure real(dp) function binomial(n, k) result(c)
        integer, intent(in) :: n, k
        integer :: i

        c = 1.0_dp
        do i = 1, k
            c = c*real(n - k + i, dp)/real(i, dp)
        end do
    end function binomial

    pure subroutine solve_dense(a, b, work, status)
        !! b <- a^(-1) b by Gaussian elimination with partial pivoting.
        real(dp), intent(inout) :: a(:, :), b(:, :)
        real(dp), intent(inout) :: work(:)
        integer, intent(out) :: status

        real(dp) :: f
        integer :: n, i, k, pivot

        n = size(a, 1)
        status = 1
        do k = 1, n
            pivot = k - 1 + maxloc(abs(a(k:, k)), 1)
            if (abs(a(pivot, k)) < 1.0e-300_dp) return
            if (pivot /= k) then
                work = a(k, :)
                a(k, :) = a(pivot, :)
                a(pivot, :) = work
                work(:size(b, 2)) = b(k, :)
                b(k, :) = b(pivot, :)
                b(pivot, :) = work(:size(b, 2))
            end if
            do i = k + 1, n
                f = a(i, k)/a(k, k)
                a(i, k:) = a(i, k:) - f*a(k, k:)
                b(i, :) = b(i, :) - f*b(k, :)
            end do
        end do
        do k = n, 1, -1
            b(k, :) = (b(k, :) - matmul(a(k, k + 1:), b(k + 1:, :)))/a(k, k)
        end do
        status = 0
    end subroutine solve_dense

end module fortfem_spline_extraction_1d
