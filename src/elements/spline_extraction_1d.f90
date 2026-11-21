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

end module fortfem_spline_extraction_1d
