program test_spline_extraction_1d
    !! Bezier extraction reproduces the Cox-de Boor B-splines of
    !! evaluate_bspline_basis (values and first derivatives) on knot vectors
    !! with graded breakpoints and interior multiplicities 1..p.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_bspline_feec, only: evaluate_bspline_basis
    use fortfem_interval_bernstein_jets, only: interval_bernstein_jet
    use fortfem_spline_extraction_1d, only: build_graded_breakpoints, &
        build_open_knot_vector, build_bezier_extraction_1d, &
        SPLINE_GRID_UNIFORM, SPLINE_GRID_SQRT_S, SPLINE_GRID_BUMP_EDGE
    implicit none

    real(dp), allocatable :: breakpoints(:), knots(:), extraction(:, :, :)
    real(dp), allocatable :: element_breaks(:)
    integer, allocatable :: first_basis(:), multiplicity(:)
    integer :: degree, status, grading, n
    real(dp) :: error
    logical :: ok

    ok = .true.
    do degree = 2, 5
        do grading = SPLINE_GRID_UNIFORM, SPLINE_GRID_BUMP_EDGE, &
                SPLINE_GRID_BUMP_EDGE - SPLINE_GRID_UNIFORM
            call build_graded_breakpoints(7, grading, breakpoints, status)
            ok = ok .and. status == 0
            n = size(breakpoints)
            allocate(multiplicity(n))
            multiplicity = 1
            multiplicity(3) = 2
            multiplicity(5) = degree
            call build_open_knot_vector(degree, breakpoints, multiplicity, &
                knots, status)
            ok = ok .and. status == 0
            call build_bezier_extraction_1d(knots, degree, extraction, &
                first_basis, element_breaks, status)
            ok = ok .and. status == 0 .and. size(extraction, 3) == n - 1
            ok = ok .and. maxval(abs(element_breaks - breakpoints)) < 1.0e-14_dp
            error = extraction_error(knots, degree, extraction, first_basis, &
                element_breaks)
            ok = ok .and. error < 1.0e-11_dp
            if (error >= 1.0e-11_dp) print *, "degree", degree, "error", error
            deallocate(multiplicity)
        end do
    end do
    call check_condition(ok, "Bezier extraction reproduces Cox-de Boor B-splines")

    call build_graded_breakpoints(4, SPLINE_GRID_SQRT_S, breakpoints, status)
    call check_condition(status == 0 .and. &
        all(breakpoints(2:) > breakpoints(:4)) .and. &
        breakpoints(5) - breakpoints(4) < breakpoints(2) - breakpoints(1), &
        "sqrt(s) grading refines towards s = 1")

    allocate(multiplicity(5))
    multiplicity = 1
    multiplicity(3) = 4
    call build_open_knot_vector(3, breakpoints, multiplicity, knots, status)
    call check_condition(status /= 0, &
        "interior multiplicity above the degree is rejected")
    call check_condition(jets_match_formula(), &
        "generated Bernstein jets match the closed form for p = 0..8")
    call check_summary("spline_extraction_1d")

contains

    real(dp) function extraction_error(knots, degree, extraction, first_basis, &
            breaks) result(error)
        real(dp), intent(in) :: knots(:), extraction(:, :, :), breaks(:)
        integer, intent(in) :: degree, first_basis(:)

        real(dp), allocatable :: values(:), derivatives(:)
        real(dp) :: bern(degree + 1), dbern(degree + 1), h, t, x
        integer :: e, q, a, status

        error = 0.0_dp
        do e = 1, size(extraction, 3)
            h = breaks(e + 1) - breaks(e)
            do q = 1, 6
                t = (real(q, dp) - 0.37_dp)/6.0_dp
                x = breaks(e) + h*t
                call evaluate_bspline_basis(knots, degree, x, values, &
                    derivatives, status)
                if (status /= 0) then
                    error = huge(1.0_dp)
                    return
                end if
                call bernstein(degree, t, bern, dbern)
                do a = 1, degree + 1
                    error = max(error, abs(values(first_basis(e) + a - 1) - &
                        dot_product(extraction(a, :, e), bern)))
                    error = max(error, h*abs(derivatives(first_basis(e) + a - 1) &
                        - dot_product(extraction(a, :, e), dbern)/h))
                end do
                values(first_basis(e):first_basis(e) + degree) = 0.0_dp
                error = max(error, maxval(abs(values)))
            end do
        end do
    end function extraction_error

    logical function jets_match_formula() result(ok)
        real(dp), allocatable :: value(:), first(:), second(:)
        real(dp), allocatable :: ref(:), dref(:), plus(:), minus(:), dummy(:)
        real(dp), parameter :: delta = 1.0e-5_dp
        real(dp) :: t
        integer :: degree, q, status

        ok = .true.
        do degree = 0, 8
            allocate(value(degree + 1), first(degree + 1), second(degree + 1))
            allocate(ref(degree + 1), dref(degree + 1), plus(degree + 1))
            allocate(minus(degree + 1), dummy(degree + 1))
            do q = 1, 5
                t = (real(q, dp) - 0.4_dp)/5.0_dp
                call interval_bernstein_jet(degree, t, value, first, second, &
                    status)
                call bernstein(degree, t, ref, dref)
                call bernstein(degree, t + delta, dummy, plus)
                call bernstein(degree, t - delta, dummy, minus)
                ok = ok .and. status == 0
                ok = ok .and. maxval(abs(value - ref)) < 1.0e-13_dp
                ok = ok .and. maxval(abs(first - dref)) < 1.0e-12_dp
                ok = ok .and. maxval(abs(second - (plus - minus)/(2*delta))) &
                    < 1.0e-5_dp*real(degree + 1, dp)**2
            end do
            deallocate(value, first, second, ref, dref, plus, minus, dummy)
        end do
    end function jets_match_formula

    pure subroutine bernstein(degree, t, value, first)
        integer, intent(in) :: degree
        real(dp), intent(in) :: t
        real(dp), intent(out) :: value(:), first(:)

        real(dp) :: binom
        integer :: i, j

        do i = 0, degree
            binom = 1.0_dp
            do j = 1, i
                binom = binom*real(degree - i + j, dp)/real(j, dp)
            end do
            value(i + 1) = binom*t**i*(1.0_dp - t)**(degree - i)
            first(i + 1) = binom*(real(i, dp)*t**max(i - 1, 0)* &
                (1.0_dp - t)**(degree - i) - real(degree - i, dp)*t**i* &
                (1.0_dp - t)**max(degree - i - 1, 0))
        end do
    end subroutine bernstein

end program test_spline_extraction_1d
