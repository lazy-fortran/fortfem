program test_fortfem_toroidal_harmonics
    ! Public-adapter oracle for Hobson toroidal P/Q harmonics.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortfem_fourier, only: toroidal_p, toroidal_q, &
        toroidal_p_derivative, toroidal_q_derivative, &
        toroidal_p_second_derivative, toroidal_q_second_derivative
    implicit none

    integer :: nfail

    nfail = 0
    call check("P_-1/2^0(2)", toroidal_p(0, 0, 2.0_dp), &
        0.90128629936044729874_dp)
    call check("Q_-1/2^0(2)", toroidal_q(0, 0, 2.0_dp), &
        1.6566381702365941664_dp)
    call check("P_3/2^0(2)", toroidal_p(2, 0, 2.0_dp), &
        3.2439396660408049155_dp)
    call check("Q_3/2^0(2)", toroidal_q(2, 0, 2.0_dp), &
        0.045158724151576976637_dp)
    call check("dP_5/2^2(2)", toroidal_p_derivative(3, 2, 2.0_dp), &
        48.525686151003148510_dp)
    call check("dQ_3/2^1(2.5)", toroidal_q_derivative(2, 1, 2.5_dp), &
        0.067949860642696846331_dp)
    call check_ode("P second derivative", .true., 3, 1, 2.2_dp)
    call check_ode("Q second derivative", .false., 3, 1, 2.2_dp)
    call check_true("invalid x", ieee_is_nan(toroidal_p(0, 0, 1.0_dp)))
    call check_true("invalid order", ieee_is_nan(toroidal_q(0, -1, 2.0_dp)))

    if (nfail /= 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) FAILED"
        error stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine check(label, got, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, expected

        if (.not. abs(got - expected) <= &
            3.0e-12_dp*(1.0_dp + abs(expected))) then
            nfail = nfail + 1
            write (error_unit, "(a,2(a,es22.14))") &
                "FAIL: "//label, " got=", got, " expected=", expected
        end if
    end subroutine check

    subroutine check_true(label, condition)
        character(*), intent(in) :: label
        logical, intent(in) :: condition

        if (.not. condition) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL: "//label
        end if
    end subroutine check_true

    subroutine check_ode(label, first_kind, degree_index, order, x)
        character(*), intent(in) :: label
        logical, intent(in) :: first_kind
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: x
        real(dp) :: value, first, second, degree, denominator, residual

        if (first_kind) then
            value = toroidal_p(degree_index, order, x)
            first = toroidal_p_derivative(degree_index, order, x)
            second = toroidal_p_second_derivative(degree_index, order, x)
        else
            value = toroidal_q(degree_index, order, x)
            first = toroidal_q_derivative(degree_index, order, x)
            second = toroidal_q_second_derivative(degree_index, order, x)
        end if
        degree = real(degree_index, dp) - 0.5_dp
        denominator = 1.0_dp - x*x
        residual = denominator*second - 2.0_dp*x*first + &
            (degree*(degree + 1.0_dp) - &
            real(order*order, dp)/denominator)*value
        call check(label, residual, 0.0_dp)
    end subroutine check_ode

end program test_fortfem_toroidal_harmonics
