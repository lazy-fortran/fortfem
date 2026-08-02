program test_axis_regular_radial_basis
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_axis_regular_radial_basis, &
        evaluate_axis_regular_radial_basis_jvp, &
        evaluate_axis_regular_radial_basis_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: term_count = 3, sample_count = 4
    integer, parameter :: powers(term_count) = [2, 4, 6]
    complex(dp), parameter :: coefficients(term_count) = [ &
        cmplx(1.2_dp, -0.3_dp, dp), cmplx(-0.7_dp, 0.5_dp, dp), &
        cmplx(0.2_dp, 0.9_dp, dp)]
    complex(dp), parameter :: coefficients_dot(term_count) = [ &
        cmplx(-0.2_dp, 0.4_dp, dp), cmplx(0.6_dp, -0.1_dp, dp), &
        cmplx(-0.3_dp, -0.5_dp, dp)]
    complex(dp), parameter :: values_bar(sample_count) = [ &
        cmplx(0.1_dp, -0.4_dp, dp), cmplx(-0.8_dp, 0.3_dp, dp), &
        cmplx(0.5_dp, 0.7_dp, dp), cmplx(-0.2_dp, -0.6_dp, dp)]
    real(dp), parameter :: rho(sample_count) = [0.0_dp, 0.2_dp, 0.6_dp, 1.1_dp]
    real(dp), parameter :: rho_fd(sample_count) = [0.1_dp, 0.25_dp, 0.6_dp, 1.1_dp]
    real(dp), parameter :: rho_dot(sample_count) = &
        [0.03_dp, -0.07_dp, 0.11_dp, -0.04_dp]
    real(dp), parameter :: epsilon = 1.0e-7_dp
    real(dp), parameter :: tolerance = 2.0e-8_dp

    type(fortsparse_status_t) :: status
    complex(dp) :: values(sample_count), expected(sample_count)
    complex(dp) :: values_dot(sample_count), plus(sample_count), minus(sample_count)
    complex(dp) :: coefficients_bar(term_count)
    real(dp) :: rho_bar(sample_count)
    real(dp) :: forward_pairing, reverse_pairing
    real(dp) :: invalid_rho(sample_count)
    complex(dp) :: invalid_coefficients(term_count)
    logical :: all_passed

    all_passed = .true.

    expected = coefficients(1)*rho**2 + coefficients(2)*rho**4 + &
        coefficients(3)*rho**6
    call evaluate_axis_regular_radial_basis( &
        2, powers, coefficients, rho, values, status)
    call record(status%code == 0 .and. maxval(abs(values - expected)) < tolerance, &
        "axis-regular polynomial matches an independent analytical oracle")
    call record(abs(values(1)) < tolerance, &
        "non-axisymmetric radial basis vanishes at the axis")

    call evaluate_axis_regular_radial_basis_jvp( &
        2, powers, coefficients, rho_fd, coefficients_dot, rho_dot, &
        values_dot, status)
    call evaluate_axis_regular_radial_basis( &
        2, powers, coefficients + epsilon*coefficients_dot, &
        rho_fd + epsilon*rho_dot, plus, status)
    call evaluate_axis_regular_radial_basis( &
        2, powers, coefficients - epsilon*coefficients_dot, &
        rho_fd - epsilon*rho_dot, minus, status)
    call record(maxval(abs(values_dot - (plus - minus)/(2.0_dp*epsilon))) < &
        tolerance, "coefficient-and-rho JVP matches central differences")

    call evaluate_axis_regular_radial_basis_vjp( &
        2, powers, coefficients, rho_fd, values_bar, coefficients_bar, &
        rho_bar, status)
    forward_pairing = real(sum(conjg(values_bar)*values_dot), dp)
    reverse_pairing = real(sum(conjg(coefficients_bar)*coefficients_dot), dp) + &
        sum(rho_bar*rho_dot)
    call record(status%code == 0 .and. &
        abs(forward_pairing - reverse_pairing) < tolerance, &
        "radial-basis VJP satisfies the complex real-part adjoint oracle")

    call evaluate_axis_regular_radial_basis( &
        2, [2, 3], coefficients(:2), rho, values, status)
    call record(status%code /= 0 .and. all(values == cmplx(0.0_dp, 0.0_dp, dp)), &
        "basis rejects a radial power with inconsistent parity")

    call evaluate_axis_regular_radial_basis( &
        2, [0, 2], coefficients(:2), rho, values, status)
    call record(status%code /= 0, "basis rejects a power below abs(m)")

    call evaluate_axis_regular_radial_basis( &
        2, [2, 2], coefficients(:2), rho, values, status)
    call record(status%code /= 0, "basis rejects duplicate powers")

    invalid_rho = rho
    invalid_rho(2) = -0.1_dp
    call evaluate_axis_regular_radial_basis( &
        2, powers, coefficients, invalid_rho, values, status)
    call record(status%code /= 0, "basis rejects negative radial samples")

    invalid_coefficients = coefficients
    invalid_coefficients(1) = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
    call evaluate_axis_regular_radial_basis( &
        2, powers, invalid_coefficients, rho, values, status)
    call record(status%code /= 0, "basis rejects non-finite coefficients")

    call evaluate_axis_regular_radial_basis( &
        2, powers, coefficients, rho, values(:sample_count - 1), status)
    call record(status%code /= 0, "basis rejects an incompatible output shape")

    call check_summary("axis-regular radial basis")
    if (.not. all_passed) error stop 1

contains

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_axis_regular_radial_basis
