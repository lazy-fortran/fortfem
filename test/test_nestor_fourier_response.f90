program test_nestor_fourier_response
    !! Independent manufactured oracle for the neutral NESTOR-like modal map.
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        apply_nestor_fourier_response_map, &
        apply_nestor_fourier_response_map_jvp, &
        apply_nestor_fourier_response_map_vjp, &
        evaluate_nestor_fourier_reciprocity_defect
    use fortfem_fourier, only: toroidal_p_derivative, toroidal_p_second_derivative
    use fortsparse, only: fortsparse_status_t, FORTSPARSE_OK
    implicit none

    integer, parameter :: dp = real64, n = 3
    real(dp), parameter :: scale = 2.3_dp, eta = 0.91_dp
    real(dp), parameter :: step = 1.0e-6_dp
    integer :: degree_indices(n), orders(n)
    type(fortsparse_status_t) :: status
    complex(dp) :: coefficients(n), coefficients_dot(n), response(n)
    complex(dp) :: response_dot(n), response_plus(n), response_minus(n)
    complex(dp) :: response_bar(n), coefficients_bar(n)
    complex(dp) :: source_one(n), source_two(n)
    real(dp) :: scale_dot, eta_dot, scale_bar, eta_bar
    real(dp) :: factor, factor_dot, radial, radial_first, radial_second
    complex(dp) :: expected(n), expected_dot(n)
    real(dp) :: defect, expected_defect
    real(dp) :: lhs, rhs, error
    integer :: mode

    degree_indices = [0, 1, 3]
    orders = [0, 2, 1]
    coefficients = [cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.2_dp, 0.6_dp, dp)]
    coefficients_dot = [cmplx(0.03_dp, 0.01_dp, dp), cmplx(-0.02_dp, 0.04_dp, dp), &
        cmplx(0.01_dp, -0.03_dp, dp)]
    scale_dot = 0.017_dp
    eta_dot = -0.013_dp

    call apply_nestor_fourier_response_map( &
        degree_indices, orders, coefficients, scale, eta, .false., response, status)
    do mode = 1, n
        radial_first = toroidal_p_derivative( &
            degree_indices(mode), orders(mode), cosh(eta))
        expected(mode) = -sinh(eta)/scale*radial_first*coefficients(mode)
    end do
    call check_condition(status%code == FORTSPARSE_OK .and. &
        maxval(abs(response-expected)) < 2.0e-13_dp, &
        "NESTOR Fourier response matches the independent radial oracle")

    call apply_nestor_fourier_response_map_jvp( &
        degree_indices, orders, coefficients, scale, eta, .false., coefficients_dot, &
        scale_dot, eta_dot, response_dot, status)
    do mode = 1, n
        radial = toroidal_p_second_derivative( &
            degree_indices(mode), orders(mode), cosh(eta))
        radial_first = toroidal_p_derivative( &
            degree_indices(mode), orders(mode), cosh(eta))
        factor_dot = (sinh(eta)/scale**2*radial_first*scale_dot) - &
            cosh(eta)/scale*radial_first*eta_dot - &
            sinh(eta)**2/scale*radial*eta_dot
        expected_dot(mode) = factor_dot*coefficients(mode) + &
            (-sinh(eta)/scale*radial_first)*coefficients_dot(mode)
    end do
    call apply_nestor_fourier_response_map( &
        degree_indices, orders, coefficients+step*coefficients_dot, &
        scale+step*scale_dot, eta+step*eta_dot, .false., response_plus, status)
    call apply_nestor_fourier_response_map( &
        degree_indices, orders, coefficients-step*coefficients_dot, &
        scale-step*scale_dot, eta-step*eta_dot, .false., response_minus, status)
    error = maxval(abs(response_dot-(response_plus-response_minus)/(2.0_dp*step)))
    call check_condition(status%code == FORTSPARSE_OK .and. &
        maxval(abs(response_dot-expected_dot)) < 2.0e-13_dp .and. &
        error < 2.0e-8_dp, &
        "NESTOR Fourier JVP matches analytic and finite-difference oracles")

    response_bar = [cmplx(0.8_dp, -0.3_dp, dp), cmplx(-0.2_dp, 0.6_dp, dp), &
        cmplx(0.1_dp, 0.4_dp, dp)]
    call apply_nestor_fourier_response_map_vjp( &
        degree_indices, orders, coefficients, scale, eta, .false., response_bar, &
        coefficients_bar, scale_bar, eta_bar, status)
    lhs = real(sum(conjg(response_bar)*response_dot), dp)
    rhs = real(sum(conjg(coefficients_bar)*coefficients_dot), dp) + &
        scale_bar*scale_dot + eta_bar*eta_dot
    call check_condition(status%code == FORTSPARSE_OK .and. &
        abs(lhs-rhs) < 2.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "NESTOR Fourier VJP satisfies the real complex adjoint identity")

    source_one = [cmplx(0.2_dp, 0.1_dp, dp), cmplx(-0.4_dp, 0.3_dp, dp), &
        cmplx(0.7_dp, -0.2_dp, dp)]
    source_two = [cmplx(-0.3_dp, 0.5_dp, dp), cmplx(0.6_dp, -0.2_dp, dp), &
        cmplx(-0.1_dp, 0.4_dp, dp)]
    do mode = 1, n
        radial_first = toroidal_p_derivative( &
            degree_indices(mode), orders(mode), cosh(eta))
        expected(mode) = -sinh(eta)/scale*radial_first
    end do
    expected_defect = abs(real(sum(conjg(source_one)*expected*source_two), dp) - &
        real(sum(conjg(source_two)*expected*source_one), dp))/max(1.0_dp, &
        abs(real(sum(conjg(source_one)*expected*source_two), dp)), &
        abs(real(sum(conjg(source_two)*expected*source_one), dp)))
    call evaluate_nestor_fourier_reciprocity_defect( &
        degree_indices, orders, source_one, source_two, scale, eta, .false., &
        defect, status)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        abs(defect-expected_defect) < 2.0e-14_dp .and. defect < 2.0e-14_dp, &
        "NESTOR Fourier response passes the independent reciprocity oracle")

    call apply_nestor_fourier_response_map( &
        degree_indices(2:), orders(2:), coefficients(2:), scale, eta, .true., &
        response(2:), status)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        all(ieee_is_finite(real(response(2:), dp))) .and. &
        all(ieee_is_finite(aimag(response(2:)))), &
        "NESTOR Fourier Q-branch response remains finite")
    call apply_nestor_fourier_response_map( &
        degree_indices, orders, coefficients, 0.0_dp, eta, .false., response, status)
    call check_condition(status%code /= FORTSPARSE_OK, &
        "NESTOR Fourier response rejects a nonpositive scale")

    call check_summary("NESTOR Fourier response")
end program test_nestor_fourier_response
