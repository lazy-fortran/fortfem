program test_toroidal_spectral_trace_ad
    use check, only: check_condition, check_summary
    use fortfem_fourier, only: &
        evaluate_toroidal_harmonic_p, evaluate_toroidal_spectral_trace, &
        evaluate_toroidal_spectral_trace_jvp, evaluate_toroidal_spectral_trace_vjp, &
        toroidal_poisson_exterior_dtn_p
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: scale = 1.7_dp, eta = 1.1_dp
    real(dp), parameter :: theta = 0.37_dp, phi = -0.61_dp
    real(dp), parameter :: step = 1.0e-5_dp
    integer :: degree_indices(3), orders(3), status
    complex(dp) :: coefficients(3), coefficients_dot(3)
    complex(dp), allocatable :: coefficients_bar(:)
    complex(dp) :: value, normal_derivative, value_dot, normal_derivative_dot
    complex(dp) :: value_plus, value_minus, normal_plus, normal_minus
    complex(dp) :: value_bar, normal_bar
    real(dp) :: scale_dot, eta_dot, theta_dot, phi_dot
    real(dp) :: scale_bar, eta_bar, theta_bar, phi_bar
    real(dp) :: value_oracle, flux_oracle, dtn_oracle
    real(dp) :: error, adjoint_error, lhs, rhs
    logical :: all_passed

    all_passed = .true.
    degree_indices = [0, 2, 3]
    orders = [0, 1, 2]
    coefficients = [ &
        cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.2_dp, 0.6_dp, dp)]
    coefficients_dot = [ &
        cmplx(0.03_dp, 0.01_dp, dp), cmplx(-0.02_dp, 0.04_dp, dp), &
        cmplx(0.01_dp, -0.03_dp, dp)]
    scale_dot = 0.017_dp
    eta_dot = -0.013_dp
    theta_dot = 0.021_dp
    phi_dot = -0.019_dp

    call evaluate_toroidal_spectral_trace( &
        degree_indices, orders, coefficients, scale, eta, theta, phi, .false., &
        value, normal_derivative, status)
    call record_condition(status == 0, &
        "toroidal spectral trace accepts a P-branch modal expansion")
    call evaluate_toroidal_harmonic_p( &
        degree_indices(1), orders(1), eta, theta, phi, value_oracle, status)
    call toroidal_poisson_exterior_dtn_p( &
        degree_indices(1), orders(1), scale, eta, theta, phi, value_oracle, &
        flux_oracle, dtn_oracle, status)
    call record_condition(status == 0, &
        "toroidal harmonic oracle evaluates the reference mode")
    call evaluate_toroidal_spectral_trace( &
        degree_indices(1:1), orders(1:1), [cmplx(1.0_dp, 0.0_dp, dp)], scale, eta, &
        theta, phi, .false., value, normal_derivative, status)
    call record_condition(status == 0 .and. abs(real(value, dp) - value_oracle) < &
        2.0e-13_dp .and. abs(real(normal_derivative, dp) - flux_oracle) < &
        2.0e-13_dp, "P-branch trace matches the analytical toroidal mode")

    call evaluate_toroidal_spectral_trace_jvp( &
        degree_indices, orders, coefficients, scale, eta, theta, phi, .false., &
        coefficients_dot, scale_dot, eta_dot, theta_dot, phi_dot, value_dot, &
        normal_derivative_dot, status)
    call evaluate_toroidal_spectral_trace( &
        degree_indices, orders, coefficients + step*coefficients_dot, &
        scale + step*scale_dot, eta + step*eta_dot, theta + step*theta_dot, &
        phi + step*phi_dot, .false., value_plus, normal_plus, status)
    call evaluate_toroidal_spectral_trace( &
        degree_indices, orders, coefficients - step*coefficients_dot, &
        scale - step*scale_dot, eta - step*eta_dot, theta - step*theta_dot, &
        phi - step*phi_dot, .false., value_minus, normal_minus, status)
    error = max(abs(value_dot - (value_plus - value_minus)/(2.0_dp*step)), &
        abs(normal_derivative_dot - &
        (normal_plus - normal_minus)/(2.0_dp*step)))
    call record_condition(status == 0 .and. error < 2.0e-8_dp, &
        "toroidal spectral trace JVP matches central reassembly")

    value_bar = cmplx(0.8_dp, -0.3_dp, dp)
    normal_bar = cmplx(-0.2_dp, 0.6_dp, dp)
    call evaluate_toroidal_spectral_trace_vjp( &
        degree_indices, orders, coefficients, scale, eta, theta, phi, .false., &
        value_bar, normal_bar, coefficients_bar, scale_bar, eta_bar, theta_bar, &
        phi_bar, status)
    lhs = real(conjg(value_bar)*value_dot + &
        conjg(normal_bar)*normal_derivative_dot, dp)
    rhs = real(sum(conjg(coefficients_bar)*coefficients_dot), dp) + &
        scale_bar*scale_dot + eta_bar*eta_dot + theta_bar*theta_dot + &
        phi_bar*phi_dot
    adjoint_error = abs(lhs - rhs)
    call record_condition(status == 0 .and. adjoint_error < &
        5.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "toroidal spectral trace VJP satisfies the real complex adjoint")

    call evaluate_toroidal_spectral_trace( &
        degree_indices(2:), orders(2:), coefficients(2:), scale, eta, theta, phi, .true., &
        value, normal_derivative, status)
    call record_condition(status == 0 .and. abs(value) > tiny(1.0_dp), &
        "Q-branch toroidal spectral trace is available for exterior data")

    call check_summary("toroidal spectral trace")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_toroidal_spectral_trace_ad
