program test_toroidal_poisson_dtn_ad
    use check, only: check_condition, check_summary
    use fortfem_fourier, only: &
        evaluate_toroidal_ampere_field_p, &
        evaluate_toroidal_ampere_field_p_jvp, &
        evaluate_toroidal_ampere_field_p_vjp, &
        evaluate_toroidal_harmonic_p, &
        evaluate_toroidal_harmonic_p_jvp, &
        evaluate_toroidal_harmonic_p_vjp, &
        toroidal_poisson_exterior_dtn_p, &
        toroidal_poisson_exterior_dtn_p_jvp, &
        toroidal_poisson_exterior_dtn_p_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: degree_index = 2, order = 1
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: scale, eta, theta, phi
    real(dp) :: scale_dot, eta_dot, theta_dot, phi_dot
    real(dp) :: scale_bar, eta_bar, theta_bar, phi_bar
    real(dp) :: value, value_dot, value_plus, value_minus, value_bar
    real(dp) :: normal_derivative, normal_derivative_dot
    real(dp) :: normal_plus, normal_minus, normal_bar
    real(dp) :: dtn_value, dtn_value_dot, dtn_plus, dtn_minus, dtn_bar
    real(dp) :: field(3), field_dot(3), field_plus(3), field_minus(3)
    real(dp) :: field_bar(3), lhs, rhs, relative_error
    integer :: status

    scale = 2.4_dp
    eta = 0.83_dp
    theta = 0.41_dp
    phi = 0.67_dp
    scale_dot = 0.031_dp
    eta_dot = -0.017_dp
    theta_dot = 0.023_dp
    phi_dot = -0.029_dp

    call evaluate_toroidal_harmonic_p( &
        degree_index, order, eta, theta, phi, value, status)
    call evaluate_toroidal_harmonic_p_jvp( &
        degree_index, order, eta, theta, phi, eta_dot, theta_dot, phi_dot, &
        value_dot, status)
    call evaluate_toroidal_harmonic_p( &
        degree_index, order, eta + step*eta_dot, theta + step*theta_dot, &
        phi + step*phi_dot, value_plus, status)
    call evaluate_toroidal_harmonic_p( &
        degree_index, order, eta - step*eta_dot, theta - step*theta_dot, &
        phi - step*phi_dot, value_minus, status)
    relative_error = abs(value_dot - (value_plus - value_minus)/(2.0_dp*step))/ &
        max(1.0_dp, abs(value_dot))
    call check_condition(relative_error < 2.0e-7_dp, &
        "Toroidal harmonic JVP matches re-evaluation")
    value_bar = 0.73_dp
    call evaluate_toroidal_harmonic_p_vjp( &
        degree_index, order, eta, theta, phi, value_bar, eta_bar, theta_bar, &
        phi_bar, status)
    lhs = value_bar*value_dot
    rhs = eta_bar*eta_dot + theta_bar*theta_dot + phi_bar*phi_dot
    call check_condition(abs(lhs - rhs) < 2.0e-9_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Toroidal harmonic VJP satisfies the adjoint identity")

    call evaluate_toroidal_ampere_field_p( &
        degree_index, order, scale, eta, theta, phi, field, status)
    call evaluate_toroidal_ampere_field_p_jvp( &
        degree_index, order, scale, eta, theta, phi, scale_dot, eta_dot, &
        theta_dot, phi_dot, field_dot, status)
    call evaluate_toroidal_ampere_field_p( &
        degree_index, order, scale + step*scale_dot, eta + step*eta_dot, &
        theta + step*theta_dot, phi + step*phi_dot, field_plus, status)
    call evaluate_toroidal_ampere_field_p( &
        degree_index, order, scale - step*scale_dot, eta - step*eta_dot, &
        theta - step*theta_dot, phi - step*phi_dot, field_minus, status)
    relative_error = maxval(abs(field_dot - (field_plus - field_minus)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(field_dot)))
    call check_condition(relative_error < 3.0e-7_dp, &
        "Toroidal Ampere JVP matches re-evaluation")
    field_bar = [0.37_dp, -0.21_dp, 0.48_dp]
    call evaluate_toroidal_ampere_field_p_vjp( &
        degree_index, order, scale, eta, theta, phi, field_bar, scale_bar, &
        eta_bar, theta_bar, phi_bar, status)
    lhs = sum(field_bar*field_dot)
    rhs = scale_bar*scale_dot + eta_bar*eta_dot + theta_bar*theta_dot + &
        phi_bar*phi_dot
    call check_condition(abs(lhs - rhs) < 4.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Toroidal Ampere VJP satisfies the adjoint identity")

    call toroidal_poisson_exterior_dtn_p( &
        degree_index, order, scale, eta, theta, phi, value, &
        normal_derivative, dtn_value, status)
    call toroidal_poisson_exterior_dtn_p_jvp( &
        degree_index, order, scale, eta, theta, phi, scale_dot, eta_dot, &
        theta_dot, phi_dot, value_dot, normal_derivative_dot, dtn_value_dot, &
        status)
    call toroidal_poisson_exterior_dtn_p( &
        degree_index, order, scale + step*scale_dot, eta + step*eta_dot, &
        theta + step*theta_dot, phi + step*phi_dot, value_plus, normal_plus, &
        dtn_plus, status)
    call toroidal_poisson_exterior_dtn_p( &
        degree_index, order, scale - step*scale_dot, eta - step*eta_dot, &
        theta - step*theta_dot, phi - step*phi_dot, value_minus, normal_minus, &
        dtn_minus, status)
    relative_error = max(abs(value_dot - (value_plus - value_minus)/(2.0_dp*step)), &
        abs(normal_derivative_dot - (normal_plus - normal_minus)/(2.0_dp*step)), &
        abs(dtn_value_dot - (dtn_plus - dtn_minus)/(2.0_dp*step)))/ &
        max(1.0_dp, abs(value_dot), abs(normal_derivative_dot), abs(dtn_value_dot))
    call check_condition(relative_error < 5.0e-7_dp, &
        "Toroidal Poisson DtN JVP matches re-evaluation")
    value_bar = 0.29_dp
    normal_bar = -0.17_dp
    dtn_bar = 0.44_dp
    call toroidal_poisson_exterior_dtn_p_vjp( &
        degree_index, order, scale, eta, theta, phi, value_bar, normal_bar, &
        dtn_bar, scale_bar, eta_bar, theta_bar, phi_bar, status)
    lhs = value_bar*value_dot + normal_bar*normal_derivative_dot + &
        dtn_bar*dtn_value_dot
    rhs = scale_bar*scale_dot + eta_bar*eta_dot + theta_bar*theta_dot + &
        phi_bar*phi_dot
    call check_condition(abs(lhs - rhs) < 7.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Toroidal Poisson DtN VJP satisfies the adjoint identity")

    call check_summary("Differentiable toroidal harmonics and Poisson DtN")
end program test_toroidal_poisson_dtn_ad
