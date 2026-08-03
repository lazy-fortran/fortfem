program test_toroidal_vector_coordinates_ad
    use check, only: check_condition, check_summary
    use fortfem_core, only: &
        toroidal_vector_to_cartesian, toroidal_vector_to_cartesian_jvp, &
        toroidal_vector_to_cartesian_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: eta, theta, phi
    real(dp) :: eta_dot, theta_dot, phi_dot
    real(dp) :: eta_bar, theta_bar, phi_bar
    real(dp) :: components(3), components_dot(3), components_bar(3)
    real(dp) :: cartesian(3), cartesian_dot(3), cartesian_bar(3)
    real(dp) :: cartesian_plus(3), cartesian_minus(3), cartesian_dot_fd(3)
    real(dp) :: lhs, rhs, relative_error
    integer :: status

    eta = 0.83_dp
    theta = 0.41_dp
    phi = 0.67_dp
    eta_dot = -0.019_dp
    theta_dot = 0.027_dp
    phi_dot = -0.031_dp
    components = [0.42_dp, -0.18_dp, 0.63_dp]
    components_dot = [-0.023_dp, 0.017_dp, 0.029_dp]

    call toroidal_vector_to_cartesian(eta, theta, phi, components, cartesian)
    call toroidal_vector_to_cartesian_jvp( &
        eta, theta, phi, components, eta_dot, theta_dot, phi_dot, &
        components_dot, cartesian_dot, status)
    call toroidal_vector_to_cartesian( &
        eta + step*eta_dot, theta + step*theta_dot, phi + step*phi_dot, &
        components + step*components_dot, cartesian_plus)
    call toroidal_vector_to_cartesian( &
        eta - step*eta_dot, theta - step*theta_dot, phi - step*phi_dot, &
        components - step*components_dot, cartesian_minus)
    cartesian_dot_fd = (cartesian_plus - cartesian_minus)/(2.0_dp*step)
    relative_error = maxval(abs(cartesian_dot - cartesian_dot_fd))/ &
        max(1.0_dp, maxval(abs(cartesian_dot)))
    call check_condition(status == 0 .and. relative_error < 3.0e-8_dp, &
        "Toroidal vector map JVP matches re-evaluation")

    cartesian_bar = [0.29_dp, -0.37_dp, 0.18_dp]
    call toroidal_vector_to_cartesian_vjp( &
        eta, theta, phi, components, cartesian_bar, eta_bar, theta_bar, &
        phi_bar, components_bar, status)
    lhs = sum(cartesian_bar*cartesian_dot)
    rhs = eta_bar*eta_dot + theta_bar*theta_dot + phi_bar*phi_dot + &
        sum(components_bar*components_dot)
    call check_condition(status == 0 .and. &
        abs(lhs - rhs) < 5.0e-9_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Toroidal vector map VJP satisfies the adjoint identity")

    call check_summary("Differentiable toroidal vector coordinate map")
end program test_toroidal_vector_coordinates_ad
