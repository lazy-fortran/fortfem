program test_toroidal_coordinates_ad
    use check, only: check_condition, check_summary
    use fortfem_core, only: &
        cartesian_to_toroidal, cartesian_to_toroidal_jvp, &
        cartesian_to_toroidal_vjp, toroidal_point_to_cartesian, &
        toroidal_point_to_cartesian_jvp, toroidal_point_to_cartesian_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: scale, scale_dot, scale_bar
    real(dp) :: eta, theta, phi, eta_dot, theta_dot, phi_dot
    real(dp) :: eta_bar, theta_bar, phi_bar
    real(dp) :: point(3), point_dot(3), point_bar(3)
    real(dp) :: point_plus(3), point_minus(3), point_dot_fd(3)
    real(dp) :: point_bar_fd(3), scale_bar_fd
    real(dp) :: eta_plus, theta_plus, phi_plus
    real(dp) :: eta_minus, theta_minus, phi_minus
    real(dp) :: lhs, rhs, relative_error
    integer :: status

    scale = 2.7_dp
    eta = 0.91_dp
    theta = 0.37_dp
    phi = 0.62_dp
    scale_dot = 0.023_dp
    eta_dot = -0.017_dp
    theta_dot = 0.031_dp
    phi_dot = -0.029_dp

    call toroidal_point_to_cartesian(scale, eta, theta, phi, point)
    call toroidal_point_to_cartesian_jvp( &
        scale, eta, theta, phi, scale_dot, eta_dot, theta_dot, phi_dot, &
        point_dot, status)
    call toroidal_point_to_cartesian( &
        scale + step*scale_dot, eta + step*eta_dot, &
        theta + step*theta_dot, phi + step*phi_dot, point_plus)
    call toroidal_point_to_cartesian( &
        scale - step*scale_dot, eta - step*eta_dot, &
        theta - step*theta_dot, phi - step*phi_dot, point_minus)
    point_dot_fd = (point_plus - point_minus)/(2.0_dp*step)
    relative_error = maxval(abs(point_dot - point_dot_fd))/ &
        max(1.0_dp, maxval(abs(point_dot)))
    call check_condition(status == 0 .and. relative_error < 2.0e-8_dp, &
        "Toroidal Cartesian map JVP matches re-evaluation")

    point_bar = [0.31_dp, -0.47_dp, 0.28_dp]
    call toroidal_point_to_cartesian_vjp( &
        scale, eta, theta, phi, point_bar, scale_bar, eta_bar, theta_bar, &
        phi_bar, status)
    lhs = sum(point_bar*point_dot)
    rhs = scale_bar*scale_dot + eta_bar*eta_dot + theta_bar*theta_dot + &
        phi_bar*phi_dot
    call check_condition(status == 0 .and. &
        abs(lhs - rhs) < 3.0e-9_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Toroidal Cartesian map VJP satisfies the adjoint identity")

    call cartesian_to_toroidal(point, scale, eta, theta, phi)
    call cartesian_to_toroidal_jvp( &
        point, scale, point_dot, scale_dot, eta_dot, theta_dot, phi_dot, status)
    call cartesian_to_toroidal( &
        point + step*point_dot, scale + step*scale_dot, eta_plus, theta_plus, &
        phi_plus)
    call cartesian_to_toroidal( &
        point - step*point_dot, scale - step*scale_dot, eta_minus, theta_minus, &
        phi_minus)
    relative_error = max(abs(eta_dot - (eta_plus - eta_minus)/(2.0_dp*step)), &
        abs(theta_dot - (theta_plus - theta_minus)/(2.0_dp*step)), &
        abs(phi_dot - (phi_plus - phi_minus)/(2.0_dp*step)))/ &
        max(1.0_dp, abs(eta_dot), abs(theta_dot), abs(phi_dot))
    call check_condition(status == 0 .and. relative_error < 4.0e-7_dp, &
        "Toroidal inverse map JVP matches re-evaluation")

    eta_bar = 0.22_dp
    theta_bar = -0.36_dp
    phi_bar = 0.41_dp
    call cartesian_to_toroidal_vjp( &
        point, scale, eta_bar, theta_bar, phi_bar, point_bar, scale_bar, status)
    lhs = eta_bar*eta_dot + theta_bar*theta_dot + phi_bar*phi_dot
    rhs = sum(point_bar*point_dot) + scale_bar*scale_dot
    call check_condition(status == 0 .and. &
        abs(lhs - rhs) < 5.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Toroidal inverse map VJP satisfies the adjoint identity")

    call check_summary("Differentiable toroidal coordinate maps")
end program test_toroidal_coordinates_ad
