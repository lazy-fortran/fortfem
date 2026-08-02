program test_fortfem_core_facade
    use check, only: check_condition, check_summary
    use fortfem_core, only: &
        cell_complex_betti_numbers, cell_complex_euler_characteristic, &
        cell_complex_t, initialize_cell_complex, &
        toroidal_point_to_cartesian, toroidal_point_to_cartesian_jvp, &
        toroidal_point_to_cartesian_vjp, validate_cell_complex
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: scale = 2.5_dp, eta = 0.9_dp
    real(dp), parameter :: theta = 0.4_dp, phi = -0.7_dp
    real(dp), parameter :: scale_dot = 0.13_dp, eta_dot = -0.08_dp
    real(dp), parameter :: theta_dot = 0.11_dp, phi_dot = -0.06_dp
    real(dp), parameter :: finite_difference_step = 1.0e-7_dp
    real(dp) :: point(3), point_dot(3), point_bar(3)
    real(dp) :: point_plus(3), point_minus(3)
    real(dp) :: scale_bar, eta_bar, theta_bar, phi_bar
    real(dp) :: expected(3), finite_difference(3), adjoint_left, adjoint_right
    real(dp) :: denominator, radial_factor, vertical_factor
    integer :: status, euler, betti(4)
    integer :: boundary_1(3, 3), boundary_2(3, 1)
    type(cell_complex_t) :: complex

    denominator = cosh(eta) - cos(theta)
    radial_factor = scale*sinh(eta)/denominator
    vertical_factor = scale*sin(theta)/denominator
    expected = [radial_factor*cos(phi), radial_factor*sin(phi), vertical_factor]
    call toroidal_point_to_cartesian(scale, eta, theta, phi, point)
    call check_condition(maxval(abs(point - expected)) < 2.0e-14_dp, &
        "core facade exposes the analytical toroidal embedding")

    call toroidal_point_to_cartesian_jvp( &
        scale, eta, theta, phi, scale_dot, eta_dot, theta_dot, phi_dot, &
        point_dot, status)
    call toroidal_point_to_cartesian( &
        scale + finite_difference_step*scale_dot, eta + finite_difference_step*eta_dot, &
        theta + finite_difference_step*theta_dot, phi + finite_difference_step*phi_dot, &
        point_plus)
    call toroidal_point_to_cartesian( &
        scale - finite_difference_step*scale_dot, eta - finite_difference_step*eta_dot, &
        theta - finite_difference_step*theta_dot, phi - finite_difference_step*phi_dot, &
        point_minus)
    finite_difference = (point_plus - point_minus)/(2.0_dp*finite_difference_step)
    call check_condition(status == 0 .and. maxval(abs(point_dot - finite_difference)) < 3.0e-8_dp, &
        "core facade toroidal JVP matches an independent finite difference")

    point_bar = [0.4_dp, -0.7_dp, 0.2_dp]
    call toroidal_point_to_cartesian_vjp( &
        scale, eta, theta, phi, point_bar, scale_bar, eta_bar, theta_bar, phi_bar, status)
    adjoint_left = dot_product(point_bar, point_dot)
    adjoint_right = scale_bar*scale_dot + eta_bar*eta_dot + &
        theta_bar*theta_dot + phi_bar*phi_dot
    call check_condition(status == 0 .and. abs(adjoint_left - adjoint_right) < 2.0e-12_dp, &
        "core facade toroidal VJP satisfies the real adjoint identity")

    boundary_1 = reshape([ &
        -1, 0, 1, &
         1,-1, 0, &
         0, 1,-1], [3, 3])
    boundary_2(:, 1) = [1, 1, 1]
    call initialize_cell_complex(complex, 3, boundary_1, boundary_2, status=status)
    call validate_cell_complex(complex, status)
    call check_condition(status == 0, &
        "core facade accepts a valid oriented triangle cell complex")
    call cell_complex_euler_characteristic(complex, euler, status)
    call check_condition(status == 0 .and. euler == 1, &
        "core facade preserves the triangle Euler characteristic")
    call cell_complex_betti_numbers(complex, betti, status)
    call check_condition(status == 0 .and. all(betti == [1, 0, 0, 0]), &
        "core facade preserves the independent triangle Betti invariant")

    call check_summary("fortfem core facade")

end program test_fortfem_core_facade
