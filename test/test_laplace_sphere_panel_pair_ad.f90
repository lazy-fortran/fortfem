program test_laplace_sphere_panel_pair_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        integrate_laplace_sphere_panel_p0_3d, &
        integrate_laplace_sphere_panel_p0_3d_jvp, &
        integrate_laplace_sphere_panel_p0_3d_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: sphere_radius = 1.0_dp
    real(dp), parameter :: step = 2.0e-6_dp
    real(dp) :: first_vertices(3, 3), first_vertices_bar(3, 3)
    real(dp) :: first_vertices_dot(3, 3)
    real(dp) :: second_vertices(3, 3), second_vertices_bar(3, 3)
    real(dp) :: second_vertices_dot(3, 3)
    real(dp) :: radius_bar, radius_dot, value, value_bar, value_dot
    real(dp) :: value_minus, value_plus, lhs, rhs
    integer :: quadrature_degree, status, status_minus, status_plus

    quadrature_degree = 2
    first_vertices(:, 1) = [1.00_dp, 0.00_dp, 0.10_dp]
    first_vertices(:, 2) = [0.75_dp, 0.35_dp, 0.05_dp]
    first_vertices(:, 3) = [0.72_dp, 0.04_dp, 0.38_dp]
    second_vertices(:, 1) = [-1.00_dp, 0.00_dp, -0.10_dp]
    second_vertices(:, 2) = [-0.75_dp, 0.32_dp, -0.04_dp]
    second_vertices(:, 3) = [-0.72_dp, 0.03_dp, -0.36_dp]
    first_vertices_dot = reshape([ &
        0.012_dp, -0.008_dp, 0.006_dp, &
        -0.017_dp, 0.011_dp, 0.009_dp, &
        0.009_dp, 0.014_dp, -0.010_dp], [3, 3])
    second_vertices_dot = reshape([ &
        -0.006_dp, 0.013_dp, -0.005_dp, &
        0.015_dp, -0.009_dp, 0.007_dp, &
        0.004_dp, 0.018_dp, -0.011_dp], [3, 3])
    radius_dot = 0.025_dp
    value_bar = -0.37_dp

    call integrate_laplace_sphere_panel_p0_3d( &
        first_vertices, second_vertices, sphere_radius, quadrature_degree, &
        value, status)
    call check_condition(status == 0, "regular sphere Laplace panel pair succeeds")
    call check_condition(value > 0.0_dp, &
        "regular sphere Laplace single-layer value is positive")

    call integrate_laplace_sphere_panel_p0_3d_jvp( &
        first_vertices, second_vertices, sphere_radius, quadrature_degree, &
        first_vertices_dot, second_vertices_dot, radius_dot, value_dot, status)
    call integrate_laplace_sphere_panel_p0_3d( &
        first_vertices + step*first_vertices_dot, &
        second_vertices + step*second_vertices_dot, &
        sphere_radius + step*radius_dot, quadrature_degree, value_plus, &
        status_plus)
    call integrate_laplace_sphere_panel_p0_3d( &
        first_vertices - step*first_vertices_dot, &
        second_vertices - step*second_vertices_dot, &
        sphere_radius - step*radius_dot, quadrature_degree, value_minus, &
        status_minus)
    call check_condition(status == 0, "regular sphere Laplace pair JVP succeeds")
    call check_condition(status_plus == 0, "positive sphere pair perturbation succeeds")
    call check_condition(status_minus == 0, "negative sphere pair perturbation succeeds")
    call check_condition( &
        abs(value_dot - (value_plus - value_minus)/(2.0_dp*step)) < 3.0e-8_dp, &
        "regular sphere Laplace pair JVP matches central differences")

    call integrate_laplace_sphere_panel_p0_3d_vjp( &
        first_vertices, second_vertices, sphere_radius, quadrature_degree, &
        value_bar, first_vertices_bar, second_vertices_bar, radius_bar, status)
    lhs = value_bar*value_dot
    rhs = sum(first_vertices_bar*first_vertices_dot) + &
        sum(second_vertices_bar*second_vertices_dot) + radius_bar*radius_dot
    call check_condition(status == 0, "regular sphere Laplace pair VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 3.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "regular sphere Laplace pair products obey the adjoint identity")

    call integrate_laplace_sphere_panel_p0_3d( &
        first_vertices, first_vertices, sphere_radius, quadrature_degree, &
        value, status)
    call check_condition(status /= 0, "sphere pair rejects coincident panels")

    call integrate_laplace_sphere_panel_p0_3d( &
        first_vertices, second_vertices, 0.0_dp, quadrature_degree, value, status)
    call check_condition(status /= 0, "sphere pair rejects nonpositive radius")

    call check_summary("Differentiable sphere Laplace panel pair")
end program test_laplace_sphere_panel_pair_ad
