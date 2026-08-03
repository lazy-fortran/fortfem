program test_laplace_torus_panel_pair_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        integrate_laplace_torus_panel_p0_3d, &
        integrate_laplace_torus_panel_p0_3d_jvp, &
        integrate_laplace_torus_panel_p0_3d_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: step = 2.0e-6_dp
    real(dp) :: first_parameters(2, 3), first_parameters_bar(2, 3)
    real(dp) :: first_parameters_dot(2, 3)
    real(dp) :: second_parameters(2, 3), second_parameters_bar(2, 3)
    real(dp) :: second_parameters_dot(2, 3)
    real(dp) :: major_radius_bar, major_radius_dot, minor_radius_bar
    real(dp) :: minor_radius_dot, value, value_bar, value_dot
    real(dp) :: value_minus, value_plus, lhs, rhs
    integer :: quadrature_degree, status, status_minus, status_plus

    quadrature_degree = 2
    first_parameters = reshape([ &
        0.10_dp, 0.20_dp, &
        0.50_dp, 0.25_dp, &
        0.20_dp, 0.65_dp], [2, 3])
    second_parameters = reshape([ &
        1.00_dp, 2.40_dp, &
        1.30_dp, 2.50_dp, &
        0.90_dp, 2.80_dp], [2, 3])
    first_parameters_dot = reshape([ &
        0.012_dp, -0.008_dp, &
        -0.017_dp, 0.011_dp, &
        0.009_dp, 0.014_dp], [2, 3])
    second_parameters_dot = reshape([ &
        -0.006_dp, 0.013_dp, &
        0.015_dp, -0.009_dp, &
        0.004_dp, 0.018_dp], [2, 3])
    major_radius_dot = 0.025_dp
    minor_radius_dot = -0.014_dp
    value_bar = -0.37_dp

    call integrate_laplace_torus_panel_p0_3d( &
        first_parameters, second_parameters, major_radius, minor_radius, &
        quadrature_degree, value, status)
    call check_condition(status == 0, "regular torus Laplace panel pair succeeds")
    call check_condition(value > 0.0_dp, &
        "regular torus Laplace single-layer value is positive")

    call integrate_laplace_torus_panel_p0_3d_jvp( &
        first_parameters, second_parameters, major_radius, minor_radius, &
        quadrature_degree, first_parameters_dot, second_parameters_dot, &
        major_radius_dot, minor_radius_dot, value_dot, status)
    call integrate_laplace_torus_panel_p0_3d( &
        first_parameters + step*first_parameters_dot, &
        second_parameters + step*second_parameters_dot, &
        major_radius + step*major_radius_dot, minor_radius + step*minor_radius_dot, &
        quadrature_degree, value_plus, status_plus)
    call integrate_laplace_torus_panel_p0_3d( &
        first_parameters - step*first_parameters_dot, &
        second_parameters - step*second_parameters_dot, &
        major_radius - step*major_radius_dot, minor_radius - step*minor_radius_dot, &
        quadrature_degree, value_minus, status_minus)
    call check_condition(status == 0, "regular torus Laplace pair JVP succeeds")
    call check_condition(status_plus == 0, "positive torus pair perturbation succeeds")
    call check_condition(status_minus == 0, "negative torus pair perturbation succeeds")
    call check_condition( &
        abs(value_dot - (value_plus - value_minus)/(2.0_dp*step)) < 2.0e-8_dp, &
        "regular torus Laplace pair JVP matches central differences")

    call integrate_laplace_torus_panel_p0_3d_vjp( &
        first_parameters, second_parameters, major_radius, minor_radius, &
        quadrature_degree, value_bar, first_parameters_bar, &
        second_parameters_bar, major_radius_bar, minor_radius_bar, status)
    lhs = value_bar*value_dot
    rhs = sum(first_parameters_bar*first_parameters_dot) + &
        sum(second_parameters_bar*second_parameters_dot) + &
        major_radius_bar*major_radius_dot + minor_radius_bar*minor_radius_dot
    call check_condition(status == 0, "regular torus Laplace pair VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 2.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "regular torus Laplace pair products obey the adjoint identity")

    call integrate_laplace_torus_panel_p0_3d( &
        first_parameters, first_parameters, major_radius, minor_radius, &
        quadrature_degree, value, status)
    call check_condition(status /= 0, "torus pair rejects coincident panels")

    call check_summary("Differentiable torus Laplace panel pair")
end program test_laplace_torus_panel_pair_ad
