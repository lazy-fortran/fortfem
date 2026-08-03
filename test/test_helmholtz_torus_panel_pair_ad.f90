program test_helmholtz_torus_panel_pair_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        integrate_helmholtz_torus_panel_p0_3d, &
        integrate_helmholtz_torus_panel_p0_3d_jvp, &
        integrate_helmholtz_torus_panel_p0_3d_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 1.3_dp, step = 2.0e-6_dp
    real(dp) :: first_parameters(2, 3), first_parameters_bar(2, 3)
    real(dp) :: first_parameters_dot(2, 3)
    real(dp) :: second_parameters(2, 3), second_parameters_bar(2, 3)
    real(dp) :: second_parameters_dot(2, 3)
    real(dp) :: major_radius_bar, major_radius_dot, minor_radius_bar
    real(dp) :: minor_radius_dot, wave_number_bar, wave_number_dot
    complex(dp) :: value, value_bar, value_dot, value_minus, value_plus
    real(dp) :: lhs, rhs
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
    wave_number_dot = 0.06_dp
    value_bar = cmplx(-0.37_dp, 0.23_dp, dp)

    call integrate_helmholtz_torus_panel_p0_3d( &
        first_parameters, second_parameters, major_radius, minor_radius, &
        wave_number, quadrature_degree, value, status)
    call check_condition(status == 0, "regular torus Helmholtz panel pair succeeds")

    call integrate_helmholtz_torus_panel_p0_3d_jvp( &
        first_parameters, second_parameters, major_radius, minor_radius, &
        wave_number, quadrature_degree, first_parameters_dot, &
        second_parameters_dot, major_radius_dot, minor_radius_dot, &
        wave_number_dot, value_dot, status)
    call integrate_helmholtz_torus_panel_p0_3d( &
        first_parameters + step*first_parameters_dot, &
        second_parameters + step*second_parameters_dot, &
        major_radius + step*major_radius_dot, minor_radius + step*minor_radius_dot, &
        wave_number + step*wave_number_dot, quadrature_degree, value_plus, &
        status_plus)
    call integrate_helmholtz_torus_panel_p0_3d( &
        first_parameters - step*first_parameters_dot, &
        second_parameters - step*second_parameters_dot, &
        major_radius - step*major_radius_dot, minor_radius - step*minor_radius_dot, &
        wave_number - step*wave_number_dot, quadrature_degree, value_minus, &
        status_minus)
    call check_condition(status == 0, "regular torus Helmholtz pair JVP succeeds")
    call check_condition(status_plus == 0, "positive torus Helmholtz perturbation succeeds")
    call check_condition(status_minus == 0, "negative torus Helmholtz perturbation succeeds")
    call check_condition( &
        abs(value_dot - (value_plus - value_minus)/(2.0_dp*step)) < 3.0e-8_dp, &
        "regular torus Helmholtz pair JVP matches central differences")

    call integrate_helmholtz_torus_panel_p0_3d_vjp( &
        first_parameters, second_parameters, major_radius, minor_radius, &
        wave_number, quadrature_degree, value_bar, first_parameters_bar, &
        second_parameters_bar, major_radius_bar, minor_radius_bar, &
        wave_number_bar, status)
    lhs = real(conjg(value_bar)*value_dot, dp)
    rhs = sum(first_parameters_bar*first_parameters_dot) + &
        sum(second_parameters_bar*second_parameters_dot) + &
        major_radius_bar*major_radius_dot + minor_radius_bar*minor_radius_dot + &
        wave_number_bar*wave_number_dot
    call check_condition(status == 0, "regular torus Helmholtz pair VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 3.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "regular torus Helmholtz pair products obey the adjoint identity")

    call integrate_helmholtz_torus_panel_p0_3d( &
        first_parameters, first_parameters, major_radius, minor_radius, &
        wave_number, quadrature_degree, value, status)
    call check_condition(status /= 0, "Helmholtz torus pair rejects coincident panels")

    call check_summary("Differentiable torus Helmholtz panel pair")
end program test_helmholtz_torus_panel_pair_ad
