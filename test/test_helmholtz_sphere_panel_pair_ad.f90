program test_helmholtz_sphere_panel_pair_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        integrate_helmholtz_sphere_panel_p0_3d, &
        integrate_helmholtz_sphere_panel_p0_3d_jvp, &
        integrate_helmholtz_sphere_panel_p0_3d_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: sphere_radius = 1.0_dp
    real(dp), parameter :: wave_number = 1.3_dp, step = 2.0e-6_dp
    real(dp) :: first_vertices(3, 3), first_vertices_bar(3, 3)
    real(dp) :: first_vertices_dot(3, 3)
    real(dp) :: second_vertices(3, 3), second_vertices_bar(3, 3)
    real(dp) :: second_vertices_dot(3, 3)
    real(dp) :: radius_bar, radius_dot, wave_number_bar, wave_number_dot
    complex(dp) :: value, value_bar, value_dot, value_minus, value_plus
    real(dp) :: lhs, rhs
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
    wave_number_dot = 0.06_dp
    value_bar = cmplx(-0.37_dp, 0.23_dp, dp)

    call integrate_helmholtz_sphere_panel_p0_3d( &
        first_vertices, second_vertices, sphere_radius, wave_number, &
        quadrature_degree, value, status)
    call check_condition(status == 0, &
        "regular sphere Helmholtz panel pair succeeds")

    call integrate_helmholtz_sphere_panel_p0_3d_jvp( &
        first_vertices, second_vertices, sphere_radius, wave_number, &
        quadrature_degree, first_vertices_dot, second_vertices_dot, &
        radius_dot, wave_number_dot, value_dot, status)
    call integrate_helmholtz_sphere_panel_p0_3d( &
        first_vertices + step*first_vertices_dot, &
        second_vertices + step*second_vertices_dot, &
        sphere_radius + step*radius_dot, wave_number + step*wave_number_dot, &
        quadrature_degree, value_plus, status_plus)
    call integrate_helmholtz_sphere_panel_p0_3d( &
        first_vertices - step*first_vertices_dot, &
        second_vertices - step*second_vertices_dot, &
        sphere_radius - step*radius_dot, wave_number - step*wave_number_dot, &
        quadrature_degree, value_minus, status_minus)
    call check_condition(status == 0, &
        "regular sphere Helmholtz pair JVP succeeds")
    call check_condition(status_plus == 0, &
        "positive sphere Helmholtz perturbation succeeds")
    call check_condition(status_minus == 0, &
        "negative sphere Helmholtz perturbation succeeds")
    call check_condition( &
        abs(value_dot - (value_plus - value_minus)/(2.0_dp*step)) < 4.0e-8_dp, &
        "regular sphere Helmholtz pair JVP matches central differences")

    call integrate_helmholtz_sphere_panel_p0_3d_vjp( &
        first_vertices, second_vertices, sphere_radius, wave_number, &
        quadrature_degree, value_bar, first_vertices_bar, second_vertices_bar, &
        radius_bar, wave_number_bar, status)
    lhs = real(conjg(value_bar)*value_dot, dp)
    rhs = sum(first_vertices_bar*first_vertices_dot) + &
        sum(second_vertices_bar*second_vertices_dot) + radius_bar*radius_dot + &
        wave_number_bar*wave_number_dot
    call check_condition(status == 0, &
        "regular sphere Helmholtz pair VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 4.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "regular sphere Helmholtz pair products obey the adjoint identity")

    call integrate_helmholtz_sphere_panel_p0_3d( &
        first_vertices, first_vertices, sphere_radius, wave_number, &
        quadrature_degree, value, status)
    call check_condition(status /= 0, &
        "sphere Helmholtz pair rejects coincident panels")

    call check_summary("Differentiable sphere Helmholtz panel pair")
end program test_helmholtz_sphere_panel_pair_ad
