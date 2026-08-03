program test_helmholtz_bem_panel_pair_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: assemble_helmholtz_single_layer_p0_3d, &
        integrate_helmholtz_single_layer_regular_panel_pair_p0_3d, &
        integrate_helmholtz_single_layer_regular_panel_pair_p0_3d_jvp, &
        integrate_helmholtz_single_layer_regular_panel_pair_p0_3d_vjp, &
        integrate_laplace_single_layer_regular_panel_pair_p0_3d
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: quadrature_degree = 4
    real(dp), parameter :: step = 2.0e-7_dp
    integer, parameter :: triangles(3, 2) = reshape([1, 2, 3, 4, 5, 6], [3, 2])
    real(dp) :: first(3, 3), first_bar(3, 3), first_dot(3, 3)
    real(dp) :: second(3, 3), second_bar(3, 3), second_dot(3, 3)
    real(dp) :: vertices(3, 6), wave_number, wave_number_bar, wave_number_dot
    complex(dp), allocatable :: assembled(:, :)
    complex(dp) :: value, value_bar, value_dot, value_minus, value_plus
    real(dp) :: laplace_value, lhs, rhs
    integer :: status, status_minus, status_plus

    first = reshape([ &
        0.1_dp, -0.2_dp, 0.05_dp, &
        1.2_dp, 0.15_dp, -0.1_dp, &
        -0.25_dp, 1.1_dp, 0.2_dp], [3, 3])
    second = reshape([ &
        0.3_dp, -0.1_dp, 1.4_dp, &
        1.15_dp, 0.4_dp, 1.65_dp, &
        -0.2_dp, 0.95_dp, 1.8_dp], [3, 3])
    first_dot = reshape([ &
        0.02_dp, -0.01_dp, 0.03_dp, &
        -0.015_dp, 0.025_dp, -0.02_dp, &
        0.01_dp, 0.035_dp, -0.005_dp], [3, 3])
    second_dot = reshape([ &
        -0.01_dp, 0.02_dp, 0.015_dp, &
        0.03_dp, -0.025_dp, 0.01_dp, &
        -0.02_dp, 0.005_dp, 0.025_dp], [3, 3])
    wave_number = 2.3_dp
    wave_number_dot = -0.17_dp
    value_bar = cmplx(-0.7_dp, 0.35_dp, dp)

    call integrate_helmholtz_single_layer_regular_panel_pair_p0_3d( &
        first, second, wave_number, quadrature_degree, value, status)
    vertices(:, 1:3) = first
    vertices(:, 4:6) = second
    call assemble_helmholtz_single_layer_p0_3d( &
        vertices, triangles, wave_number, quadrature_degree, assembled, &
        status_plus)
    call check_condition( &
        status == 0 .and. status_plus == 0, &
        "regular Helmholtz panel-pair integration succeeds")
    call check_condition( &
        abs(value - assembled(1, 2)) < 5.0e-15_dp, &
        "panel-pair primal matches the established Galerkin assembler")
    call integrate_laplace_single_layer_regular_panel_pair_p0_3d( &
        first, second, quadrature_degree, laplace_value, status)
    call integrate_helmholtz_single_layer_regular_panel_pair_p0_3d( &
        first, second, 0.0_dp, quadrature_degree, value, status_plus)
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. &
        abs(value - cmplx(laplace_value, 0.0_dp, dp)) < 5.0e-15_dp, &
        "zero-wave-number Helmholtz kernel reduces to Laplace")

    call integrate_helmholtz_single_layer_regular_panel_pair_p0_3d_jvp( &
        first, second, wave_number, quadrature_degree, first_dot, second_dot, &
        wave_number_dot, value_dot, status)
    call integrate_helmholtz_single_layer_regular_panel_pair_p0_3d( &
        first + step*first_dot, second + step*second_dot, &
        wave_number + step*wave_number_dot, quadrature_degree, value_plus, &
        status_plus)
    call integrate_helmholtz_single_layer_regular_panel_pair_p0_3d( &
        first - step*first_dot, second - step*second_dot, &
        wave_number - step*wave_number_dot, quadrature_degree, value_minus, &
        status_minus)
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "regular Helmholtz panel-pair JVP succeeds")
    call check_condition( &
        abs(value_dot - (value_plus - value_minus)/(2.0_dp*step)) < &
        3.0e-10_dp, "panel-pair JVP matches central differences")

    call integrate_helmholtz_single_layer_regular_panel_pair_p0_3d_vjp( &
        first, second, wave_number, quadrature_degree, value_bar, first_bar, &
        second_bar, wave_number_bar, status)
    lhs = real(conjg(value_bar)*value_dot, dp)
    rhs = sum(first_bar*first_dot) + sum(second_bar*second_dot) + &
        wave_number_bar*wave_number_dot
    call check_condition( &
        status == 0, "regular Helmholtz panel-pair VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 3.0e-13_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "panel-pair products obey the real-complex adjoint identity")

    call check_summary("Helmholtz BEM panel-pair derivatives")
end program test_helmholtz_bem_panel_pair_ad
