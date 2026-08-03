program test_helmholtz_bem_state_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: solve_helmholtz_dirichlet_p0_3d, &
        solve_helmholtz_dirichlet_p0_3d_jvp, &
        solve_helmholtz_dirichlet_p0_3d_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: quadrature_degree = 4
    real(dp), parameter :: step = 2.0e-7_dp
    integer, parameter :: triangles(3, 4) = reshape([ &
        1, 3, 2, 1, 2, 4, 2, 3, 4, 3, 1, 4], [3, 4])
    real(dp) :: vertices(3, 4), vertices_bar(3, 4), vertices_dot(3, 4)
    real(dp) :: wave_number, wave_number_bar, wave_number_dot
    complex(dp), allocatable :: density(:), density_dot(:)
    complex(dp), allocatable :: density_minus(:), density_plus(:)
    complex(dp) :: density_bar(4), boundary_value, boundary_value_bar
    complex(dp) :: boundary_value_dot
    real(dp) :: lhs, rhs
    integer :: status, status_minus, status_plus

    vertices = reshape([ &
        0.1_dp, -0.2_dp, 0.05_dp, &
        1.2_dp, 0.15_dp, -0.1_dp, &
        -0.25_dp, 1.1_dp, 0.2_dp, &
        0.3_dp, -0.1_dp, 1.4_dp], [3, 4])
    vertices_dot = reshape([ &
        0.02_dp, -0.01_dp, 0.03_dp, &
        -0.015_dp, 0.025_dp, -0.02_dp, &
        0.01_dp, 0.035_dp, -0.005_dp, &
        -0.01_dp, 0.02_dp, 0.015_dp], [3, 4])
    boundary_value = cmplx(1.3_dp, -0.25_dp, dp)
    boundary_value_dot = cmplx(-0.08_dp, 0.04_dp, dp)
    wave_number = 1.7_dp
    wave_number_dot = -0.13_dp
    density_bar = [ &
        cmplx(0.2_dp, -0.1_dp, dp), cmplx(-0.3_dp, 0.25_dp, dp), &
        cmplx(0.15_dp, 0.2_dp, dp), cmplx(0.4_dp, -0.15_dp, dp)]

    call solve_helmholtz_dirichlet_p0_3d( &
        vertices, triangles, boundary_value, wave_number, quadrature_degree, &
        density, status)
    call solve_helmholtz_dirichlet_p0_3d_jvp( &
        vertices, triangles, boundary_value, wave_number, quadrature_degree, &
        vertices_dot, boundary_value_dot, wave_number_dot, density_dot, status)
    call solve_helmholtz_dirichlet_p0_3d( &
        vertices + step*vertices_dot, triangles, &
        boundary_value + step*boundary_value_dot, &
        wave_number + step*wave_number_dot, quadrature_degree, density_plus, &
        status_plus)
    call solve_helmholtz_dirichlet_p0_3d( &
        vertices - step*vertices_dot, triangles, &
        boundary_value - step*boundary_value_dot, &
        wave_number - step*wave_number_dot, quadrature_degree, density_minus, &
        status_minus)
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "Helmholtz BEM state JVP succeeds")
    call check_condition(maxval(abs( &
        density_dot - (density_plus - density_minus)/(2.0_dp*step))) < &
        3.0e-8_dp, "Helmholtz BEM density JVP matches central differences")

    call solve_helmholtz_dirichlet_p0_3d_vjp( &
        vertices, triangles, boundary_value, wave_number, quadrature_degree, &
        density_bar, vertices_bar, boundary_value_bar, wave_number_bar, status)
    lhs = real(sum(conjg(density_bar)*density_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + &
        real(conjg(boundary_value_bar)*boundary_value_dot, dp) + &
        wave_number_bar*wave_number_dot
    call check_condition(status == 0, "Helmholtz BEM state VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Helmholtz BEM state products obey the adjoint identity")

    call check_summary("Helmholtz BEM state derivatives")
end program test_helmholtz_bem_state_ad
