program test_laplace_bem_state_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: solve_laplace_dirichlet_p0_3d, &
        solve_laplace_dirichlet_p0_3d_jvp, &
        solve_laplace_dirichlet_p0_3d_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: quadrature_degree = 4
    real(dp), parameter :: step = 2.0e-7_dp
    integer, parameter :: triangles(3, 4) = reshape([ &
        1, 3, 2, 1, 2, 4, 2, 3, 4, 3, 1, 4], [3, 4])
    real(dp) :: vertices(3, 4), vertices_bar(3, 4), vertices_dot(3, 4)
    real(dp), allocatable :: density(:), density_dot(:)
    real(dp), allocatable :: density_minus(:), density_plus(:)
    real(dp) :: density_bar(4), boundary_value, boundary_value_bar
    real(dp) :: boundary_value_dot, capacity, capacity_bar, capacity_dot
    real(dp) :: capacity_minus, capacity_plus, lhs, rhs
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
    boundary_value = 1.3_dp
    boundary_value_dot = -0.08_dp
    density_bar = [0.2_dp, -0.3_dp, 0.15_dp, 0.4_dp]
    capacity_bar = -0.6_dp

    call solve_laplace_dirichlet_p0_3d( &
        vertices, triangles, boundary_value, quadrature_degree, density, &
        capacity, status)
    call solve_laplace_dirichlet_p0_3d_jvp( &
        vertices, triangles, boundary_value, quadrature_degree, vertices_dot, &
        boundary_value_dot, density_dot, capacity_dot, status)
    call solve_laplace_dirichlet_p0_3d( &
        vertices + step*vertices_dot, triangles, &
        boundary_value + step*boundary_value_dot, quadrature_degree, &
        density_plus, capacity_plus, status_plus)
    call solve_laplace_dirichlet_p0_3d( &
        vertices - step*vertices_dot, triangles, &
        boundary_value - step*boundary_value_dot, quadrature_degree, &
        density_minus, capacity_minus, status_minus)
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "Laplace BEM state JVP succeeds")
    call check_condition(maxval(abs( &
        density_dot - (density_plus - density_minus)/(2.0_dp*step))) < &
        2.0e-8_dp, "Laplace BEM density JVP matches central differences")
    call check_condition(abs( &
        capacity_dot - (capacity_plus - capacity_minus)/(2.0_dp*step)) < &
        2.0e-8_dp, "Laplace BEM capacity JVP matches central differences")

    call solve_laplace_dirichlet_p0_3d_vjp( &
        vertices, triangles, boundary_value, quadrature_degree, density_bar, &
        capacity_bar, vertices_bar, boundary_value_bar, status)
    lhs = dot_product(density_bar, density_dot) + capacity_bar*capacity_dot
    rhs = sum(vertices_bar*vertices_dot) + &
        boundary_value_bar*boundary_value_dot
    call check_condition(status == 0, "Laplace BEM state VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 2.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Laplace BEM state products obey the adjoint identity")

    call check_summary("Laplace BEM state derivatives")
end program test_laplace_bem_state_ad
