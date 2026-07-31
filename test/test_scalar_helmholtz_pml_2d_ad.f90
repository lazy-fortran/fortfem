program test_scalar_helmholtz_pml_2d_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        solve_scalar_helmholtz_pml_p1_2d, &
        solve_scalar_helmholtz_pml_p1_2d_jvp, &
        solve_scalar_helmholtz_pml_p1_2d_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-6_dp
    real(dp) :: vertices(2, 9), vertices_dot(2, 9), vertices_bar(2, 9)
    complex(dp) :: stretch(2, 8), stretch_dot(2, 8), stretch_bar(2, 8)
    complex(dp) :: volume_load(9), volume_load_dot(9), volume_load_bar(9)
    complex(dp) :: dirichlet_values(8), dirichlet_values_dot(8)
    complex(dp) :: dirichlet_values_bar(8), solution_bar(9)
    complex(dp), allocatable :: solution(:), solution_dot(:)
    complex(dp), allocatable :: solution_plus(:), solution_minus(:)
    integer :: triangles(3, 8), dirichlet_nodes(8)
    real(dp) :: wave_number, wave_number_dot, wave_number_bar
    real(dp) :: lhs, rhs
    integer :: status, status_minus, status_plus

    vertices = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, &
        0.0_dp, 0.5_dp, 1.0_dp, 0.5_dp, 2.0_dp, 0.5_dp, &
        0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp], [2, 9])
    vertices_dot = reshape([ &
        0.002_dp, -0.001_dp, -0.001_dp, 0.003_dp, 0.001_dp, -0.002_dp, &
        -0.002_dp, 0.001_dp, 0.002_dp, -0.003_dp, -0.001_dp, 0.002_dp, &
        0.001_dp, -0.002_dp, -0.002_dp, 0.001_dp, 0.002_dp, -0.001_dp], &
        [2, 9])
    triangles = reshape([ &
        1, 2, 5, 1, 5, 4, 2, 3, 6, 2, 6, 5, &
        4, 5, 8, 4, 8, 7, 5, 6, 9, 5, 9, 8], [3, 8])
    dirichlet_nodes = [1, 2, 3, 4, 6, 7, 8, 9]
    stretch = reshape([ &
        cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
        cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
        cmplx(1.1_dp, 0.1_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
        cmplx(1.2_dp, 0.2_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
        cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
        cmplx(1.1_dp, 0.15_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
        cmplx(1.3_dp, 0.25_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
        cmplx(1.2_dp, 0.05_dp, dp), cmplx(1.0_dp, 0.0_dp, dp)], [2, 8])
    stretch_dot = reshape([ &
        cmplx(0.01_dp, -0.02_dp, dp), cmplx(-0.01_dp, 0.01_dp, dp), &
        cmplx(0.02_dp, 0.01_dp, dp), cmplx(-0.02_dp, 0.02_dp, dp), &
        cmplx(-0.01_dp, 0.03_dp, dp), cmplx(0.01_dp, -0.01_dp, dp), &
        cmplx(0.02_dp, -0.02_dp, dp), cmplx(-0.01_dp, 0.01_dp, dp), &
        cmplx(0.01_dp, 0.02_dp, dp), cmplx(-0.02_dp, -0.01_dp, dp), &
        cmplx(-0.01_dp, 0.01_dp, dp), cmplx(0.02_dp, 0.03_dp, dp), &
        cmplx(0.01_dp, -0.02_dp, dp), cmplx(-0.02_dp, 0.01_dp, dp), &
        cmplx(0.02_dp, 0.01_dp, dp), cmplx(-0.01_dp, -0.02_dp, dp)], [2, 8])
    volume_load = cmplx(0.0_dp, 0.0_dp, dp)
    volume_load(5) = cmplx(0.7_dp, -0.2_dp, dp)
    volume_load_dot = cmplx(0.0_dp, 0.0_dp, dp)
    volume_load_dot(5) = cmplx(-0.03_dp, 0.04_dp, dp)
    dirichlet_values = [ &
        cmplx(0.1_dp, 0.0_dp, dp), cmplx(0.2_dp, 0.1_dp, dp), &
        cmplx(0.3_dp, -0.1_dp, dp), cmplx(0.0_dp, 0.2_dp, dp), &
        cmplx(0.4_dp, 0.0_dp, dp), cmplx(-0.1_dp, 0.1_dp, dp), &
        cmplx(0.2_dp, -0.2_dp, dp), cmplx(0.0_dp, 0.0_dp, dp)]
    dirichlet_values_dot = [ &
        cmplx(0.01_dp, -0.01_dp, dp), cmplx(-0.02_dp, 0.02_dp, dp), &
        cmplx(0.01_dp, 0.00_dp, dp), cmplx(0.00_dp, 0.01_dp, dp), &
        cmplx(-0.01_dp, 0.02_dp, dp), cmplx(0.02_dp, -0.01_dp, dp), &
        cmplx(-0.02_dp, 0.01_dp, dp), cmplx(0.01_dp, 0.01_dp, dp)]
    wave_number = 1.7_dp
    wave_number_dot = -0.06_dp
    solution_bar = [ &
        cmplx(0.1_dp, -0.1_dp, dp), cmplx(-0.2_dp, 0.1_dp, dp), &
        cmplx(0.05_dp, 0.2_dp, dp), cmplx(0.1_dp, 0.0_dp, dp), &
        cmplx(0.3_dp, -0.2_dp, dp), cmplx(-0.1_dp, 0.05_dp, dp), &
        cmplx(0.2_dp, 0.1_dp, dp), cmplx(-0.05_dp, -0.1_dp, dp), &
        cmplx(0.1_dp, 0.2_dp, dp)]

    call solve_scalar_helmholtz_pml_p1_2d( &
        vertices, triangles, stretch, wave_number, volume_load, &
        dirichlet_nodes, dirichlet_values, solution, status)
    call check_condition(status == 0, "2D scalar PML primal solve succeeds")
    call solve_scalar_helmholtz_pml_p1_2d_jvp( &
        vertices, triangles, stretch, wave_number, volume_load, &
        dirichlet_nodes, dirichlet_values, vertices_dot, stretch_dot, &
        wave_number_dot, volume_load_dot, dirichlet_values_dot, solution_dot, &
        status)
    call solve_scalar_helmholtz_pml_p1_2d( &
        vertices + step*vertices_dot, triangles, stretch + step*stretch_dot, &
        wave_number + step*wave_number_dot, volume_load + step*volume_load_dot, &
        dirichlet_nodes, dirichlet_values + step*dirichlet_values_dot, &
        solution_plus, status_plus)
    call solve_scalar_helmholtz_pml_p1_2d( &
        vertices - step*vertices_dot, triangles, stretch - step*stretch_dot, &
        wave_number - step*wave_number_dot, volume_load - step*volume_load_dot, &
        dirichlet_nodes, dirichlet_values - step*dirichlet_values_dot, &
        solution_minus, status_minus)
    call check_condition(status == 0, "2D scalar PML JVP succeeds")
    call check_condition(status_plus == 0 .and. status_minus == 0, &
        "2D scalar PML perturbed solves succeed")
    call check_condition(maxval(abs(solution_dot - &
        (solution_plus - solution_minus)/(2.0_dp*step))) < 2.0e-7_dp, &
        "2D scalar PML JVP matches central differences")

    call solve_scalar_helmholtz_pml_p1_2d_vjp( &
        vertices, triangles, stretch, wave_number, volume_load, &
        dirichlet_nodes, dirichlet_values, solution, solution_bar, &
        vertices_bar, stretch_bar, wave_number_bar, volume_load_bar, &
        dirichlet_values_bar, status)
    lhs = real(sum(conjg(solution_bar)*solution_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + &
        real(sum(conjg(stretch_bar)*stretch_dot), dp) + &
        wave_number_bar*wave_number_dot + &
        real(sum(conjg(volume_load_bar)*volume_load_dot), dp) + &
        real(sum(conjg(dirichlet_values_bar)*dirichlet_values_dot), dp)
    call check_condition(status == 0, "2D scalar PML VJP succeeds")
    call check_condition(abs(lhs - rhs) < 5.0e-9_dp*max(1.0_dp, abs(lhs)), &
        "2D scalar PML products obey the real adjoint identity")
    call check_summary("Differentiable scalar Helmholtz PML in 2D")
end program test_scalar_helmholtz_pml_2d_ad
