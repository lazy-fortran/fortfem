program test_scalar_helmholtz_pml_3d_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        solve_scalar_helmholtz_pml_p1_3d, &
        solve_scalar_helmholtz_pml_p1_3d_jvp, &
        solve_scalar_helmholtz_pml_p1_3d_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-6_dp
    real(dp) :: vertices(3, 5), vertices_dot(3, 5), vertices_bar(3, 5)
    complex(dp) :: stretch(3, 2), stretch_dot(3, 2), stretch_bar(3, 2)
    complex(dp) :: volume_load(5), volume_load_dot(5), volume_load_bar(5)
    complex(dp) :: dirichlet_values(4), dirichlet_values_dot(4)
    complex(dp) :: dirichlet_values_bar(4), solution_bar(5)
    complex(dp), allocatable :: solution(:), solution_dot(:)
    complex(dp), allocatable :: solution_plus(:), solution_minus(:)
    integer :: tetrahedra(4, 2), dirichlet_nodes(4)
    real(dp) :: wave_number, wave_number_dot, wave_number_bar
    real(dp) :: lhs, rhs
    integer :: status, status_minus, status_plus

    vertices = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, &
        0.2_dp, 0.2_dp, 0.2_dp], [3, 5])
    vertices_dot = reshape([ &
        0.001_dp, -0.002_dp, 0.001_dp, &
        -0.002_dp, 0.001_dp, 0.002_dp, &
        0.002_dp, 0.001_dp, -0.001_dp, &
        -0.001_dp, 0.002_dp, 0.001_dp, &
        0.003_dp, -0.002_dp, 0.001_dp], [3, 5])
    tetrahedra = reshape([1, 2, 3, 5, 1, 3, 4, 5], [4, 2])
    dirichlet_nodes = [1, 2, 3, 4]
    stretch = reshape([ &
        cmplx(1.1_dp, 0.1_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
        cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.2_dp, 0.2_dp, dp), &
        cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp)], [3, 2])
    stretch_dot = reshape([ &
        cmplx(0.01_dp, -0.02_dp, dp), cmplx(0.02_dp, 0.01_dp, dp), &
        cmplx(-0.01_dp, 0.02_dp, dp), cmplx(0.01_dp, 0.01_dp, dp), &
        cmplx(0.02_dp, 0.00_dp, dp), cmplx(-0.01_dp, 0.02_dp, dp)], [3, 2])
    volume_load = cmplx(0.0_dp, 0.0_dp, dp)
    volume_load(5) = cmplx(0.6_dp, -0.15_dp, dp)
    volume_load_dot = cmplx(0.0_dp, 0.0_dp, dp)
    volume_load_dot(5) = cmplx(-0.02_dp, 0.03_dp, dp)
    dirichlet_values = [ &
        cmplx(0.1_dp, 0.0_dp, dp), cmplx(0.2_dp, 0.1_dp, dp), &
        cmplx(0.0_dp, -0.1_dp, dp), cmplx(0.3_dp, 0.2_dp, dp)]
    dirichlet_values_dot = [ &
        cmplx(0.01_dp, -0.01_dp, dp), cmplx(-0.02_dp, 0.01_dp, dp), &
        cmplx(0.00_dp, 0.02_dp, dp), cmplx(0.01_dp, 0.01_dp, dp)]
    wave_number = 1.4_dp
    wave_number_dot = 0.04_dp
    solution_bar = [ &
        cmplx(0.1_dp, -0.1_dp, dp), cmplx(-0.2_dp, 0.1_dp, dp), &
        cmplx(0.05_dp, 0.2_dp, dp), cmplx(0.1_dp, 0.0_dp, dp), &
        cmplx(0.3_dp, -0.2_dp, dp)]

    call solve_scalar_helmholtz_pml_p1_3d( &
        vertices, tetrahedra, stretch, wave_number, volume_load, &
        dirichlet_nodes, dirichlet_values, solution, status)
    call check_condition(status == 0, "3D scalar PML primal solve succeeds")
    call solve_scalar_helmholtz_pml_p1_3d_jvp( &
        vertices, tetrahedra, stretch, wave_number, volume_load, &
        dirichlet_nodes, dirichlet_values, vertices_dot, stretch_dot, &
        wave_number_dot, volume_load_dot, dirichlet_values_dot, solution_dot, &
        status)
    call solve_scalar_helmholtz_pml_p1_3d( &
        vertices + step*vertices_dot, tetrahedra, stretch + step*stretch_dot, &
        wave_number + step*wave_number_dot, volume_load + step*volume_load_dot, &
        dirichlet_nodes, dirichlet_values + step*dirichlet_values_dot, &
        solution_plus, status_plus)
    call solve_scalar_helmholtz_pml_p1_3d( &
        vertices - step*vertices_dot, tetrahedra, stretch - step*stretch_dot, &
        wave_number - step*wave_number_dot, volume_load - step*volume_load_dot, &
        dirichlet_nodes, dirichlet_values - step*dirichlet_values_dot, &
        solution_minus, status_minus)
    call check_condition(status == 0, "3D scalar PML JVP succeeds")
    call check_condition(status_plus == 0 .and. status_minus == 0, &
        "3D scalar PML perturbed solves succeed")
    call check_condition(maxval(abs(solution_dot - &
        (solution_plus - solution_minus)/(2.0_dp*step))) < 3.0e-7_dp, &
        "3D scalar PML JVP matches central differences")

    call solve_scalar_helmholtz_pml_p1_3d_vjp( &
        vertices, tetrahedra, stretch, wave_number, volume_load, &
        dirichlet_nodes, dirichlet_values, solution, solution_bar, &
        vertices_bar, stretch_bar, wave_number_bar, volume_load_bar, &
        dirichlet_values_bar, status)
    lhs = real(sum(conjg(solution_bar)*solution_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + &
        real(sum(conjg(stretch_bar)*stretch_dot), dp) + &
        wave_number_bar*wave_number_dot + &
        real(sum(conjg(volume_load_bar)*volume_load_dot), dp) + &
        real(sum(conjg(dirichlet_values_bar)*dirichlet_values_dot), dp)
    call check_condition(status == 0, "3D scalar PML VJP succeeds")
    call check_condition(abs(lhs - rhs) < 5.0e-9_dp*max(1.0_dp, abs(lhs)), &
        "3D scalar PML products obey the real adjoint identity")
    call check_summary("Differentiable scalar Helmholtz PML in 3D")
end program test_scalar_helmholtz_pml_3d_ad
