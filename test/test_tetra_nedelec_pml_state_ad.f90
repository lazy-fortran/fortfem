program test_tetra_nedelec_pml_state_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_tetra_nedelec_pml_csc, &
        solve_tetra_nedelec_pml, solve_tetra_nedelec_pml_jvp, &
        solve_tetra_nedelec_pml_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_z_t, fortsparse_status_t
    implicit none

    integer, parameter :: order = 2
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(3, 5), vertices_dot(3, 5), vertices_bar(3, 5)
    integer :: tetrahedra(4, 2)
    integer, allocatable :: dirichlet_dofs(:), boundary_operator_dofs(:)
    complex(dp) :: stretch(3, 2), stretch_dot(3, 2), stretch_bar(3, 2)
    complex(dp) :: boundary_operator(3, 3), boundary_operator_dot(3, 3)
    complex(dp) :: boundary_operator_bar(3, 3)
    complex(dp), allocatable :: exact(:), volume_load(:), volume_load_dot(:)
    complex(dp), allocatable :: dirichlet_values(:), dirichlet_values_dot(:)
    complex(dp), allocatable :: solution(:), solution_dot(:), plus(:), minus(:)
    complex(dp), allocatable :: solution_bar(:)
    complex(dp), allocatable :: volume_load_bar(:), dirichlet_values_bar(:)
    complex(dp), allocatable :: stretch_plus(:, :), stretch_minus(:, :)
    complex(dp), allocatable :: load_plus(:), load_minus(:)
    complex(dp) :: dirichlet_plus(2), dirichlet_minus(2)
    complex(dp) :: operator_plus(3, 3), operator_minus(3, 3)
    type(csc_z_t) :: matrix
    type(fortsparse_status_t) :: sparse_status
    real(dp) :: wave_number, wave_number_dot, wave_bar
    real(dp) :: lhs, rhs, relative_error
    integer :: dof, n, status

    vertices = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, &
        0.0_dp, 0.0_dp, -1.0_dp], [3, 5])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.015_dp, &
        -0.03_dp, 0.01_dp, 0.02_dp, &
        0.025_dp, -0.015_dp, 0.005_dp, &
        -0.01_dp, 0.03_dp, -0.02_dp, &
        0.02_dp, -0.01_dp, 0.025_dp], [3, 5])
    tetrahedra = reshape([1, 2, 3, 4, 1, 3, 2, 5], [4, 2])
    stretch(:, 1) = [ &
        cmplx(1.1_dp, 0.35_dp, dp), cmplx(0.95_dp, 0.12_dp, dp), &
        cmplx(1.05_dp, 0.27_dp, dp)]
    stretch(:, 2) = [ &
        cmplx(1.0_dp, 0.22_dp, dp), cmplx(1.08_dp, 0.18_dp, dp), &
        cmplx(0.97_dp, 0.31_dp, dp)]
    stretch_dot(:, 1) = [ &
        cmplx(0.04_dp, -0.03_dp, dp), cmplx(-0.02_dp, 0.01_dp, dp), &
        cmplx(0.03_dp, 0.02_dp, dp)]
    stretch_dot(:, 2) = [ &
        cmplx(-0.01_dp, 0.02_dp, dp), cmplx(0.025_dp, -0.015_dp, dp), &
        cmplx(0.02_dp, 0.01_dp, dp)]
    wave_number = 1.1_dp
    wave_number_dot = 0.13_dp
    boundary_operator_dofs = [2, 3, 4]
    boundary_operator = reshape([ &
        cmplx(0.5_dp, -0.2_dp, dp), cmplx(0.1_dp, 0.05_dp, dp), &
        cmplx(-0.03_dp, 0.02_dp, dp), cmplx(0.1_dp, 0.05_dp, dp), &
        cmplx(0.7_dp, -0.1_dp, dp), cmplx(0.08_dp, -0.04_dp, dp), &
        cmplx(-0.03_dp, 0.02_dp, dp), cmplx(0.08_dp, -0.04_dp, dp), &
        cmplx(0.6_dp, -0.3_dp, dp)], [3, 3])
    boundary_operator_dot = reshape([ &
        cmplx(0.02_dp, -0.01_dp, dp), cmplx(-0.01_dp, 0.015_dp, dp), &
        cmplx(0.004_dp, 0.003_dp, dp), cmplx(0.005_dp, 0.01_dp, dp), &
        cmplx(-0.015_dp, 0.006_dp, dp), cmplx(0.012_dp, -0.008_dp, dp), &
        cmplx(0.003_dp, 0.004_dp, dp), cmplx(-0.008_dp, 0.002_dp, dp), &
        cmplx(0.01_dp, -0.012_dp, dp)], [3, 3])

    call assemble_tetra_nedelec_pml_csc( &
        vertices, tetrahedra, order, stretch, wave_number, matrix, sparse_status)
    call check_condition(sparse_status%code == 0, "Vector PML state matrix assembles")
    n = matrix%nrow
    allocate(exact(n), volume_load(n), volume_load_dot(n))
    do dof = 1, n
        exact(dof) = cmplx(sin(0.3_dp*real(dof, dp)), &
            cos(0.2_dp*real(dof, dp)), dp)
        volume_load_dot(dof) = cmplx(0.002_dp*dof, -0.001_dp*dof, dp)
    end do
    volume_load = csc_matvec(matrix, exact)
    volume_load(boundary_operator_dofs) = volume_load(boundary_operator_dofs) + &
        matmul(boundary_operator, exact(boundary_operator_dofs))
    allocate(dirichlet_dofs(2), dirichlet_values(2), dirichlet_values_dot(2))
    dirichlet_dofs = [1, n]
    dirichlet_values = exact(dirichlet_dofs)
    dirichlet_values_dot = [cmplx(0.003_dp, -0.004_dp, dp), &
        cmplx(-0.002_dp, 0.001_dp, dp)]

    call solve_tetra_nedelec_pml( &
        vertices, tetrahedra, order, stretch, wave_number, volume_load, &
        dirichlet_dofs, dirichlet_values, solution, status, &
        boundary_operator_dofs=boundary_operator_dofs, &
        boundary_operator=boundary_operator)
    call check_condition(status == 0, "Vector PML state solve succeeds")
    call check_condition(maxval(abs(solution - exact)) < 2.0e-10_dp, &
        "Vector PML state solve recovers a manufactured field")

    call solve_tetra_nedelec_pml_jvp( &
        vertices, tetrahedra, order, stretch, wave_number, volume_load, &
        dirichlet_dofs, dirichlet_values, vertices_dot, stretch_dot, &
        wave_number_dot, volume_load_dot, dirichlet_values_dot, solution_dot, &
        status, boundary_operator_dofs=boundary_operator_dofs, &
        boundary_operator=boundary_operator, &
        boundary_operator_dot=boundary_operator_dot)
    call check_condition(status == 0, "Vector PML state JVP succeeds")

    allocate(stretch_plus(3, 2), stretch_minus(3, 2))
    stretch_plus = stretch + step*stretch_dot
    stretch_minus = stretch - step*stretch_dot
    allocate(load_plus(n), load_minus(n))
    load_plus = volume_load + step*volume_load_dot
    load_minus = volume_load - step*volume_load_dot
    dirichlet_plus = dirichlet_values + step*dirichlet_values_dot
    dirichlet_minus = dirichlet_values - step*dirichlet_values_dot
    operator_plus = boundary_operator + step*boundary_operator_dot
    operator_minus = boundary_operator - step*boundary_operator_dot
    call solve_tetra_nedelec_pml( &
        vertices + step*vertices_dot, tetrahedra, order, stretch_plus, &
        wave_number + step*wave_number_dot, load_plus, dirichlet_dofs, &
        dirichlet_plus, plus, status, &
        boundary_operator_dofs=boundary_operator_dofs, &
        boundary_operator=operator_plus)
    call solve_tetra_nedelec_pml( &
        vertices - step*vertices_dot, tetrahedra, order, stretch_minus, &
        wave_number - step*wave_number_dot, load_minus, dirichlet_dofs, &
        dirichlet_minus, minus, status, &
        boundary_operator_dofs=boundary_operator_dofs, &
        boundary_operator=operator_minus)
    relative_error = maxval(abs(solution_dot - (plus - minus)/(2.0_dp*step)))/ &
        max(1.0_dp, maxval(abs(solution_dot)))
    call check_condition(relative_error < 2.0e-6_dp, &
        "Vector PML state JVP matches independent re-solves")

    allocate(volume_load_bar(n), dirichlet_values_bar(2))
    allocate(solution_bar(n))
    do dof = 1, n
        solution_bar(dof) = cmplx(0.01_dp*dof, -0.02_dp*dof, dp)
    end do
    call solve_tetra_nedelec_pml_vjp( &
        vertices, tetrahedra, order, stretch, wave_number, volume_load, &
        dirichlet_dofs, dirichlet_values, solution, &
        solution_bar, vertices_bar, &
        stretch_bar, wave_bar, volume_load_bar, dirichlet_values_bar, status, &
        boundary_operator_dofs=boundary_operator_dofs, &
        boundary_operator=boundary_operator, &
        boundary_operator_bar=boundary_operator_bar)
    lhs = real(sum(conjg(solution_bar)*solution_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + &
        real(sum(conjg(stretch_bar)*stretch_dot), dp) + wave_bar*wave_number_dot + &
        real(sum(conjg(volume_load_bar)*volume_load_dot), dp) + &
        real(sum(conjg(dirichlet_values_bar)*dirichlet_values_dot), dp) + &
        real(sum(conjg(boundary_operator_bar)*boundary_operator_dot), dp)
    call check_condition(status == 0, "Vector PML state VJP succeeds")
    call check_condition(abs(lhs - rhs) < 3.0e-7_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Vector PML state VJP satisfies the real-complex adjoint identity")

    call check_summary("Differentiable tetrahedral Nedelec PML state")
end program test_tetra_nedelec_pml_state_ad
