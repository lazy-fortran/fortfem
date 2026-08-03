program test_tetra_nedelec_curvilinear_pml_state_ad
    use check, only: check_condition, check_summary
    use fortfem_assembly_tetra_nedelec_3d, only: &
        assemble_tetra_nedelec_curvilinear_pml_csc, &
        assemble_tetra_nedelec_curvilinear_pml_csc_jvp, &
        assemble_tetra_nedelec_curvilinear_pml_csc_vjp
    use fortfem_tetra_nedelec_solver_3d, only: &
        solve_tetra_nedelec_curvilinear_pml
    use fortfem_tetra_nedelec_pml_state_3d, only: &
        solve_tetra_nedelec_curvilinear_pml_jvp, &
        solve_tetra_nedelec_curvilinear_pml_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_z_t, fortsparse_status_t
    implicit none

    integer, parameter :: order = 1
    real(dp), parameter :: step = 2.0e-7_dp
    integer :: tetrahedra(4, 2)
    real(dp) :: vertices(3, 5), vertices_dot(3, 5), vertices_bar(3, 5)
    complex(dp) :: stretch(3, 3, 2), stretch_dot(3, 3, 2)
    complex(dp) :: stretch_bar(3, 3, 2)
    integer, allocatable :: dirichlet_dofs(:)
    complex(dp), allocatable :: dirichlet_values(:), dirichlet_values_dot(:)
    complex(dp), allocatable :: volume_load(:), volume_load_dot(:)
    complex(dp), allocatable :: solution(:), solution_dot(:), plus(:), minus(:)
    complex(dp), allocatable :: solution_bar(:), volume_load_bar(:)
    complex(dp), allocatable :: dirichlet_values_bar(:)
    type(csc_z_t) :: matrix
    type(fortsparse_status_t) :: sparse_status
    real(dp) :: wave_number, wave_number_dot, wave_number_bar
    real(dp) :: lhs, rhs, relative_error
    integer :: dof, n, status
    logical :: all_passed

    all_passed = .true.
    tetrahedra = reshape([1, 2, 3, 4, 1, 3, 2, 5], [4, 2])
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
    stretch(:, :, 1) = reshape([ &
        cmplx(1.1_dp, 0.2_dp, dp), cmplx(0.08_dp, -0.03_dp, dp), &
        cmplx(-0.02_dp, 0.01_dp, dp), cmplx(0.04_dp, 0.06_dp, dp), &
        cmplx(0.95_dp, 0.15_dp, dp), cmplx(0.05_dp, -0.02_dp, dp), &
        cmplx(0.03_dp, 0.01_dp, dp), cmplx(-0.04_dp, 0.02_dp, dp), &
        cmplx(1.05_dp, 0.18_dp, dp)], [3, 3])
    stretch(:, :, 2) = reshape([ &
        cmplx(1.0_dp, 0.16_dp, dp), cmplx(-0.04_dp, 0.02_dp, dp), &
        cmplx(0.03_dp, -0.01_dp, dp), cmplx(0.06_dp, 0.04_dp, dp), &
        cmplx(1.08_dp, 0.12_dp, dp), cmplx(-0.02_dp, 0.03_dp, dp), &
        cmplx(0.01_dp, 0.02_dp, dp), cmplx(0.05_dp, -0.02_dp, dp), &
        cmplx(0.97_dp, 0.22_dp, dp)], [3, 3])
    stretch_dot(:, :, 1) = reshape([ &
        cmplx(0.03_dp, -0.01_dp, dp), cmplx(-0.02_dp, 0.04_dp, dp), &
        cmplx(0.01_dp, 0.02_dp, dp), cmplx(0.05_dp, 0.01_dp, dp), &
        cmplx(-0.03_dp, 0.02_dp, dp), cmplx(0.02_dp, -0.03_dp, dp), &
        cmplx(-0.01_dp, 0.03_dp, dp), cmplx(0.04_dp, -0.02_dp, dp), &
        cmplx(0.01_dp, 0.01_dp, dp)], [3, 3])
    stretch_dot(:, :, 2) = reshape([ &
        cmplx(-0.01_dp, 0.02_dp, dp), cmplx(0.025_dp, -0.015_dp, dp), &
        cmplx(0.02_dp, 0.01_dp, dp), cmplx(0.01_dp, 0.03_dp, dp), &
        cmplx(-0.02_dp, 0.01_dp, dp), cmplx(0.03_dp, -0.02_dp, dp), &
        cmplx(0.02_dp, 0.02_dp, dp), cmplx(-0.01_dp, 0.04_dp, dp), &
        cmplx(0.015_dp, -0.01_dp, dp)], [3, 3])
    wave_number = 1.1_dp
    wave_number_dot = 0.13_dp

    call assemble_tetra_nedelec_curvilinear_pml_csc( &
        vertices, tetrahedra, order, stretch, wave_number, matrix, sparse_status)
    call record_condition(sparse_status%code == 0, &
        "Curvilinear PML state matrix assembles")
    n = matrix%nrow
    allocate(volume_load(n), volume_load_dot(n), solution_bar(n))
    do dof = 1, n
        volume_load(dof) = cmplx(0.01_dp*dof, -0.02_dp*dof, dp)
        volume_load_dot(dof) = cmplx(0.002_dp*dof, -0.001_dp*dof, dp)
        solution_bar(dof) = cmplx(0.01_dp*dof, -0.02_dp*dof, dp)
    end do
    allocate(dirichlet_dofs(2), dirichlet_values(2), dirichlet_values_dot(2))
    dirichlet_dofs = [1, n]
    dirichlet_values = [cmplx(0.1_dp, -0.2_dp, dp), &
        cmplx(-0.15_dp, 0.05_dp, dp)]
    dirichlet_values_dot = [cmplx(0.003_dp, -0.004_dp, dp), &
        cmplx(-0.002_dp, 0.001_dp, dp)]

    call solve_tetra_nedelec_curvilinear_pml( &
        vertices, tetrahedra, order, stretch, wave_number, volume_load, &
        dirichlet_dofs, dirichlet_values, solution, status)
    call record_condition(status == 0, "Curvilinear PML state solve succeeds")

    call solve_tetra_nedelec_curvilinear_pml_jvp( &
        vertices, tetrahedra, order, stretch, wave_number, volume_load, &
        dirichlet_dofs, dirichlet_values, vertices_dot, stretch_dot, &
        wave_number_dot, volume_load_dot, dirichlet_values_dot, solution_dot, &
        status)
    call record_condition(status == 0, "Curvilinear PML state JVP succeeds")

    call solve_tetra_nedelec_curvilinear_pml( &
        vertices + step*vertices_dot, tetrahedra, order, &
        stretch + step*stretch_dot, wave_number + step*wave_number_dot, &
        volume_load + step*volume_load_dot, dirichlet_dofs, &
        dirichlet_values + step*dirichlet_values_dot, plus, status)
    call solve_tetra_nedelec_curvilinear_pml( &
        vertices - step*vertices_dot, tetrahedra, order, &
        stretch - step*stretch_dot, wave_number - step*wave_number_dot, &
        volume_load - step*volume_load_dot, dirichlet_dofs, &
        dirichlet_values - step*dirichlet_values_dot, minus, status)
    relative_error = maxval(abs(solution_dot - &
        (plus - minus)/(2.0_dp*step)))/max(1.0_dp, maxval(abs(solution_dot)))
    call record_condition(relative_error < 3.0e-6_dp, &
        "Curvilinear PML state JVP matches independent re-solves")

    allocate(volume_load_bar(n), dirichlet_values_bar(2))
    call solve_tetra_nedelec_curvilinear_pml_vjp( &
        vertices, tetrahedra, order, stretch, wave_number, volume_load, &
        dirichlet_dofs, dirichlet_values, solution, solution_bar, vertices_bar, &
        stretch_bar, wave_number_bar, volume_load_bar, dirichlet_values_bar, &
        status)
    lhs = real(sum(conjg(solution_bar)*solution_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + &
        real(sum(conjg(stretch_bar)*stretch_dot), dp) + &
        wave_number_bar*wave_number_dot + &
        real(sum(conjg(volume_load_bar)*volume_load_dot), dp) + &
        real(sum(conjg(dirichlet_values_bar)*dirichlet_values_dot), dp)
    call record_condition(status == 0, "Curvilinear PML state VJP succeeds")
    call record_condition( &
        abs(lhs - rhs) < 5.0e-7_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Curvilinear PML state VJP satisfies the real-complex adjoint identity")

    call check_summary("Curvilinear tetrahedral Nedelec PML state")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_nedelec_curvilinear_pml_state_ad
