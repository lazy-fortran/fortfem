program test_maxwell_torus_curved_rbc_pairing_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_torus_curved_rwg_rbc_pairing, &
        assemble_maxwell_torus_curved_rwg_rbc_pairing_jvp, &
        assemble_maxwell_torus_curved_rwg_rbc_pairing_vjp, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: step = 2.0e-6_dp
    real(dp) :: lhs, rhs, major_dot, minor_dot
    real(dp) :: major_radius_bar, minor_radius_bar
    real(dp), allocatable :: vertices(:, :), parameters(:, :)
    real(dp), allocatable :: vertices_dot(:, :), parameters_dot(:, :)
    real(dp), allocatable :: vertices_bar(:, :), parameters_bar(:, :)
    real(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    real(dp), allocatable :: matrix_minus(:, :), matrix_plus(:, :)
    real(dp), allocatable :: matrix_bar(:, :)
    integer, allocatable :: triangles(:, :)
    integer :: status, status_minus, status_plus, column, vertex
    logical :: all_passed

    all_passed = .true.
    major_dot = 0.003_dp
    minor_dot = -0.002_dp
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 4, 6, vertices, triangles, parameters)
    allocate( &
        vertices_dot(size(vertices, 1), size(vertices, 2)), &
        parameters_dot(size(parameters, 1), size(parameters, 2)))
    do vertex = 1, size(vertices, 2)
        vertices_dot(:, vertex) = [ &
            0.002_dp*sin(real(vertex, dp)), &
            -0.001_dp*cos(real(2*vertex, dp)), &
            0.0015_dp*sin(real(3*vertex, dp))]
        parameters_dot(:, vertex) = [ &
            0.002_dp*cos(real(vertex, dp)), &
            -0.001_dp*sin(real(2*vertex, dp))]
    end do
    call assemble_maxwell_torus_curved_rwg_rbc_pairing_jvp( &
        vertices, triangles, parameters, major_radius, minor_radius, 6, &
        vertices_dot, parameters_dot, major_dot, minor_dot, matrix, matrix_dot, &
        status)
    call assemble_maxwell_torus_curved_rwg_rbc_pairing( &
        vertices + step*vertices_dot, triangles, parameters + step*parameters_dot, &
        major_radius + step*major_dot, minor_radius + step*minor_dot, 6, &
        matrix_plus, status_plus)
    call assemble_maxwell_torus_curved_rwg_rbc_pairing( &
        vertices - step*vertices_dot, triangles, parameters - step*parameters_dot, &
        major_radius - step*major_dot, minor_radius - step*minor_dot, 6, &
        matrix_minus, status_minus)
    call record_condition(status == 0, "torus RWG-RBC pairing JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed torus RWG-RBC pairings succeed")
    call record_condition(maxval(abs(matrix_dot - &
        (matrix_plus - matrix_minus)/(2.0_dp*step))) < 3.0e-7_dp, &
        "torus RWG-RBC pairing JVP matches reassembly")

    allocate( &
        matrix_bar(size(matrix, 1), size(matrix, 2)), &
        vertices_bar(size(vertices, 1), size(vertices, 2)), &
        parameters_bar(size(parameters, 1), size(parameters, 2)))
    do column = 1, size(matrix_bar, 2)
        matrix_bar(:, column) = 0.004_dp*sin(real(column, dp)) + &
            0.003_dp*cos(real([(vertex, vertex=1, size(matrix_bar, 1))], dp))
    end do
    call assemble_maxwell_torus_curved_rwg_rbc_pairing_vjp( &
        vertices, triangles, parameters, major_radius, minor_radius, 6, &
        matrix_bar, vertices_bar, parameters_bar, major_radius_bar, &
        minor_radius_bar, status)
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(vertices_bar*vertices_dot) + sum(parameters_bar*parameters_dot) + &
        major_dot*major_radius_bar + minor_dot*minor_radius_bar
    call record_condition(status == 0, "torus RWG-RBC pairing VJP succeeds")
    call record_condition(abs(lhs - rhs) < 2.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "torus RWG-RBC pairing products obey the adjoint identity")

    call check_summary("Differentiable exact-curved torus RWG-RBC pairing")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_rbc_pairing_ad
