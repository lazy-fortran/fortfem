program test_maxwell_sphere_curved_rbc_pairing_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_maxwell_sphere_curved_rwg_rbc_pairing, &
        assemble_maxwell_sphere_curved_rwg_rbc_pairing_jvp, &
        assemble_maxwell_sphere_curved_rwg_rbc_pairing_vjp
    use fortfem_core, only: generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: radius = 1.3_dp
    real(dp), parameter :: step = 2.0e-6_dp
    real(dp) :: radius_dot, lhs, rhs
    real(dp) :: radius_bar
    real(dp), allocatable :: vertices(:, :), vertices_dot(:, :)
    real(dp), allocatable :: vertices_bar(:, :)
    real(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    real(dp), allocatable :: matrix_minus(:, :), matrix_plus(:, :)
    real(dp), allocatable :: matrix_bar(:, :)
    integer, allocatable :: triangles(:, :)
    integer :: status, status_minus, status_plus, column, vertex
    logical :: all_passed

    all_passed = .true.
    radius_dot = -0.004_dp
    call generate_sphere_surface_mesh(radius, 0, vertices, triangles)
    allocate(vertices_dot(size(vertices, 1), size(vertices, 2)))
    do vertex = 1, size(vertices, 2)
        vertices_dot(:, vertex) = [ &
            0.001_dp*sin(real(vertex, dp)), &
            -0.0015_dp*cos(real(2*vertex, dp)), &
            0.0008_dp*sin(real(3*vertex, dp))]
    end do
    call assemble_maxwell_sphere_curved_rwg_rbc_pairing_jvp( &
        vertices, triangles, radius, 6, vertices_dot, radius_dot, matrix, &
        matrix_dot, status)
    call assemble_maxwell_sphere_curved_rwg_rbc_pairing( &
        vertices + step*vertices_dot, triangles, radius + step*radius_dot, 6, &
        matrix_plus, status_plus)
    call assemble_maxwell_sphere_curved_rwg_rbc_pairing( &
        vertices - step*vertices_dot, triangles, radius - step*radius_dot, 6, &
        matrix_minus, status_minus)
    call record_condition(status == 0, &
        "sphere RWG-RBC pairing JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed sphere RWG-RBC pairings succeed")
    call record_condition(maxval(abs(matrix_dot - &
        (matrix_plus - matrix_minus)/(2.0_dp*step))) < 4.0e-7_dp, &
        "sphere RWG-RBC pairing JVP matches reassembly")

    allocate( &
        matrix_bar(size(matrix, 1), size(matrix, 2)), &
        vertices_bar(size(vertices, 1), size(vertices, 2)))
    do column = 1, size(matrix_bar, 2)
        matrix_bar(:, column) = 0.004_dp*sin(real(column, dp)) + &
            0.003_dp*cos(real([(vertex, vertex=1, size(matrix_bar, 1))], dp))
    end do
    call assemble_maxwell_sphere_curved_rwg_rbc_pairing_vjp( &
        vertices, triangles, radius, 6, matrix_bar, vertices_bar, radius_bar, &
        status)
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(vertices_bar*vertices_dot) + radius_dot*radius_bar
    call record_condition(status == 0, &
        "sphere RWG-RBC pairing VJP succeeds")
    call record_condition(abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "sphere RWG-RBC pairing products obey the adjoint identity")

    call check_summary("Differentiable exact-curved sphere RWG-RBC pairing")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_rbc_pairing_ad
