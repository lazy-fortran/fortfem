program test_maxwell_sphere_curved_rwg_mass_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_sphere_curved_rwg_mass_matrix, &
        assemble_maxwell_sphere_curved_rwg_mass_matrix_jvp, &
        assemble_maxwell_sphere_curved_rwg_mass_matrix_vjp, &
        generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: radius = 1.0_dp, radius_dot = 0.017_dp
    real(dp), parameter :: step = 2.0e-6_dp
    real(dp), allocatable :: vertices(:, :), vertices_dot(:, :)
    real(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    real(dp), allocatable :: matrix_plus(:, :), matrix_minus(:, :)
    real(dp), allocatable :: matrix_bar(:, :), vertices_bar(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp) :: radius_bar, lhs, rhs, jvp_error
    integer :: i, j, status, status_minus, status_plus
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(radius, 0, vertices, triangles)
    allocate(vertices_dot(size(vertices, 1), size(vertices, 2)))
    do j = 1, size(vertices_dot, 2)
        vertices_dot(:, j) = [ &
            0.013_dp*sin(real(2*j, dp)), &
            -0.009_dp*cos(real(3*j, dp)), &
            0.011_dp*sin(real(5*j, dp))]
    end do

    call assemble_maxwell_sphere_curved_rwg_mass_matrix_jvp( &
        vertices, triangles, radius, 3, vertices_dot, radius_dot, matrix, &
        matrix_dot, status)
    call assemble_maxwell_sphere_curved_rwg_mass_matrix( &
        vertices + step*vertices_dot, triangles, radius + step*radius_dot, 3, &
        matrix_plus, status_plus)
    call assemble_maxwell_sphere_curved_rwg_mass_matrix( &
        vertices - step*vertices_dot, triangles, radius - step*radius_dot, 3, &
        matrix_minus, status_minus)
    if (status == 0) then
        if (status_plus == 0 .and. status_minus == 0) then
            jvp_error = maxval(abs(matrix_dot - &
                (matrix_plus - matrix_minus)/(2.0_dp*step)))
        else
            jvp_error = huge(1.0_dp)
        end if
    else
        jvp_error = huge(1.0_dp)
    end if
    call record_condition(status == 0, &
        "spherical RWG mass geometry JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed spherical RWG masses assemble")
    call record_condition(jvp_error < 2.0e-7_dp, &
        "spherical RWG mass geometry JVP matches reassembly")

    allocate(matrix_bar(size(matrix, 1), size(matrix, 2)))
    do j = 1, size(matrix_bar, 2)
        do i = 1, size(matrix_bar, 1)
            matrix_bar(i, j) = sin(real(2*i + 3*j, dp))
        end do
    end do
    allocate(vertices_bar(size(vertices, 1), size(vertices, 2)))
    call assemble_maxwell_sphere_curved_rwg_mass_matrix_vjp( &
        vertices, triangles, radius, 3, matrix_bar, vertices_bar, radius_bar, &
        status)
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(vertices_bar*vertices_dot) + radius_bar*radius_dot
    call record_condition(status == 0, "spherical RWG mass geometry VJP succeeds")
    call record_condition(abs(lhs - rhs) < &
        5.0e-9_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "spherical RWG mass geometry VJP satisfies the adjoint identity")

    call check_summary("Differentiable spherical RWG mass geometry")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_rwg_mass_ad
