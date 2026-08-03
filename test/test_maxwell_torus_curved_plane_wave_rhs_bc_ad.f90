program test_maxwell_torus_curved_plane_wave_rhs_bc_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d, &
        assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d_jvp, &
        assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d_vjp
    use fortfem_core, only: generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 0.43_dp, step = 2.0e-6_dp
    real(dp) :: major_radius_dot, minor_radius_dot, wave_number_dot
    real(dp) :: major_radius_bar, minor_radius_bar, wave_number_bar
    real(dp) :: direction(3), direction_dot(3), direction_bar(3)
    complex(dp) :: polarization(3), polarization_dot(3)
    complex(dp), allocatable :: polarization_bar(:)
    complex(dp), allocatable :: right_hand_side(:), right_hand_side_dot(:)
    complex(dp), allocatable :: right_hand_side_minus(:), right_hand_side_plus(:)
    complex(dp), allocatable :: right_hand_side_bar(:)
    real(dp), allocatable :: vertices(:, :), parameters(:, :)
    real(dp), allocatable :: vertices_dot(:, :), parameters_dot(:, :)
    real(dp), allocatable :: vertices_plus(:, :), vertices_minus(:, :)
    real(dp), allocatable :: parameters_plus(:, :), parameters_minus(:, :)
    real(dp), allocatable :: vertices_bar(:, :), parameters_bar(:, :)
    integer, allocatable :: triangles(:, :)
    integer :: basis, status, status_minus, status_plus, vertex
    real(dp) :: lhs, rhs, jvp_error
    logical :: all_passed

    all_passed = .true.
    direction = [1.0_dp, 2.0_dp, -1.0_dp]/sqrt(6.0_dp)
    direction_dot = 0.0_dp
    polarization = cmplx([2.0_dp, -1.0_dp, 0.0_dp]/sqrt(5.0_dp), &
        [0.0_dp, 0.0_dp, 0.0_dp], dp)
    polarization_dot = cmplx([0.004_dp, 0.0_dp, 0.004_dp], &
        [0.003_dp, 0.0_dp, 0.003_dp], dp)
    major_radius_dot = 0.017_dp
    minor_radius_dot = -0.009_dp
    wave_number_dot = 0.031_dp
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 4, 6, vertices, triangles, parameters)
    allocate( &
        vertices_dot(size(vertices, 1), size(vertices, 2)), &
        parameters_dot(size(parameters, 1), size(parameters, 2)), &
        vertices_plus(size(vertices, 1), size(vertices, 2)), &
        vertices_minus(size(vertices, 1), size(vertices, 2)), &
        parameters_plus(size(parameters, 1), size(parameters, 2)), &
        parameters_minus(size(parameters, 1), size(parameters, 2)))
    do vertex = 1, size(vertices, 2)
        vertices_dot(:, vertex) = [ &
            0.002_dp*sin(real(vertex, dp)), &
            -0.001_dp*cos(real(2*vertex, dp)), &
            0.0015_dp*sin(real(3*vertex, dp))]
        parameters_dot(:, vertex) = [ &
            0.002_dp*cos(real(vertex, dp)), &
            -0.001_dp*sin(real(2*vertex, dp))]
    end do
    vertices_plus = vertices + step*vertices_dot
    vertices_minus = vertices - step*vertices_dot
    parameters_plus = parameters + step*parameters_dot
    parameters_minus = parameters - step*parameters_dot

    call assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d_jvp( &
        vertices, triangles, parameters, major_radius, minor_radius, direction, &
        polarization, wave_number, 6, vertices_dot, parameters_dot, &
        major_radius_dot, minor_radius_dot, direction_dot, polarization_dot, &
        wave_number_dot, right_hand_side, right_hand_side_dot, status)
    call assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d( &
        vertices_plus, triangles, parameters_plus, &
        major_radius + step*major_radius_dot, minor_radius + step*minor_radius_dot, &
        direction, polarization + step*polarization_dot, &
        wave_number + step*wave_number_dot, 6, right_hand_side_plus, status_plus)
    call assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d( &
        vertices_minus, triangles, parameters_minus, &
        major_radius - step*major_radius_dot, minor_radius - step*minor_radius_dot, &
        direction, polarization - step*polarization_dot, &
        wave_number - step*wave_number_dot, 6, right_hand_side_minus, status_minus)
    jvp_error = maxval(abs(right_hand_side_dot - &
        (right_hand_side_plus - right_hand_side_minus)/(2.0_dp*step)))
    call record_condition(status == 0, &
        "exact-torus BC plane-wave RHS JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed exact-torus BC plane-wave RHS assembles")
    call record_condition(jvp_error < 3.0e-7_dp, &
        "exact-torus BC plane-wave RHS JVP matches reassembly")

    allocate( &
        right_hand_side_bar(size(right_hand_side)), &
        vertices_bar(size(vertices, 1), size(vertices, 2)), &
        parameters_bar(size(parameters, 1), size(parameters, 2)))
    do basis = 1, size(right_hand_side_bar)
        right_hand_side_bar(basis) = cmplx( &
            sin(real(2*basis, dp)), cos(real(3*basis, dp)), dp)
    end do
    call assemble_maxwell_torus_curved_plane_wave_rhs_bc_3d_vjp( &
        vertices, triangles, parameters, major_radius, minor_radius, direction, &
        polarization, wave_number, 6, right_hand_side_bar, vertices_bar, &
        parameters_bar, major_radius_bar, minor_radius_bar, direction_bar, &
        polarization_bar, wave_number_bar, status)
    lhs = real(sum(conjg(right_hand_side_bar)*right_hand_side_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + sum(parameters_bar*parameters_dot) + &
        major_radius_bar*major_radius_dot + minor_radius_bar*minor_radius_dot + &
        dot_product(direction_bar, direction_dot) + &
        real(sum(conjg(polarization_bar)*polarization_dot), dp) + &
        wave_number_bar*wave_number_dot
    call record_condition(status == 0, &
        "exact-torus BC plane-wave RHS VJP succeeds")
    call record_condition(abs(lhs - rhs) < &
        5.0e-9_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "exact-torus BC plane-wave RHS products obey the adjoint identity")

    call check_summary("Differentiable exact-curved torus BC plane-wave RHS")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_plane_wave_rhs_bc_ad
