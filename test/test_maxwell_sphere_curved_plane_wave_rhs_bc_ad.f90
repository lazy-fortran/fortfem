program test_maxwell_sphere_curved_plane_wave_rhs_bc_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d, &
        assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d_jvp, &
        assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d_vjp, &
        generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: radius = 1.3_dp, wave_number = 0.43_dp
    real(dp), parameter :: step = 2.0e-6_dp
    real(dp) :: radius_dot, wave_number_dot, radius_bar, wave_number_bar
    real(dp) :: direction(3), direction_dot(3), direction_bar(3)
    complex(dp) :: polarization(3), polarization_dot(3)
    complex(dp), allocatable :: polarization_bar(:)
    complex(dp), allocatable :: right_hand_side(:), right_hand_side_dot(:)
    complex(dp), allocatable :: right_hand_side_minus(:), right_hand_side_plus(:)
    complex(dp), allocatable :: right_hand_side_bar(:)
    real(dp), allocatable :: vertices(:, :), vertices_dot(:, :), vertices_bar(:, :)
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
    radius_dot = -0.004_dp
    wave_number_dot = 0.031_dp
    call generate_sphere_surface_mesh(radius, 0, vertices, triangles)
    allocate(vertices_dot(size(vertices, 1), size(vertices, 2)))
    do vertex = 1, size(vertices, 2)
        vertices_dot(:, vertex) = [ &
            0.001_dp*sin(real(vertex, dp)), &
            -0.0015_dp*cos(real(2*vertex, dp)), &
            0.0008_dp*sin(real(3*vertex, dp))]
    end do

    call assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d_jvp( &
        vertices, triangles, radius, direction, polarization, wave_number, 6, &
        vertices_dot, radius_dot, direction_dot, polarization_dot, wave_number_dot, &
        right_hand_side, right_hand_side_dot, status)
    call assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d( &
        vertices + step*vertices_dot, triangles, radius + step*radius_dot, &
        direction, polarization + step*polarization_dot, &
        wave_number + step*wave_number_dot, 6, right_hand_side_plus, status_plus)
    call assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d( &
        vertices - step*vertices_dot, triangles, radius - step*radius_dot, &
        direction, polarization - step*polarization_dot, &
        wave_number - step*wave_number_dot, 6, right_hand_side_minus, status_minus)
    jvp_error = maxval(abs(right_hand_side_dot - &
        (right_hand_side_plus - right_hand_side_minus)/(2.0_dp*step)))
    call record_condition(status == 0, &
        "exact-sphere BC plane-wave RHS JVP succeeds")
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed exact-sphere BC plane-wave RHS assembles")
    call record_condition(jvp_error < 3.0e-7_dp, &
        "exact-sphere BC plane-wave RHS JVP matches reassembly")

    allocate( &
        right_hand_side_bar(size(right_hand_side)), &
        vertices_bar(size(vertices, 1), size(vertices, 2)))
    do basis = 1, size(right_hand_side_bar)
        right_hand_side_bar(basis) = cmplx( &
            sin(real(2*basis, dp)), cos(real(3*basis, dp)), dp)
    end do
    call assemble_maxwell_sphere_curved_plane_wave_rhs_bc_3d_vjp( &
        vertices, triangles, radius, direction, polarization, wave_number, 6, &
        right_hand_side_bar, vertices_bar, radius_bar, direction_bar, &
        polarization_bar, wave_number_bar, status)
    lhs = real(sum(conjg(right_hand_side_bar)*right_hand_side_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + radius_bar*radius_dot + &
        dot_product(direction_bar, direction_dot) + &
        real(sum(conjg(polarization_bar)*polarization_dot), dp) + &
        wave_number_bar*wave_number_dot
    call record_condition(status == 0, &
        "exact-sphere BC plane-wave RHS VJP succeeds")
    call record_condition(abs(lhs - rhs) < &
        5.0e-9_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "exact-sphere BC plane-wave RHS products obey the adjoint identity")

    call check_summary("Differentiable exact-curved sphere BC plane-wave RHS")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_plane_wave_rhs_bc_ad
