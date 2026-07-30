program test_maxwell_torus_curved_wave
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d, &
        evaluate_maxwell_torus_curved_far_field_rwg_3d, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: impedance = 1.7_dp, wave_number = 0.45_dp
    complex(dp), allocatable :: coefficients(:), right_hand_side(:)
    complex(dp), allocatable :: scaled_right_hand_side(:)
    complex(dp) :: far_field(3), polarization(3), received, transmitted
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    real(dp), allocatable :: scaled_vertices(:, :)
    real(dp) :: direction(3), relative_scaling_error
    integer :: basis, status
    logical :: all_passed

    all_passed = .true.
    direction = [1.0_dp, 2.0_dp, -1.0_dp]/sqrt(6.0_dp)
    polarization = cmplx([2.0_dp, -1.0_dp, 0.0_dp]/sqrt(5.0_dp), 0.0_dp, dp)
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 5, 6, vertices, triangles, parameters)
    call assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        -direction, polarization, wave_number, 12, right_hand_side, status)
    call record_condition(status == 0, &
        "exact-torus plane-wave RWG trace load assembles")
    allocate(coefficients(size(right_hand_side)))
    do basis = 1, size(coefficients)
        coefficients(basis) = cmplx( &
            cos(real(2*basis, dp)), sin(real(3*basis, dp)), dp)
    end do
    call evaluate_maxwell_torus_curved_far_field_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        coefficients, direction, wave_number, impedance, 12, far_field, status)
    received = sum(polarization*far_field)
    transmitted = -cmplx( &
        0.0_dp, wave_number*impedance/(4.0_dp*acos(-1.0_dp)), dp)* &
        sum(coefficients*right_hand_side)
    call record_condition(status == 0 .and. abs(received - transmitted) < &
        2.0e-13_dp*max(1.0_dp, abs(transmitted)), &
        "exact-torus RWG radiation obeys receive-transmit reciprocity")

    call assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        -direction, polarization, 0.0_dp, 12, right_hand_side, status)
    scaled_vertices = 3.0_dp*vertices
    call assemble_maxwell_torus_curved_plane_wave_rhs_rwg_3d( &
        scaled_vertices, triangles, parameters, 3.0_dp*major_radius, &
        3.0_dp*minor_radius, -direction, polarization, 0.0_dp, 12, &
        scaled_right_hand_side, status)
    relative_scaling_error = &
        maxval(abs(scaled_right_hand_side - 9.0_dp*right_hand_side))/ &
        maxval(abs(scaled_right_hand_side))
    call record_condition(status == 0 .and. relative_scaling_error < &
        5.0e-14_dp, "constant-field trace load obeys analytical area scaling")

    call check_summary("Exact-curved torus Maxwell wave traces")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_wave
