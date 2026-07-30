program test_maxwell_sphere_curved_wave
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d, &
        evaluate_maxwell_sphere_curved_far_field_rwg_3d, &
        generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: coefficients(:), right_hand_side(:)
    complex(dp), allocatable :: scaled_right_hand_side(:)
    complex(dp) :: far_field(3), polarization(3), scaled_far_field(3)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: scaled_vertices(:, :), vertices(:, :)
    real(dp) :: direction(3), error
    integer :: status
    logical :: all_passed

    all_passed = .true.
    direction = [0.0_dp, 0.0_dp, 1.0_dp]
    polarization = cmplx([1.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, dp)
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d( &
        vertices, triangles, 1.0_dp, direction, polarization, 1.3_dp, 12, &
        right_hand_side, status)
    call record_condition(status == 0, "curved incident trace assembles")
    allocate(coefficients(size(right_hand_side)))
    coefficients = right_hand_side*cmplx(1.0_dp, 0.5_dp, dp)
    call evaluate_maxwell_sphere_curved_far_field_rwg_3d( &
        vertices, triangles, 1.0_dp, coefficients, direction, 1.3_dp, &
        2.0_dp, 12, far_field, status)
    call record_condition(status == 0, "curved far field evaluates")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d( &
        scaled_vertices, triangles, 2.0_dp, direction, polarization, 0.65_dp, &
        12, scaled_right_hand_side, status)
    error = maxval(abs(scaled_right_hand_side - 4.0_dp*right_hand_side))/ &
        maxval(abs(scaled_right_hand_side))
    call record_condition(status == 0 .and. error < 5.0e-14_dp, &
        "curved incident trace obeys electromagnetic similarity scaling")
    call evaluate_maxwell_sphere_curved_far_field_rwg_3d( &
        scaled_vertices, triangles, 2.0_dp, coefficients, direction, 0.65_dp, &
        2.0_dp, 12, scaled_far_field, status)
    error = maxval(abs(scaled_far_field - 2.0_dp*far_field))/ &
        maxval(abs(scaled_far_field))
    call record_condition(status == 0 .and. error < 5.0e-14_dp, &
        "curved far field obeys electromagnetic similarity scaling")

    call check_summary("Curved-sphere incident and far fields")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_wave
