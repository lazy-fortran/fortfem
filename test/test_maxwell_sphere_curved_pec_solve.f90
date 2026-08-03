program test_maxwell_sphere_curved_pec_solve
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_maxwell_sphere_curved_efie_rwg_3d, &
        assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d, &
        solve_maxwell_pec_sphere_curved_efie_rwg_3d
    use fortfem_core, only: generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: density(:), matrix(:, :), right_hand_side(:)
    complex(dp), allocatable :: scaled_density(:)
    complex(dp) :: polarization(3)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: scaled_vertices(:, :), vertices(:, :)
    real(dp) :: direction(3), error, residual
    integer :: status
    logical :: all_passed

    all_passed = .true.
    direction = [0.0_dp, 0.0_dp, 1.0_dp]
    polarization = cmplx([1.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, dp)
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call solve_maxwell_pec_sphere_curved_efie_rwg_3d( &
        vertices, triangles, 1.0_dp, direction, polarization, 0.65_dp, 1.8_dp, &
        4, 1.0e-5_dp, 2, density, status)
    call record_condition(status == 0 .and. size(density) > 0, &
        "curved-sphere PEC EFIE solve succeeds")
    call assemble_maxwell_sphere_curved_efie_rwg_3d( &
        vertices, triangles, 1.0_dp, 0.65_dp, 1.8_dp, 4, 1.0e-5_dp, 2, matrix, &
        status)
    call assemble_maxwell_sphere_curved_plane_wave_rhs_rwg_3d( &
        vertices, triangles, 1.0_dp, direction, polarization, 0.65_dp, 4, &
        right_hand_side, status)
    residual = sqrt(sum(abs(matmul(matrix, density) - right_hand_side)**2))/ &
        sqrt(sum(abs(right_hand_side)**2))
    call record_condition(status == 0 .and. residual < 3.0e-12_dp, &
        "curved PEC density satisfies the independently assembled EFIE")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call solve_maxwell_pec_sphere_curved_efie_rwg_3d( &
        scaled_vertices, triangles, 2.0_dp, direction, polarization, 0.325_dp, &
        1.8_dp, 4, 1.0e-5_dp, 2, scaled_density, status)
    error = maxval(abs(scaled_density - density))/maxval(abs(density))
    call record_condition(status == 0 .and. error < 8.0e-11_dp, &
        "curved PEC RWG coefficients obey electromagnetic similarity")

    call check_summary("Curved-sphere PEC Maxwell solve")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_pec_solve
