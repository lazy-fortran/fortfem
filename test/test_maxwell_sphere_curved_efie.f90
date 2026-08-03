program test_maxwell_sphere_curved_efie
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_maxwell_sphere_curved_efie_rwg_3d, &
        assemble_maxwell_sphere_curved_potential_operators_rwg_3d
    use fortfem_core, only: generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: scaled_vertices(:, :), vertices(:, :)
    complex(dp), allocatable :: direct(:, :), matrix(:, :), scalar(:, :)
    complex(dp), allocatable :: scaled_matrix(:, :), vector(:, :)
    real(dp) :: error, impedance, wave_number
    integer :: status
    logical :: all_passed

    all_passed = .true.
    wave_number = 0.55_dp
    impedance = 2.4_dp
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call assemble_maxwell_sphere_curved_efie_rwg_3d( &
        vertices, triangles, 1.0_dp, wave_number, impedance, 4, 1.0e-5_dp, 2, &
        matrix, status)
    call record_condition(status == 0 .and. maxval(abs(matrix)) > 0.0_dp, &
        "curved-sphere EFIE matrix assembles")
    call assemble_maxwell_sphere_curved_potential_operators_rwg_3d( &
        vertices, triangles, 1.0_dp, wave_number, 4, 1.0e-5_dp, 2, vector, &
        scalar, status)
    direct = cmplx(0.0_dp, wave_number*impedance, dp)*vector - &
        cmplx(0.0_dp, impedance/wave_number, dp)*scalar
    call record_condition(maxval(abs(matrix - direct)) < &
        3.0e-13_dp*max(1.0_dp, maxval(abs(direct))), &
        "curved EFIE has the analytical vector-scalar decomposition")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call assemble_maxwell_sphere_curved_efie_rwg_3d( &
        scaled_vertices, triangles, 2.0_dp, wave_number/2.0_dp, impedance, 4, &
        1.0e-5_dp, 2, scaled_matrix, status)
    error = maxval(abs(scaled_matrix - 4.0_dp*matrix))/ &
        maxval(abs(scaled_matrix))
    call record_condition(status == 0 .and. error < 5.0e-11_dp, &
        "curved EFIE obeys analytical electromagnetic similarity scaling")

    call check_summary("Curved-sphere Maxwell EFIE")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_efie
