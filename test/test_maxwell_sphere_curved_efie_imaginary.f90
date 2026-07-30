program test_maxwell_sphere_curved_efie_imaginary
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d, &
        assemble_maxwell_sphere_curved_potential_operators_rwg_3d, &
        generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: scaled_vertices(:, :), vertices(:, :)
    complex(dp), allocatable :: direct(:, :), matrix(:, :), scalar(:, :)
    complex(dp), allocatable :: scaled(:, :), vector(:, :)
    real(dp) :: decay_rate, error, impedance
    integer :: status
    logical :: all_passed

    all_passed = .true.
    decay_rate = 0.8_dp
    impedance = 1.7_dp
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d( &
        vertices, triangles, 1.0_dp, decay_rate, impedance, 4, 1.0e-5_dp, 2, &
        matrix, status)
    call record_condition(status == 0 .and. &
        maxval(abs(aimag(matrix))) < 3.0e-14_dp, &
        "curved imaginary-wave EFIE is real")
    call assemble_maxwell_sphere_curved_potential_operators_rwg_3d( &
        vertices, triangles, 1.0_dp, decay_rate, 4, 1.0e-5_dp, 2, vector, &
        scalar, status, .true.)
    direct = -impedance*(decay_rate*vector + scalar/decay_rate)
    call record_condition(maxval(abs(matrix - direct)) < &
        3.0e-13_dp*max(1.0_dp, maxval(abs(direct))), &
        "curved imaginary EFIE has the decaying-kernel decomposition")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call assemble_maxwell_sphere_curved_efie_imaginary_rwg_3d( &
        scaled_vertices, triangles, 2.0_dp, decay_rate/2.0_dp, impedance, 4, &
        1.0e-5_dp, 2, scaled, status)
    error = maxval(abs(scaled - 4.0_dp*matrix))/maxval(abs(scaled))
    call record_condition(status == 0 .and. error < 5.0e-11_dp, &
        "curved imaginary EFIE obeys analytical similarity scaling")

    call check_summary("Curved-sphere imaginary-wave Maxwell EFIE")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_efie_imaginary
