program test_maxwell_torus_curved_reciprocity_slow
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_maxwell_torus_curved_far_field_rwg_3d, &
        generate_torus_surface_mesh, &
        solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 0.45_dp, impedance = 1.7_dp
    complex(dp), allocatable :: currents(:, :)
    complex(dp) :: first_far_field(3), second_far_field(3)
    complex(dp), parameter :: z_polarization(3) = [ &
        cmplx(0.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp), &
        cmplx(1.0_dp, 0.0_dp, dp)]
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    real(dp) :: amplitude, reciprocity_error
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 3, vertices, triangles, parameters)
    call solve_maxwell_pec_torus_curved_regularized_cfie_rwg_multiple_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp], [3, 2]), &
        reshape([z_polarization, z_polarization], [3, 2]), wave_number, &
        impedance, 6, 3.0e-4_dp, 1, 0.12_dp, currents, status)
    if (status /= 0) error stop "exact-torus batched Maxwell CFIE solve failed"
    call evaluate_maxwell_torus_curved_far_field_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        currents(:, 1), [0.0_dp, 1.0_dp, 0.0_dp], wave_number, impedance, 8, &
        first_far_field, status)
    if (status /= 0) error stop "first exact-torus far field failed"
    call evaluate_maxwell_torus_curved_far_field_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        currents(:, 2), [-1.0_dp, 0.0_dp, 0.0_dp], wave_number, impedance, 8, &
        second_far_field, status)
    if (status /= 0) error stop "second exact-torus far field failed"
    amplitude = max(abs(first_far_field(3)), abs(second_far_field(3)))
    reciprocity_error = abs(first_far_field(3) - second_far_field(3))/amplitude
    call record_condition(amplitude > 1.0e-5_dp, &
        "exact-torus scattered far field is nontrivial")
    call record_condition(reciprocity_error < 0.2_dp, &
        "exact-torus CFIE scattering obeys Lorentz reciprocity")

    call check_summary("Exact-curved torus Maxwell CFIE reciprocity")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_reciprocity_slow
