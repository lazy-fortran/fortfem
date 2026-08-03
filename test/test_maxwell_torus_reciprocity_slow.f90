program test_maxwell_torus_reciprocity_slow
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: evaluate_maxwell_efie_far_field_rwg_3d, &
        solve_maxwell_pec_regularized_cfie_rwg_multiple_3d
    use fortfem_core, only: generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: wave_number = 0.45_dp, impedance = 1.7_dp
    complex(dp), allocatable :: currents(:, :)
    complex(dp) :: first_far_field(3), second_far_field(3)
    complex(dp), parameter :: z_polarization(3) = [ &
        cmplx(0.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp), &
        cmplx(1.0_dp, 0.0_dp, dp)]
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    real(dp) :: reciprocity_error
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        2.0_dp, 0.6_dp, 4, 6, vertices, triangles, parameters)
    call solve_maxwell_pec_regularized_cfie_rwg_multiple_3d( &
        vertices, triangles, reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp], [3, 2]), &
        reshape([z_polarization, z_polarization], [3, 2]), &
        wave_number, impedance, 3, 2.0e-3_dp, 1, currents, status)
    if (status /= 0) error stop "batched toroidal Maxwell CFIE solve failed"
    call evaluate_maxwell_efie_far_field_rwg_3d( &
        vertices, triangles, currents(:, 1), [0.0_dp, 1.0_dp, 0.0_dp], &
        wave_number, impedance, 4, first_far_field, status)
    if (status /= 0) error stop "first toroidal Maxwell far field failed"
    call evaluate_maxwell_efie_far_field_rwg_3d( &
        vertices, triangles, currents(:, 2), [-1.0_dp, 0.0_dp, 0.0_dp], &
        wave_number, impedance, 4, second_far_field, status)
    if (status /= 0) error stop "second toroidal Maxwell far field failed"

    reciprocity_error = abs(first_far_field(3) - second_far_field(3))/ &
        max(abs(first_far_field(3)), abs(second_far_field(3)))
    write (*, '(a,es12.4)') "toroidal Maxwell reciprocity error: ", &
        reciprocity_error
    call record_condition( &
        abs(first_far_field(3)) > 1.0e-5_dp .and. &
        abs(second_far_field(3)) > 1.0e-5_dp, &
        "toroidal PEC produces nontrivial reciprocal scattering amplitudes")
    call record_condition(reciprocity_error < 8.0e-2_dp, &
        "toroidal Maxwell CFIE satisfies Lorentz far-field reciprocity")
    call check_summary("Toroidal Maxwell CFIE reciprocity")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_reciprocity_slow
