program test_maxwell_torus_curved_efie
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_torus_curved_efie_rwg_3d, &
        assemble_maxwell_torus_curved_potential_operators_rwg_3d, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: impedance = 1.7_dp, wave_number = 0.6_dp
    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    complex(dp), allocatable :: efie(:, :), scalar(:, :), vector(:, :)
    complex(dp), allocatable :: scaled_efie(:, :), scaled_scalar(:, :)
    complex(dp), allocatable :: scaled_vector(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    real(dp), allocatable :: scaled_vertices(:, :)
    real(dp) :: decomposition_error, scalar_scaling_error, scaling_error
    real(dp) :: symmetry_error, vector_scaling_error
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 4, vertices, triangles, parameters)
    call assemble_maxwell_torus_curved_potential_operators_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number, 3, 1.0e-5_dp, 1, vector, scalar, status)
    call assemble_maxwell_torus_curved_efie_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        wave_number, impedance, 3, 1.0e-5_dp, 1, efie, status)
    symmetry_error = max( &
        maxval(abs(vector - transpose(vector))), &
        maxval(abs(scalar - transpose(scalar))))
    decomposition_error = maxval(abs(efie - &
        cmplx(0.0_dp, wave_number*impedance, dp)*vector + &
        cmplx(0.0_dp, impedance/wave_number, dp)*scalar))
    call record_condition(status == 0 .and. symmetry_error < 3.0e-13_dp, &
        "exact-torus Maxwell potential operators are reciprocal")
    call record_condition(decomposition_error < &
        3.0e-13_dp*max(1.0_dp, maxval(abs(efie))), &
        "exact-torus EFIE has the Maxwell vector/scalar decomposition")

    scaled_vertices = 2.0_dp*vertices
    call assemble_maxwell_torus_curved_potential_operators_rwg_3d( &
        scaled_vertices, triangles, parameters, 2.0_dp*major_radius, &
        2.0_dp*minor_radius, wave_number/2.0_dp, 3, 1.0e-5_dp, 1, &
        scaled_vector, scaled_scalar, status)
    call assemble_maxwell_torus_curved_efie_rwg_3d( &
        scaled_vertices, triangles, parameters, 2.0_dp*major_radius, &
        2.0_dp*minor_radius, wave_number/2.0_dp, impedance, 3, 1.0e-5_dp, 1, &
        scaled_efie, status)
    vector_scaling_error = maxval(abs(scaled_vector - 8.0_dp*vector))/ &
        maxval(abs(scaled_vector))
    scalar_scaling_error = maxval(abs(scaled_scalar - 2.0_dp*scalar))/ &
        maxval(abs(scaled_scalar))
    scaling_error = maxval(abs(scaled_efie - 4.0_dp*efie))/ &
        maxval(abs(scaled_efie))
    call record_condition(status == 0 .and. vector_scaling_error < &
        2.0e-13_dp .and. scalar_scaling_error < 2.0e-13_dp, &
        "exact-torus potentials obey analytical similarity scaling")
    call record_condition(scaling_error < 2.0e-13_dp, &
        "exact-torus EFIE obeys analytical area scaling")

    call check_summary("Exact-curved torus Maxwell EFIE")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_efie
