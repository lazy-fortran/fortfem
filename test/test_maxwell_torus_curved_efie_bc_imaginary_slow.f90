program test_maxwell_torus_curved_efie_bc_imaginary_slow
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_torus_curved_efie_bc_imaginary_3d, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: decay_rate = 0.5_dp, impedance = 1.4_dp
    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    complex(dp), allocatable :: matrix(:, :), scaled(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    real(dp), allocatable :: scaled_vertices(:, :)
    real(dp) :: scaling_error, symmetry_error
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 3, vertices, triangles, parameters)
    call assemble_maxwell_torus_curved_efie_bc_imaginary_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, &
        decay_rate, impedance, 1, 3.0e-4_dp, 1, matrix, status)
    symmetry_error = maxval(abs(matrix - transpose(matrix)))
    call record_condition(status == 0 .and. symmetry_error < 3.0e-13_dp .and. &
        maxval(abs(aimag(matrix))) < 3.0e-14_dp .and. &
        maxval(abs(matrix)) > 0.0_dp, &
        "exact-torus BC Yukawa EFIE is real and reciprocal")

    scaled_vertices = 2.0_dp*vertices
    call assemble_maxwell_torus_curved_efie_bc_imaginary_3d( &
        scaled_vertices, triangles, parameters, 2.0_dp*major_radius, &
        2.0_dp*minor_radius, decay_rate/2.0_dp, impedance, 1, 3.0e-4_dp, 1, &
        scaled, status)
    scaling_error = maxval(abs(scaled - matrix))/maxval(abs(scaled))
    call record_condition(status == 0 .and. scaling_error < 3.0e-12_dp, &
        "exact-torus BC Yukawa regularizer is similarity invariant")

    call check_summary("Exact-curved torus BC Yukawa EFIE")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_efie_bc_imaginary_slow
