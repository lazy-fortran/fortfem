program test_maxwell_torus_curved_mfie_limit
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_maxwell_torus_curved_mfie_rwg_rbc_3d
    use fortfem_core, only: &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    complex(dp), allocatable :: matrix(:, :), scaled(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    real(dp), allocatable :: scaled_vertices(:, :)
    real(dp) :: scaling_error
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 3, vertices, triangles, parameters)
    call assemble_maxwell_torus_curved_mfie_rwg_rbc_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, 0.4_dp, &
        4, 0.16_dp, matrix, status)
    scaled_vertices = 2.5_dp*vertices
    call assemble_maxwell_torus_curved_mfie_rwg_rbc_3d( &
        scaled_vertices, triangles, parameters, 2.5_dp*major_radius, &
        2.5_dp*minor_radius, 0.4_dp/2.5_dp, 4, 0.16_dp, scaled, status)
    scaling_error = maxval(abs(scaled - 2.5_dp*matrix))/maxval(abs(scaled))
    call record_condition(status == 0 .and. maxval(abs(matrix)) > 0.0_dp .and. &
        scaling_error < 8.0e-13_dp, &
        "extrapolated torus MFIE obeys analytical BC length scaling")

    call check_summary("Exact-curved torus Maxwell MFIE limit")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_mfie_limit
