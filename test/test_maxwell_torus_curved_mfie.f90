program test_maxwell_torus_curved_mfie
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d, &
        assemble_maxwell_torus_curved_rwg_rbc_pairing, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    complex(dp), allocatable :: exterior(:, :), interior(:, :)
    complex(dp), allocatable :: exterior_fine(:, :), interior_fine(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: pairing(:, :), parameters(:, :), vertices(:, :)
    real(dp) :: coarse_error, fine_error
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 3, vertices, triangles, parameters)
    call assemble_maxwell_torus_curved_rwg_rbc_pairing( &
        vertices, triangles, parameters, major_radius, minor_radius, 5, &
        pairing, status)
    call assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, 0.4_dp, &
        5, 0.12_dp, exterior, status)
    call assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, 0.4_dp, &
        5, -0.12_dp, interior, status)
    call assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, 0.4_dp, &
        12, 0.06_dp, exterior_fine, status)
    call assemble_maxwell_torus_curved_mfie_offset_trace_rwg_rbc_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, 0.4_dp, &
        12, -0.06_dp, interior_fine, status)
    coarse_error = jump_error(exterior, interior, pairing)
    fine_error = jump_error(exterior_fine, interior_fine, pairing)
    call record_condition(status == 0 .and. fine_error < coarse_error .and. &
        fine_error < 0.35_dp, &
        "curved torus magnetic traces converge to the Maxwell jump relation")

    call check_summary("Exact-curved torus Maxwell MFIE traces")
    if (.not. all_passed) error stop 1

contains

    pure function jump_error(exterior_trace, interior_trace, mass) result(error)
        complex(dp), intent(in) :: exterior_trace(:, :), interior_trace(:, :)
        real(dp), intent(in) :: mass(:, :)
        real(dp) :: error

        error = min( &
            maxval(abs( &
            exterior_trace - interior_trace - cmplx(mass, 0.0_dp, dp))), &
            maxval(abs( &
            exterior_trace - interior_trace + cmplx(mass, 0.0_dp, dp))))/ &
            maxval(abs(mass))
    end function jump_error

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_mfie
