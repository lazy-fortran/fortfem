program test_maxwell_sphere_curved_mfie_trace
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d, &
        generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: coarse(:, :), fine(:, :), medium(:, :)
    complex(dp), allocatable :: scaled(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: scaled_vertices(:, :), vertices(:, :)
    real(dp) :: coarse_change, error, fine_change
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d( &
        vertices, triangles, 1.0_dp, 0.75_dp, 8, 0.20_dp, coarse, status)
    call assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d( &
        vertices, triangles, 1.0_dp, 0.75_dp, 8, 0.10_dp, medium, status)
    call assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d( &
        vertices, triangles, 1.0_dp, 0.75_dp, 8, 0.05_dp, fine, status)
    call record_condition(status == 0 .and. maxval(abs(fine)) > 0.0_dp, &
        "curved RBC-tested exterior magnetic trace assembles")
    coarse_change = maxval(abs(medium - coarse))
    fine_change = maxval(abs(fine - medium))
    call record_condition(fine_change < coarse_change, &
        "curved exterior magnetic trace converges toward the boundary")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d( &
        scaled_vertices, triangles, 2.0_dp, 0.375_dp, 8, 0.05_dp, scaled, &
        status)
    error = maxval(abs(scaled - 2.0_dp*fine))/maxval(abs(scaled))
    call record_condition(status == 0 .and. error < 3.0e-12_dp, &
        "curved magnetic trace obeys analytical BC length scaling")

    call check_summary("Curved-sphere exterior Maxwell MFIE trace")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_mfie_trace
