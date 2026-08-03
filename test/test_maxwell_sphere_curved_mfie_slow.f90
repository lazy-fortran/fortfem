program test_maxwell_sphere_curved_mfie_slow
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d, &
        assemble_maxwell_sphere_curved_mfie_rwg_rbc_3d
    use fortfem_core, only: generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: matrix(:, :), refined(:, :), scaled(:, :)
    complex(dp), allocatable :: trace_full(:, :), trace_half(:, :)
    complex(dp), allocatable :: trace_quarter(:, :)
    complex(dp), allocatable :: direct(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: scaled_vertices(:, :), vertices(:, :)
    real(dp) :: error
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call assemble_maxwell_sphere_curved_mfie_rwg_rbc_3d( &
        vertices, triangles, 1.0_dp, 0.7_dp, 12, 0.16_dp, matrix, status)
    call record_condition(status == 0 .and. maxval(abs(matrix)) > 0.0_dp, &
        "quadratically extrapolated curved MFIE assembles")
    call assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d( &
        vertices, triangles, 1.0_dp, 0.7_dp, 12, 0.16_dp, trace_full, status)
    call assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d( &
        vertices, triangles, 1.0_dp, 0.7_dp, 12, 0.08_dp, trace_half, status)
    call assemble_maxwell_sphere_curved_mfie_exterior_trace_rwg_rbc_3d( &
        vertices, triangles, 1.0_dp, 0.7_dp, 12, 0.04_dp, trace_quarter, status)
    direct = trace_full/3.0_dp - 2.0_dp*trace_half + &
        8.0_dp*trace_quarter/3.0_dp
    call record_condition(maxval(abs(matrix - direct)) < &
        3.0e-13_dp*max(1.0_dp, maxval(abs(direct))), &
        "curved MFIE matches independent zero-offset extrapolation")
    call assemble_maxwell_sphere_curved_mfie_rwg_rbc_3d( &
        vertices, triangles, 1.0_dp, 0.7_dp, 12, 0.08_dp, refined, status)
    call record_condition(maxval(abs(refined - matrix)) < &
        maxval(abs(trace_half - trace_full)), &
        "zero-offset extrapolation improves the raw curved trace sequence")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call assemble_maxwell_sphere_curved_mfie_rwg_rbc_3d( &
        scaled_vertices, triangles, 2.0_dp, 0.35_dp, 12, 0.16_dp, scaled, &
        status)
    error = maxval(abs(scaled - 2.0_dp*matrix))/maxval(abs(scaled))
    call record_condition(status == 0 .and. error < 4.0e-12_dp, &
        "extrapolated curved MFIE obeys analytical BC length scaling")

    call check_summary("Curved-sphere Maxwell MFIE")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_mfie_slow
