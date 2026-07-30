program test_maxwell_sphere_curved_coincident
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_maxwell_rwg_surface_space, generate_sphere_surface_mesh, &
        integrate_maxwell_sphere_curved_coincident_rwg_pair_3d
    use fortfem_kinds, only: dp
    implicit none

    complex(dp) :: reciprocal, reference, scaled, values(3)
    integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: scaled_vertices(:, :), vertices(:, :)
    real(dp) :: errors(2), scaling_error
    integer :: basis, candidate, other_basis, panel, status
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call build_maxwell_rwg_surface_space( &
        vertices, triangles, edge_vertices, edge_triangles, status)
    basis = 1
    panel = edge_triangles(1, basis)
    other_basis = 0
    do candidate = 1, size(edge_triangles, 2)
        if (candidate /= basis .and. &
            any(edge_triangles(:, candidate) == panel)) then
            other_basis = candidate
            exit
        end if
    end do
    call integrate_maxwell_sphere_curved_coincident_rwg_pair_3d( &
        vertices, triangles, edge_vertices, edge_triangles, basis, panel, &
        basis, 1.0_dp, 0.9_dp, 48, reference, status)
    call record_condition(status == 0 .and. abs(reference) > 0.0_dp, &
        "curved coincident RWG product integral is finite")
    call integrate_maxwell_sphere_curved_coincident_rwg_pair_3d( &
        vertices, triangles, edge_vertices, edge_triangles, basis, panel, &
        other_basis, 1.0_dp, 0.9_dp, 24, values(1), status)
    call integrate_maxwell_sphere_curved_coincident_rwg_pair_3d( &
        vertices, triangles, edge_vertices, edge_triangles, other_basis, panel, &
        basis, 1.0_dp, 0.9_dp, 24, reciprocal, status)
    call record_condition(other_basis > 0 .and. status == 0 .and. &
        abs(values(1) - reciprocal) < 2.0e-12_dp*max(1.0_dp, &
        abs(reciprocal)), "curved coincident kernel obeys RWG reciprocity")
    call integrate_maxwell_sphere_curved_coincident_rwg_pair_3d( &
        vertices, triangles, edge_vertices, edge_triangles, basis, panel, &
        basis, 1.0_dp, 0.9_dp, 8, values(1), status)
    call integrate_maxwell_sphere_curved_coincident_rwg_pair_3d( &
        vertices, triangles, edge_vertices, edge_triangles, basis, panel, &
        basis, 1.0_dp, 0.9_dp, 16, values(2), status)
    call integrate_maxwell_sphere_curved_coincident_rwg_pair_3d( &
        vertices, triangles, edge_vertices, edge_triangles, basis, panel, &
        basis, 1.0_dp, 0.9_dp, 32, values(3), status)
    errors = abs(values(1:2) - reference)
    if (errors(2) >= 0.03_dp*errors(1) .or. &
        abs(values(3) - reference) >= 0.7_dp*errors(2)) &
        write (*, *) "curved self errors", errors, abs(values(3) - reference)
    call record_condition(errors(2) < 0.03_dp*errors(1) .and. &
        abs(values(3) - reference) < 0.7_dp*errors(2), &
        "radial-Duffy curved self integral converges under refinement")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call integrate_maxwell_sphere_curved_coincident_rwg_pair_3d( &
        scaled_vertices, triangles, edge_vertices, edge_triangles, basis, &
        panel, basis, 2.0_dp, 0.45_dp, 32, scaled, status)
    scaling_error = abs(scaled - 8.0_dp*values(3))/abs(scaled)
    call record_condition(status == 0 .and. scaling_error < 8.0e-13_dp, &
        "curved coincident Helmholtz product obeys analytical R-cubed scaling")

    call check_summary("Curved-sphere coincident RWG product")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_coincident
