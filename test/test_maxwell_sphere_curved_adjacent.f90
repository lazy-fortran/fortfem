program test_maxwell_sphere_curved_adjacent
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_maxwell_rwg_surface_space, generate_sphere_surface_mesh, &
        integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d
    use fortfem_kinds, only: dp
    implicit none

    complex(dp) :: reciprocal, reference, scaled, values(2)
    integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: scaled_vertices(:, :), vertices(:, :)
    real(dp) :: errors(2), scaling_error
    integer :: basis, first_panel, second_panel, status
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call build_maxwell_rwg_surface_space( &
        vertices, triangles, edge_vertices, edge_triangles, status)
    basis = 1
    first_panel = edge_triangles(1, basis)
    second_panel = edge_triangles(2, basis)
    call integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d( &
        vertices, triangles, edge_vertices, edge_triangles, basis, &
        first_panel, basis, second_panel, 1.0_dp, 0.8_dp, 6, 1.0e-10_dp, 5, &
        reference, status)
    call record_condition(status == 0 .and. abs(reference) > 0.0_dp, &
        "adaptive curved adjacent RWG integral is finite")
    call integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d( &
        vertices, triangles, edge_vertices, edge_triangles, basis, &
        first_panel, basis, second_panel, 1.0_dp, 0.8_dp, 6, 1.0e-3_dp, 1, &
        values(1), status)
    call integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d( &
        vertices, triangles, edge_vertices, edge_triangles, basis, &
        first_panel, basis, second_panel, 1.0_dp, 0.8_dp, 6, 1.0e-7_dp, 4, &
        values(2), status)
    errors = abs(values - reference)
    call record_condition(errors(2) < 0.1_dp*errors(1), &
        "curved adjacent integral converges with adaptive refinement")

    call integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d( &
        vertices, triangles, edge_vertices, edge_triangles, basis, &
        second_panel, basis, first_panel, 1.0_dp, 0.8_dp, 6, 1.0e-10_dp, 5, &
        reciprocal, status)
    call record_condition(status == 0 .and. &
        abs(reciprocal - reference) < 2.0e-11_dp*max(1.0_dp, abs(reference)), &
        "adaptive curved adjacent kernel obeys reciprocity")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d( &
        scaled_vertices, triangles, edge_vertices, edge_triangles, basis, &
        first_panel, basis, second_panel, 2.0_dp, 0.4_dp, 6, 1.0e-10_dp, 5, &
        scaled, status)
    scaling_error = abs(scaled - 8.0_dp*reference)/abs(scaled)
    call record_condition(status == 0 .and. scaling_error < 2.0e-11_dp, &
        "curved adjacent Helmholtz product obeys analytical R-cubed scaling")

    call check_summary("Curved-sphere adjacent RWG product")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_adjacent
