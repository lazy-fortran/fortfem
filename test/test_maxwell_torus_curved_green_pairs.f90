program test_maxwell_torus_curved_green_pairs
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        integrate_maxwell_torus_curved_adjacent_rwg_pair_3d, &
        integrate_maxwell_torus_curved_coincident_rwg_pair_3d
    use fortfem_core, only: generate_torus_surface_mesh
    use fortfem_feec, only: build_maxwell_rwg_surface_space
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    complex(dp) :: adjacent, adjacent_transpose, coincident
    complex(dp) :: reference, scalar_coincident, scalar_scaled, scaled
    complex(dp) :: refined(2)
    integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), scaled_vertices(:, :)
    real(dp), allocatable :: vertices(:, :)
    real(dp) :: errors(2), scalar_scaling_error, scaling_error
    integer :: basis, panel, status
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 4, 6, vertices, triangles, parameters)
    call build_maxwell_rwg_surface_space( &
        vertices, triangles, edge_vertices, edge_triangles, status)
    basis = 1
    panel = edge_triangles(1, basis)
    call integrate_maxwell_torus_curved_coincident_rwg_pair_3d( &
        vertices, triangles, parameters, edge_vertices, edge_triangles, &
        basis, panel, basis, major_radius, minor_radius, 0.0_dp, 8, &
        coincident, status, scalar_coincident)
    call integrate_maxwell_torus_curved_coincident_rwg_pair_3d( &
        vertices, triangles, parameters, edge_vertices, edge_triangles, &
        basis, panel, basis, major_radius, minor_radius, 0.0_dp, 16, &
        refined(1), status)
    call integrate_maxwell_torus_curved_coincident_rwg_pair_3d( &
        vertices, triangles, parameters, edge_vertices, edge_triangles, &
        basis, panel, basis, major_radius, minor_radius, 0.0_dp, 32, &
        refined(2), status)
    call integrate_maxwell_torus_curved_coincident_rwg_pair_3d( &
        vertices, triangles, parameters, edge_vertices, edge_triangles, &
        basis, panel, basis, major_radius, minor_radius, 0.0_dp, 48, &
        reference, status)
    errors = abs([coincident, refined(1)] - reference)
    call record_condition(status == 0 .and. errors(2) < 0.2_dp*errors(1) .and. &
        abs(refined(2) - reference) < 0.7_dp*errors(2), &
        "radial-Duffy torus self-panel integral converges")
    call record_condition(real(reference, dp) > 0.0_dp .and. &
        abs(aimag(reference)) < 1.0e-14_dp, &
        "Laplace torus self-panel energy is real and positive")

    scaled_vertices = 3.0_dp*vertices
    call integrate_maxwell_torus_curved_coincident_rwg_pair_3d( &
        scaled_vertices, triangles, parameters, edge_vertices, edge_triangles, &
        basis, panel, basis, 3.0_dp*major_radius, 3.0_dp*minor_radius, &
        0.0_dp, 8, scaled, status, scalar_scaled)
    scaling_error = abs(scaled - 27.0_dp*coincident)/abs(scaled)
    scalar_scaling_error = &
        abs(scalar_scaled - scalar_coincident/3.0_dp)/abs(scalar_scaled)
    call record_condition(scaling_error < 8.0e-14_dp, &
        "vector Green pair obeys its analytical cubic scale law")
    call record_condition(scalar_scaling_error < 8.0e-14_dp, &
        "reference scalar Green pair obeys its inverse scale law")

    call integrate_maxwell_torus_curved_adjacent_rwg_pair_3d( &
        vertices, triangles, parameters, edge_vertices, edge_triangles, &
        basis, edge_triangles(1, basis), basis, edge_triangles(2, basis), &
        major_radius, minor_radius, 0.7_dp, 6, 1.0e-7_dp, 3, adjacent, status)
    call integrate_maxwell_torus_curved_adjacent_rwg_pair_3d( &
        vertices, triangles, parameters, edge_vertices, edge_triangles, &
        basis, edge_triangles(2, basis), basis, edge_triangles(1, basis), &
        major_radius, minor_radius, 0.7_dp, 6, 1.0e-7_dp, 3, &
        adjacent_transpose, status)
    call record_condition(status == 0 .and. abs(adjacent - &
        adjacent_transpose) < 3.0e-13_dp*max(1.0_dp, abs(adjacent)), &
        "adaptive adjacent-panel Green pairing is reciprocal")

    call check_summary("Exact-curved torus Maxwell Green pairs")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_green_pairs
