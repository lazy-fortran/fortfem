program test_maxwell_sphere_curved_vector_potential
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_sphere_curved_vector_potential_rwg_3d, &
        build_maxwell_rwg_surface_space, generate_sphere_surface_mesh, &
        integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d, &
        integrate_maxwell_sphere_curved_coincident_rwg_pair_3d
    use fortfem_kinds, only: dp
    implicit none

    integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: scaled_vertices(:, :), vertices(:, :)
    complex(dp), allocatable :: matrix(:, :), scaled_matrix(:, :)
    complex(dp) :: direct_entry, panel_entry
    real(dp) :: error, symmetry_error
    integer :: basis, first_slot, panel, second_slot, status
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call assemble_maxwell_sphere_curved_vector_potential_rwg_3d( &
        vertices, triangles, 1.0_dp, 0.7_dp, 4, 1.0e-5_dp, 2, matrix, status)
    call record_condition(status == 0 .and. maxval(abs(matrix)) > 0.0_dp, &
        "global curved RWG vector potential assembles")
    symmetry_error = maxval(abs(matrix - transpose(matrix)))
    call record_condition(symmetry_error < 3.0e-12_dp, &
        "global curved RWG vector potential is complex symmetric")

    call build_maxwell_rwg_surface_space( &
        vertices, triangles, edge_vertices, edge_triangles, status)
    basis = 1
    direct_entry = cmplx(0.0_dp, 0.0_dp, dp)
    do first_slot = 1, 2
        do second_slot = 1, 2
            panel = edge_triangles(first_slot, basis)
            if (panel == edge_triangles(second_slot, basis)) then
                call integrate_maxwell_sphere_curved_coincident_rwg_pair_3d( &
                    vertices, triangles, edge_vertices, edge_triangles, basis, &
                    panel, basis, 1.0_dp, 0.7_dp, 4, panel_entry, status)
            else
                call integrate_maxwell_sphere_curved_adjacent_rwg_pair_3d( &
                    vertices, triangles, edge_vertices, edge_triangles, basis, &
                    panel, basis, edge_triangles(second_slot, basis), 1.0_dp, &
                    0.7_dp, 4, 1.0e-5_dp, 2, panel_entry, status)
            end if
            direct_entry = direct_entry + panel_entry
        end do
    end do
    call record_condition(abs(matrix(basis, basis) - direct_entry) < &
        3.0e-13_dp*max(1.0_dp, abs(direct_entry)), &
        "global diagonal matches direct singular-panel composition")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call assemble_maxwell_sphere_curved_vector_potential_rwg_3d( &
        scaled_vertices, triangles, 2.0_dp, 0.35_dp, 4, 1.0e-5_dp, 2, &
        scaled_matrix, status)
    error = maxval(abs(scaled_matrix - 8.0_dp*matrix))/ &
        maxval(abs(scaled_matrix))
    call record_condition(status == 0 .and. error < 3.0e-11_dp, &
        "global curved vector potential obeys analytical R-cubed scaling")

    call check_summary("Global curved-sphere RWG vector potential")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_vector_potential
