program test_maxwell_fem_bem_coupling_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_efie_rwg_3d, &
        assemble_maxwell_fem_bem_boundary_matrix_3d, &
        assemble_maxwell_fem_bem_system_3d, &
        assemble_maxwell_rwg_mass_matrix, &
        assemble_tetra_nedelec_curl_mass_element, &
        build_maxwell_rwg_surface_space, build_tetra_edge_dof_map, &
        map_maxwell_rwg_to_tetra_nedelec_edges
    use fortfem_api, only: solve_maxwell_fem_bem_system_3d
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: boundary_matrix(:, :), efie(:, :)
    complex(dp), allocatable :: scaled_matrix(:, :)
    complex(dp), allocatable :: system_matrix(:, :)
    complex(dp), allocatable :: right_hand_side(:), solved_current(:)
    complex(dp), allocatable :: solved_field(:)
    integer, allocatable :: edge_dofs(:, :), edge_orientations(:, :)
    integer, allocatable :: edges(:, :), rwg_dofs(:), rwg_triangles(:, :)
    integer, allocatable :: rwg_vertices(:, :)
    real(dp), allocatable :: rwg_scales(:), trace_coefficients(:)
    real(dp), allocatable :: element_matrix(:, :), rwg_mass(:, :)
    real(dp) :: electric_field(3), vertices(3, 4)
    real(dp), allocatable :: volume_coefficients(:)
    integer :: boundary_triangles(3, 4), edge, status, tetrahedra(4, 1)
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    tetrahedra(:, 1) = [1, 2, 3, 4]
    boundary_triangles(:, 1) = [1, 3, 2]
    boundary_triangles(:, 2) = [1, 2, 4]
    boundary_triangles(:, 3) = [1, 4, 3]
    boundary_triangles(:, 4) = [2, 3, 4]

    call assemble_maxwell_fem_bem_boundary_matrix_3d( &
        vertices, tetrahedra, boundary_triangles, 0.8_dp, 1.7_dp, 8, &
        1.0e-3_dp, 1, boundary_matrix, status)
    call record_condition(status == 0 .and. &
        maxval(abs(boundary_matrix - transpose(boundary_matrix))) < 2.0e-13_dp, &
        "Pulled-back Maxwell FEM-BEM boundary matrix is complex symmetric")

    call build_tetra_edge_dof_map( &
        tetrahedra, edges, edge_dofs, edge_orientations, status)
    call build_maxwell_rwg_surface_space( &
        vertices, boundary_triangles, rwg_vertices, rwg_triangles, status)
    call map_maxwell_rwg_to_tetra_nedelec_edges( &
        vertices, tetrahedra, rwg_vertices, rwg_dofs, rwg_scales, status)
    allocate(volume_coefficients(size(edges, 2)))
    electric_field = [0.7_dp, -0.4_dp, 0.2_dp]
    do edge = 1, size(edges, 2)
        volume_coefficients(edge) = dot_product( &
            electric_field, vertices(:, edges(2, edge)) - &
            vertices(:, edges(1, edge)))
    end do
    allocate(trace_coefficients(size(rwg_dofs)))
    trace_coefficients = rwg_scales*volume_coefficients(rwg_dofs)
    call assemble_maxwell_efie_rwg_3d( &
        vertices, boundary_triangles, 0.8_dp, 1.7_dp, 8, 1.0e-3_dp, 1, &
        efie, status)
    call record_condition(abs(dot_product( &
        cmplx(volume_coefficients, 0.0_dp, dp), &
        matmul(boundary_matrix, volume_coefficients)) - dot_product( &
        cmplx(trace_coefficients, 0.0_dp, dp), &
        matmul(efie, trace_coefficients))) < 2.0e-12_dp, &
        "FEM-BEM pullback preserves the exact RWG trace energy")

    call assemble_maxwell_fem_bem_boundary_matrix_3d( &
        2.0_dp*vertices, tetrahedra, boundary_triangles, 0.4_dp, 1.7_dp, 8, &
        1.0e-3_dp, 1, scaled_matrix, status)
    call record_condition(maxval(abs(scaled_matrix - boundary_matrix)) < &
        2.0e-12_dp, &
        "Nedelec Maxwell boundary pullback is wave-scaled geometry invariant")

    call assemble_maxwell_fem_bem_system_3d( &
        vertices, tetrahedra, boundary_triangles, 1.2_dp, -0.64_dp, &
        0.8_dp, 1.7_dp, 8, 1.0e-3_dp, 1, system_matrix, status)
    call assemble_tetra_nedelec_curl_mass_element( &
        vertices, 1, 8, element_matrix, status, &
        curl_coefficient=1.2_dp, mass_coefficient=-0.64_dp)
    call assemble_maxwell_rwg_mass_matrix( &
        vertices, boundary_triangles, 4, rwg_mass, status)
    call record_condition(status == 0 .and. &
        size(system_matrix, 1) == 12 .and. &
        maxval(abs(system_matrix - transpose(system_matrix))) < 2.0e-12_dp .and. &
        maxval(abs(system_matrix(:6, :6) - element_matrix)) < 2.0e-12_dp .and. &
        maxval(abs(system_matrix(7:, 7:) + efie)) < 2.0e-12_dp, &
        "Mixed Maxwell FEM-BEM system contains the verified volume and EFIE blocks")
    call record_condition(sqrt(sum(abs( &
        matmul(system_matrix(7:, :6), volume_coefficients) - &
        matmul(rwg_mass, trace_coefficients))**2)) < 2.0e-12_dp, &
        "Mixed Maxwell system uses the exact RWG trace duality block")
    allocate(right_hand_side(12))
    right_hand_side(:6) = matmul( &
        cmplx(element_matrix, 0.0_dp, dp), volume_coefficients)
    right_hand_side(7:) = matmul(rwg_mass, trace_coefficients)
    call solve_maxwell_fem_bem_system_3d( &
        vertices, tetrahedra, boundary_triangles, 1.2_dp, -0.64_dp, &
        0.8_dp, 1.7_dp, 8, 1.0e-3_dp, 1, right_hand_side, solved_field, &
        solved_current, status)
    call record_condition(status == 0 .and. &
        sqrt(sum(abs(solved_field - volume_coefficients)**2)) < 2.0e-11_dp .and. &
        sqrt(sum(abs(solved_current)**2)) < 2.0e-11_dp, &
        "Mixed FEM-BEM solve recovers an analytical constant Nedelec field")
    call check_summary("Three-dimensional Maxwell FEM-BEM trace coupling")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_fem_bem_coupling_3d
