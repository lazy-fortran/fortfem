program test_maxwell_torus_fem_bem_nonzero_state
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_maxwell_fem_bem_torus_curved_system_3d, &
        solve_maxwell_fem_bem_torus_curved_system_3d
    use fortfem_core, only: generate_solid_torus_tetra_mesh
    use fortfem_feec, only: build_tetra_nedelec_dof_map
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp
    real(dp), parameter :: minor_radius = 0.65_dp
    real(dp), parameter :: electric_field(3) = [0.7_dp, -0.4_dp, 0.2_dp]
    integer, allocatable :: boundary_triangles(:, :), edges(:, :), faces(:, :)
    integer, allocatable :: edge_orientations(:, :), global_dofs(:, :)
    integer, allocatable :: face_permutations(:, :, :), tetrahedra(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    complex(dp), allocatable :: matrix(:, :), right_hand_side(:)
    complex(dp), allocatable :: exact_state(:), field(:), current(:)
    real(dp) :: field_error, current_error
    integer :: edge, index, field_count, status

    call generate_solid_torus_tetra_mesh( &
        major_radius, minor_radius, 1, 4, 4, vertices, tetrahedra, &
        boundary_triangles, parameters)
    call build_tetra_nedelec_dof_map( &
        1, tetrahedra, edges, faces, global_dofs, edge_orientations, &
        face_permutations, status)
    call check_condition(status == 0, "nonzero torus state obtains Nedelec map")
    if (status /= 0) error stop 1
    field_count = maxval(global_dofs)
    call assemble_maxwell_fem_bem_torus_curved_system_3d( &
        vertices, tetrahedra, parameters, boundary_triangles, major_radius, &
        minor_radius, 1.1_dp, -0.7_dp, 0.8_dp, 1.4_dp, 3, 1.0e-5_dp, 1, &
        matrix, status)
    call check_condition(status == 0, "nonzero torus FEM-BEM system assembles")
    if (status /= 0) error stop 1

    allocate(exact_state(size(matrix, 1)), right_hand_side(size(matrix, 1)))
    exact_state = cmplx(0.0_dp, 0.0_dp, dp)
    do edge = 1, size(edges, 2)
        exact_state(edge) = cmplx(dot_product( &
            electric_field, vertices(:, edges(2, edge)) - &
            vertices(:, edges(1, edge))), 0.0_dp, dp)
    end do
    do index = field_count + 1, size(exact_state)
        exact_state(index) = cmplx( &
            0.04_dp*cos(0.17_dp*real(index, dp)), &
            -0.03_dp*sin(0.23_dp*real(index, dp)), dp)
    end do
    right_hand_side = matmul(matrix, exact_state)
    call solve_maxwell_fem_bem_torus_curved_system_3d( &
        vertices, tetrahedra, parameters, boundary_triangles, major_radius, &
        minor_radius, 1.1_dp, -0.7_dp, 0.8_dp, 1.4_dp, 3, 1.0e-5_dp, 1, &
        right_hand_side, field, current, status)
    call check_condition(status == 0, "nonzero torus FEM-BEM state solves")
    if (status /= 0) error stop 1
    field_error = maxval(abs(field - exact_state(:field_count)))
    current_error = maxval(abs(current - exact_state(field_count + 1:)))
    call check_condition(field_error < 5.0e-9_dp, &
        "nonzero torus state recovers the volume field")
    call check_condition(current_error < 5.0e-9_dp .and. maxval(abs(current)) > 1.0e-3_dp, &
        "nonzero torus state recovers a nonzero surface current")
    call check_summary("Nonzero torus Maxwell FEM-BEM state")
end program test_maxwell_torus_fem_bem_nonzero_state
