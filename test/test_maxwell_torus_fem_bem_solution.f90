program test_maxwell_torus_fem_bem_solution
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_fem_bem_torus_curved_system_3d, &
        build_tetra_nedelec_dof_map, generate_solid_torus_tetra_mesh, &
        solve_maxwell_fem_bem_torus_curved_system_3d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp
    real(dp), parameter :: minor_radius = 0.65_dp
    real(dp), parameter :: electric_field(3) = [0.5_dp, -0.3_dp, 0.2_dp]
    real(dp), parameter :: electric_rotation(3) = [0.0_dp, 0.0_dp, 0.35_dp]
    integer, allocatable :: boundary_triangles(:, :), faces(:, :)
    integer, allocatable :: edge_orientations(:, :), edges(:, :)
    integer, allocatable :: face_permutations(:, :, :), global_dofs(:, :)
    integer, allocatable :: tetrahedra(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    complex(dp), allocatable :: matrix(:, :), right_hand_side(:)
    complex(dp), allocatable :: exact_state(:), field(:), current(:)
    real(dp) :: edge_delta(3), edge_midpoint(3)
    integer :: edge, field_count, status

    call generate_solid_torus_tetra_mesh( &
        major_radius, minor_radius, 1, 4, 4, vertices, tetrahedra, &
        boundary_triangles, parameters)
    call build_tetra_nedelec_dof_map( &
        1, tetrahedra, edges, faces, global_dofs, edge_orientations, &
        face_permutations, status)
    call check_condition(status == 0, &
        "torus FEM-BEM solution obtains a compatible Nedelec map")
    if (status /= 0) error stop 1
    field_count = maxval(global_dofs)
    call assemble_maxwell_fem_bem_torus_curved_system_3d( &
        vertices, tetrahedra, parameters, boundary_triangles, major_radius, &
        minor_radius, 1.1_dp, -0.7_dp, 0.8_dp, 1.4_dp, 3, 1.0e-5_dp, 1, &
        matrix, status)
    call check_condition(status == 0, &
        "torus FEM-BEM manufactured system assembles")
    if (status /= 0) error stop 1

    allocate(exact_state(size(matrix, 1)), right_hand_side(size(matrix, 1)))
    exact_state = cmplx(0.0_dp, 0.0_dp, dp)
    do edge = 1, size(edges, 2)
        edge_delta = vertices(:, edges(2, edge)) - vertices(:, edges(1, edge))
        edge_midpoint = 0.5_dp*(vertices(:, edges(2, edge)) + &
            vertices(:, edges(1, edge)))
        exact_state(edge) = cmplx(dot_product( &
            manufactured_vector(edge_midpoint), edge_delta), 0.0_dp, dp)
    end do
    right_hand_side = matmul(matrix, exact_state)
    call solve_maxwell_fem_bem_torus_curved_system_3d( &
        vertices, tetrahedra, parameters, boundary_triangles, major_radius, &
        minor_radius, 1.1_dp, -0.7_dp, 0.8_dp, 1.4_dp, 3, 1.0e-5_dp, 1, &
        right_hand_side, field, current, status)
    call check_condition(status == 0, &
        "torus FEM-BEM manufactured system solves")
    if (status /= 0) error stop 1
    call check_condition(size(field) == field_count .and. size(current) > 0, &
        "torus FEM-BEM solve returns volume and surface unknowns")
    call check_condition(maxval(abs(field - exact_state(:field_count))) < 3.0e-9_dp, &
        "torus FEM-BEM solve recovers an affine Nedelec field")
    call check_condition(maxval(abs(current)) < 3.0e-9_dp, &
        "torus FEM-BEM manufactured surface current is zero")
    call check_condition(maxval(abs(matrix - transpose(matrix))) < 3.0e-10_dp, &
        "torus FEM-BEM manufactured system is reciprocal")
    call check_summary("Torus Maxwell FEM-BEM manufactured solution")

contains

    pure function manufactured_vector(point) result(vector)
        real(dp), intent(in) :: point(3)
        real(dp) :: vector(3)

        vector = electric_field + [ &
            electric_rotation(2)*point(3) - electric_rotation(3)*point(2), &
            electric_rotation(3)*point(1) - electric_rotation(1)*point(3), &
            electric_rotation(1)*point(2) - electric_rotation(2)*point(1)]
    end function manufactured_vector
end program test_maxwell_torus_fem_bem_solution
