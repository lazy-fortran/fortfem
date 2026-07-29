program test_nedelec_global_assembly
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_assembly_nedelec_2d, only: assemble_nedelec_curl_mass
    use check, only: check_condition, check_summary
    implicit none

    type(mesh_2d_t) :: mesh
    real(dp), allocatable :: matrix(:, :), x_dofs(:), y_dofs(:)
    real(dp) :: x_energy, y_energy, cross_energy
    integer :: edge_id, dof_id, vertex_a, vertex_b
    logical :: all_passed

    all_passed = .true.
    mesh%n_vertices = 4
    mesh%n_triangles = 2
    mesh%has_triangles = .true.
    allocate(mesh%vertices(2, 4), mesh%triangles(3, 2))
    mesh%vertices(:, 1) = [0.0_dp, 0.0_dp]
    mesh%vertices(:, 2) = [1.0_dp, 0.0_dp]
    mesh%vertices(:, 3) = [1.0_dp, 1.0_dp]
    mesh%vertices(:, 4) = [0.0_dp, 1.0_dp]
    mesh%triangles(:, 1) = [1, 2, 3]
    mesh%triangles(:, 2) = [1, 3, 4]
    call mesh%build_edge_connectivity()

    allocate(matrix(mesh%n_edges, mesh%n_edges))
    allocate(x_dofs(mesh%n_edges), y_dofs(mesh%n_edges))
    call assemble_nedelec_curl_mass(mesh, matrix)

    do edge_id = 1, mesh%n_edges
        vertex_a = mesh%edges(1, edge_id)
        vertex_b = mesh%edges(2, edge_id)
        dof_id = mesh%edge_to_dof(edge_id) + 1
        x_dofs(dof_id) = mesh%vertices(1, vertex_b) - &
            mesh%vertices(1, vertex_a)
        y_dofs(dof_id) = mesh%vertices(2, vertex_b) - &
            mesh%vertices(2, vertex_a)
    end do

    x_energy = dot_product(x_dofs, matmul(matrix, x_dofs))
    y_energy = dot_product(y_dofs, matmul(matrix, y_dofs))
    cross_energy = dot_product(x_dofs, matmul(matrix, y_dofs))

    call record_condition(abs(x_energy - 1.0_dp) < 1.0e-13_dp, &
        "Constant x field has exact unit curl-mass energy")
    call record_condition(abs(y_energy - 1.0_dp) < 1.0e-13_dp, &
        "Constant y field has exact unit curl-mass energy")
    call record_condition(abs(cross_energy) < 1.0e-13_dp, &
        "Orthogonal constant fields have zero cross energy")

    call check_summary("Nedelec global assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_nedelec_global_assembly
