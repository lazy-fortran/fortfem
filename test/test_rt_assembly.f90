program test_rt_assembly
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_assembly_rt_2d, only: assemble_rt_div_mass_element, &
        assemble_rt_div_mass
    use check, only: check_condition, check_summary
    implicit none

    logical :: all_passed

    all_passed = .true.
    call test_reference_element_matrix()
    call test_global_constant_fields()
    call check_summary("Raviart-Thomas assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine test_reference_element_matrix()
        type(mesh_2d_t) :: mesh
        real(dp) :: element_matrix(3, 3)
        real(dp), parameter :: expected(3, 3) = reshape([ &
            7.0_dp / 3.0_dp, 2.0_dp, 11.0_dp / 6.0_dp, &
            2.0_dp, 13.0_dp / 6.0_dp, 2.0_dp, &
            11.0_dp / 6.0_dp, 2.0_dp, 7.0_dp / 3.0_dp], [3, 3])

        mesh%n_vertices = 3
        mesh%n_triangles = 1
        allocate(mesh%vertices(2, 3), mesh%triangles(3, 1))
        mesh%vertices(:, 1) = [0.0_dp, 0.0_dp]
        mesh%vertices(:, 2) = [1.0_dp, 0.0_dp]
        mesh%vertices(:, 3) = [0.0_dp, 1.0_dp]
        mesh%triangles(:, 1) = [1, 2, 3]

        call assemble_rt_div_mass_element(mesh, 1, element_matrix)

        call record_condition( &
            maxval(abs(element_matrix - expected)) < 1.0e-13_dp, &
            "Reference div-div plus mass matrix matches exact integration")
    end subroutine test_reference_element_matrix

    subroutine test_global_constant_fields()
        type(mesh_2d_t) :: mesh
        real(dp), allocatable :: matrix(:, :), x_dofs(:), y_dofs(:)
        real(dp) :: x_energy, y_energy, cross_energy
        real(dp) :: edge_vector(2)
        integer :: edge_id, dof_id, vertex_a, vertex_b

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
        call assemble_rt_div_mass(mesh, matrix)

        do edge_id = 1, mesh%n_edges
            vertex_a = mesh%edges(1, edge_id)
            vertex_b = mesh%edges(2, edge_id)
            edge_vector = mesh%vertices(:, vertex_b) - &
                mesh%vertices(:, vertex_a)
            dof_id = mesh%edge_to_dof(edge_id) + 1
            x_dofs(dof_id) = edge_vector(2)
            y_dofs(dof_id) = -edge_vector(1)
        end do

        x_energy = dot_product(x_dofs, matmul(matrix, x_dofs))
        y_energy = dot_product(y_dofs, matmul(matrix, y_dofs))
        cross_energy = dot_product(x_dofs, matmul(matrix, y_dofs))

        call record_condition(abs(x_energy - 1.0_dp) < 1.0e-13_dp, &
            "Constant x flux has exact unit div-mass energy")
        call record_condition(abs(y_energy - 1.0_dp) < 1.0e-13_dp, &
            "Constant y flux has exact unit div-mass energy")
        call record_condition(abs(cross_energy) < 1.0e-13_dp, &
            "Orthogonal constant fluxes have zero cross energy")
    end subroutine test_global_constant_fields

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_rt_assembly
