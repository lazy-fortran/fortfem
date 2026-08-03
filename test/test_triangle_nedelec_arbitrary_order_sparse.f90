program test_triangle_nedelec_arbitrary_order_sparse
    use check, only: check_condition, check_summary
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    use fortfem_feec, only: &
        assemble_triangle_nedelec_curl_mass_csc, &
        build_triangle_discrete_gradient, build_triangle_trimmed_dof_map, &
        initialize_triangle_lagrange_basis, triangle_lagrange_basis_t, &
        triangle_lagrange_nodes
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    implicit none

    type(mesh_2d_t) :: mesh
    type(triangle_lagrange_basis_t) :: lagrange_basis
    type(csc_t) :: matrix
    type(fortsparse_status_t) :: sparse_status
    real(dp), allocatable :: dofs(:), gradient_matrix(:, :)
    real(dp), allocatable :: local_dofs(:), matrix_times_dofs(:)
    real(dp), allocatable :: nodal_values(:), nodes(:, :)
    integer, allocatable :: global_dofs(:, :), transforms(:, :)
    logical, allocatable :: assigned(:)
    real(dp) :: exact_energy, physical_x, physical_y, tensor(2, 2)
    real(dp) :: vertices(2, 3)
    integer :: dof, global_dof_count, local_status, node, order, triangle
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
    tensor = reshape([2.0_dp, 0.5_dp, 0.5_dp, 3.0_dp], [2, 2])

    do order = 1, 4
        call build_triangle_trimmed_dof_map( &
            mesh, order, global_dofs, transforms, global_dof_count, &
            local_status)
        call initialize_triangle_lagrange_basis( &
            order, lagrange_basis, local_status)
        call triangle_lagrange_nodes(lagrange_basis, nodes)
        call build_triangle_discrete_gradient( &
            order, gradient_matrix, local_status)
        allocate(dofs(global_dof_count), assigned(global_dof_count))
        allocate(local_dofs(size(gradient_matrix, 1)))
        allocate(nodal_values(size(nodes, 2)))
        dofs = 0.0_dp
        assigned = .false.

        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            do node = 1, size(nodes, 2)
                physical_x = vertices(1, 1) + &
                    (vertices(1, 2) - vertices(1, 1)) * nodes(1, node) + &
                    (vertices(1, 3) - vertices(1, 1)) * nodes(2, node)
                physical_y = vertices(2, 1) + &
                    (vertices(2, 2) - vertices(2, 1)) * nodes(1, node) + &
                    (vertices(2, 3) - vertices(2, 1)) * nodes(2, node)
                nodal_values(node) = physical_x**order + &
                    2.0_dp * physical_y**order
            end do
            local_dofs = matmul(gradient_matrix, nodal_values)
            do dof = 1, size(local_dofs)
                if (assigned(global_dofs(dof, triangle))) then
                    call record_condition(abs( &
                        dofs(global_dofs(dof, triangle)) - &
                        real(transforms(dof, triangle), dp) * &
                        local_dofs(dof)) < 3.0e-10_dp, &
                        "Neighbouring elements agree on each shared edge moment")
                else
                    dofs(global_dofs(dof, triangle)) = &
                        real(transforms(dof, triangle), dp) * local_dofs(dof)
                    assigned(global_dofs(dof, triangle)) = .true.
                end if
            end do
        end do

        call assemble_triangle_nedelec_curl_mass_csc( &
            mesh, order, 2 * order, matrix, sparse_status, &
            curl_coefficient=0.0_dp, mass_tensor=tensor)
        allocate(matrix_times_dofs(global_dof_count))
        matrix_times_dofs = csc_matvec(matrix, dofs)
        exact_energy = 14.0_dp * real(order * order, dp) / &
            real(2 * order - 1, dp) + 2.0_dp
        call record_condition(sparse_status%code == 0 .and. &
            matrix%nrow == global_dof_count .and. &
            matrix%ncol == global_dof_count, &
            "Sparse arbitrary-order operator has the exact global dimension")
        call record_condition(abs( &
            dot_product(dofs, matrix_times_dofs) - exact_energy) < &
            8.0e-8_dp, &
            "Sparse operator reproduces exact anisotropic gradient energy")

        deallocate( &
            assigned, dofs, global_dofs, gradient_matrix, local_dofs, &
            matrix_times_dofs, nodal_values, nodes, transforms)
    end do

    call assemble_triangle_nedelec_curl_mass_csc( &
        mesh, 0, 2, matrix, sparse_status)
    call record_condition(sparse_status%code /= 0, &
        "Sparse arbitrary-order assembly rejects order zero")

    call check_summary("Arbitrary-order triangle Nedelec sparse assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_nedelec_arbitrary_order_sparse
