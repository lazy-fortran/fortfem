program test_triangle_nedelec_discontinuous_curl
    use check, only: check_condition, check_summary
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    use fortfem_api, only: &
        assemble_triangle_nedelec_curl_csc, &
        build_triangle_discontinuous_dof_map, &
        build_triangle_trimmed_dof_map, initialize_triangle_lagrange_basis, &
        interpolate_triangle_nedelec, triangle_lagrange_basis_t, &
        triangle_lagrange_nodes
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    implicit none

    type(mesh_2d_t) :: mesh
    type(triangle_lagrange_basis_t) :: scalar_basis
    type(csc_t) :: curl_matrix
    type(fortsparse_status_t) :: sparse_status
    real(dp), allocatable :: curl_action(:), local_nedelec_dofs(:)
    real(dp), allocatable :: nedelec_dofs(:), nodes(:, :), scalar_dofs(:)
    integer, allocatable :: nedelec_global_dofs(:, :)
    integer, allocatable :: nedelec_transforms(:, :)
    integer, allocatable :: scalar_global_dofs(:, :)
    logical, allocatable :: assigned(:)
    real(dp) :: exact_pairing, physical_point(2), vertices(2, 3)
    integer :: active_order, dof, local_status, node, order
    integer :: nedelec_dof_count, scalar_dof_count, triangle
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

    do order = 1, 4
        active_order = order
        call build_triangle_trimmed_dof_map( &
            mesh, order, nedelec_global_dofs, nedelec_transforms, &
            nedelec_dof_count, local_status)
        call build_triangle_discontinuous_dof_map( &
            mesh, order - 1, scalar_global_dofs, scalar_dof_count, &
            local_status)
        call initialize_triangle_lagrange_basis( &
            order - 1, scalar_basis, local_status)
        call triangle_lagrange_nodes(scalar_basis, nodes)
        allocate(nedelec_dofs(nedelec_dof_count))
        allocate(assigned(nedelec_dof_count))
        allocate(scalar_dofs(scalar_dof_count))
        allocate(local_nedelec_dofs(size(nedelec_global_dofs, 1)))
        nedelec_dofs = 0.0_dp
        scalar_dofs = 0.0_dp
        assigned = .false.

        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            call interpolate_triangle_nedelec( &
                vertices, order, 2 * order + 2, manufactured_field, &
                local_nedelec_dofs, local_status)
            do dof = 1, size(local_nedelec_dofs)
                if (.not. assigned(nedelec_global_dofs(dof, triangle))) then
                    nedelec_dofs(nedelec_global_dofs(dof, triangle)) = &
                        real(nedelec_transforms(dof, triangle), dp) * &
                        local_nedelec_dofs(dof)
                    assigned(nedelec_global_dofs(dof, triangle)) = .true.
                end if
            end do
            do node = 1, size(nodes, 2)
                physical_point = vertices(:, 1) + &
                    (vertices(:, 2) - vertices(:, 1)) * nodes(1, node) + &
                    (vertices(:, 3) - vertices(:, 1)) * nodes(2, node)
                scalar_dofs(scalar_global_dofs(node, triangle)) = &
                    physical_point(1)**(order - 1)
            end do
        end do

        call assemble_triangle_nedelec_curl_csc( &
            mesh, order, 2 * order, curl_matrix, sparse_status)
        allocate(curl_action(scalar_dof_count))
        curl_action = csc_matvec(curl_matrix, nedelec_dofs)
        exact_pairing = real(order + 1, dp) / real(2 * order - 1, dp)
        call record_condition(sparse_status%code == 0 .and. &
            curl_matrix%nrow == scalar_dof_count .and. &
            curl_matrix%ncol == nedelec_dof_count, &
            "Mixed Nedelec-DG curl form has exact rectangular dimensions")
        call record_condition(abs( &
            dot_product(scalar_dofs, curl_action) - exact_pairing) < &
            2.0e-9_dp, &
            "Mixed Nedelec-DG curl form reproduces an exact polynomial pair")

        deallocate( &
            assigned, curl_action, local_nedelec_dofs, nedelec_dofs, &
            nedelec_global_dofs, nedelec_transforms, nodes, scalar_dofs, &
            scalar_global_dofs)
    end do

    call assemble_triangle_nedelec_curl_csc( &
        mesh, 0, 2, curl_matrix, sparse_status)
    call record_condition(sparse_status%code /= 0, &
        "Mixed Nedelec-DG curl form rejects order zero")

    call check_summary("Arbitrary-order Nedelec-DG curl form")
    if (.not. all_passed) error stop 1

contains

    subroutine manufactured_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value(1) = -y * x**(active_order - 1)
        value(2) = x**active_order
    end subroutine manufactured_field

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_nedelec_discontinuous_curl
