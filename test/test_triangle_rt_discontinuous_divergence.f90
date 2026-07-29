program test_triangle_rt_discontinuous_divergence
    use check, only: check_condition, check_summary
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    use fortfem_api, only: &
        assemble_triangle_rt_divergence_csc, &
        build_triangle_discontinuous_dof_map, &
        build_triangle_trimmed_dof_map, initialize_triangle_lagrange_basis, &
        interpolate_triangle_rt, triangle_lagrange_basis_t, &
        triangle_lagrange_nodes
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    implicit none

    type(mesh_2d_t) :: mesh
    type(triangle_lagrange_basis_t) :: scalar_basis
    type(csc_t) :: divergence_matrix
    type(fortsparse_status_t) :: sparse_status
    real(dp), allocatable :: divergence_action(:), local_rt_dofs(:)
    real(dp), allocatable :: nodes(:, :), rt_dofs(:), scalar_dofs(:)
    integer, allocatable :: rt_global_dofs(:, :), rt_transforms(:, :)
    integer, allocatable :: scalar_global_dofs(:, :)
    logical, allocatable :: assigned(:)
    real(dp) :: exact_pairing, physical_point(2), vertices(2, 3)
    integer :: active_degree, degree, dof, local_status, node
    integer :: rt_dof_count, scalar_dof_count, triangle
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

    do degree = 0, 3
        active_degree = degree
        call build_triangle_trimmed_dof_map( &
            mesh, degree + 1, rt_global_dofs, rt_transforms, rt_dof_count, &
            local_status)
        call build_triangle_discontinuous_dof_map( &
            mesh, degree, scalar_global_dofs, scalar_dof_count, local_status)
        call initialize_triangle_lagrange_basis( &
            degree, scalar_basis, local_status)
        call triangle_lagrange_nodes(scalar_basis, nodes)
        allocate(rt_dofs(rt_dof_count), assigned(rt_dof_count))
        allocate(scalar_dofs(scalar_dof_count))
        allocate(local_rt_dofs(size(rt_global_dofs, 1)))
        rt_dofs = 0.0_dp
        scalar_dofs = 0.0_dp
        assigned = .false.

        do triangle = 1, mesh%n_triangles
            vertices = mesh%vertices(:, mesh%triangles(:, triangle))
            call interpolate_triangle_rt( &
                vertices, degree, 2 * degree + 4, manufactured_flux, &
                local_rt_dofs, local_status)
            do dof = 1, size(local_rt_dofs)
                if (.not. assigned(rt_global_dofs(dof, triangle))) then
                    rt_dofs(rt_global_dofs(dof, triangle)) = &
                        real(rt_transforms(dof, triangle), dp) * &
                        local_rt_dofs(dof)
                    assigned(rt_global_dofs(dof, triangle)) = .true.
                end if
            end do
            do node = 1, size(nodes, 2)
                physical_point = vertices(:, 1) + &
                    (vertices(:, 2) - vertices(:, 1)) * nodes(1, node) + &
                    (vertices(:, 3) - vertices(:, 1)) * nodes(2, node)
                scalar_dofs(scalar_global_dofs(node, triangle)) = &
                    physical_point(1)**degree
            end do
        end do

        call assemble_triangle_rt_divergence_csc( &
            mesh, degree, 2 * degree + 2, divergence_matrix, sparse_status)
        allocate(divergence_action(scalar_dof_count))
        divergence_action = csc_matvec(divergence_matrix, rt_dofs)
        exact_pairing = real(degree + 2, dp) / real(2 * degree + 1, dp)
        call record_condition(sparse_status%code == 0 .and. &
            divergence_matrix%nrow == scalar_dof_count .and. &
            divergence_matrix%ncol == rt_dof_count, &
            "Mixed RT-DG divergence form has exact rectangular dimensions")
        call record_condition(abs( &
            dot_product(scalar_dofs, divergence_action) - exact_pairing) < &
            2.0e-9_dp, &
            "Mixed RT-DG divergence form reproduces an exact polynomial pair")

        deallocate( &
            assigned, divergence_action, local_rt_dofs, nodes, rt_dofs, &
            rt_global_dofs, rt_transforms, scalar_dofs, scalar_global_dofs)
    end do

    call build_triangle_discontinuous_dof_map( &
        mesh, -1, scalar_global_dofs, scalar_dof_count, local_status)
    call record_condition(local_status /= 0, &
        "Discontinuous triangle map rejects negative degree")

    call check_summary("Arbitrary-order RT-DG divergence form")
    if (.not. all_passed) error stop 1

contains

    subroutine manufactured_flux(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value(1) = x**(active_degree + 1)
        value(2) = x**active_degree * y
    end subroutine manufactured_flux

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_rt_discontinuous_divergence
