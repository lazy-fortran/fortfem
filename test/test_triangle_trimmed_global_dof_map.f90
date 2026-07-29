program test_triangle_trimmed_global_dof_map
    use check, only: check_condition, check_summary
    use fortfem_api, only: build_triangle_trimmed_dof_map
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    implicit none

    type(mesh_2d_t) :: mesh
    integer, allocatable :: global_dofs(:, :), transforms(:, :)
    integer :: edge_dofs(3), edge_signs(3), global_dof_count
    integer :: local_dof, moment, order, status
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
    call mesh%build_edge_dof_numbering()

    call build_triangle_trimmed_dof_map( &
        mesh, 1, global_dofs, transforms, global_dof_count, status)
    call mesh%get_triangle_edge_dofs(1, edge_dofs, edge_signs)
    call record_condition(status == 0 .and. global_dof_count == mesh%n_edges, &
        "Order-one trimmed map preserves the existing global dimension")
    call record_condition(all(global_dofs(:, 1) == edge_dofs + 1) .and. &
        all(transforms(:, 1) == edge_signs), &
        "Order-one trimmed map preserves existing edge numbering and signs")
    deallocate(global_dofs, transforms)

    order = 4
    call build_triangle_trimmed_dof_map( &
        mesh, order, global_dofs, transforms, global_dof_count, status)
    call record_condition(status == 0 .and. &
        size(global_dofs, 1) == order * (order + 2), &
        "Higher-order trimmed map has the exact local dimension")
    call record_condition(global_dof_count == &
        mesh%n_edges * order + mesh%n_triangles * order * (order - 1), &
        "Higher-order trimmed map has exact edge and cell global dimensions")

    do moment = 1, order
        call record_condition(global_dofs(2 * order + moment, 1) == &
            global_dofs(moment, 2), &
            "Both cells map a shared-edge moment to one global degree of freedom")
        if (mod(moment, 2) == 1) then
            call record_condition(transforms(2 * order + moment, 1) == -1, &
                "Reversed shared edge applies odd-moment negative parity")
        else
            call record_condition(transforms(2 * order + moment, 1) == 1, &
                "Reversed shared edge applies even-moment positive parity")
        end if
        call record_condition(transforms(moment, 2) == 1, &
            "Aligned shared edge leaves each higher moment unchanged")
    end do

    do local_dof = 3 * order + 1, size(global_dofs, 1)
        call record_condition(global_dofs(local_dof, 1) > mesh%n_edges * order, &
            "Cell moments follow all shared edge blocks")
        call record_condition(global_dofs(local_dof, 1) /= &
            global_dofs(local_dof, 2), &
            "Cell moments remain local to their triangle")
        call record_condition(transforms(local_dof, 1) == 1 .and. &
            transforms(local_dof, 2) == 1, &
            "Cell moments need no edge orientation transform")
    end do

    call build_triangle_trimmed_dof_map( &
        mesh, 0, global_dofs, transforms, global_dof_count, status)
    call record_condition(status /= 0, &
        "Trimmed global map rejects order zero")

    call check_summary("Arbitrary-order triangle trimmed global DOF map")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_trimmed_global_dof_map
