program test_triangle_full_vector_global_dof_map
    use check, only: check_condition, check_summary
    use fortfem_api, only: build_triangle_full_vector_dof_map
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    implicit none

    type(mesh_2d_t) :: mesh
    integer, allocatable :: global_dofs(:, :), transforms(:, :)
    integer :: degree, global_dof_count, local_dof, moment, status
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

    call build_triangle_full_vector_dof_map( &
        mesh, 1, global_dofs, transforms, global_dof_count, status)
    call record_condition(status == 0 .and. size(global_dofs, 1) == 6, &
        "Degree-one full vector map consists of six edge moments")
    call record_condition(global_dof_count == 2 * mesh%n_edges, &
        "Degree-one full vector map has no cell moments")
    deallocate(global_dofs, transforms)

    degree = 4
    call build_triangle_full_vector_dof_map( &
        mesh, degree, global_dofs, transforms, global_dof_count, status)
    call record_condition(size(global_dofs, 1) == &
        (degree + 1) * (degree + 2), &
        "Full vector map has exact local [P_k]^2 dimension")
    call record_condition(global_dof_count == &
        mesh%n_edges * (degree + 1) + &
        mesh%n_triangles * (degree * degree - 1), &
        "Full vector map has exact global edge and cell dimensions")

    do moment = 1, degree + 1
        call record_condition(global_dofs(2 * (degree + 1) + moment, 1) == &
            global_dofs(moment, 2), &
            "Both cells share each full-family diagonal edge moment")
        if (mod(moment, 2) == 1) then
            call record_condition( &
                transforms(2 * (degree + 1) + moment, 1) == -1, &
                "Reversed full-family edge applies negative odd parity")
        else
            call record_condition( &
                transforms(2 * (degree + 1) + moment, 1) == 1, &
                "Reversed full-family edge applies positive even parity")
        end if
    end do
    do local_dof = 3 * (degree + 1) + 1, size(global_dofs, 1)
        call record_condition(global_dofs(local_dof, 1) /= &
            global_dofs(local_dof, 2), &
            "Full-family cell moments remain local to each triangle")
        call record_condition(transforms(local_dof, 1) == 1 .and. &
            transforms(local_dof, 2) == 1, &
            "Full-family cell moments need no orientation transform")
    end do

    call build_triangle_full_vector_dof_map( &
        mesh, 0, global_dofs, transforms, global_dof_count, status)
    call record_condition(status /= 0, &
        "Full vector global map rejects degree zero")

    call check_summary("Full polynomial vector global DOF map")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_full_vector_global_dof_map
