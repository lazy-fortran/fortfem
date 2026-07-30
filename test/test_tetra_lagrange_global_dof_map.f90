program test_tetra_lagrange_global_dof_map
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_tetra_lagrange_dof_map, initialize_tetra_lagrange, &
        tetra_lagrange_barycentric_indices, tetra_lagrange_t
    implicit none

    type(tetra_lagrange_t) :: basis
    integer, allocatable :: barycentric_indices(:, :), global_dofs(:, :)
    integer :: degree, global_count, local_dof, shared_count, status
    integer :: tetrahedra(4, 2)
    logical :: all_passed

    all_passed = .true.
    tetrahedra(:, 1) = [1, 2, 3, 4]
    tetrahedra(:, 2) = [3, 2, 1, 5]
    do degree = 1, 4
        call initialize_tetra_lagrange(degree, basis, status)
        call tetra_lagrange_barycentric_indices( &
            basis, barycentric_indices)
        call build_tetra_lagrange_dof_map( &
            degree, tetrahedra, global_dofs, global_count, status)
        shared_count = 0
        do local_dof = 1, size(global_dofs, 1)
            if (barycentric_indices(4, local_dof) /= 0) cycle
            shared_count = shared_count + 1
            call record_condition(any( &
                global_dofs(:, 2) == global_dofs(local_dof, 1)), &
                "Tetrahedral H1 face nodes are shared globally")
        end do
        call record_condition(shared_count == (degree + 1)*(degree + 2)/2, &
            "Tetrahedral H1 shared face has the analytic node count")
        call record_condition(global_count == &
            2*(degree + 1)*(degree + 2)*(degree + 3)/6 - shared_count, &
            "Tetrahedral H1 global dimension deduplicates the shared face")
        deallocate(barycentric_indices, global_dofs)
    end do

    call build_tetra_lagrange_dof_map( &
        0, tetrahedra, global_dofs, global_count, status)
    call record_condition(status /= 0, &
        "Conforming tetrahedral H1 topology rejects degree zero")
    tetrahedra(:, 2) = [1, 1, 2, 5]
    call build_tetra_lagrange_dof_map( &
        2, tetrahedra, global_dofs, global_count, status)
    call record_condition(status /= 0, &
        "Tetrahedral H1 topology rejects repeated vertices")

    call check_summary("Arbitrary-order tetrahedral H1 topology")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_lagrange_global_dof_map
