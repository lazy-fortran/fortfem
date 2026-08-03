program test_tetra_edge_dof_map
    use check, only: check_condition, check_summary
    use fortfem_feec, only: build_tetra_edge_dof_map
    implicit none

    integer, allocatable :: edges(:, :), global_dofs(:, :), orientations(:, :)
    integer :: status
    integer :: tetrahedra(4, 2)
    logical :: all_passed

    all_passed = .true.
    tetrahedra(:, 1) = [1, 2, 3, 4]
    tetrahedra(:, 2) = [1, 3, 2, 5]
    call build_tetra_edge_dof_map( &
        tetrahedra, edges, global_dofs, orientations, status)

    call record_condition(status == 0, &
        "Two-tetrahedron edge topology builds")
    call record_condition(size(edges, 2) == 9, &
        "Shared tetrahedral face reuses its three global edges")
    call record_condition(global_dofs(1, 1) == global_dofs(2, 2), &
        "Shared edge 1-2 has one global degree of freedom")
    call record_condition(global_dofs(2, 1) == global_dofs(1, 2), &
        "Shared edge 1-3 has one global degree of freedom")
    call record_condition(global_dofs(4, 1) == global_dofs(4, 2), &
        "Shared edge 2-3 has one global degree of freedom")
    call record_condition( &
        orientations(4, 1) == -orientations(4, 2), &
        "Reversed shared edge exposes opposite local orientations")
    call record_condition(all(global_dofs >= 1) .and. &
        all(global_dofs <= size(edges, 2)), &
        "Every tetrahedral local edge maps to a valid global degree")

    call check_summary("Tetrahedral edge degree-of-freedom map")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_edge_dof_map
