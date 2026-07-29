program test_edge_orientations
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use check, only: check_condition, check_summary
    implicit none

    type(mesh_2d_t) :: mesh
    integer :: first_dofs(3), second_dofs(3)
    integer :: first_signs(3), second_signs(3)
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
    call mesh%get_triangle_edge_dofs(1, first_dofs, first_signs)
    call mesh%get_triangle_edge_dofs(2, second_dofs, second_signs)

    call record_condition(first_dofs(3) == second_dofs(1), &
        "Adjacent cells share one global edge degree of freedom")
    call record_condition(first_signs(3) == -1, &
        "First cell reverses the shared global edge orientation")
    call record_condition(second_signs(1) == 1, &
        "Second cell follows the shared global edge orientation")
    call record_condition( &
        first_signs(3) == -second_signs(1), &
        "Interior edge traces receive opposite local orientation signs")

    call check_summary("Triangle edge orientations")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_edge_orientations
