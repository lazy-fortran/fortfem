program test_broken_skeleton_spaces
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        BROKEN_SPACE_H1, BROKEN_SPACE_HCURL, BROKEN_SPACE_HDIV, &
        BROKEN_SPACE_L2, SKELETON_SPACE_NORMAL, SKELETON_SPACE_SCALAR, &
        SKELETON_SPACE_TANGENTIAL, broken_space_layout_global_count, &
        broken_space_layout_maps, broken_space_layout_t, &
        initialize_broken_space_layout, initialize_skeleton_space_layout, &
        skeleton_space_layout_global_count, skeleton_space_layout_maps, &
        skeleton_space_layout_t
    implicit none

    integer, parameter :: cell_count = 3, local_dof_count = 2
    integer, parameter :: facet_count = 3
    integer, parameter :: cell_sign(2, 3) = reshape([ &
        1, -1, -1, 1, 1, 1], [2, 3])
    integer, parameter :: facet_cells(2, 3) = reshape([ &
        1, 2, 2, 0, 3, 3], [2, 3])
    integer, parameter :: facet_sign(2, 3) = reshape([ &
        1, -1, -1, 0, 1, 1], [2, 3])
    integer, parameter :: bad_facet_sign(2, 3) = reshape([ &
        1, 1, -1, 1, 1, 1], [2, 3])
    integer, parameter :: expected_broken_map(2, 3) = reshape([ &
        1, 2, 3, 4, 5, 6], [2, 3])
    integer, parameter :: expected_skeleton_map(2, 3) = reshape([ &
        1, 2, 3, 4, 5, 6], [2, 3])
    type(broken_space_layout_t) :: broken
    type(skeleton_space_layout_t) :: skeleton
    integer, allocatable :: broken_map(:, :), broken_sign(:, :)
    integer, allocatable :: skeleton_map(:, :), cells(:, :), signs(:, :)
    integer :: status
    logical :: all_passed

    all_passed = .true.

    call initialize_broken_space_layout( &
        broken, BROKEN_SPACE_HCURL, cell_count, local_dof_count, &
        local_sign=cell_sign, status=status)
    call record_condition(status == 0, &
        "broken H(curl) layout accepts signed local maps")
    call record_condition(broken_space_layout_global_count(broken) == 6, &
        "broken layout owns one copy of every local degree of freedom")
    call broken_space_layout_maps(broken, broken_map, broken_sign, status)
    call record_condition(status == 0, "broken layout exports its maps")
    call record_condition(all(broken_map == expected_broken_map), &
        "broken layout uses deterministic cell-major numbering")
    call record_condition(all(broken_sign == cell_sign), &
        "broken layout preserves H(curl) orientations")

    call initialize_broken_space_layout( &
        broken, BROKEN_SPACE_H1, 2, local_dof_count, status=status)
    call record_condition(status == 0, &
        "scalar broken layout supplies positive signs")
    call broken_space_layout_maps(broken, broken_map, broken_sign, status)
    call record_condition(status == 0 .and. all(broken_sign == 1), &
        "H1 broken layout has no orientation flips")
    call initialize_broken_space_layout( &
        broken, BROKEN_SPACE_HDIV, cell_count, local_dof_count, &
        local_sign=cell_sign, status=status)
    call record_condition(status == 0, "H(div) layout accepts signed face maps")
    call initialize_broken_space_layout( &
        broken, BROKEN_SPACE_L2, cell_count, local_dof_count, status=status)
    call record_condition(status == 0, "L2 layout accepts independent cell blocks")

    call initialize_skeleton_space_layout( &
        skeleton, SKELETON_SPACE_TANGENTIAL, cell_count, facet_cells, &
        facet_sign, local_dof_count, status)
    call record_condition(status == 0, &
        "tangential skeleton accepts interior and boundary facets")
    call record_condition(skeleton_space_layout_global_count(skeleton) == 6, &
        "skeleton layout owns one trace block per facet")
    call skeleton_space_layout_maps( &
        skeleton, skeleton_map, cells, signs, status)
    call record_condition(status == 0, &
        "skeleton layout exports maps and adjacency")
    call record_condition(all(skeleton_map == expected_skeleton_map), &
        "skeleton layout uses deterministic facet-major numbering")
    call record_condition(all(cells == facet_cells) .and. all(signs == facet_sign), &
        "skeleton layout preserves side incidence and signs")
    call record_condition(signs(2, 2) == 0 .and. cells(2, 2) == 0, &
        "skeleton layout marks a boundary side explicitly")

    call initialize_broken_space_layout( &
        broken, BROKEN_SPACE_L2, 2, 0, status=status)
    call record_condition(status /= 0, &
        "broken layout rejects zero local dimension")
    call initialize_skeleton_space_layout( &
        skeleton, SKELETON_SPACE_NORMAL, cell_count, facet_cells, &
        bad_facet_sign, local_dof_count, status)
    call record_condition(status /= 0, &
        "skeleton layout rejects a nonzero boundary sign")

    if (allocated(broken_map)) deallocate(broken_map)
    if (allocated(broken_sign)) deallocate(broken_sign)
    if (allocated(skeleton_map)) deallocate(skeleton_map)
    if (allocated(cells)) deallocate(cells)
    if (allocated(signs)) deallocate(signs)
    call check_summary("broken and skeleton space layouts")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_broken_skeleton_spaces
