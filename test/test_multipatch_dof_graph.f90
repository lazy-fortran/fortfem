program test_multipatch_dof_graph
    use check, only: check_condition, check_summary
    use fortfem_feec, only: build_multipatch_signed_dof_map
    implicit none

    integer :: local_to_global(9), status, global_count
    integer :: patch_offsets(4), left_patch(4), left_local(4)
    integer :: right_patch(4), right_local(4), interface_sign(4)
    integer :: relation, left, right, distinct_count
    integer :: seen(5)
    integer :: bad_offsets(3), bad_left_patch(2), bad_left_local(2)
    integer :: bad_right_patch(2), bad_right_local(2), bad_sign(2)
    logical :: all_passed

    all_passed = .true.
    patch_offsets = [1, 4, 7, 10]
    left_patch = [1, 2, 1, 2]
    left_local = [3, 3, 2, 2]
    right_patch = [2, 3, 2, 3]
    right_local = [1, 1, 2, 2]
    interface_sign = [1, -1, 1, 1]
    call build_multipatch_signed_dof_map( &
        patch_offsets, left_patch, left_local, right_patch, right_local, &
        interface_sign, local_to_global, global_count, status)

    distinct_count = 0
    seen = 0
    do relation = 1, size(local_to_global)
        if (all(seen(1:distinct_count) /= abs(local_to_global(relation)))) then
            distinct_count = distinct_count + 1
            seen(distinct_count) = abs(local_to_global(relation))
        end if
    end do
    call record_condition(status == 0 .and. global_count == 5 .and. &
        all(local_to_global /= 0) .and. distinct_count == global_count, &
        "arbitrary multipatch map creates one signed ID per equivalence class")

    do relation = 1, size(left_patch)
        left = patch_offsets(left_patch(relation)) + left_local(relation) - 1
        right = patch_offsets(right_patch(relation)) + right_local(relation) - 1
        call record_condition(local_to_global(left) == &
            interface_sign(relation)*local_to_global(right), &
            "multipatch map preserves every interface orientation relation")
    end do

    bad_offsets = [1, 3, 5]
    bad_left_patch = [1, 1]
    bad_left_local = [1, 1]
    bad_right_patch = [2, 2]
    bad_right_local = [1, 1]
    bad_sign = [1, -1]
    call build_multipatch_signed_dof_map( &
        bad_offsets, bad_left_patch, bad_left_local, bad_right_patch, &
        bad_right_local, bad_sign, local_to_global(1:4), global_count, status)
    call record_condition(status == 2, &
        "multipatch map rejects an inconsistent signed cycle")

    bad_sign = [1, 0]
    call build_multipatch_signed_dof_map( &
        bad_offsets, bad_left_patch, bad_left_local, bad_right_patch, &
        bad_right_local, bad_sign, local_to_global(1:4), global_count, status)
    call record_condition(status == 1, "multipatch map rejects invalid orientation signs")
    call check_summary("multipatch signed DOF graph")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_multipatch_dof_graph
