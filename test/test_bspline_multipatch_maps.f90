program test_bspline_multipatch_maps
    use check, only: check_condition, check_summary
    use fortfem_api, only: BSPLINE_FACE_X_MAX, BSPLINE_FACE_X_MIN, &
        build_bspline_feec_2d_multipatch_maps
    implicit none

    integer, parameter :: patch_count = 3, nx = 3, ny = 3
    integer, parameter :: h1_local = nx*ny
    integer, parameter :: hcurl_local = (nx - 1)*ny + nx*(ny - 1)
    integer, parameter :: l2_local = (nx - 1)*(ny - 1)
    integer :: nx_all(patch_count), ny_all(patch_count)
    integer :: left_patch(3), right_patch(3), face_left(3), face_right(3)
    logical :: reversed(3)
    integer, allocatable :: h1_offsets(:), h1_map(:)
    integer, allocatable :: hcurl_offsets(:), hcurl_map(:)
    integer, allocatable :: l2_offsets(:), l2_map(:)
    integer :: h1_count, hcurl_count, l2_count, status
    integer :: interface, trace, left, right, right_trace
    integer :: h1_distinct, hcurl_distinct, l2_distinct
    logical :: all_passed

    all_passed = .true.
    nx_all = nx
    ny_all = ny
    left_patch = [1, 2, 3]
    right_patch = [2, 3, 1]
    face_left = [BSPLINE_FACE_X_MAX, BSPLINE_FACE_X_MAX, BSPLINE_FACE_X_MAX]
    face_right = [BSPLINE_FACE_X_MIN, BSPLINE_FACE_X_MIN, BSPLINE_FACE_X_MIN]
    reversed = [.false., .true., .false.]

    call build_bspline_feec_2d_multipatch_maps( &
        nx_all, ny_all, left_patch, right_patch, face_left, face_right, &
        reversed, h1_offsets, h1_map, h1_count, hcurl_offsets, hcurl_map, &
        hcurl_count, l2_offsets, l2_map, l2_count, status)

    call record_condition(status == 0, "arbitrary 2D B-spline patch graph builds")
    call record_condition(all(h1_offsets == [1, 10, 19, 28]) .and. &
        all(hcurl_offsets == [1, 13, 25, 37]) .and. &
        all(l2_offsets == [1, 5, 9, 13]), &
        "packed multipatch offsets use the documented local ordering")
    call record_condition(h1_count == 18 .and. hcurl_count == 30 .and. &
        l2_count == 12 .and. all(h1_map /= 0) .and. all(hcurl_map /= 0) .and. &
        all(l2_map /= 0), "interface identifications reduce independent DOFs")

    h1_distinct = count_unique_abs(h1_map)
    hcurl_distinct = count_unique_abs(hcurl_map)
    l2_distinct = count_unique_abs(l2_map)
    call record_condition(h1_distinct == h1_count .and. &
        hcurl_distinct == hcurl_count .and. l2_distinct == l2_count, &
        "all packed maps have one global ID per signed equivalence class")

    do interface = 1, size(left_patch)
        do trace = 1, ny
            left = h1_offsets(left_patch(interface)) + nx + &
                (trace - 1)*nx - 1
            right_trace = ny + 1 - trace
            if (.not. reversed(interface)) right_trace = trace
            right = h1_offsets(right_patch(interface)) + 1 + &
                (right_trace - 1)*nx - 1
            call record_condition(h1_map(left) == h1_map(right), &
                "H1 traces are identified across every patch edge")
        end do

        do trace = 1, ny - 1
            right_trace = ny - trace
            if (.not. reversed(interface)) right_trace = trace
            left = hcurl_offsets(left_patch(interface)) + (nx - 1)*ny + &
                nx + (trace - 1)*nx - 1
            right = hcurl_offsets(right_patch(interface)) + (nx - 1)*ny + &
                1 + (right_trace - 1)*nx - 1
            call record_condition(hcurl_map(left) == merge(1, -1, &
                .not. reversed(interface))*hcurl_map(right), &
                "H(curl) traces preserve reversal orientation signs")
        end do
    end do

    call check_summary("arbitrary 2D B-spline multipatch maps")
    if (.not. all_passed) error stop 1

contains

    integer function count_unique_abs(values) result(count)
        integer, intent(in) :: values(:)
        integer, allocatable :: representatives(:)
        integer :: value, index
        logical :: already_seen

        allocate(representatives(size(values)))
        representatives = 0
        count = 0
        do value = 1, size(values)
            already_seen = .false.
            do index = 1, count
                if (representatives(index) == abs(values(value))) then
                    already_seen = .true.
                    exit
                end if
            end do
            if (.not. already_seen) then
                count = count + 1
                representatives(count) = abs(values(value))
            end if
        end do
    end function count_unique_abs

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_bspline_multipatch_maps
