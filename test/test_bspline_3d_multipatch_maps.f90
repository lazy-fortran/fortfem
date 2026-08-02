program test_bspline_3d_multipatch_maps
    use check, only: check_condition, check_summary
    use fortfem_api, only: BSPLINE_FACE_Z_MAX, BSPLINE_FACE_Z_MIN, &
        build_bspline_feec_3d_interface_dofs, &
        build_bspline_feec_3d_multipatch_maps
    implicit none

    integer, parameter :: patch_count = 3, nx = 3, ny = 3, nz = 3
    integer, parameter :: h1_local = nx*ny*nz
    integer, parameter :: hcurl_local = (nx - 1)*ny*nz + &
        nx*(ny - 1)*nz + nx*ny*(nz - 1)
    integer, parameter :: hdiv_local = nx*(ny - 1)*(nz - 1) + &
        (nx - 1)*ny*(nz - 1) + (nx - 1)*(ny - 1)*nz
    integer, parameter :: l2_local = (nx - 1)*(ny - 1)*(nz - 1)
    integer :: nx_all(patch_count), ny_all(patch_count), nz_all(patch_count)
    integer :: left_patch(3), right_patch(3), face_left(3), face_right(3)
    logical :: swap_axes(3), reverse_u(3), reverse_v(3)
    integer, allocatable :: h1_offsets(:), h1_map(:)
    integer, allocatable :: hcurl_offsets(:), hcurl_map(:)
    integer, allocatable :: hdiv_offsets(:), hdiv_map(:)
    integer, allocatable :: l2_offsets(:), l2_map(:)
    integer, allocatable :: h1_left(:), h1_right(:)
    integer, allocatable :: hcurl_left(:), hcurl_right(:), hcurl_sign(:)
    integer, allocatable :: hdiv_left(:), hdiv_right(:), hdiv_sign(:)
    integer :: h1_count, hcurl_count, hdiv_count, l2_count, status
    integer :: interface, trace, left, right, trace_status
    integer :: h1_distinct, hcurl_distinct, hdiv_distinct, l2_distinct
    logical :: all_passed

    all_passed = .true.
    nx_all = nx
    ny_all = ny
    nz_all = nz
    left_patch = [1, 2, 3]
    right_patch = [2, 3, 1]
    face_left = [BSPLINE_FACE_Z_MAX, BSPLINE_FACE_Z_MAX, BSPLINE_FACE_Z_MAX]
    face_right = [BSPLINE_FACE_Z_MIN, BSPLINE_FACE_Z_MIN, BSPLINE_FACE_Z_MIN]
    swap_axes = [.false., .true., .false.]
    reverse_u = [.false., .false., .true.]
    reverse_v = [.false., .true., .false.]

    call build_bspline_feec_3d_multipatch_maps( &
        nx_all, ny_all, nz_all, left_patch, right_patch, face_left, face_right, &
        swap_axes, reverse_u, reverse_v, h1_offsets, h1_map, h1_count, &
        hcurl_offsets, hcurl_map, hcurl_count, hdiv_offsets, hdiv_map, &
        hdiv_count, l2_offsets, l2_map, l2_count, status)

    call record_condition(status == 0, "arbitrary 3D B-spline patch graph builds")
    call record_condition(all(h1_offsets == [1, 28, 55, 82]) .and. &
        all(hcurl_offsets == [1, 55, 109, 163]) .and. &
        all(hdiv_offsets == [1, 37, 73, 109]) .and. &
        all(l2_offsets == [1, 9, 17, 25]), &
        "3D packed offsets use the tensor-product FEEC local ordering")
    call record_condition(h1_count == 54 .and. hcurl_count == 126 .and. &
        hdiv_count == 96 .and. l2_count == 24 .and. all(h1_map /= 0) .and. &
        all(hcurl_map /= 0) .and. all(hdiv_map /= 0) .and. all(l2_map /= 0), &
        "3D face identifications reduce independent H1/H(curl)/H(div) DOFs")

    h1_distinct = count_unique_abs(h1_map)
    hcurl_distinct = count_unique_abs(hcurl_map)
    hdiv_distinct = count_unique_abs(hdiv_map)
    l2_distinct = count_unique_abs(l2_map)
    call record_condition(h1_distinct == h1_count .and. &
        hcurl_distinct == hcurl_count .and. hdiv_distinct == hdiv_count .and. &
        l2_distinct == l2_count, &
        "3D packed maps have one global ID per signed equivalence class")

    do interface = 1, size(left_patch)
        call build_bspline_feec_3d_interface_dofs( &
            nx, ny, nz, nx, ny, nz, face_left(interface), face_right(interface), &
            swap_axes(interface), reverse_u(interface), reverse_v(interface), &
            h1_left, h1_right, hcurl_left, hcurl_right, hcurl_sign, &
            hdiv_left, hdiv_right, hdiv_sign, trace_status)
        call record_condition(trace_status == 0, &
            "3D interface extractor accepts every graph relation")
        do trace = 1, size(h1_left)
            left = h1_offsets(left_patch(interface)) + h1_left(trace) - 1
            right = h1_offsets(right_patch(interface)) + h1_right(trace) - 1
            call record_condition(h1_map(left) == h1_map(right), &
                "3D H1 face traces are identified")
        end do
        do trace = 1, size(hcurl_left)
            left = hcurl_offsets(left_patch(interface)) + hcurl_left(trace) - 1
            right = hcurl_offsets(right_patch(interface)) + hcurl_right(trace) - 1
            call record_condition(hcurl_map(left) == hcurl_sign(trace)*hcurl_map(right), &
                "3D H(curl) face orientations are preserved")
        end do
        do trace = 1, size(hdiv_left)
            left = hdiv_offsets(left_patch(interface)) + hdiv_left(trace) - 1
            right = hdiv_offsets(right_patch(interface)) + hdiv_right(trace) - 1
            call record_condition(hdiv_map(left) == hdiv_sign(trace)*hdiv_map(right), &
                "3D H(div) face orientations are preserved")
        end do
    end do

    call check_summary("arbitrary 3D B-spline multipatch maps")
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

end program test_bspline_3d_multipatch_maps
