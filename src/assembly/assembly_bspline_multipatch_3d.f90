module fortfem_assembly_bspline_multipatch_3d
    !! Sparse quotient-complex assembly for two conforming spline volumes.
    use fortfem_bspline_feec, only: build_bspline_feec_3d_operators
    use fortfem_bspline_multipatch, only: &
        BSPLINE_FACE_X_MAX, BSPLINE_FACE_X_MIN, BSPLINE_FACE_Y_MAX, &
        BSPLINE_FACE_Y_MIN, build_bspline_feec_3d_two_patch_maps
    use fortfem_kinds, only: dp
    use fortsparse, only: &
        csc_from_triplet, csc_t, FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: build_bspline_feec_3d_two_patch_operators_csc

contains

    subroutine build_bspline_feec_3d_two_patch_operators_csc( &
            knots_x_left, knots_y_left, knots_z_left, &
            degree_x_left, degree_y_left, degree_z_left, &
            knots_x_right, knots_y_right, knots_z_right, &
            degree_x_right, degree_y_right, degree_z_right, face_left, &
            face_right, swap_axes, reverse_u, reverse_v, gradient, curl, &
            divergence, status)
        real(dp), intent(in) :: knots_x_left(:), knots_y_left(:)
        real(dp), intent(in) :: knots_z_left(:)
        real(dp), intent(in) :: knots_x_right(:), knots_y_right(:)
        real(dp), intent(in) :: knots_z_right(:)
        integer, intent(in) :: degree_x_left, degree_y_left, degree_z_left
        integer, intent(in) :: degree_x_right, degree_y_right, degree_z_right
        integer, intent(in) :: face_left, face_right
        logical, intent(in) :: swap_axes, reverse_u, reverse_v
        type(csc_t), intent(out) :: gradient, curl, divergence
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), h1_left(:), h1_right(:)
        integer, allocatable :: hcurl_left(:), hcurl_right(:)
        integer, allocatable :: hcurl_left_sign(:), hcurl_right_sign(:)
        integer, allocatable :: hdiv_left(:), hdiv_right(:)
        integer, allocatable :: hdiv_left_sign(:), hdiv_right_sign(:)
        integer, allocatable :: l2_left(:), l2_right(:), rows(:)
        real(dp), allocatable :: curl_left(:, :), curl_right(:, :)
        real(dp), allocatable :: divergence_left(:, :), divergence_right(:, :)
        real(dp), allocatable :: gradient_left(:, :), gradient_right(:, :)
        real(dp), allocatable :: values(:)
        integer :: entry, local_status, nx_left, nx_right, ny_left, ny_right
        integer :: nz_left, nz_right

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse two-patch 3D isogeometric FEEC assembly failed")
        nx_left = size(knots_x_left) - degree_x_left - 1
        ny_left = size(knots_y_left) - degree_y_left - 1
        nz_left = size(knots_z_left) - degree_z_left - 1
        nx_right = size(knots_x_right) - degree_x_right - 1
        ny_right = size(knots_y_right) - degree_y_right - 1
        nz_right = size(knots_z_right) - degree_z_right - 1
        if (.not. compatible_trace_spaces()) return
        call build_bspline_feec_3d_two_patch_maps( &
            nx_left, ny_left, nz_left, nx_right, ny_right, nz_right, &
            face_left, face_right, swap_axes, reverse_u, reverse_v, h1_left, &
            h1_right, hcurl_left, hcurl_right, hcurl_left_sign, &
            hcurl_right_sign, hdiv_left, hdiv_right, hdiv_left_sign, &
            hdiv_right_sign, l2_left, l2_right, local_status)
        if (local_status /= 0) return
        call build_bspline_feec_3d_operators( &
            knots_x_left, knots_y_left, knots_z_left, degree_x_left, &
            degree_y_left, degree_z_left, gradient_left, curl_left, &
            divergence_left, local_status)
        if (local_status /= 0) return
        call build_bspline_feec_3d_operators( &
            knots_x_right, knots_y_right, knots_z_right, degree_x_right, &
            degree_y_right, degree_z_right, gradient_right, curl_right, &
            divergence_right, local_status)
        if (local_status /= 0) return

        allocate(rows(size(gradient_left) + size(gradient_right)))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        call append_mapped_matrix( &
            gradient_left, hcurl_left, h1_left, hcurl_left_sign, rows=rows, &
            columns=columns, values=values, entry=entry)
        call append_mapped_matrix( &
            gradient_right, hcurl_right, h1_right, hcurl_right_sign, &
            rows=rows, columns=columns, values=values, entry=entry, &
            skip_mapped_rows_through=size(hcurl_left))
        call csc_from_triplet( &
            max(maxval(hcurl_left), maxval(hcurl_right)), &
            max(maxval(h1_left), maxval(h1_right)), rows(:entry), &
            columns(:entry), values(:entry), gradient, status)
        if (status%code /= 0) return

        deallocate(rows, columns, values)
        allocate(rows(size(curl_left) + size(curl_right)))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        call append_mapped_matrix( &
            curl_left, hdiv_left, hcurl_left, hdiv_left_sign, &
            hcurl_left_sign, rows, columns, values, entry)
        call append_mapped_matrix( &
            curl_right, hdiv_right, hcurl_right, hdiv_right_sign, &
            hcurl_right_sign, rows, columns, values, entry, &
            skip_mapped_rows_through=size(hdiv_left))
        call csc_from_triplet( &
            max(maxval(hdiv_left), maxval(hdiv_right)), &
            max(maxval(hcurl_left), maxval(hcurl_right)), rows(:entry), &
            columns(:entry), values(:entry), curl, status)
        if (status%code /= 0) return

        deallocate(rows, columns, values)
        allocate(rows(size(divergence_left) + size(divergence_right)))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        call append_mapped_matrix( &
            divergence_left, l2_left, hdiv_left, column_signs=hdiv_left_sign, &
            rows=rows, columns=columns, values=values, entry=entry)
        call append_mapped_matrix( &
            divergence_right, l2_right, hdiv_right, &
            column_signs=hdiv_right_sign, rows=rows, columns=columns, &
            values=values, entry=entry)
        call csc_from_triplet( &
            max(maxval(l2_left), maxval(l2_right)), &
            max(maxval(hdiv_left), maxval(hdiv_right)), rows(:entry), &
            columns(:entry), values(:entry), divergence, status)

    contains

        logical function compatible_trace_spaces() result(compatible)
            real(dp), allocatable :: left_u(:), left_v(:)
            real(dp), allocatable :: right_u(:), right_v(:)
            integer :: degree_left_u, degree_left_v
            integer :: degree_right_u, degree_right_v

            call trace_spaces( &
                .true., face_left, left_u, left_v, degree_left_u, degree_left_v)
            call trace_spaces( &
                .false., face_right, right_u, right_v, degree_right_u, &
                degree_right_v)
            if (swap_axes) then
                compatible = matching_space( &
                    left_u, degree_left_u, right_v, degree_right_v, &
                    reverse_v) .and. matching_space( &
                    left_v, degree_left_v, right_u, degree_right_u, reverse_u)
            else
                compatible = matching_space( &
                    left_u, degree_left_u, right_u, degree_right_u, &
                    reverse_u) .and. matching_space( &
                    left_v, degree_left_v, right_v, degree_right_v, reverse_v)
            end if
        end function compatible_trace_spaces

        subroutine trace_spaces( &
                left_patch, face, knots_u, knots_v, degree_u, degree_v)
            logical, intent(in) :: left_patch
            integer, intent(in) :: face
            real(dp), allocatable, intent(out) :: knots_u(:), knots_v(:)
            integer, intent(out) :: degree_u, degree_v

            if (left_patch) then
                select case (face)
                case (BSPLINE_FACE_X_MIN, BSPLINE_FACE_X_MAX)
                    knots_u = knots_y_left
                    knots_v = knots_z_left
                    degree_u = degree_y_left
                    degree_v = degree_z_left
                case (BSPLINE_FACE_Y_MIN, BSPLINE_FACE_Y_MAX)
                    knots_u = knots_x_left
                    knots_v = knots_z_left
                    degree_u = degree_x_left
                    degree_v = degree_z_left
                case default
                    knots_u = knots_x_left
                    knots_v = knots_y_left
                    degree_u = degree_x_left
                    degree_v = degree_y_left
                end select
            else
                select case (face)
                case (BSPLINE_FACE_X_MIN, BSPLINE_FACE_X_MAX)
                    knots_u = knots_y_right
                    knots_v = knots_z_right
                    degree_u = degree_y_right
                    degree_v = degree_z_right
                case (BSPLINE_FACE_Y_MIN, BSPLINE_FACE_Y_MAX)
                    knots_u = knots_x_right
                    knots_v = knots_z_right
                    degree_u = degree_x_right
                    degree_v = degree_z_right
                case default
                    knots_u = knots_x_right
                    knots_v = knots_y_right
                    degree_u = degree_x_right
                    degree_v = degree_y_right
                end select
            end if
        end subroutine trace_spaces

        logical function matching_space( &
                left_knots, left_degree, right_knots, right_degree, reversed) &
                result(matches)
            real(dp), intent(in) :: left_knots(:), right_knots(:)
            integer, intent(in) :: left_degree, right_degree
            logical, intent(in) :: reversed
            real(dp), allocatable :: normalized_left(:), normalized_right(:)

            matches = left_degree == right_degree .and. &
                size(left_knots) == size(right_knots)
            if (.not. matches) return
            normalized_left = (left_knots - left_knots(1))/ &
                (left_knots(size(left_knots)) - left_knots(1))
            normalized_right = (right_knots - right_knots(1))/ &
                (right_knots(size(right_knots)) - right_knots(1))
            if (reversed) then
                normalized_right = &
                    1.0_dp - normalized_right(size(normalized_right):1:-1)
            end if
            matches = maxval(abs(normalized_left - normalized_right)) <= &
                128.0_dp*epsilon(1.0_dp)
        end function matching_space
    end subroutine build_bspline_feec_3d_two_patch_operators_csc

    pure subroutine append_mapped_matrix( &
            local_matrix, row_map, column_map, row_signs, column_signs, rows, &
            columns, values, entry, skip_mapped_rows_through)
        real(dp), intent(in) :: local_matrix(:, :)
        integer, intent(in) :: row_map(:), column_map(:)
        integer, intent(in), optional :: row_signs(:), column_signs(:)
        integer, intent(inout) :: rows(:), columns(:), entry
        real(dp), intent(inout) :: values(:)
        integer, intent(in), optional :: skip_mapped_rows_through

        integer :: column, row, sign

        do column = 1, size(local_matrix, 2)
            do row = 1, size(local_matrix, 1)
                if (local_matrix(row, column) == 0.0_dp) cycle
                if (present(skip_mapped_rows_through)) then
                    if (row_map(row) <= skip_mapped_rows_through) cycle
                end if
                entry = entry + 1
                rows(entry) = row_map(row)
                columns(entry) = column_map(column)
                sign = 1
                if (present(row_signs)) sign = sign*row_signs(row)
                if (present(column_signs)) sign = sign*column_signs(column)
                values(entry) = real(sign, dp)*local_matrix(row, column)
            end do
        end do
    end subroutine append_mapped_matrix

end module fortfem_assembly_bspline_multipatch_3d
