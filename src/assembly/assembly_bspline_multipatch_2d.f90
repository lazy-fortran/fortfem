module fortfem_assembly_bspline_multipatch_2d
    !! Sparse quotient-complex assembly for two conforming spline patches.
    use fortfem_bspline_feec, only: build_bspline_feec_2d_operators
    use fortfem_bspline_multipatch, only: &
        BSPLINE_FACE_X_MAX, BSPLINE_FACE_X_MIN, &
        build_bspline_feec_2d_two_patch_maps
    use fortfem_kinds, only: dp
    use fortsparse, only: &
        csc_from_triplet, csc_t, FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: build_bspline_feec_2d_two_patch_operators_csc

contains

    subroutine build_bspline_feec_2d_two_patch_operators_csc( &
            knots_x_left, knots_y_left, degree_x_left, degree_y_left, &
            knots_x_right, knots_y_right, degree_x_right, degree_y_right, &
            face_left, face_right, reversed, gradient, curl, status)
        real(dp), intent(in) :: knots_x_left(:), knots_y_left(:)
        real(dp), intent(in) :: knots_x_right(:), knots_y_right(:)
        integer, intent(in) :: degree_x_left, degree_y_left
        integer, intent(in) :: degree_x_right, degree_y_right
        integer, intent(in) :: face_left, face_right
        logical, intent(in) :: reversed
        type(csc_t), intent(out) :: gradient, curl
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), h1_left(:), h1_right(:)
        integer, allocatable :: hcurl_left(:), hcurl_right(:)
        integer, allocatable :: hcurl_left_sign(:), hcurl_right_sign(:)
        integer, allocatable :: l2_left(:), l2_right(:), rows(:)
        real(dp), allocatable :: curl_left(:, :), curl_right(:, :)
        real(dp), allocatable :: gradient_left(:, :), gradient_right(:, :)
        real(dp), allocatable :: values(:)
        integer :: entry, local_status, nx_left, nx_right, ny_left, ny_right

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse two-patch isogeometric FEEC assembly failed")
        nx_left = size(knots_x_left) - degree_x_left - 1
        ny_left = size(knots_y_left) - degree_y_left - 1
        nx_right = size(knots_x_right) - degree_x_right - 1
        ny_right = size(knots_y_right) - degree_y_right - 1
        if (.not. compatible_trace_spaces()) return
        call build_bspline_feec_2d_two_patch_maps( &
            nx_left, ny_left, nx_right, ny_right, face_left, face_right, &
            reversed, h1_left, h1_right, hcurl_left, hcurl_right, &
            hcurl_left_sign, hcurl_right_sign, l2_left, l2_right, local_status)
        if (local_status /= 0) return
        call build_bspline_feec_2d_operators( &
            knots_x_left, knots_y_left, degree_x_left, degree_y_left, &
            gradient_left, curl_left, local_status)
        if (local_status /= 0) return
        call build_bspline_feec_2d_operators( &
            knots_x_right, knots_y_right, degree_x_right, degree_y_right, &
            gradient_right, curl_right, local_status)
        if (local_status /= 0) return

        allocate(rows(size(gradient_left) + size(gradient_right)))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        call append_mapped_matrix( &
            gradient_left, hcurl_left, h1_left, hcurl_left_sign, rows, &
            columns, values, entry)
        call append_mapped_matrix( &
            gradient_right, hcurl_right, h1_right, hcurl_right_sign, rows, &
            columns, values, entry, skip_mapped_rows_through= &
            size(hcurl_left))
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
            curl_left, l2_left, hcurl_left, hcurl_left_sign, rows, columns, &
            values, entry, sign_columns=.true.)
        call append_mapped_matrix( &
            curl_right, l2_right, hcurl_right, hcurl_right_sign, rows, &
            columns, values, entry, sign_columns=.true.)
        call csc_from_triplet( &
            max(maxval(l2_left), maxval(l2_right)), &
            max(maxval(hcurl_left), maxval(hcurl_right)), rows(:entry), &
            columns(:entry), values(:entry), curl, status)

    contains

        logical function compatible_trace_spaces() result(compatible)
            real(dp), allocatable :: left(:), right(:)
            integer :: left_degree, right_degree

            if (face_left == BSPLINE_FACE_X_MIN .or. &
                face_left == BSPLINE_FACE_X_MAX) then
                left = knots_y_left
                left_degree = degree_y_left
            else
                left = knots_x_left
                left_degree = degree_x_left
            end if
            if (face_right == BSPLINE_FACE_X_MIN .or. &
                face_right == BSPLINE_FACE_X_MAX) then
                right = knots_y_right
                right_degree = degree_y_right
            else
                right = knots_x_right
                right_degree = degree_x_right
            end if
            compatible = left_degree == right_degree .and. &
                size(left) == size(right)
            if (.not. compatible) return
            left = (left - left(1))/(left(size(left)) - left(1))
            right = (right - right(1))/(right(size(right)) - right(1))
            if (reversed) right = 1.0_dp - right(size(right):1:-1)
            compatible = maxval(abs(left - right)) <= &
                128.0_dp*epsilon(1.0_dp)
        end function compatible_trace_spaces
    end subroutine build_bspline_feec_2d_two_patch_operators_csc

    pure subroutine append_mapped_matrix( &
            local_matrix, row_map, column_map, signs, rows, columns, values, &
            entry, sign_columns, skip_mapped_rows_through)
        real(dp), intent(in) :: local_matrix(:, :)
        integer, intent(in) :: row_map(:), column_map(:), signs(:)
        integer, intent(inout) :: rows(:), columns(:), entry
        real(dp), intent(inout) :: values(:)
        logical, intent(in), optional :: sign_columns
        integer, intent(in), optional :: skip_mapped_rows_through

        integer :: column, row
        logical :: use_column_sign

        use_column_sign = .false.
        if (present(sign_columns)) use_column_sign = sign_columns
        do column = 1, size(local_matrix, 2)
            do row = 1, size(local_matrix, 1)
                if (local_matrix(row, column) == 0.0_dp) cycle
                if (present(skip_mapped_rows_through)) then
                    if (row_map(row) <= skip_mapped_rows_through) cycle
                end if
                entry = entry + 1
                rows(entry) = row_map(row)
                columns(entry) = column_map(column)
                if (use_column_sign) then
                    values(entry) = real(signs(column), dp)* &
                        local_matrix(row, column)
                else
                    values(entry) = real(signs(row), dp)* &
                        local_matrix(row, column)
                end if
            end do
        end do
    end subroutine append_mapped_matrix

end module fortfem_assembly_bspline_multipatch_2d
