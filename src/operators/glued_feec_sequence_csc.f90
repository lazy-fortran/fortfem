module fortfem_glued_feec_sequence_csc
    !! CSC reference assembly for signed local-to-global FEEC maps.
    !!
    !! This is the sparse companion to the dense glued FEEC contract.  The
    !! integer local-to-global maps are fixed topology; duplicate entries are
    !! compressed by FortSparse, while orientation signs are retained in the
    !! triplet values.  The path owns no mesh or constitutive data.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_is_valid, csc_t, &
        fortsparse_status_t, status_set, FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_glued_feec_sequence_csc
    public :: assemble_glued_feec_sequence_csc_jvp
    public :: assemble_glued_feec_sequence_csc_vjp

contains

    subroutine assemble_glued_feec_sequence_csc( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, &
            gradient, curl, divergence, status)
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        integer, intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        integer, intent(in) :: l2_map(:, :)
        integer, intent(in) :: scalar_count, hcurl_count, hdiv_count, l2_count
        type(csc_t), intent(out) :: gradient, curl, divergence
        type(fortsparse_status_t), intent(out) :: status

        call validate_sequence_inputs( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, status)
        if (status%code /= FORTSPARSE_OK) return
        call assemble_one_csc( &
            local_gradient, hcurl_map, scalar_map, hcurl_count, scalar_count, &
            gradient, status)
        if (status%code /= FORTSPARSE_OK) return
        call assemble_one_csc( &
            local_curl, hdiv_map, hcurl_map, hdiv_count, hcurl_count, curl, status)
        if (status%code /= FORTSPARSE_OK) return
        call assemble_one_csc( &
            local_divergence, l2_map, hdiv_map, l2_count, hdiv_count, divergence, status)
    end subroutine assemble_glued_feec_sequence_csc

    subroutine assemble_glued_feec_sequence_csc_jvp( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, &
            local_gradient_dot, local_curl_dot, local_divergence_dot, gradient_dot, &
            curl_dot, divergence_dot, status)
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        integer, intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        integer, intent(in) :: l2_map(:, :)
        integer, intent(in) :: scalar_count, hcurl_count, hdiv_count, l2_count
        real(dp), intent(in) :: local_gradient_dot(:, :, :), local_curl_dot(:, :, :)
        real(dp), intent(in) :: local_divergence_dot(:, :, :)
        type(csc_t), intent(out) :: gradient_dot, curl_dot, divergence_dot
        type(fortsparse_status_t), intent(out) :: status

        call validate_sequence_inputs( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_direction( &
            local_gradient, local_curl, local_divergence, local_gradient_dot, &
            local_curl_dot, local_divergence_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "glued FEEC CSC JVP has incompatible local increments")
            return
        end if
        call assemble_one_csc( &
            local_gradient_dot, hcurl_map, scalar_map, hcurl_count, scalar_count, &
            gradient_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        call assemble_one_csc( &
            local_curl_dot, hdiv_map, hcurl_map, hdiv_count, hcurl_count, &
            curl_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        call assemble_one_csc( &
            local_divergence_dot, l2_map, hdiv_map, l2_count, hdiv_count, &
            divergence_dot, status)
    end subroutine assemble_glued_feec_sequence_csc_jvp

    subroutine assemble_glued_feec_sequence_csc_vjp( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, &
            gradient_bar, curl_bar, divergence_bar, local_gradient_bar, &
            local_curl_bar, local_divergence_bar, status)
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        integer, intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        integer, intent(in) :: l2_map(:, :)
        integer, intent(in) :: scalar_count, hcurl_count, hdiv_count, l2_count
        type(csc_t), intent(in) :: gradient_bar, curl_bar, divergence_bar
        real(dp), intent(out) :: local_gradient_bar(:, :, :), local_curl_bar(:, :, :)
        real(dp), intent(out) :: local_divergence_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        local_gradient_bar = 0.0_dp
        local_curl_bar = 0.0_dp
        local_divergence_bar = 0.0_dp
        call validate_sequence_inputs( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_cotangent( &
            gradient_bar, curl_bar, divergence_bar, scalar_count, hcurl_count, &
            hdiv_count, l2_count)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "glued FEEC CSC VJP has invalid matrix cotangents")
            return
        end if
        if (any(shape(local_gradient_bar) /= shape(local_gradient)) .or. &
            any(shape(local_curl_bar) /= shape(local_curl)) .or. &
            any(shape(local_divergence_bar) /= shape(local_divergence))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "glued FEEC CSC VJP has incompatible local cotangents")
            return
        end if
        call scatter_one_bar( &
            gradient_bar, hcurl_map, scalar_map, local_gradient_bar)
        call scatter_one_bar(curl_bar, hdiv_map, hcurl_map, local_curl_bar)
        call scatter_one_bar( &
            divergence_bar, l2_map, hdiv_map, local_divergence_bar)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_glued_feec_sequence_csc_vjp

    subroutine assemble_one_csc( &
            local_matrix, row_map, column_map, row_count, column_count, matrix, status)
        real(dp), intent(in) :: local_matrix(:, :, :)
        integer, intent(in) :: row_map(:, :), column_map(:, :)
        integer, intent(in) :: row_count, column_count
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: rows(:), columns(:)
        real(dp), allocatable :: values(:)
        integer :: local_rows, local_columns, cells, triplet_count
        integer :: cell, local_row, local_column, entry

        local_rows = size(local_matrix, 1)
        local_columns = size(local_matrix, 2)
        cells = size(local_matrix, 3)
        triplet_count = local_rows*local_columns*cells
        allocate(rows(triplet_count), columns(triplet_count), values(triplet_count))
        entry = 0
        do cell = 1, cells
            do local_row = 1, local_rows
                do local_column = 1, local_columns
                    entry = entry + 1
                    rows(entry) = abs(row_map(local_row, cell))
                    columns(entry) = abs(column_map(local_column, cell))
                    values(entry) = signed_product( &
                        row_map(local_row, cell), column_map(local_column, cell)) * &
                        local_matrix(local_row, local_column, cell)
                end do
            end do
        end do
        call csc_from_triplet( &
            row_count, column_count, rows, columns, values, matrix, status)
    end subroutine assemble_one_csc

    subroutine scatter_one_bar(matrix_bar, row_map, column_map, local_bar)
        type(csc_t), intent(in) :: matrix_bar
        integer, intent(in) :: row_map(:, :), column_map(:, :)
        real(dp), intent(out) :: local_bar(:, :, :)
        integer :: cell, local_row, local_column

        local_bar = 0.0_dp
        do cell = 1, size(local_bar, 3)
            do local_row = 1, size(local_bar, 1)
                do local_column = 1, size(local_bar, 2)
                    local_bar(local_row, local_column, cell) = &
                        signed_product(row_map(local_row, cell), column_map(local_column, cell)) * &
                        csc_entry(matrix_bar, abs(row_map(local_row, cell)), &
                        abs(column_map(local_column, cell)))
                end do
            end do
        end do
    end subroutine scatter_one_bar

    subroutine validate_sequence_inputs( &
            local_gradient, local_curl, local_divergence, scalar_map, hcurl_map, &
            hdiv_map, l2_map, scalar_count, hcurl_count, hdiv_count, l2_count, status)
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        integer, intent(in) :: scalar_map(:, :), hcurl_map(:, :), hdiv_map(:, :)
        integer, intent(in) :: l2_map(:, :)
        integer, intent(in) :: scalar_count, hcurl_count, hdiv_count, l2_count
        type(fortsparse_status_t), intent(out) :: status
        integer :: cells

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "glued FEEC CSC maps have incompatible arrays")
        cells = size(local_gradient, 3)
        if (scalar_count < 1 .or. hcurl_count < 1 .or. hdiv_count < 1 .or. &
            l2_count < 1 .or. cells < 1 .or. size(local_gradient, 1) < 1 .or. &
            size(local_gradient, 2) < 1 .or. size(local_curl, 1) < 1 .or. &
            size(local_divergence, 1) < 1) return
        if (size(local_curl, 2) /= size(local_gradient, 1) .or. &
            size(local_curl, 3) /= cells .or. size(local_divergence, 2) /= size(local_curl, 1) .or. &
            size(local_divergence, 3) /= cells .or. size(scalar_map, 1) /= size(local_gradient, 2) .or. &
            size(scalar_map, 2) /= cells .or. size(hcurl_map, 1) /= size(local_gradient, 1) .or. &
            size(hcurl_map, 2) /= cells .or. size(hdiv_map, 1) /= size(local_curl, 1) .or. &
            size(hdiv_map, 2) /= cells .or. size(l2_map, 1) /= size(local_divergence, 1) .or. &
            size(l2_map, 2) /= cells) return
        if (any(abs(scalar_map) < 1) .or. any(abs(scalar_map) > scalar_count) .or. &
            any(abs(hcurl_map) < 1) .or. any(abs(hcurl_map) > hcurl_count) .or. &
            any(abs(hdiv_map) < 1) .or. any(abs(hdiv_map) > hdiv_count) .or. &
            any(abs(l2_map) < 1) .or. any(abs(l2_map) > l2_count)) return
        if (.not. all(ieee_is_finite(local_gradient)) .or. &
            .not. all(ieee_is_finite(local_curl)) .or. &
            .not. all(ieee_is_finite(local_divergence))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_sequence_inputs

    logical function valid_direction( &
            local_gradient, local_curl, local_divergence, local_gradient_dot, &
            local_curl_dot, local_divergence_dot) result(valid)
        real(dp), intent(in) :: local_gradient(:, :, :), local_curl(:, :, :)
        real(dp), intent(in) :: local_divergence(:, :, :)
        real(dp), intent(in) :: local_gradient_dot(:, :, :), local_curl_dot(:, :, :)
        real(dp), intent(in) :: local_divergence_dot(:, :, :)

        valid = all(shape(local_gradient_dot) == shape(local_gradient)) .and. &
            all(shape(local_curl_dot) == shape(local_curl)) .and. &
            all(shape(local_divergence_dot) == shape(local_divergence)) .and. &
            all(ieee_is_finite(local_gradient_dot)) .and. &
            all(ieee_is_finite(local_curl_dot)) .and. &
            all(ieee_is_finite(local_divergence_dot))
    end function valid_direction

    logical function valid_cotangent( &
            gradient_bar, curl_bar, divergence_bar, scalar_count, hcurl_count, &
            hdiv_count, l2_count) result(valid)
        type(csc_t), intent(in) :: gradient_bar, curl_bar, divergence_bar
        integer, intent(in) :: scalar_count, hcurl_count, hdiv_count, l2_count

        valid = csc_is_valid(gradient_bar) .and. csc_is_valid(curl_bar) .and. &
            csc_is_valid(divergence_bar) .and. gradient_bar%nrow == hcurl_count .and. &
            gradient_bar%ncol == scalar_count .and. curl_bar%nrow == hdiv_count .and. &
            curl_bar%ncol == hcurl_count .and. divergence_bar%nrow == l2_count .and. &
            divergence_bar%ncol == hdiv_count .and. finite_csc(gradient_bar) .and. &
            finite_csc(curl_bar) .and. finite_csc(divergence_bar)
    end function valid_cotangent

    logical function finite_csc(matrix) result(valid)
        type(csc_t), intent(in) :: matrix

        valid = .true.
        if (matrix%nnz > 0) valid = all(ieee_is_finite(matrix%val))
    end function finite_csc

    pure real(dp) function signed_product(first, second) result(value)
        integer, intent(in) :: first, second

        value = real(sign(1, first)*sign(1, second), dp)
    end function signed_product

    pure real(dp) function csc_entry(matrix, row, column) result(value)
        type(csc_t), intent(in) :: matrix
        integer, intent(in) :: row, column
        integer :: entry

        value = 0.0_dp
        do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
            if (matrix%row_idx(entry) == row) then
                value = matrix%val(entry)
                return
            end if
        end do
    end function csc_entry

end module fortfem_glued_feec_sequence_csc
