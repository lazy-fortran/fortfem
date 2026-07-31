module fortfem_hdg_global_skeleton_csc
    !! Sparse CSC assembly for oriented local HDG skeleton blocks.
    !!
    !! The local-to-global map and signs are frozen topology.  Duplicate
    !! entries are deliberately handed to FortSparse's triplet compressor so
    !! the sparse path has the same signed-map contract as the dense reference
    !! assembler, while retaining a canonical CSC pattern for direct solves.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_is_valid, csc_t, &
        fortsparse_status_t, status_set, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_OK
    implicit none
    private

    public :: assemble_hdg_global_skeleton_csc
    public :: assemble_hdg_global_skeleton_csc_jvp
    public :: assemble_hdg_global_skeleton_csc_vjp

contains

    subroutine assemble_hdg_global_skeleton_csc( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix, global_rhs, status)
        !! Assemble signed local blocks into a duplicate-compressed CSC matrix.
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        type(csc_t), intent(out) :: global_matrix
        real(dp), intent(out) :: global_rhs(:)
        type(fortsparse_status_t), intent(out) :: status

        call assemble_csc_data(local_matrix, local_rhs, local_to_global, &
            local_sign, global_count, global_matrix, global_rhs, status)
    end subroutine assemble_hdg_global_skeleton_csc

    subroutine assemble_hdg_global_skeleton_csc_jvp( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            local_matrix_dot, local_rhs_dot, global_matrix_dot, global_rhs_dot, &
            status)
        !! Apply the fixed-topology JVP of sparse skeleton assembly.
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        real(dp), intent(in) :: local_matrix_dot(:, :, :), local_rhs_dot(:, :)
        type(csc_t), intent(out) :: global_matrix_dot
        real(dp), intent(out) :: global_rhs_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: local_count, element_count

        global_rhs_dot = 0.0_dp
        call validate_skeleton_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_rhs_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        local_count = size(local_matrix, 1)
        element_count = size(local_matrix, 3)
        if (.not. valid_direction(local_matrix_dot, local_rhs_dot, local_count, &
                element_count)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "HDG sparse skeleton JVP has incompatible increments")
            return
        end if
        call assemble_csc_data(local_matrix_dot, local_rhs_dot, local_to_global, &
            local_sign, global_count, global_matrix_dot, global_rhs_dot, status)
    end subroutine assemble_hdg_global_skeleton_csc_jvp

    subroutine assemble_hdg_global_skeleton_csc_vjp( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix_bar, global_rhs_bar, local_matrix_bar, local_rhs_bar, &
            status)
        !! Apply the real reverse product of sparse skeleton assembly.
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        type(csc_t), intent(in) :: global_matrix_bar
        real(dp), intent(in) :: global_rhs_bar(:)
        real(dp), intent(out) :: local_matrix_bar(:, :, :), local_rhs_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: local_count, element_count
        integer :: element, local_row, local_column, row, column

        local_matrix_bar = 0.0_dp
        local_rhs_bar = 0.0_dp
        call validate_skeleton_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_rhs_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. csc_is_valid(global_matrix_bar) .or. &
            global_matrix_bar%nrow /= global_count .or. &
            global_matrix_bar%ncol /= global_count .or. &
            (global_matrix_bar%nnz > 0 .and. &
            .not. all(ieee_is_finite(global_matrix_bar%val)))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "HDG sparse skeleton VJP has an invalid matrix cotangent")
            return
        end if
        local_count = size(local_matrix, 1)
        element_count = size(local_matrix, 3)
        if (size(local_matrix_bar, 1) /= local_count .or. &
            size(local_matrix_bar, 2) /= local_count .or. &
            size(local_matrix_bar, 3) /= element_count .or. &
            size(local_rhs_bar, 1) /= local_count .or. &
            size(local_rhs_bar, 2) /= element_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "HDG sparse skeleton VJP has incompatible cotangents")
            return
        end if
        do element = 1, element_count
            do local_row = 1, local_count
                row = local_to_global(local_row, element)
                local_rhs_bar(local_row, element) = &
                    local_sign(local_row, element)*global_rhs_bar(row)
                do local_column = 1, local_count
                    column = local_to_global(local_column, element)
                    local_matrix_bar(local_row, local_column, element) = &
                        local_sign(local_row, element)*local_sign(local_column, element)* &
                        csc_entry(global_matrix_bar, row, column)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_hdg_global_skeleton_csc_vjp

    subroutine assemble_csc_data( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix, global_rhs, status)
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        type(csc_t), intent(out) :: global_matrix
        real(dp), intent(out) :: global_rhs(:)
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: rows(:), columns(:)
        real(dp), allocatable :: values(:)
        integer :: local_count, element_count, triplet_count
        integer :: element, local_row, local_column, row, column, entry

        global_rhs = 0.0_dp
        call validate_skeleton_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_rhs, status)
        if (status%code /= FORTSPARSE_OK) return
        local_count = size(local_matrix, 1)
        element_count = size(local_matrix, 3)
        triplet_count = local_count*local_count*element_count
        allocate(rows(triplet_count), columns(triplet_count), values(triplet_count))
        entry = 0
        do element = 1, element_count
            do local_row = 1, local_count
                row = local_to_global(local_row, element)
                global_rhs(row) = global_rhs(row) + &
                    local_sign(local_row, element)*local_rhs(local_row, element)
                do local_column = 1, local_count
                    column = local_to_global(local_column, element)
                    entry = entry + 1
                    rows(entry) = row
                    columns(entry) = column
                    values(entry) = local_sign(local_row, element)* &
                        local_sign(local_column, element)* &
                        local_matrix(local_row, local_column, element)
                end do
            end do
        end do
        call csc_from_triplet(global_count, global_count, rows, columns, values, &
            global_matrix, status)
    end subroutine assemble_csc_data

    subroutine validate_skeleton_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_rhs, status)
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        real(dp), intent(in) :: global_rhs(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: local_count, element_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "HDG sparse skeleton assembly received incompatible arrays")
        local_count = size(local_matrix, 1)
        element_count = size(local_matrix, 3)
        if (global_count < 1 .or. local_count < 1 .or. element_count < 1 .or. &
            size(local_matrix, 2) /= local_count .or. &
            size(local_rhs, 1) /= local_count .or. &
            size(local_rhs, 2) /= element_count .or. &
            size(local_to_global, 1) /= local_count .or. &
            size(local_to_global, 2) /= element_count .or. &
            size(local_sign, 1) /= local_count .or. &
            size(local_sign, 2) /= element_count .or. &
            size(global_rhs) /= global_count) return
        if (any(local_to_global < 1) .or. any(local_to_global > global_count) .or. &
            any(abs(local_sign) /= 1) .or. &
            any(.not. ieee_is_finite(local_matrix)) .or. &
            any(.not. ieee_is_finite(local_rhs)) .or. &
            any(.not. ieee_is_finite(global_rhs))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_skeleton_inputs

    pure logical function valid_direction( &
            local_matrix_dot, local_rhs_dot, local_count, element_count) result(valid)
        real(dp), intent(in) :: local_matrix_dot(:, :, :), local_rhs_dot(:, :)
        integer, intent(in) :: local_count, element_count

        valid = size(local_matrix_dot, 1) == local_count .and. &
            size(local_matrix_dot, 2) == local_count .and. &
            size(local_matrix_dot, 3) == element_count .and. &
            size(local_rhs_dot, 1) == local_count .and. &
            size(local_rhs_dot, 2) == element_count .and. &
            all(ieee_is_finite(local_matrix_dot)) .and. &
            all(ieee_is_finite(local_rhs_dot))
    end function valid_direction

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

end module fortfem_hdg_global_skeleton_csc
