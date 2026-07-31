module fortfem_hdg_global_skeleton
    !! Dense reference assembly for oriented local HDG skeleton blocks.
    !!
    !! The integer local-to-global map and signs are frozen topology.  Each
    !! local block contributes with one sign on each row/column, which makes
    !! the reference path applicable to scalar, H(curl), H(div), and IGA
    !! skeleton traces before a sparse backend is selected.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_hdg_global_skeleton
    public :: assemble_hdg_global_skeleton_jvp
    public :: assemble_hdg_global_skeleton_vjp

contains

    subroutine assemble_hdg_global_skeleton( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix, global_rhs, status)
        !! Assemble all signed local trace blocks into a dense global matrix.
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        real(dp), intent(out) :: global_matrix(:, :), global_rhs(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: local_count, element_count, element, local_row, local_column
        integer :: row, column

        global_matrix = 0.0_dp
        global_rhs = 0.0_dp
        call validate_skeleton_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix, global_rhs, status)
        if (status%code /= FORTSPARSE_OK) return
        local_count = size(local_matrix, 1)
        element_count = size(local_matrix, 3)
        do element = 1, element_count
            do local_row = 1, local_count
                row = local_to_global(local_row, element)
                global_rhs(row) = global_rhs(row) + &
                    local_sign(local_row, element)*local_rhs(local_row, element)
                do local_column = 1, local_count
                    column = local_to_global(local_column, element)
                    global_matrix(row, column) = global_matrix(row, column) + &
                        local_sign(local_row, element)*local_sign(local_column, element)* &
                        local_matrix(local_row, local_column, element)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_hdg_global_skeleton

    subroutine assemble_hdg_global_skeleton_jvp( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            local_matrix_dot, local_rhs_dot, global_matrix_dot, global_rhs_dot, status)
        !! Apply the frozen-topology JVP of global skeleton assembly.
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        real(dp), intent(in) :: local_matrix_dot(:, :, :), local_rhs_dot(:, :)
        real(dp), intent(out) :: global_matrix_dot(:, :), global_rhs_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: local_count, element_count, element, local_row, local_column
        integer :: row, column

        global_matrix_dot = 0.0_dp
        global_rhs_dot = 0.0_dp
        call validate_skeleton_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix_dot, global_rhs_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        local_count = size(local_matrix, 1)
        element_count = size(local_matrix, 3)
        if (.not. valid_skeleton_direction( &
            local_matrix_dot, local_rhs_dot, local_count, element_count)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "HDG skeleton assembly JVP has incompatible increments")
            return
        end if
        do element = 1, element_count
            do local_row = 1, local_count
                row = local_to_global(local_row, element)
                global_rhs_dot(row) = global_rhs_dot(row) + &
                    local_sign(local_row, element)*local_rhs_dot(local_row, element)
                do local_column = 1, local_count
                    column = local_to_global(local_column, element)
                    global_matrix_dot(row, column) = &
                        global_matrix_dot(row, column) + &
                        local_sign(local_row, element)*local_sign(local_column, element)* &
                        local_matrix_dot(local_row, local_column, element)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_hdg_global_skeleton_jvp

    subroutine assemble_hdg_global_skeleton_vjp( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix_bar, global_rhs_bar, local_matrix_bar, local_rhs_bar, status)
        !! Apply the real reverse product of global skeleton assembly.
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        real(dp), intent(in) :: global_matrix_bar(:, :), global_rhs_bar(:)
        real(dp), intent(out) :: local_matrix_bar(:, :, :), local_rhs_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: local_count, element_count, element, local_row, local_column
        integer :: row, column

        local_matrix_bar = 0.0_dp
        local_rhs_bar = 0.0_dp
        call validate_skeleton_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix_bar, global_rhs_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        local_count = size(local_matrix, 1)
        element_count = size(local_matrix, 3)
        if (size(local_matrix_bar, 1) /= local_count .or. &
            size(local_matrix_bar, 2) /= local_count .or. &
            size(local_matrix_bar, 3) /= element_count .or. &
            size(local_rhs_bar, 1) /= local_count .or. &
            size(local_rhs_bar, 2) /= element_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "HDG skeleton assembly VJP has incompatible cotangents")
            return
        end if
        do element = 1, element_count
            do local_row = 1, local_count
                row = local_to_global(local_row, element)
                local_rhs_bar(local_row, element) = local_rhs_bar(local_row, element) + &
                    local_sign(local_row, element)*global_rhs_bar(row)
                do local_column = 1, local_count
                    column = local_to_global(local_column, element)
                    local_matrix_bar(local_row, local_column, element) = &
                        local_sign(local_row, element)*local_sign(local_column, element)* &
                        global_matrix_bar(row, column)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_hdg_global_skeleton_vjp

    pure logical function valid_skeleton_direction( &
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
    end function valid_skeleton_direction

    subroutine validate_skeleton_inputs( &
            local_matrix, local_rhs, local_to_global, local_sign, global_count, &
            global_matrix, global_rhs, status)
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer, intent(in) :: local_to_global(:, :), local_sign(:, :)
        integer, intent(in) :: global_count
        real(dp), intent(in) :: global_matrix(:, :), global_rhs(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: local_count, element_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "HDG skeleton assembly received incompatible arrays")
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
            size(global_matrix, 1) /= global_count .or. &
            size(global_matrix, 2) /= global_count .or. &
            size(global_rhs) /= global_count) return
        if (any(local_to_global < 1) .or. any(local_to_global > global_count) .or. &
            any(abs(local_sign) /= 1) .or. &
            any(.not. ieee_is_finite(local_matrix)) .or. &
            any(.not. ieee_is_finite(local_rhs)) .or. &
            any(.not. ieee_is_finite(global_matrix)) .or. &
            any(.not. ieee_is_finite(global_rhs))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_skeleton_inputs

end module fortfem_hdg_global_skeleton
