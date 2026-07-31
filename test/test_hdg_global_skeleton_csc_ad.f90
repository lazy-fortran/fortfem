program test_hdg_global_skeleton_csc_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_hdg_global_skeleton_csc, &
        assemble_hdg_global_skeleton_csc_jvp, assemble_hdg_global_skeleton_csc_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: local_count = 2, element_count = 2, global_count = 3
    real(dp), parameter :: local_matrix(local_count, local_count, element_count) = &
        reshape([2.0_dp, 5.0_dp, 3.0_dp, 7.0_dp, &
        11.0_dp, 17.0_dp, 13.0_dp, 19.0_dp], &
        [local_count, local_count, element_count])
    real(dp), parameter :: local_rhs(local_count, element_count) = &
        reshape([23.0_dp, 29.0_dp, 31.0_dp, 37.0_dp], [local_count, element_count])
    integer, parameter :: local_to_global(local_count, element_count) = &
        reshape([1, 2, 2, 3], [local_count, element_count])
    integer, parameter :: local_sign(local_count, element_count) = &
        reshape([1, -1, -1, 1], [local_count, element_count])
    real(dp), parameter :: local_matrix_dot(local_count, local_count, element_count) = &
        reshape([0.2_dp, -0.4_dp, 0.6_dp, -0.8_dp, &
        1.0_dp, -1.2_dp, 1.4_dp, -1.6_dp], &
        [local_count, local_count, element_count])
    real(dp), parameter :: local_rhs_dot(local_count, element_count) = &
        reshape([0.3_dp, -0.5_dp, 0.7_dp, -0.9_dp], [local_count, element_count])
    integer, parameter :: bar_rows(7) = [1, 2, 1, 2, 3, 2, 3]
    integer, parameter :: bar_columns(7) = [1, 1, 2, 2, 2, 3, 3]
    real(dp), parameter :: bar_values(7) = [0.4_dp, -0.2_dp, 0.1_dp, &
        0.3_dp, -0.5_dp, 0.7_dp, -0.1_dp]
    real(dp), parameter :: global_rhs_bar(global_count) = [0.2_dp, -0.4_dp, 0.5_dp]
    type(csc_t) :: global_matrix, global_matrix_dot, global_matrix_bar
    type(fortsparse_status_t) :: status
    real(dp) :: global_rhs(global_count), global_rhs_dot(global_count)
    real(dp) :: local_matrix_bar(local_count, local_count, element_count)
    real(dp) :: local_rhs_bar(local_count, element_count)
    real(dp) :: global_matrix_dense(global_count, global_count)
    real(dp) :: global_matrix_dot_dense(global_count, global_count)
    real(dp) :: global_matrix_expected(global_count, global_count)
    real(dp) :: global_matrix_dot_expected(global_count, global_count)
    real(dp) :: global_rhs_expected(global_count), global_rhs_dot_expected(global_count)
    real(dp) :: lhs, rhs

    call assemble_hdg_global_skeleton_csc( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        global_matrix, global_rhs, status)
    call oracle_assembly(global_matrix_expected, global_rhs_expected, &
        local_matrix, local_rhs)
    call csc_to_dense(global_matrix, global_matrix_dense)
    call check_condition(status%code == 0 .and. &
        maxval(abs(global_matrix_dense - global_matrix_expected)) < 1.0e-14_dp .and. &
        maxval(abs(global_rhs - global_rhs_expected)) < 1.0e-14_dp, &
        "sparse HDG skeleton assembly matches the independent triplet oracle")

    call assemble_hdg_global_skeleton_csc_jvp( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        local_matrix_dot, local_rhs_dot, global_matrix_dot, global_rhs_dot, status)
    call oracle_assembly(global_matrix_dot_expected, global_rhs_dot_expected, &
        local_matrix_dot, local_rhs_dot)
    call csc_to_dense(global_matrix_dot, global_matrix_dot_dense)
    call check_condition(status%code == 0 .and. &
        maxval(abs(global_matrix_dot_dense - global_matrix_dot_expected)) < 1.0e-14_dp .and. &
        maxval(abs(global_rhs_dot - global_rhs_dot_expected)) < 1.0e-14_dp, &
        "sparse HDG skeleton JVP matches an independent signed-map oracle")

    call csc_from_triplet(global_count, global_count, bar_rows, bar_columns, &
        bar_values, global_matrix_bar, status)
    call assemble_hdg_global_skeleton_csc_vjp( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        global_matrix_bar, global_rhs_bar, local_matrix_bar, local_rhs_bar, status)
    lhs = sparse_dot(global_matrix_bar, global_matrix_dot) + &
        dot_product(global_rhs_bar, global_rhs_dot)
    rhs = sum(local_matrix_bar*local_matrix_dot) + &
        sum(local_rhs_bar*local_rhs_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "sparse HDG skeleton VJP satisfies the real dot-product identity")

    call assemble_hdg_global_skeleton_csc( &
        local_matrix, local_rhs, reshape([1, 4, 2, 3], [local_count, element_count]), &
        local_sign, global_count, global_matrix, global_rhs, status)
    call check_condition(status%code /= 0, &
        "sparse HDG skeleton assembly rejects an out-of-range trace index")
    call check_summary("HDG global skeleton CSC AD")

contains

    subroutine oracle_assembly(matrix, vector, blocks, loads)
        real(dp), intent(out) :: matrix(:, :), vector(:)
        real(dp), intent(in) :: blocks(:, :, :), loads(:, :)
        integer :: element, local_row, local_column, row, column

        matrix = 0.0_dp
        vector = 0.0_dp
        do element = 1, element_count
            do local_row = 1, local_count
                row = local_to_global(local_row, element)
                vector(row) = vector(row) + local_sign(local_row, element)* &
                    loads(local_row, element)
                do local_column = 1, local_count
                    column = local_to_global(local_column, element)
                    matrix(row, column) = matrix(row, column) + &
                        local_sign(local_row, element)*local_sign(local_column, element)* &
                        blocks(local_row, local_column, element)
                end do
            end do
        end do
    end subroutine oracle_assembly

    subroutine csc_to_dense(matrix, dense)
        type(csc_t), intent(in) :: matrix
        real(dp), intent(out) :: dense(:, :)
        integer :: column, entry

        dense = 0.0_dp
        do column = 1, matrix%ncol
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                dense(matrix%row_idx(entry), column) = matrix%val(entry)
            end do
        end do
    end subroutine csc_to_dense

    pure real(dp) function sparse_dot(left, right) result(value)
        type(csc_t), intent(in) :: left, right
        integer :: column, entry, position

        value = 0.0_dp
        do column = 1, left%ncol
            do entry = left%col_ptr(column), left%col_ptr(column + 1) - 1
                do position = right%col_ptr(column), right%col_ptr(column + 1) - 1
                    if (right%row_idx(position) == left%row_idx(entry)) then
                        value = value + left%val(entry)*right%val(position)
                        exit
                    end if
                end do
            end do
        end do
    end function sparse_dot

end program test_hdg_global_skeleton_csc_ad
