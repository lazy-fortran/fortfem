program test_multipatch_signed_trace_assembly_csc
    use check, only: check_condition, check_summary
    use fortfem_multipatch_signed_trace_assembly, only: &
        assemble_multipatch_signed_trace_assembly_csc, &
        assemble_multipatch_signed_trace_assembly_csc_jvp, &
        assemble_multipatch_signed_trace_assembly_csc_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: local_count = 2, patch_count = 2, global_count = 3
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp) :: local_matrix(local_count, local_count, patch_count)
    real(dp) :: local_rhs(local_count, patch_count)
    integer :: local_to_global(local_count, patch_count), local_sign(local_count, patch_count)
    real(dp) :: local_matrix_dot(local_count, local_count, patch_count)
    real(dp) :: local_rhs_dot(local_count, patch_count)
    type(csc_t) :: global_matrix, global_matrix_dot, global_matrix_plus
    type(csc_t) :: global_matrix_minus, global_matrix_bar
    real(dp) :: global_rhs(global_count), global_rhs_dot(global_count)
    real(dp) :: global_rhs_plus(global_count), global_rhs_minus(global_count)
    real(dp) :: global_rhs_bar(global_count)
    real(dp) :: local_matrix_bar(local_count, local_count, patch_count)
    real(dp) :: local_rhs_bar(local_count, patch_count)
    real(dp) :: expected_matrix(global_count, global_count), expected_rhs(global_count)
    real(dp) :: matrix_dense(global_count, global_count), matrix_dot_dense(global_count, global_count)
    real(dp) :: matrix_plus_dense(global_count, global_count), matrix_minus_dense(global_count, global_count)
    real(dp) :: expected_dot(global_count, global_count), expected_rhs_dot(global_count)
    real(dp) :: lhs, rhs
    integer, parameter :: bar_rows(7) = [1, 2, 1, 2, 2, 3, 3]
    integer, parameter :: bar_columns(7) = [1, 1, 2, 2, 3, 2, 3]
    real(dp), parameter :: bar_values(7) = [0.4_dp, -0.2_dp, 0.1_dp, 0.3_dp, -0.5_dp, 0.7_dp, -0.1_dp]
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    local_matrix = reshape([2.0_dp, 0.3_dp, 0.3_dp, 1.4_dp, &
        1.7_dp, -0.2_dp, -0.2_dp, 2.1_dp], shape(local_matrix))
    local_rhs = reshape([0.4_dp, -0.7_dp, 1.1_dp, 0.2_dp], shape(local_rhs))
    local_to_global = reshape([1, 2, 2, 3], shape(local_to_global))
    local_sign = reshape([1, -1, -1, 1], shape(local_sign))
    local_matrix_dot = reshape([0.1_dp, -0.2_dp, 0.05_dp, 0.3_dp, &
        -0.15_dp, 0.08_dp, 0.12_dp, -0.04_dp], shape(local_matrix_dot))
    local_rhs_dot = reshape([-0.1_dp, 0.2_dp, 0.03_dp, 0.15_dp], shape(local_rhs_dot))
    global_rhs_bar = [0.2_dp, -0.4_dp, 0.5_dp]

    call assemble_multipatch_signed_trace_assembly_csc( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        global_matrix, global_rhs, status)
    call oracle(expected_matrix, expected_rhs, local_matrix, local_rhs)
    call csc_to_dense(global_matrix, matrix_dense)
    call record_condition(status%code == 0 .and. &
        maxval(abs(matrix_dense - expected_matrix)) < 1.0e-14_dp .and. &
        maxval(abs(global_rhs - expected_rhs)) < 1.0e-14_dp, &
        "multipatch CSC assembly matches an independent nested-loop oracle")

    call assemble_multipatch_signed_trace_assembly_csc_jvp( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        local_matrix_dot, local_rhs_dot, global_matrix_dot, global_rhs_dot, status)
    call oracle(expected_dot, expected_rhs_dot, local_matrix_dot, local_rhs_dot)
    call csc_to_dense(global_matrix_dot, matrix_dot_dense)
    call record_condition(status%code == 0 .and. &
        maxval(abs(matrix_dot_dense - expected_dot)) < 1.0e-14_dp .and. &
        maxval(abs(global_rhs_dot - expected_rhs_dot)) < 1.0e-14_dp, &
        "multipatch CSC JVP matches an independent signed-map oracle")
    call assemble_multipatch_signed_trace_assembly_csc( &
        local_matrix + step*local_matrix_dot, local_rhs + step*local_rhs_dot, &
        local_to_global, local_sign, global_count, global_matrix_plus, global_rhs_plus, status)
    call assemble_multipatch_signed_trace_assembly_csc( &
        local_matrix - step*local_matrix_dot, local_rhs - step*local_rhs_dot, &
        local_to_global, local_sign, global_count, global_matrix_minus, global_rhs_minus, status)
    call csc_to_dense(global_matrix_plus, matrix_plus_dense)
    call csc_to_dense(global_matrix_minus, matrix_minus_dense)
    call record_condition(status%code == 0 .and. &
        maxval(abs(matrix_dot_dense - (matrix_plus_dense - matrix_minus_dense)/(2.0_dp*step))) < 2.0e-8_dp .and. &
        maxval(abs(global_rhs_dot - (global_rhs_plus - global_rhs_minus)/(2.0_dp*step))) < 2.0e-8_dp, &
        "multipatch CSC JVP matches an independent central difference")

    call csc_from_triplet(global_count, global_count, bar_rows, bar_columns, &
        bar_values, global_matrix_bar, status)
    call assemble_multipatch_signed_trace_assembly_csc_vjp( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        global_matrix_bar, global_rhs_bar, local_matrix_bar, local_rhs_bar, status)
    lhs = sparse_dot(global_matrix_bar, global_matrix_dot) + dot_product(global_rhs_bar, global_rhs_dot)
    rhs = sum(local_matrix_bar*local_matrix_dot) + sum(local_rhs_bar*local_rhs_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "multipatch CSC VJP satisfies the real dot-product identity")

    local_to_global(1, 1) = local_to_global(2, 1)
    call assemble_multipatch_signed_trace_assembly_csc( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        global_matrix, global_rhs, status)
    call record_condition(status%code /= 0, "multipatch CSC rejects duplicate IDs within a patch")
    local_to_global = reshape([1, 2, 2, 3], shape(local_to_global))
    local_sign(1, 1) = 0
    call assemble_multipatch_signed_trace_assembly_csc( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        global_matrix, global_rhs, status)
    call record_condition(status%code /= 0, "multipatch CSC rejects an invalid orientation sign")

    call check_summary("multipatch signed trace CSC assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine oracle(matrix, vector, blocks, loads)
        real(dp), intent(out) :: matrix(:, :), vector(:)
        real(dp), intent(in) :: blocks(:, :, :), loads(:, :)
        integer :: patch, local_row, local_column, row, column

        matrix = 0.0_dp
        vector = 0.0_dp
        do patch = 1, patch_count
            do local_row = 1, local_count
                row = local_to_global(local_row, patch)
                vector(row) = vector(row) + local_sign(local_row, patch)*loads(local_row, patch)
                do local_column = 1, local_count
                    column = local_to_global(local_column, patch)
                    matrix(row, column) = matrix(row, column) + &
                        local_sign(local_row, patch)*local_sign(local_column, patch)* &
                        blocks(local_row, local_column, patch)
                end do
            end do
        end do
    end subroutine oracle

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

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_multipatch_signed_trace_assembly_csc
