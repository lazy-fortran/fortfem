program test_hdg_global_skeleton_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_hdg_global_skeleton, &
        assemble_hdg_global_skeleton_jvp, assemble_hdg_global_skeleton_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: local_count = 2, element_count = 2, global_count = 3
    real(dp), parameter :: local_matrix(local_count, local_count, element_count) = &
        reshape([2.0_dp, 0.3_dp, 0.3_dp, 1.4_dp, &
        1.7_dp, -0.2_dp, -0.2_dp, 2.1_dp], [local_count, local_count, element_count])
    real(dp), parameter :: local_rhs(local_count, element_count) = &
        reshape([0.4_dp, -0.7_dp, 1.1_dp, 0.2_dp], [local_count, element_count])
    integer, parameter :: local_to_global(local_count, element_count) = &
        reshape([1, 2, 2, 3], [local_count, element_count])
    integer, parameter :: local_sign(local_count, element_count) = &
        reshape([1, 1, -1, 1], [local_count, element_count])
    real(dp), parameter :: local_matrix_dot(local_count, local_count, element_count) = &
        reshape([0.1_dp, -0.2_dp, 0.05_dp, 0.3_dp, &
        -0.15_dp, 0.08_dp, 0.12_dp, -0.04_dp], &
        [local_count, local_count, element_count])
    real(dp), parameter :: local_rhs_dot(local_count, element_count) = &
        reshape([-0.1_dp, 0.2_dp, 0.03_dp, 0.15_dp], [local_count, element_count])
    real(dp), parameter :: global_matrix_bar(global_count, global_count) = reshape([ &
        0.4_dp, -0.2_dp, 0.1_dp, 0.3_dp, -0.5_dp, 0.7_dp, &
        -0.1_dp, 0.6_dp, 0.2_dp], [global_count, global_count])
    real(dp), parameter :: global_rhs_bar(global_count) = [0.2_dp, -0.4_dp, 0.5_dp]
    real(dp) :: global_matrix(global_count, global_count), global_rhs(global_count)
    real(dp) :: global_matrix_dot(global_count, global_count), global_rhs_dot(global_count)
    real(dp) :: local_matrix_bar(local_count, local_count, element_count)
    real(dp) :: local_rhs_bar(local_count, element_count)
    real(dp) :: global_matrix_expected(global_count, global_count)
    real(dp) :: global_rhs_expected(global_count)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_hdg_global_skeleton( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        global_matrix, global_rhs, status)
    call oracle_assembly(global_matrix_expected, global_rhs_expected)
    call check_condition(status%code == 0 .and. &
        maxval(abs(global_matrix - global_matrix_expected)) < 1.0e-14_dp .and. &
        maxval(abs(global_rhs - global_rhs_expected)) < 1.0e-14_dp, &
        "HDG skeleton assembly matches the independent signed-map oracle")

    call assemble_hdg_global_skeleton_jvp( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        local_matrix_dot, local_rhs_dot, global_matrix_dot, global_rhs_dot, status)
    call check_condition(status%code == 0, &
        "HDG skeleton assembly JVP accepts local block directions")

    call assemble_hdg_global_skeleton_vjp( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        global_matrix_bar, global_rhs_bar, local_matrix_bar, local_rhs_bar, status)
    lhs = sum(global_matrix_bar*global_matrix_dot) + &
        dot_product(global_rhs_bar, global_rhs_dot)
    rhs = sum(local_matrix_bar*local_matrix_dot) + sum(local_rhs_bar*local_rhs_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "HDG skeleton assembly VJP satisfies the real dot-product identity")

    call assemble_hdg_global_skeleton( &
        local_matrix, local_rhs, reshape([1, 4, 2, 3], [local_count, element_count]), &
        local_sign, global_count, global_matrix, global_rhs, status)
    call check_condition(status%code /= 0, &
        "HDG skeleton assembly rejects an out-of-range trace index")
    call check_summary("HDG global skeleton AD")

contains

    subroutine oracle_assembly(matrix, vector)
        real(dp), intent(out) :: matrix(:, :), vector(:)
        integer :: element, local_row, local_column, row, column

        matrix = 0.0_dp
        vector = 0.0_dp
        do element = 1, element_count
            do local_row = 1, local_count
                row = local_to_global(local_row, element)
                vector(row) = vector(row) + local_sign(local_row, element)* &
                    local_rhs(local_row, element)
                do local_column = 1, local_count
                    column = local_to_global(local_column, element)
                    matrix(row, column) = matrix(row, column) + &
                        local_sign(local_row, element)*local_sign(local_column, element)* &
                        local_matrix(local_row, local_column, element)
                end do
            end do
        end do
    end subroutine oracle_assembly

end program test_hdg_global_skeleton_ad
