program test_multipatch_signed_trace_assembly
    use check, only: check_condition, check_summary
    use fortfem_multipatch_signed_trace_assembly, only: &
        assemble_multipatch_signed_trace_assembly, &
        assemble_multipatch_signed_trace_assembly_jvp, &
        assemble_multipatch_signed_trace_assembly_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: local_count = 2, patch_count = 2, global_count = 3
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp) :: local_matrix(local_count, local_count, patch_count)
    real(dp) :: local_rhs(local_count, patch_count)
    integer :: local_to_global(local_count, patch_count)
    integer :: local_sign(local_count, patch_count)
    real(dp) :: local_matrix_dot(local_count, local_count, patch_count)
    real(dp) :: local_rhs_dot(local_count, patch_count)
    real(dp) :: global_matrix(global_count, global_count), global_rhs(global_count)
    real(dp) :: global_matrix_dot(global_count, global_count), global_rhs_dot(global_count)
    real(dp) :: global_matrix_plus(global_count, global_count), global_rhs_plus(global_count)
    real(dp) :: global_matrix_minus(global_count, global_count), global_rhs_minus(global_count)
    real(dp) :: global_matrix_bar(global_count, global_count), global_rhs_bar(global_count)
    real(dp) :: local_matrix_bar(local_count, local_count, patch_count)
    real(dp) :: local_rhs_bar(local_count, patch_count)
    real(dp) :: expected_matrix(global_count, global_count), expected_rhs(global_count)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status
    logical :: all_passed
    integer :: patch, local_row, local_column, row, column

    all_passed = .true.
    local_matrix = reshape([2.0_dp, 0.3_dp, 0.3_dp, 1.4_dp, &
        1.7_dp, -0.2_dp, -0.2_dp, 2.1_dp], shape(local_matrix))
    local_rhs = reshape([0.4_dp, -0.7_dp, 1.1_dp, 0.2_dp], shape(local_rhs))
    local_to_global = reshape([1, 2, 2, 3], shape(local_to_global))
    local_sign = reshape([1, 1, -1, 1], shape(local_sign))
    local_matrix_dot = reshape([0.1_dp, -0.2_dp, 0.05_dp, 0.3_dp, &
        -0.15_dp, 0.08_dp, 0.12_dp, -0.04_dp], shape(local_matrix_dot))
    local_rhs_dot = reshape([-0.1_dp, 0.2_dp, 0.03_dp, 0.15_dp], shape(local_rhs_dot))
    global_matrix_bar = reshape([0.4_dp, -0.2_dp, 0.1_dp, 0.3_dp, -0.5_dp, 0.7_dp, &
        -0.1_dp, 0.6_dp, 0.2_dp], shape(global_matrix_bar))
    global_rhs_bar = [0.2_dp, -0.4_dp, 0.5_dp]

    call assemble_multipatch_signed_trace_assembly( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        global_matrix, global_rhs, status)
    call oracle(expected_matrix, expected_rhs, local_matrix, local_rhs)
    call record_condition(status%code == 0 .and. &
        maxval(abs(global_matrix - expected_matrix)) < 1.0e-14_dp .and. &
        maxval(abs(global_rhs - expected_rhs)) < 1.0e-14_dp, &
        "signed multipatch assembly matches an independent nested-loop oracle")

    call assemble_multipatch_signed_trace_assembly_jvp( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        local_matrix_dot, local_rhs_dot, global_matrix_dot, global_rhs_dot, status)
    call assemble_multipatch_signed_trace_assembly( &
        local_matrix + step*local_matrix_dot, local_rhs + step*local_rhs_dot, &
        local_to_global, local_sign, global_count, global_matrix_plus, global_rhs_plus, status)
    call assemble_multipatch_signed_trace_assembly( &
        local_matrix - step*local_matrix_dot, local_rhs - step*local_rhs_dot, &
        local_to_global, local_sign, global_count, global_matrix_minus, global_rhs_minus, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(global_matrix_dot - (global_matrix_plus - global_matrix_minus)/(2.0_dp*step))) < 2.0e-8_dp .and. &
        maxval(abs(global_rhs_dot - (global_rhs_plus - global_rhs_minus)/(2.0_dp*step))) < 2.0e-8_dp, &
        "signed multipatch JVP matches an independent central difference")

    call assemble_multipatch_signed_trace_assembly_vjp( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        global_matrix_bar, global_rhs_bar, local_matrix_bar, local_rhs_bar, status)
    lhs = sum(global_matrix_bar*global_matrix_dot) + dot_product(global_rhs_bar, global_rhs_dot)
    rhs = sum(local_matrix_bar*local_matrix_dot) + sum(local_rhs_bar*local_rhs_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "signed multipatch VJP satisfies the real dot-product identity")

    local_to_global(1, 1) = local_to_global(2, 1)
    call assemble_multipatch_signed_trace_assembly( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        global_matrix, global_rhs, status)
    call record_condition(status%code /= 0, &
        "signed multipatch assembly rejects duplicate IDs within a patch")
    local_to_global = reshape([1, 2, 2, 3], shape(local_to_global))
    local_sign(1, 1) = 0
    call assemble_multipatch_signed_trace_assembly( &
        local_matrix, local_rhs, local_to_global, local_sign, global_count, &
        global_matrix, global_rhs, status)
    call record_condition(status%code /= 0, &
        "signed multipatch assembly rejects an invalid orientation sign")

    call check_summary("multipatch signed trace assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine oracle(matrix, vector, local_matrix, local_rhs)
        real(dp), intent(out) :: matrix(:, :), vector(:)
        real(dp), intent(in) :: local_matrix(:, :, :), local_rhs(:, :)
        integer :: element, local_row, local_column, row, column

        matrix = 0.0_dp
        vector = 0.0_dp
        do element = 1, patch_count
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
    end subroutine oracle

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_multipatch_signed_trace_assembly
