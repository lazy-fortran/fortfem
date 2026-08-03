program test_block_graph_residual
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_block_graph_residual, &
        assemble_block_graph_residual_jvp, &
        assemble_block_graph_residual_vjp
    use fortfem_generated_block_graph, only: generated_block_graph_product
    use fortfem_generated_block_graph_product_jvp, only: &
        generated_block_graph_product_jvp
    use fortfem_generated_block_graph_product_vjp, only: &
        generated_block_graph_product_vjp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: field_count = 3, edge_count = 5
    integer, parameter :: total_state_count = 6, total_value_count = 26
    integer :: field_sizes(field_count) = [2, 3, 1]
    integer :: edge_rows(edge_count) = [1, 2, 3, 1, 2]
    integer :: edge_columns(edge_count) = [1, 2, 3, 2, 1]
    integer :: block_offsets(edge_count + 1) = [1, 5, 14, 15, 21, 27]
    real(dp) :: block_values(total_value_count), block_values_dot(total_value_count)
    real(dp) :: state(total_state_count), rhs(total_state_count)
    real(dp) :: state_dot(total_state_count), rhs_dot(total_state_count)
    real(dp) :: residual(total_state_count), residual_dot(total_state_count)
    real(dp) :: residual_plus(total_state_count), residual_bar(total_state_count)
    real(dp) :: block_values_bar(total_value_count), state_bar(total_state_count)
    real(dp) :: rhs_bar(total_state_count), oracle(total_state_count)
    real(dp) :: oracle_dot(total_state_count), epsilon, fd_error, lhs, rhs_inner
    real(dp) :: matrix(total_state_count, total_state_count)
    real(dp) :: product, product_dot, product_bar, state_value_bar
    type(fortsparse_status_t) :: status
    logical :: all_passed
    integer :: edge, row, column, row_offset, column_offset, value_index

    all_passed = .true.
    call generated_block_graph_product(2.5_dp, -0.4_dp, product)
    call record_condition(abs(product + 1.0_dp) < 1.0e-14_dp, &
        "FortSym block graph product reproduces the scalar action")
    call generated_block_graph_product_jvp(2.5_dp, -0.4_dp, 0.3_dp, 0.2_dp, product_dot)
    call record_condition(abs(product_dot - 0.38_dp) < 1.0e-14_dp, &
        "FortSym block graph product reproduces the product-rule JVP")
    call generated_block_graph_product_vjp(2.5_dp, -0.4_dp, 0.7_dp, product_bar, state_value_bar)
    call record_condition(abs(product_bar + 0.28_dp) < 1.0e-14_dp .and. &
        abs(state_value_bar - 1.75_dp) < 1.0e-14_dp, &
        "FortSym block graph product reproduces the scalar VJP")
    block_values = [(0.03_dp*real(value_index, dp), value_index=1,total_value_count)]
    block_values_dot = [ &
        0.01_dp, -0.02_dp, 0.03_dp, -0.01_dp, 0.02_dp, 0.04_dp, -0.03_dp, 0.01_dp, &
        0.02_dp, -0.04_dp, 0.01_dp, 0.03_dp, 0.02_dp, -0.01_dp, 0.05_dp, -0.02_dp, &
        0.01_dp, 0.04_dp, -0.03_dp, 0.02_dp, 0.03_dp, -0.02_dp, 0.01_dp, &
        0.02_dp, -0.01_dp, 0.03_dp]
    state = [0.4_dp, -0.2_dp, 0.1_dp, 0.5_dp, -0.3_dp, 0.7_dp]
    rhs = [-0.1_dp, 0.2_dp, -0.3_dp, 0.4_dp, 0.1_dp, -0.2_dp]
    state_dot = [0.02_dp, -0.01_dp, 0.03_dp, 0.04_dp, -0.02_dp, 0.01_dp]
    rhs_dot = [-0.03_dp, 0.02_dp, -0.01_dp, 0.04_dp, 0.01_dp, -0.02_dp]
    call build_matrix(field_sizes, edge_rows, edge_columns, block_offsets, block_values, matrix)
    oracle = matmul(matrix, state) - rhs
    call assemble_block_graph_residual(field_sizes, edge_rows, edge_columns, block_offsets, &
        block_values, state, rhs, residual, status)
    call record_condition(status%code == 0, "packed block graph residual assembles")
    call record_condition(maxval(abs(residual - oracle)) < 1.0e-14_dp, &
        "packed block graph residual matches the independent matrix oracle")

    call build_matrix(field_sizes, edge_rows, edge_columns, block_offsets, block_values_dot, &
        matrix)
    oracle_dot = matmul(matrix, state)
    call build_matrix(field_sizes, edge_rows, edge_columns, block_offsets, block_values, matrix)
    oracle_dot = oracle_dot + matmul(matrix, state_dot) - rhs_dot
    call assemble_block_graph_residual_jvp(field_sizes, edge_rows, edge_columns, block_offsets, &
        block_values, state, rhs, block_values_dot, state_dot, rhs_dot, residual_dot, status)
    call record_condition(status%code == 0, "packed block graph residual JVP assembles")
    epsilon = 1.0e-7_dp
    call assemble_block_graph_residual(field_sizes, edge_rows, edge_columns, block_offsets, &
        block_values + epsilon*block_values_dot, state + epsilon*state_dot, &
        rhs + epsilon*rhs_dot, residual_plus, status)
    fd_error = maxval(abs(residual_dot - (residual_plus - residual)/epsilon))
    call record_condition(fd_error < 3.0e-8_dp .and. &
        maxval(abs(residual_dot - oracle_dot)) < 1.0e-14_dp, &
        "packed block graph residual JVP matches both oracles")

    residual_bar = [0.2_dp, -0.3_dp, 0.1_dp, 0.4_dp, -0.2_dp, 0.5_dp]
    call assemble_block_graph_residual_vjp(field_sizes, edge_rows, edge_columns, block_offsets, &
        block_values, state, rhs, residual_bar, block_values_bar, state_bar, rhs_bar, status)
    call record_condition(status%code == 0, "packed block graph residual VJP assembles")
    lhs = dot_product(residual_bar, residual_dot)
    rhs_inner = dot_product(block_values_bar, block_values_dot) + &
        dot_product(state_bar, state_dot) + dot_product(rhs_bar, rhs_dot)
    call record_condition(abs(lhs - rhs_inner) < 3.0e-13_dp, &
        "packed block graph residual VJP satisfies the adjoint identity")

    call assemble_block_graph_residual(field_sizes, edge_rows, edge_columns, &
        [1, 5, 14, 15, 21, 26], block_values, state, rhs, residual, status)
    call record_condition(status%code /= 0, "inconsistent block offsets are rejected")

    call check_summary("packed block graph residual")
    if (.not. all_passed) error stop 1

contains

    subroutine build_matrix(field_sizes, edge_rows, edge_columns, block_offsets, values, matrix)
        integer, intent(in) :: field_sizes(:), edge_rows(:), edge_columns(:)
        integer, intent(in) :: block_offsets(:)
        real(dp), intent(in) :: values(:)
        real(dp), intent(out) :: matrix(:, :)
        integer :: edge, row, column, row_offset, column_offset, value_index

        matrix = 0.0_dp
        do edge = 1, size(edge_rows)
            row_offset = field_offset(field_sizes, edge_rows(edge))
            column_offset = field_offset(field_sizes, edge_columns(edge))
            do column = 1, field_sizes(edge_columns(edge))
                do row = 1, field_sizes(edge_rows(edge))
                    value_index = block_offsets(edge) + row - 1 + &
                        (column - 1)*field_sizes(edge_rows(edge))
                    matrix(row_offset + row - 1, column_offset + column - 1) = &
                        matrix(row_offset + row - 1, column_offset + column - 1) + values(value_index)
                end do
            end do
        end do
    end subroutine build_matrix

    integer function field_offset(field_sizes, field_index) result(offset)
        integer, intent(in) :: field_sizes(:), field_index
        integer :: index

        offset = 1
        do index = 1, field_index - 1
            offset = offset + field_sizes(index)
        end do
    end function field_offset

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_block_graph_residual
