program test_complex_block_graph_residual
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_complex_block_graph_residual, &
        assemble_complex_block_graph_residual_jvp, &
        assemble_complex_block_graph_residual_vjp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: field_count = 3, edge_count = 4
    integer, parameter :: state_count = 6, value_count = 20
    integer :: field_sizes(field_count) = [2, 3, 1]
    integer :: edge_rows(edge_count) = [1, 2, 3, 1]
    integer :: edge_columns(edge_count) = [1, 2, 3, 2]
    integer :: block_offsets(edge_count + 1) = [1, 5, 14, 15, 21]
    complex(dp) :: block_values(value_count), block_values_dot(value_count)
    complex(dp) :: state(state_count), rhs(state_count)
    complex(dp) :: state_dot(state_count), rhs_dot(state_count)
    complex(dp) :: residual(state_count), residual_dot(state_count)
    complex(dp) :: residual_plus(state_count), residual_bar(state_count)
    complex(dp) :: block_values_bar(value_count), state_bar(state_count)
    complex(dp) :: rhs_bar(state_count), oracle(state_count), oracle_dot(state_count)
    complex(dp) :: matrix(state_count, state_count)
    real(dp) :: epsilon, fd_error, lhs, rhs_inner
    type(fortsparse_status_t) :: status
    logical :: all_passed
    integer :: value_index

    all_passed = .true.
    block_values = [(cmplx(0.02_dp*real(value_index, dp), &
        -0.01_dp*real(value_index, dp), dp), value_index=1,value_count)]
    block_values_dot = [(cmplx(0.01_dp*real(value_index, dp), &
        0.003_dp*real(value_index, dp), dp), value_index=1,value_count)]
    state = [cmplx(0.4_dp, -0.1_dp, dp), cmplx(-0.2_dp, 0.3_dp, dp), &
        cmplx(0.1_dp, 0.2_dp, dp), cmplx(0.5_dp, -0.2_dp, dp), &
        cmplx(-0.3_dp, 0.1_dp, dp), cmplx(0.7_dp, -0.4_dp, dp)]
    rhs = [cmplx(-0.1_dp, 0.2_dp, dp), cmplx(0.2_dp, -0.1_dp, dp), &
        cmplx(-0.3_dp, 0.3_dp, dp), cmplx(0.4_dp, -0.2_dp, dp), &
        cmplx(0.1_dp, 0.1_dp, dp), cmplx(-0.2_dp, 0.4_dp, dp)]
    state_dot = [cmplx(0.02_dp, -0.01_dp, dp), cmplx(-0.01_dp, 0.02_dp, dp), &
        cmplx(0.03_dp, 0.01_dp, dp), cmplx(0.04_dp, -0.02_dp, dp), &
        cmplx(-0.02_dp, 0.03_dp, dp), cmplx(0.01_dp, -0.03_dp, dp)]
    rhs_dot = [cmplx(-0.03_dp, 0.01_dp, dp), cmplx(0.02_dp, -0.03_dp, dp), &
        cmplx(-0.01_dp, 0.02_dp, dp), cmplx(0.04_dp, -0.01_dp, dp), &
        cmplx(0.01_dp, 0.03_dp, dp), cmplx(-0.02_dp, 0.01_dp, dp)]
    call build_matrix(field_sizes, edge_rows, edge_columns, block_offsets, block_values, matrix)
    oracle = matmul(matrix, state) - rhs
    call assemble_complex_block_graph_residual(field_sizes, edge_rows, edge_columns, &
        block_offsets, block_values, state, rhs, residual, status)
    call record_condition(status%code == 0, "complex packed block graph residual assembles")
    call record_condition(maxval(abs(residual - oracle)) < 2.0e-14_dp, &
        "complex packed block graph residual matches the independent oracle")

    call build_matrix(field_sizes, edge_rows, edge_columns, block_offsets, block_values_dot, matrix)
    oracle_dot = matmul(matrix, state)
    call build_matrix(field_sizes, edge_rows, edge_columns, block_offsets, block_values, matrix)
    oracle_dot = oracle_dot + matmul(matrix, state_dot) - rhs_dot
    call assemble_complex_block_graph_residual_jvp(field_sizes, edge_rows, edge_columns, &
        block_offsets, block_values, state, rhs, block_values_dot, state_dot, rhs_dot, &
        residual_dot, status)
    call record_condition(status%code == 0, "complex packed block graph residual JVP assembles")
    epsilon = 1.0e-7_dp
    call assemble_complex_block_graph_residual(field_sizes, edge_rows, edge_columns, &
        block_offsets, block_values + epsilon*block_values_dot, state + epsilon*state_dot, &
        rhs + epsilon*rhs_dot, residual_plus, status)
    fd_error = maxval(abs(residual_dot - (residual_plus - residual)/epsilon))
    call record_condition(fd_error < 4.0e-8_dp .and. maxval(abs(residual_dot - oracle_dot)) < 2.0e-14_dp, &
        "complex packed block graph residual JVP matches both oracles")

    residual_bar = [cmplx(0.2_dp, -0.1_dp, dp), cmplx(-0.3_dp, 0.2_dp, dp), &
        cmplx(0.1_dp, 0.3_dp, dp), cmplx(0.4_dp, -0.2_dp, dp), &
        cmplx(-0.2_dp, 0.1_dp, dp), cmplx(0.5_dp, -0.3_dp, dp)]
    call assemble_complex_block_graph_residual_vjp(field_sizes, edge_rows, edge_columns, &
        block_offsets, block_values, state, rhs, residual_bar, block_values_bar, state_bar, &
        rhs_bar, status)
    call record_condition(status%code == 0, "complex packed block graph residual VJP assembles")
    lhs = real(sum(conjg(residual_bar)*residual_dot), dp)
    rhs_inner = real(sum(conjg(block_values_bar)*block_values_dot) + &
        sum(conjg(state_bar)*state_dot) + sum(conjg(rhs_bar)*rhs_dot), dp)
    call record_condition(abs(lhs - rhs_inner) < 4.0e-13_dp, &
        "complex packed block graph residual VJP satisfies the real-part adjoint")

    call assemble_complex_block_graph_residual(field_sizes, edge_rows, edge_columns, &
        [1, 5, 14, 15, 20], block_values, state, rhs, residual, status)
    call record_condition(status%code /= 0, "inconsistent complex block offsets are rejected")

    call check_summary("complex packed block graph residual")
    if (.not. all_passed) error stop 1

contains

    subroutine build_matrix(field_sizes, edge_rows, edge_columns, block_offsets, values, matrix)
        integer, intent(in) :: field_sizes(:), edge_rows(:), edge_columns(:), block_offsets(:)
        complex(dp), intent(in) :: values(:)
        complex(dp), intent(out) :: matrix(:, :)
        integer :: edge, row, column, row_offset, column_offset, value_index

        matrix = cmplx(0.0_dp, 0.0_dp, dp)
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

end program test_complex_block_graph_residual
