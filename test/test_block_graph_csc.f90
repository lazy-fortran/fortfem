program test_block_graph_csc
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_block_graph_csc, sparse_direct_factor_adjoint_csc, &
        sparse_direct_factor_csc, sparse_direct_factor_t, sparse_direct_factor_transpose_csc, &
        sparse_direct_free, sparse_direct_solve_factored, sparse_direct_solve_factored_jvp, &
        sparse_direct_solve_factored_vjp
    use fortsparse, only: csc_is_valid, csc_matvec, csc_t, csc_z_t, fortsparse_status_t
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: field_sizes(2) = [2, 1]
    integer, parameter :: edge_rows(4) = [1, 2, 1, 1]
    integer, parameter :: edge_columns(4) = [1, 2, 2, 2]
    integer, parameter :: block_offsets(5) = [1, 5, 6, 8, 10]
    real(dp), parameter :: block_values(9) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
        2.0_dp, 0.5_dp, -0.25_dp, -0.1_dp, 0.3_dp]
    real(dp), parameter :: state(3) = [0.5_dp, -0.2_dp, 0.7_dp]
    real(dp), parameter :: complex_state(3) = [0.5_dp, -0.2_dp, 0.7_dp]
    real(dp), parameter :: block_values_dot(9) = [0.02_dp, -0.01_dp, 0.03_dp, 0.04_dp, &
        -0.02_dp, 0.01_dp, -0.03_dp, 0.04_dp, 0.02_dp]
    complex(dp), parameter :: complex_values_dot(9) = [ &
        cmplx(0.02_dp, -0.01_dp, dp), cmplx(-0.01_dp, 0.03_dp, dp), &
        cmplx(0.03_dp, 0.02_dp, dp), cmplx(0.04_dp, -0.02_dp, dp), &
        cmplx(-0.02_dp, 0.01_dp, dp), cmplx(0.01_dp, -0.03_dp, dp), &
        cmplx(-0.03_dp, 0.04_dp, dp), cmplx(0.04_dp, 0.01_dp, dp), &
        cmplx(0.02_dp, -0.04_dp, dp)]
    real(dp) :: expected(3), dense_matrix(3, 3), applied(3)
    complex(dp) :: complex_values(9), complex_applied(3), complex_expected(3)
    real(dp) :: real_rhs(3), real_rhs_dot(3), real_solution(3), real_solution_dot(3)
    real(dp) :: real_solution_bar(3), real_rhs_bar(3), real_values_bar(7)
    real(dp) :: real_lhs, real_rhs_inner
    complex(dp) :: complex_rhs(3), complex_rhs_dot(3), complex_solution(3)
    complex(dp) :: complex_solution_dot(3), complex_solution_bar(3), complex_rhs_bar(3)
    complex(dp) :: complex_values_bar(7)
    real(dp) :: complex_lhs, complex_rhs_inner
    type(csc_t) :: dot_matrix
    type(csc_z_t) :: complex_dot_matrix
    type(sparse_direct_factor_t) :: factor, transpose_factor
    type(sparse_direct_factor_t) :: complex_factor, complex_adjoint_factor
    integer :: direct_status
    type(csc_t) :: matrix
    type(csc_z_t) :: complex_matrix
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    dense_matrix = 0.0_dp
    dense_matrix(1:2, 1:2) = reshape(block_values(1:4), [2, 2])
    dense_matrix(3, 3) = block_values(5)
    dense_matrix(1, 3) = dense_matrix(1, 3) + block_values(6)
    dense_matrix(2, 3) = dense_matrix(2, 3) + block_values(7)
    dense_matrix(1, 3) = dense_matrix(1, 3) + block_values(8)
    dense_matrix(2, 3) = dense_matrix(2, 3) + block_values(9)
    expected = matmul(dense_matrix, state)
    call assemble_block_graph_csc( &
        field_sizes, edge_rows, edge_columns, block_offsets, block_values, matrix, status)
    call record_condition(status%code == 0 .and. csc_is_valid(matrix), &
        "real packed graph converts to a valid CSC matrix")
    applied = csc_matvec(matrix, state)
    call record_condition(matrix%nnz == 7 .and. maxval(abs(applied - expected)) < 1.0e-14_dp, &
        "real CSC graph matvec matches the independent dense oracle")

    complex_values = cmplx(block_values, [0.1_dp, -0.2_dp, 0.3_dp, -0.4_dp, &
        0.2_dp, 0.05_dp, -0.15_dp, 0.04_dp, -0.06_dp], dp)
    complex_expected(1:2) = [ &
        complex_values(1)*complex_state(1) + complex_values(3)*complex_state(2) + &
        (complex_values(6) + complex_values(8))*complex_state(3), &
        complex_values(2)*complex_state(1) + complex_values(4)*complex_state(2) + &
        (complex_values(7) + complex_values(9))*complex_state(3)]
    complex_expected(3) = complex_values(5)*complex_state(3)
    call assemble_block_graph_csc( &
        field_sizes, edge_rows, edge_columns, block_offsets, complex_values, complex_matrix, status)
    call record_condition(status%code == 0 .and. csc_is_valid(complex_matrix), &
        "complex packed graph converts to a valid CSC matrix")
    complex_applied = csc_matvec(complex_matrix, cmplx(complex_state, 0.0_dp, dp))
    call record_condition(complex_matrix%nnz == 7 .and. &
        maxval(abs(complex_applied - complex_expected)) < 1.0e-14_dp, &
        "complex CSC graph matvec matches the independent oracle")

    call assemble_block_graph_csc( &
        field_sizes, edge_rows, edge_columns, block_offsets, block_values_dot, dot_matrix, status)
    call sparse_direct_factor_csc( &
        factor, matrix%nrow, matrix%col_ptr, matrix%row_idx, matrix%val, direct_status)
    call record_condition(status%code == 0 .and. direct_status == 0, &
        "real packed graph CSC can be retained and factored")
    real_rhs = csc_matvec(matrix, state)
    real_rhs_dot = [0.03_dp, -0.02_dp, 0.01_dp]
    call sparse_direct_solve_factored(factor, real_rhs, real_solution, direct_status)
    call sparse_direct_solve_factored_jvp( &
        factor, matrix%nrow, matrix%col_ptr, matrix%row_idx, dot_matrix%val, &
        real_solution, real_rhs_dot, real_solution_dot, direct_status)
    call sparse_direct_factor_transpose_csc( &
        transpose_factor, matrix%nrow, matrix%col_ptr, matrix%row_idx, matrix%val, direct_status)
    real_solution_bar = [0.4_dp, -0.3_dp, 0.2_dp]
    call sparse_direct_solve_factored_vjp( &
        transpose_factor, matrix%nrow, matrix%col_ptr, matrix%row_idx, real_solution, &
        real_solution_bar, real_rhs_bar, real_values_bar, direct_status)
    real_lhs = dot_product(real_solution_bar, real_solution_dot)
    real_rhs_inner = dot_product(real_rhs_bar, real_rhs_dot) + &
        dot_product(real_values_bar, dot_matrix%val)
    call record_condition(direct_status == 0 .and. maxval(abs(real_solution - state)) < 1.0e-13_dp .and. &
        abs(real_lhs - real_rhs_inner) < 1.0e-12_dp, &
        "real retained graph factor exposes solve JVP/VJP composition")
    call sparse_direct_free(factor)
    call sparse_direct_free(transpose_factor)

    call assemble_block_graph_csc( &
        field_sizes, edge_rows, edge_columns, block_offsets, complex_values_dot, &
        complex_dot_matrix, status)
    call sparse_direct_factor_csc( &
        complex_factor, complex_matrix%nrow, complex_matrix%col_ptr, complex_matrix%row_idx, &
        complex_matrix%val, direct_status)
    call record_condition(status%code == 0 .and. direct_status == 0, &
        "complex packed graph CSC can be retained and factored")
    complex_rhs = csc_matvec(complex_matrix, cmplx(complex_state, 0.0_dp, dp))
    complex_rhs_dot = [cmplx(0.03_dp, -0.01_dp, dp), cmplx(-0.02_dp, 0.02_dp, dp), &
        cmplx(0.01_dp, 0.03_dp, dp)]
    call sparse_direct_solve_factored(complex_factor, complex_rhs, complex_solution, direct_status)
    call sparse_direct_solve_factored_jvp( &
        complex_factor, complex_matrix%nrow, complex_matrix%col_ptr, complex_matrix%row_idx, &
        complex_dot_matrix%val, complex_solution, complex_rhs_dot, complex_solution_dot, direct_status)
    call sparse_direct_factor_adjoint_csc( &
        complex_adjoint_factor, complex_matrix%nrow, complex_matrix%col_ptr, &
        complex_matrix%row_idx, complex_matrix%val, direct_status)
    complex_solution_bar = [cmplx(0.4_dp, 0.1_dp, dp), cmplx(-0.3_dp, 0.2_dp, dp), &
        cmplx(0.2_dp, -0.1_dp, dp)]
    call sparse_direct_solve_factored_vjp( &
        complex_adjoint_factor, complex_matrix%nrow, complex_matrix%col_ptr, &
        complex_matrix%row_idx, complex_solution, complex_solution_bar, complex_rhs_bar, &
        complex_values_bar, direct_status)
    complex_lhs = real(sum(conjg(complex_solution_bar)*complex_solution_dot), dp)
    complex_rhs_inner = real(sum(conjg(complex_rhs_bar)*complex_rhs_dot) + &
        sum(conjg(complex_values_bar)*complex_dot_matrix%val), dp)
    call record_condition(direct_status == 0 .and. maxval(abs(complex_solution - &
        cmplx(complex_state, 0.0_dp, dp))) < 1.0e-13_dp .and. &
        abs(complex_lhs - complex_rhs_inner) < 1.0e-12_dp, &
        "complex retained graph factor exposes solve JVP/VJP composition")
    call sparse_direct_free(complex_factor)
    call sparse_direct_free(complex_adjoint_factor)

    call assemble_block_graph_csc( &
        field_sizes, edge_rows, edge_columns, [1, 5, 6, 8, 9], block_values, matrix, status)
    call record_condition(status%code /= 0 .and. matrix%nrow == 0, &
        "inconsistent packed graph offsets are rejected")
    call check_summary("packed block graph CSC adapter")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_block_graph_csc
