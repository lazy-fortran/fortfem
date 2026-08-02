program test_block_graph_csc
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_block_graph_csc
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
    real(dp) :: expected(3), dense_matrix(3, 3), applied(3)
    complex(dp) :: complex_values(9), complex_applied(3), complex_expected(3)
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
