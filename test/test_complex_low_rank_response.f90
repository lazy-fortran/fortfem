program test_complex_low_rank_response
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        apply_complex_low_rank_matrix, &
        apply_complex_low_rank_matrix_jvp, &
        apply_complex_low_rank_matrix_vjp, &
        compress_complex_matrix_cross, &
        complex_low_rank_matrix_t, &
        materialize_complex_low_rank_matrix, &
        validate_complex_low_rank_matrix
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: row_count = 4, column_count = 3, rank = 2
    type(complex_low_rank_matrix_t) :: low_rank
    complex(dp) :: left(row_count, rank), right(column_count, rank)
    complex(dp) :: matrix(row_count, column_count), materialized(row_count, column_count)
    complex(dp) :: left_dot(row_count, rank), right_dot(column_count, rank)
    complex(dp) :: x(column_count), x_dot(column_count), y(row_count), y_dot(row_count)
    complex(dp) :: y_plus(row_count), y_minus(row_count), y_bar(row_count), x_bar(column_count)
    complex(dp), allocatable :: left_bar(:, :), right_bar(:, :)
    real(dp) :: epsilon, finite_difference_error, lhs, rhs
    integer :: status

    left = reshape([ &
        cmplx(1.0_dp, 0.1_dp, dp), cmplx(-0.2_dp, 0.3_dp, dp), &
        cmplx(0.4_dp, -0.5_dp, dp), cmplx(0.7_dp, 0.2_dp, dp), &
        cmplx(-0.1_dp, 0.6_dp, dp), cmplx(0.3_dp, -0.4_dp, dp), &
        cmplx(0.8_dp, 0.0_dp, dp), cmplx(-0.5_dp, 0.2_dp, dp)], shape(left))
    right = reshape([ &
        cmplx(0.6_dp, -0.2_dp, dp), cmplx(0.1_dp, 0.5_dp, dp), &
        cmplx(-0.4_dp, 0.3_dp, dp), cmplx(0.7_dp, -0.1_dp, dp), &
        cmplx(0.2_dp, 0.8_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp)], shape(right))
    matrix = matmul(left, transpose(right))

    call compress_complex_matrix_cross(matrix, 1.0e-12_dp, rank, low_rank, status)
    call check_condition(status == 0, "complex cross compression succeeds")
    call check_condition(validate_complex_low_rank_matrix(low_rank, status) .and. status == 0, &
        "compressed response validates")
    call materialize_complex_low_rank_matrix(low_rank, materialized, status)
    call check_condition(status == 0 .and. maxval(abs(materialized - matrix)) < 1.0e-11_dp, &
        "cross approximation matches the dense oracle")

    x = [ &
        cmplx(0.2_dp, -0.4_dp, dp), cmplx(-0.5_dp, 0.1_dp, dp), &
        cmplx(0.7_dp, 0.3_dp, dp)]
    call apply_complex_low_rank_matrix(low_rank, x, y, status)
    call check_condition(status == 0 .and. maxval(abs(y - matmul(matrix, x))) < 1.0e-11_dp, &
        "matrix-free low-rank action matches dense multiplication")

    left_dot = reshape([ &
        cmplx(0.03_dp, -0.02_dp, dp), cmplx(-0.01_dp, 0.04_dp, dp), &
        cmplx(0.02_dp, 0.01_dp, dp), cmplx(-0.05_dp, 0.03_dp, dp), &
        cmplx(0.04_dp, -0.01_dp, dp), cmplx(0.02_dp, -0.03_dp, dp), &
        cmplx(-0.02_dp, 0.05_dp, dp), cmplx(0.01_dp, 0.02_dp, dp)], shape(left_dot))
    right_dot = reshape([ &
        cmplx(-0.02_dp, 0.01_dp, dp), cmplx(0.03_dp, -0.04_dp, dp), &
        cmplx(0.01_dp, 0.02_dp, dp), cmplx(-0.02_dp, 0.03_dp, dp), &
        cmplx(0.04_dp, -0.01_dp, dp), cmplx(0.02_dp, 0.05_dp, dp)], shape(right_dot))
    x_dot = [ &
        cmplx(-0.03_dp, 0.02_dp, dp), cmplx(0.01_dp, 0.04_dp, dp), &
        cmplx(0.02_dp, -0.01_dp, dp)]
    call apply_complex_low_rank_matrix_jvp( &
        low_rank, left_dot, right_dot, x, x_dot, y_dot, status)
    epsilon = 1.0e-7_dp
    low_rank%left = low_rank%left + epsilon*left_dot
    low_rank%right = low_rank%right + epsilon*right_dot
    call apply_complex_low_rank_matrix(low_rank, x + epsilon*x_dot, y_plus, status)
    low_rank%left = low_rank%left - 2.0_dp*epsilon*left_dot
    low_rank%right = low_rank%right - 2.0_dp*epsilon*right_dot
    call apply_complex_low_rank_matrix(low_rank, x - epsilon*x_dot, y_minus, status)
    finite_difference_error = maxval(abs(y_dot - (y_plus - y_minus)/(2.0_dp*epsilon)))
    call check_condition(finite_difference_error < 2.0e-8_dp, &
        "low-rank JVP matches central differences")
    low_rank%left = low_rank%left + epsilon*left_dot
    low_rank%right = low_rank%right + epsilon*right_dot

    y_bar = [ &
        cmplx(0.4_dp, -0.1_dp, dp), cmplx(-0.2_dp, 0.3_dp, dp), &
        cmplx(0.1_dp, 0.5_dp, dp), cmplx(-0.3_dp, 0.2_dp, dp)]
    call apply_complex_low_rank_matrix_vjp( &
        low_rank, x, y_bar, left_bar, right_bar, x_bar, status)
    lhs = real(sum(conjg(y_bar)*y_dot), dp)
    rhs = real(sum(conjg(left_bar)*left_dot) + sum(conjg(right_bar)*right_dot) + &
        sum(conjg(x_bar)*x_dot), dp)
    call check_condition(status == 0 .and. abs(lhs - rhs) < 2.0e-12_dp, &
        "low-rank VJP satisfies the real complex adjoint identity")
    call check_summary("complex low-rank response")
end program test_complex_low_rank_response
