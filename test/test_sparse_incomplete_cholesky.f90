program test_sparse_incomplete_cholesky
    use check, only: check_condition, check_summary
    use fortfem_sparse_incomplete_cholesky, only: apply_sparse_incomplete_cholesky, &
        apply_sparse_incomplete_cholesky_jvp, &
        apply_sparse_incomplete_cholesky_vjp, build_sparse_incomplete_cholesky, &
        sparse_incomplete_cholesky_factor_t
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: n = 3
    integer, parameter :: rows(7) = [1, 2, 1, 2, 3, 2, 3]
    integer, parameter :: columns(7) = [1, 1, 2, 2, 2, 3, 3]
    real(dp), parameter :: values(7) = [4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, &
        1.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: nonsymmetric_values(7) = [4.0_dp, 1.0_dp, 2.0_dp, &
        3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: right_hand_side(n) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: right_hand_side_dot(n) = [0.3_dp, -0.4_dp, 0.5_dp]
    real(dp), parameter :: right_hand_side_bar(n) = [-0.2_dp, 0.6_dp, 0.4_dp]
    real(dp), parameter :: expected_solution(n) = [2.0_dp/9.0_dp, &
        1.0_dp/9.0_dp, 13.0_dp/9.0_dp]
    type(csc_t) :: matrix, nonsymmetric, indefinite
    type(sparse_incomplete_cholesky_factor_t) :: factor
    type(fortsparse_status_t) :: status
    real(dp), allocatable :: solution(:), solution_dot(:), solution_bar(:)
    real(dp) :: lhs, rhs

    call csc_from_triplet(n, n, rows, columns, values, matrix, status)
    call check_condition(status%code == 0, &
        "sparse IC(0) matrix fixture has a valid CSC structure")
    call build_sparse_incomplete_cholesky(matrix, factor, status)
    call check_condition(status%code == 0, &
        "sparse IC(0) accepts a symmetric positive-definite CSC matrix")
    call apply_sparse_incomplete_cholesky(factor, right_hand_side, solution, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(solution - expected_solution)) < 1.0e-13_dp, &
        "sparse IC(0) apply matches the independent tridiagonal solve")

    call apply_sparse_incomplete_cholesky_jvp( &
        factor, right_hand_side_dot, solution_dot, status)
    call apply_sparse_incomplete_cholesky_vjp( &
        factor, right_hand_side_bar, solution_bar, status)
    lhs = dot_product(solution_bar, right_hand_side_dot)
    rhs = dot_product(right_hand_side_bar, solution_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "fixed sparse IC(0) apply has the symmetric VJP identity")

    call csc_from_triplet(n, n, rows, columns, nonsymmetric_values, nonsymmetric, &
        status)
    call build_sparse_incomplete_cholesky(nonsymmetric, factor, status)
    call check_condition(status%code /= 0, &
        "sparse IC(0) rejects a nonsymmetric matrix")

    call csc_from_triplet(n, n, rows, columns, &
        [4.0_dp, 1.0_dp, 1.0_dp, -3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], &
        indefinite, status)
    call build_sparse_incomplete_cholesky(indefinite, factor, status)
    call check_condition(status%code /= 0, &
        "sparse IC(0) rejects a non-positive pivot")
    call check_summary("sparse incomplete Cholesky")
end program test_sparse_incomplete_cholesky
