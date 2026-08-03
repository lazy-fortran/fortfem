program test_sparse_incomplete_lu
    use check, only: check_condition, check_summary
    use fortfem_sparse_incomplete_lu, only: apply_sparse_incomplete_lu, &
        apply_sparse_incomplete_lu_jvp, apply_sparse_incomplete_lu_vjp, &
        build_sparse_incomplete_lu, sparse_incomplete_lu_factor_t
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: n = 3
    integer, parameter :: rows(7) = [1, 2, 1, 2, 3, 2, 3]
    integer, parameter :: columns(7) = [1, 1, 2, 2, 2, 3, 3]
    real(dp), parameter :: values(7) = [4.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, &
        1.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: right_hand_side(n) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: right_hand_side_dot(n) = [0.3_dp, -0.4_dp, 0.5_dp]
    real(dp), parameter :: solution_bar(n) = [-0.2_dp, 0.6_dp, 0.4_dp]
    real(dp), parameter :: expected_solution(n) = [0.25_dp, 0.0_dp, 1.5_dp]
    type(csc_t) :: matrix, zero_pivot
    type(sparse_incomplete_lu_factor_t) :: factor
    type(fortsparse_status_t) :: status
    real(dp), allocatable :: solution(:), solution_dot(:), right_hand_side_bar(:)
    real(dp) :: lhs, rhs

    call csc_from_triplet(n, n, rows, columns, values, matrix, status)
    call check_condition(status%code == 0, &
        "sparse ILU(0) matrix fixture has a valid CSC structure")
    call build_sparse_incomplete_lu(matrix, factor, status)
    call check_condition(status%code == 0, &
        "sparse ILU(0) accepts a nonsymmetric CSC matrix")
    call apply_sparse_incomplete_lu(factor, right_hand_side, solution, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(solution - expected_solution)) < 1.0e-13_dp, &
        "sparse ILU(0) apply matches the independent triangular solve")

    call apply_sparse_incomplete_lu_jvp( &
        factor, right_hand_side_dot, solution_dot, status)
    call apply_sparse_incomplete_lu_vjp( &
        factor, solution_bar, right_hand_side_bar, status)
    lhs = dot_product(solution_bar, solution_dot)
    rhs = dot_product(right_hand_side_bar, right_hand_side_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "fixed sparse ILU(0) VJP satisfies the real dot-product identity")

    call csc_from_triplet(n, n, rows, columns, &
        [0.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], &
        zero_pivot, status)
    call build_sparse_incomplete_lu(zero_pivot, factor, status)
    call check_condition(status%code /= 0, &
        "sparse ILU(0) rejects a zero diagonal pivot")
    call check_summary("sparse incomplete LU")
end program test_sparse_incomplete_lu
