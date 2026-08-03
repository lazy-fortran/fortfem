program test_sparse_tree_cotree
    use check, only: check_condition, check_summary
    use fortfem_feec, only: build_tree_cotree_gauge, &
        sparse_direct_solve_tree_cotree, sparse_direct_solve_tree_cotree_jvp, &
        sparse_direct_solve_tree_cotree_vjp, tree_cotree_gauge_t
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_t, csc_z_t
    implicit none

    integer, parameter :: incidence(3, 3) = reshape([ &
        -1, 1, 0, 0, -1, 1, -1, 0, 1], [3, 3])
    integer, parameter :: col_ptr(4) = [1, 4, 7, 10]
    integer, parameter :: row_idx(9) = [1, 2, 3, 1, 2, 3, 1, 2, 3]
    real(dp), parameter :: values(9) = [ &
        2.0_dp, 0.1_dp, 0.2_dp, 0.1_dp, 3.0_dp, 0.3_dp, &
        0.2_dp, 0.3_dp, 4.0_dp]
    real(dp), parameter :: rhs(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: values_dot(9) = [ &
        0.2_dp, -0.1_dp, 0.05_dp, 0.03_dp, 0.1_dp, -0.02_dp, &
        0.04_dp, -0.03_dp, 0.2_dp]
    real(dp), parameter :: rhs_dot(3) = [0.2_dp, -0.1_dp, 0.3_dp]
    type(tree_cotree_gauge_t) :: gauge
    type(csc_t) :: matrix
    type(csc_z_t) :: complex_matrix
    real(dp) :: solution(3), solution_dot(3), solution_bar(3)
    real(dp) :: matrix_values_bar(9), rhs_bar(3)
    real(dp) :: lhs, rhs_identity
    complex(dp) :: complex_solution(3), complex_rhs(3)
    integer :: status

    matrix%nrow = 3
    matrix%ncol = 3
    matrix%nnz = 9
    matrix%col_ptr = col_ptr
    matrix%row_idx = row_idx
    matrix%val = values
    call build_tree_cotree_gauge(incidence, gauge, status)
    call check_condition(status == 0, "sparse tree-cotree gauge initializes")

    call sparse_direct_solve_tree_cotree( &
        gauge, matrix, rhs, solution, status)
    call check_condition(status == 0 .and. &
        maxval(abs(solution - [0.0_dp, 0.0_dp, 0.75_dp])) < 1.0e-13_dp, &
        "sparse tree-cotree solve fixes tree unknowns and solves cotree block")

    call sparse_direct_solve_tree_cotree_jvp( &
        gauge, matrix, rhs, values_dot, rhs_dot, solution_dot, status)
    call check_condition(status == 0 .and. solution_dot(1) == 0.0_dp .and. &
        solution_dot(2) == 0.0_dp .and. &
        abs(solution_dot(3) - 0.0375_dp) < 1.0e-13_dp, &
        "sparse tree-cotree JVP keeps the fixed tree subspace")

    solution_bar = [0.4_dp, -0.3_dp, 0.8_dp]
    call sparse_direct_solve_tree_cotree_vjp( &
        gauge, matrix, rhs, solution, solution_bar, matrix_values_bar, &
        rhs_bar, status)
    lhs = dot_product(solution_bar, solution_dot)
    rhs_identity = sum(matrix_values_bar*values_dot) + &
        dot_product(rhs_bar, rhs_dot)
    call check_condition(status == 0 .and. &
        abs(lhs - rhs_identity) < 1.0e-12_dp, &
        "sparse tree-cotree derivatives satisfy the real adjoint identity")

    complex_matrix%nrow = 3
    complex_matrix%ncol = 3
    complex_matrix%nnz = 9
    complex_matrix%col_ptr = col_ptr
    complex_matrix%row_idx = row_idx
    complex_matrix%val = cmplx(values, 0.02_dp*values, dp)
    complex_rhs = cmplx(rhs, [-0.1_dp, 0.2_dp, 0.3_dp], dp)
    call sparse_direct_solve_tree_cotree( &
        gauge, complex_matrix, complex_rhs, complex_solution, status)
    call check_condition(status == 0 .and. &
        maxval(abs(complex_solution(1:2))) < 1.0e-13_dp .and. &
        abs(complex_solution(3) - complex_rhs(3)/complex_matrix%val(9)) < &
        1.0e-13_dp, "complex sparse tree-cotree solve uses the same selector")

    call check_summary("sparse tree-cotree direct solve")
end program test_sparse_tree_cotree
