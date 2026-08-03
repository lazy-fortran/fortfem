program test_krylov_solver_ad
    use check, only: check_condition, check_summary
    use fortfem_advanced_solvers, only: &
        bicgstab_solve, bicgstab_solve_jvp, bicgstab_solve_vjp, &
        gmres_solve, gmres_solve_jvp, gmres_solve_vjp, solver_options_t, &
        solver_stats_t
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: matrix(3, 3), matrix_dot(3, 3), rhs(3), rhs_dot(3)
    real(dp) :: solution(3), solution_dot(3), solution_plus(3), solution_minus(3)
    real(dp) :: solution_dot_fd(3), solution_bar(3)
    real(dp) :: matrix_bar(3, 3), rhs_bar(3), lhs, rhs_identity
    real(dp) :: relative_error
    type(solver_options_t) :: options
    type(solver_stats_t) :: stats
    integer :: status

    matrix = reshape([ &
        3.8_dp, 0.7_dp, -0.2_dp, &
        -0.1_dp, 2.9_dp, 0.4_dp, &
        0.25_dp, -0.35_dp, 2.4_dp], shape(matrix))
    matrix_dot = reshape([ &
        0.06_dp, -0.02_dp, 0.03_dp, &
        0.01_dp, 0.05_dp, -0.04_dp, &
        -0.03_dp, 0.02_dp, 0.07_dp], shape(matrix_dot))
    rhs = [0.9_dp, -0.3_dp, 0.55_dp]
    rhs_dot = [-0.05_dp, 0.08_dp, 0.06_dp]
    solution_bar = [0.24_dp, -0.29_dp, 0.37_dp]
    options%tolerance = 1.0e-11_dp
    options%tolerance_type = "absolute"
    options%max_iterations = 80
    options%restart = 8
    options%preconditioner = "none"

    solution = 0.0_dp
    call gmres_solve(matrix, rhs, solution, options, stats)
    call check_condition(stats%converged, "GMRES primal solve converges")
    call gmres_solve_jvp( &
        matrix, rhs, solution, options, matrix_dot, rhs_dot, solution_dot, &
        status)
    solution_plus = 0.0_dp
    call gmres_solve( &
        matrix + step*matrix_dot, rhs + step*rhs_dot, solution_plus, options, &
        stats)
    solution_minus = 0.0_dp
    call gmres_solve( &
        matrix - step*matrix_dot, rhs - step*rhs_dot, solution_minus, options, &
        stats)
    solution_dot_fd = (solution_plus - solution_minus)/(2.0_dp*step)
    relative_error = maxval(abs(solution_dot - solution_dot_fd))/ &
        max(1.0_dp, maxval(abs(solution_dot)))
    call check_condition(status == 0 .and. relative_error < 5.0e-7_dp, &
        "GMRES JVP matches re-evaluation of the converged state")
    call gmres_solve_vjp( &
        matrix, rhs, solution, options, solution_bar, matrix_bar, rhs_bar, &
        status)
    lhs = dot_product(solution_bar, solution_dot)
    rhs_identity = sum(matrix_bar*matrix_dot) + dot_product(rhs_bar, rhs_dot)
    call check_condition(status == 0 .and. abs(lhs - rhs_identity) < &
        5.0e-7_dp*max(1.0_dp, abs(lhs), abs(rhs_identity)), &
        "GMRES VJP satisfies the implicit-solve adjoint identity")

    solution = 0.0_dp
    call bicgstab_solve(matrix, rhs, solution, options, stats)
    call check_condition(stats%converged, "BiCGSTAB primal solve converges")
    call bicgstab_solve_jvp( &
        matrix, rhs, solution, options, matrix_dot, rhs_dot, solution_dot, &
        status)
    solution_plus = 0.0_dp
    call bicgstab_solve( &
        matrix + step*matrix_dot, rhs + step*rhs_dot, solution_plus, options, &
        stats)
    solution_minus = 0.0_dp
    call bicgstab_solve( &
        matrix - step*matrix_dot, rhs - step*rhs_dot, solution_minus, options, &
        stats)
    solution_dot_fd = (solution_plus - solution_minus)/(2.0_dp*step)
    relative_error = maxval(abs(solution_dot - solution_dot_fd))/ &
        max(1.0_dp, maxval(abs(solution_dot)))
    call check_condition(status == 0 .and. relative_error < 5.0e-7_dp, &
        "BiCGSTAB JVP matches re-evaluation of the converged state")
    call bicgstab_solve_vjp( &
        matrix, rhs, solution, options, solution_bar, matrix_bar, rhs_bar, &
        status)
    lhs = dot_product(solution_bar, solution_dot)
    rhs_identity = sum(matrix_bar*matrix_dot) + dot_product(rhs_bar, rhs_dot)
    call check_condition(status == 0 .and. abs(lhs - rhs_identity) < &
        5.0e-7_dp*max(1.0_dp, abs(lhs), abs(rhs_identity)), &
        "BiCGSTAB VJP satisfies the implicit-solve adjoint identity")

    call check_summary("Differentiable nonsymmetric Krylov states")
end program test_krylov_solver_ad
