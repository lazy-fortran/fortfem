program test_pcg_solver_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: pcg_solve, pcg_solve_jvp, pcg_solve_vjp, &
        solver_options_t, solver_stats_t
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
        5.0_dp, 1.0_dp, 0.0_dp, &
        1.0_dp, 4.0_dp, 0.3_dp, &
        0.0_dp, 0.3_dp, 2.5_dp], shape(matrix))
    matrix_dot = reshape([ &
        0.09_dp, -0.02_dp, 0.01_dp, &
        -0.02_dp, 0.05_dp, 0.03_dp, &
        0.01_dp, 0.03_dp, -0.06_dp], shape(matrix_dot))
    rhs = [0.8_dp, -0.25_dp, 0.55_dp]
    rhs_dot = [-0.06_dp, 0.11_dp, 0.04_dp]
    solution_bar = [0.22_dp, -0.35_dp, 0.41_dp]
    options%tolerance = 1.0e-12_dp
    options%tolerance_type = "absolute"
    options%max_iterations = 40
    options%preconditioner = "jacobi"

    solution = 0.0_dp
    call pcg_solve(matrix, rhs, solution, options, stats)
    call check_condition(stats%converged, "PCG primal solve converges")
    call pcg_solve_jvp( &
        matrix, rhs, solution, options, matrix_dot, rhs_dot, solution_dot, &
        status)
    solution_plus = 0.0_dp
    call pcg_solve( &
        matrix + step*matrix_dot, rhs + step*rhs_dot, solution_plus, options, &
        stats)
    solution_minus = 0.0_dp
    call pcg_solve( &
        matrix - step*matrix_dot, rhs - step*rhs_dot, solution_minus, options, &
        stats)
    solution_dot_fd = (solution_plus - solution_minus)/(2.0_dp*step)
    relative_error = maxval(abs(solution_dot - solution_dot_fd))/ &
        max(1.0_dp, maxval(abs(solution_dot)))
    call check_condition(status == 0 .and. relative_error < 2.0e-7_dp, &
        "PCG JVP matches re-evaluation of the converged state")

    call pcg_solve_vjp( &
        matrix, rhs, solution, options, solution_bar, matrix_bar, rhs_bar, &
        status)
    lhs = dot_product(solution_bar, solution_dot)
    rhs_identity = sum(matrix_bar*matrix_dot) + dot_product(rhs_bar, rhs_dot)
    call check_condition(status == 0 .and. abs(lhs - rhs_identity) < &
        3.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs_identity)), &
        "PCG VJP satisfies the implicit-solve adjoint identity")

    call check_summary("Differentiable preconditioned CG state")
end program test_pcg_solver_ad
