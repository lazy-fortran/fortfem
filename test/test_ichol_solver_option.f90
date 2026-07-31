program test_ichol_solver_option
    use check, only: check_condition, check_summary
    use fortfem_api, only: ichol_preconditioner, pcg_solve, pcg_solve_jvp, &
        solver_options, solver_options_t, solver_stats_t
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: matrix(4, 4) = reshape([ &
        4.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 4.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 4.0_dp, 1.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, 3.0_dp], [4, 4])
    real(dp), parameter :: exact_solution(4) = [1.0_dp, -2.0_dp, &
        0.5_dp, 3.0_dp]
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: right_hand_side(4), solution(4), solution_dot(4)
    real(dp) :: solution_plus(4), solution_minus(4), solution_dot_fd(4)
    real(dp) :: matrix_dot(4, 4), right_hand_side_dot(4), relative_error
    type(solver_options_t) :: options
    type(solver_stats_t) :: stats
    integer :: status

    right_hand_side = matmul(matrix, exact_solution)
    solution = 0.0_dp
    options = solver_options(method="pcg", &
        preconditioner=ichol_preconditioner(), tolerance=1.0e-12_dp, &
        max_iterations=20)

    call pcg_solve(matrix, right_hand_side, solution, options, stats)
    call check_condition(stats%converged, &
        "PCG accepts the differentiable ICHOL option")
    call check_condition(maxval(abs(solution - exact_solution)) < 1.0e-10_dp, &
        "PCG+ICHOL satisfies the independent exact-solution oracle")
    call check_condition(maxval(abs(matmul(matrix, solution) - &
        right_hand_side)) < 1.0e-10_dp, &
        "PCG+ICHOL satisfies the residual oracle")

    matrix_dot = 0.0_dp
    matrix_dot(1, 1) = 0.2_dp
    matrix_dot(2, 2) = -0.1_dp
    matrix_dot(3, 3) = 0.15_dp
    matrix_dot(4, 4) = 0.05_dp
    right_hand_side_dot = [0.1_dp, -0.2_dp, 0.05_dp, 0.3_dp]
    call pcg_solve_jvp(matrix, right_hand_side, solution, options, &
        matrix_dot, right_hand_side_dot, solution_dot, status)
    solution_plus = 0.0_dp
    call pcg_solve(matrix + step*matrix_dot, &
        right_hand_side + step*right_hand_side_dot, solution_plus, options, &
        stats)
    solution_minus = 0.0_dp
    call pcg_solve(matrix - step*matrix_dot, &
        right_hand_side - step*right_hand_side_dot, solution_minus, options, &
        stats)
    solution_dot_fd = (solution_plus - solution_minus)/(2.0_dp*step)
    relative_error = maxval(abs(solution_dot - solution_dot_fd))/ &
        max(1.0_dp, maxval(abs(solution_dot)))
    call check_condition(status == 0 .and. relative_error < 2.0e-7_dp, &
        "PCG+ICHOL JVP matches the independent re-evaluation oracle")

    call check_summary("ICHOL solver option")
end program test_ichol_solver_option
