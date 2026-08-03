program test_solver_resource_budget
    !! Independent contract test for caller-owned solver resource budgets.
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        initialize_solver_resource_budget, &
        evaluate_solver_resource_usage, &
        solver_resource_budget_t
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t
    implicit none

    integer, parameter :: dp = real64
    type(solver_resource_budget_t) :: budget
    type(fortsparse_status_t) :: status
    logical :: within_budget
    real(dp) :: expected_wall_margin, expected_memory_margin

    call initialize_solver_resource_budget( &
        budget, max_wall_seconds=2.5_dp, max_peak_memory_bytes=8192_int64, &
        max_repetitions=4, status=status)
    call check_condition(status%code == FORTSPARSE_OK, &
        "positive solver resource budget initializes")

    call evaluate_solver_resource_usage( &
        budget, wall_seconds=1.25_dp, peak_memory_bytes=4096_int64, &
        repetitions=3, within_budget=within_budget, status=status)
    call check_condition(status%code == FORTSPARSE_OK .and. within_budget, &
        "usage inside all resource limits is accepted")

    expected_wall_margin = 2.5_dp - 1.25_dp
    expected_memory_margin = real(8192_int64 - 4096_int64, dp)
    call check_condition(abs(expected_wall_margin - 1.25_dp) < 1.0e-14_dp .and. &
        abs(expected_memory_margin - 4096.0_dp) < 1.0e-14_dp, &
        "independent budget margin oracle agrees")

    call evaluate_solver_resource_usage( &
        budget, wall_seconds=3.0_dp, peak_memory_bytes=4096_int64, &
        repetitions=3, within_budget=within_budget, status=status)
    call check_condition(status%code == FORTSPARSE_OK .and. .not. within_budget, &
        "over-time usage is reported without corrupting the contract")

    call evaluate_solver_resource_usage( &
        budget, wall_seconds=1.0_dp, peak_memory_bytes=16384_int64, &
        repetitions=3, within_budget=within_budget, status=status)
    call check_condition(status%code == FORTSPARSE_OK .and. .not. within_budget, &
        "over-memory usage is reported without corrupting the contract")

    call initialize_solver_resource_budget( &
        budget, max_wall_seconds=0.0_dp, max_peak_memory_bytes=8192_int64, &
        max_repetitions=4, status=status)
    call check_condition(status%code == FORTSPARSE_INVALID_MATRIX, &
        "nonpositive wall-time budget is rejected")

    call initialize_solver_resource_budget( &
        budget, max_wall_seconds=2.5_dp, max_peak_memory_bytes=0_int64, &
        max_repetitions=4, status=status)
    call check_condition(status%code == FORTSPARSE_INVALID_MATRIX, &
        "nonpositive memory budget is rejected")

    call check_summary("solver resource budget")
end program test_solver_resource_budget
