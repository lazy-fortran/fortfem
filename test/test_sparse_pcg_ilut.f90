program test_sparse_pcg_ilut
    use check, only: check_condition, check_summary
    use fortfem_advanced_solvers, only: solve_sparse, solver_options, solver_options_t, &
        solver_stats_t
    use fortfem_sparse_matrix, only: sparse_matrix_t
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, fortsparse_status_t
    implicit none

    integer, parameter :: n = 3
    integer, parameter :: rows(7) = [1, 2, 1, 2, 3, 2, 3]
    integer, parameter :: columns(7) = [1, 1, 2, 2, 2, 3, 3]
    real(dp), parameter :: values(7) = [4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, &
        1.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: right_hand_side(n) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: expected_solution(n) = [2.0_dp/9.0_dp, &
        1.0_dp/9.0_dp, 13.0_dp/9.0_dp]
    type(sparse_matrix_t) :: matrix
    type(solver_options_t) :: options
    type(solver_stats_t) :: stats
    type(fortsparse_status_t) :: status
    real(dp) :: solution(n)
    logical :: all_passed

    all_passed = .true.
    call csc_from_triplet(n, n, rows, columns, values, matrix, status)
    call record_condition(status%code == 0, &
        "sparse PCG ILUT matrix has a valid CSC structure")
    options = solver_options(method="pcg", tolerance=1.0e-12_dp, &
        max_iterations=1, preconditioner="sparse_ilut", &
        drop_tolerance=0.0_dp, fill_level=n)
    solution = 0.0_dp
    call solve_sparse(matrix, right_hand_side, solution, options, stats)
    call record_condition(stats%converged .and. stats%iterations <= 1 .and. &
        maxval(abs(solution - expected_solution)) < 1.0e-12_dp, &
        "sparse PCG ILUT matches the exact row-factor oracle")
    options%preconditioner = "ilu"
    solution = 0.0_dp
    call solve_sparse(matrix, right_hand_side, solution, options, stats)
    call record_condition(stats%converged .and. stats%iterations <= 1 .and. &
        maxval(abs(solution - expected_solution)) < 1.0e-12_dp, &
        "sparse PCG ILU alias matches the exact row-factor oracle")
    options%preconditioner = "ichol"
    solution = 0.0_dp
    call solve_sparse(matrix, right_hand_side, solution, options, stats)
    call record_condition(stats%converged .and. stats%iterations <= 1 .and. &
        maxval(abs(solution - expected_solution)) < 1.0e-12_dp, &
        "sparse PCG ICHOL alias matches the exact sparse-factor oracle")
    call check_summary("sparse PCG ILUT integration")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_sparse_pcg_ilut
