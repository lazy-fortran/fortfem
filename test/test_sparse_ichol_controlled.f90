program test_sparse_ichol_controlled
    use check, only: check_condition, check_summary
    use fortfem_api, only: apply_sparse_incomplete_cholesky, &
        apply_sparse_incomplete_cholesky_jvp, apply_sparse_incomplete_cholesky_vjp, &
        build_sparse_ichol, sparse_incomplete_cholesky_factor_t
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
    real(dp), parameter :: expected_solution(n) = [2.0_dp/9.0_dp, &
        1.0_dp/9.0_dp, 13.0_dp/9.0_dp]
    real(dp), parameter :: expected_pivot_energy(n) = [4.0_dp, 11.0_dp/4.0_dp, &
        18.0_dp/11.0_dp]
    real(dp), parameter :: right_hand_side_dot(n) = [0.3_dp, -0.4_dp, 0.5_dp]
    real(dp), parameter :: solution_bar(n) = [-0.2_dp, 0.6_dp, 0.4_dp]
    type(csc_t) :: matrix, nonsymmetric, indefinite
    type(sparse_incomplete_cholesky_factor_t) :: factor
    type(fortsparse_status_t) :: status
    real(dp), allocatable :: solution(:), solution_dot(:), solution_bar_rhs(:)
    real(dp) :: lhs, rhs
    logical :: all_passed

    all_passed = .true.
    call csc_from_triplet(n, n, rows, columns, values, matrix, status)
    call record_condition(status%code == 0, &
        "controlled ICHOL matrix fixture has a valid CSC structure")

    call build_sparse_ichol(matrix, 0.0_dp, n, factor, status)
    call record_condition(status%code == 0, &
        "controlled ICHOL accepts a full-fill SPD matrix")
    call apply_sparse_incomplete_cholesky(factor, right_hand_side, solution, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(solution - expected_solution)) < 1.0e-13_dp, &
        "full-fill ICHOL solve matches the independent tridiagonal solution")

    call build_sparse_ichol(matrix, 0.0_dp, 0, factor, status)
    call apply_sparse_incomplete_cholesky(factor, right_hand_side, solution, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(solution - right_hand_side/expected_pivot_energy)) < 1.0e-13_dp, &
        "zero-fill ICHOL retains the hand-derived Cholesky diagonal")

    call build_sparse_ichol(matrix, 0.0_dp, 1, factor, status)
    call apply_sparse_incomplete_cholesky(factor, right_hand_side, solution, status)
    call record_condition(status%code == 0 .and. all(solution == solution), &
        "fill-controlled ICHOL returns a finite fixed-pattern apply")

    call apply_sparse_incomplete_cholesky_jvp( &
        factor, right_hand_side_dot, solution_dot, status)
    call apply_sparse_incomplete_cholesky_vjp( &
        factor, solution_bar, solution_bar_rhs, status)
    lhs = dot_product(solution_bar, solution_dot)
    rhs = dot_product(solution_bar_rhs, right_hand_side_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "fixed-factor ICHOL JVP and VJP satisfy the symmetric dot identity")

    call build_sparse_ichol(matrix, -1.0_dp, 1, factor, status)
    call record_condition(status%code /= 0, &
        "controlled ICHOL rejects a negative drop tolerance")
    call build_sparse_ichol(matrix, 0.0_dp, -1, factor, status)
    call record_condition(status%code /= 0, &
        "controlled ICHOL rejects a negative fill limit")

    call csc_from_triplet(n, n, rows, columns, nonsymmetric_values, &
        nonsymmetric, status)
    call build_sparse_ichol(nonsymmetric, 0.0_dp, n, factor, status)
    call record_condition(status%code /= 0, &
        "controlled ICHOL rejects a nonsymmetric matrix")

    call csc_from_triplet(n, n, rows, columns, &
        [4.0_dp, 1.0_dp, 1.0_dp, -3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], &
        indefinite, status)
    call build_sparse_ichol(indefinite, 0.0_dp, n, factor, status)
    call record_condition(status%code /= 0, &
        "controlled ICHOL rejects a non-positive pivot")

    call check_summary("controlled sparse ICHOL")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_sparse_ichol_controlled
