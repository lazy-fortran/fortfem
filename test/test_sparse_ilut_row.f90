program test_sparse_ilut_row
    use check, only: check_condition, check_summary
    use fortfem_api, only: apply_sparse_ilut, apply_sparse_ilut_jvp, &
        apply_sparse_ilut_vjp, build_sparse_ilut, build_sparse_ilut_row, &
        sparse_ilut_factor_t
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: n = 4
    integer, parameter :: rows(10) = [1, 2, 1, 2, 3, 2, 3, 4, 3, 4]
    integer, parameter :: columns(10) = [1, 1, 2, 2, 2, 3, 3, 3, 4, 4]
    real(dp), parameter :: values(10) = [4.0_dp, 2.0_dp, 1.0_dp, 5.0_dp, &
        3.0_dp, 1.0_dp, 6.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: right_hand_side(n) = [3.0_dp, -1.0_dp, 10.0_dp, &
        3.5_dp]
    real(dp), parameter :: expected_solution(n) = [1.0_dp, -1.0_dp, 2.0_dp, &
        0.5_dp]
    real(dp), parameter :: right_hand_side_dot(n) = [0.3_dp, -0.4_dp, 0.5_dp, &
        -0.2_dp]
    real(dp), parameter :: solution_bar(n) = [-0.2_dp, 0.6_dp, 0.4_dp, 0.7_dp]
    ! With zero retained off-diagonals the row ILUT sweep drops U entries as
    ! each row is finalized, so the independent diagonal oracle is the
    ! diagonal of the resulting no-fill row factor.
    real(dp), parameter :: expected_row_pivot(n) = [4.0_dp, 5.0_dp, 6.0_dp, &
        3.0_dp]
    type(csc_t) :: matrix
    type(sparse_ilut_factor_t) :: factor
    type(fortsparse_status_t) :: status
    real(dp), allocatable :: solution(:), solution_dot(:), right_hand_side_bar(:)
    real(dp) :: lhs, rhs
    logical :: all_passed

    all_passed = .true.
    call csc_from_triplet(n, n, rows, columns, values, matrix, status)
    call record_condition(status%code == 0, &
        "row ILUT matrix fixture has a valid CSC structure")

    call build_sparse_ilut_row(matrix, 0.0_dp, n, factor, status)
    call record_condition(status%code == 0, &
        "row ILUT accepts a full-fill, zero-drop factorization")
    call apply_sparse_ilut(factor, right_hand_side, solution, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(solution - expected_solution)) < 1.0e-13_dp, &
        "row ILUT solve matches the independent hand-LU solution")

    call build_sparse_ilut_row(matrix, 0.0_dp, 0, factor, status)
    call apply_sparse_ilut(factor, right_hand_side, solution, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(solution - right_hand_side/expected_row_pivot)) < 1.0e-13_dp, &
        "row ILUT zero-fill matches the independent no-fill diagonal oracle")

    call build_sparse_ilut_row(matrix, 0.0_dp, 1, factor, status)
    call apply_sparse_ilut(factor, right_hand_side, solution, status)
    call record_condition(status%code == 0 .and. all(solution == solution), &
        "row ILUT fill control returns a finite fixed-pattern apply")
    call apply_sparse_ilut_jvp(factor, right_hand_side_dot, solution_dot, status)
    call apply_sparse_ilut_vjp(factor, solution_bar, right_hand_side_bar, status)
    lhs = dot_product(solution_bar, solution_dot)
    rhs = dot_product(right_hand_side_bar, right_hand_side_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "row ILUT JVP and VJP satisfy the real dot-product identity")

    call build_sparse_ilut_row(matrix, -1.0_dp, 1, factor, status)
    call record_condition(status%code /= 0, &
        "row ILUT rejects a negative drop tolerance")
    call build_sparse_ilut_row(matrix, 0.0_dp, -1, factor, status)
    call record_condition(status%code /= 0, &
        "row ILUT rejects a negative row fill limit")

    call check_summary("row-oriented sparse ILUT")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_sparse_ilut_row
