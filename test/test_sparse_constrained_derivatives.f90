program test_sparse_constrained_derivatives
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: sparse_direct_solve_constrained, &
        sparse_direct_solve_constrained_jvp, &
        sparse_direct_solve_constrained_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: n = 4
    integer, parameter :: rows(10) = [1, 2, 1, 2, 3, 2, 3, 4, 3, 4]
    integer, parameter :: columns(10) = [1, 1, 2, 2, 2, 3, 3, 3, 4, 4]
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp), parameter :: values(10) = [ &
        4.0_dp, -1.0_dp, -1.0_dp, 3.5_dp, -0.6_dp, &
        -0.6_dp, 3.0_dp, -0.8_dp, -0.8_dp, 2.5_dp]
    real(dp), parameter :: values_dot(10) = [ &
        0.12_dp, -0.03_dp, 0.02_dp, -0.08_dp, 0.04_dp, &
        -0.01_dp, 0.07_dp, 0.03_dp, -0.02_dp, -0.05_dp]
    type(csc_t) :: matrix, minus_matrix, plus_matrix
    type(fortsparse_status_t) :: sparse_status
    logical :: constrained(n)
    real(dp) :: constrained_values(n), constrained_values_bar(n)
    real(dp) :: constrained_values_dot(n)
    real(dp) :: matrix_values_bar(10), rhs(n), rhs_bar(n), rhs_dot(n)
    real(dp) :: solution(n), solution_bar(n), solution_dot(n)
    real(dp) :: solution_minus(n), solution_plus(n)
    real(dp) :: lhs, rhs_product
    integer :: failures, status

    failures = 0
    constrained = [.true., .false., .false., .true.]
    constrained_values = [0.4_dp, 0.0_dp, 0.0_dp, -0.3_dp]
    constrained_values_dot = [-0.06_dp, 0.0_dp, 0.0_dp, 0.09_dp]
    rhs = [0.2_dp, 1.1_dp, -0.7_dp, 0.5_dp]
    rhs_dot = [-0.03_dp, 0.08_dp, 0.04_dp, -0.02_dp]
    solution_bar = [0.3_dp, -0.4_dp, 0.25_dp, -0.15_dp]
    call make_matrix(values, matrix)

    call sparse_direct_solve_constrained( &
        matrix, rhs, constrained, constrained_values, solution, status)
    call check(status == 0, "constrained primal solve succeeds")
    call sparse_direct_solve_constrained_jvp( &
        matrix, rhs, constrained, constrained_values, values_dot, rhs_dot, &
        constrained_values_dot, solution_dot, status)
    call check(status == 0, "constrained solve JVP succeeds")
    call make_matrix(values + step*values_dot, plus_matrix)
    call make_matrix(values - step*values_dot, minus_matrix)
    call sparse_direct_solve_constrained( &
        plus_matrix, rhs + step*rhs_dot, constrained, &
        constrained_values + step*constrained_values_dot, solution_plus, &
        status)
    call sparse_direct_solve_constrained( &
        minus_matrix, rhs - step*rhs_dot, constrained, &
        constrained_values - step*constrained_values_dot, solution_minus, &
        status)
    call check(maxval(abs( &
        solution_dot - (solution_plus - solution_minus)/(2.0_dp*step))) < &
        2.0e-9_dp, &
        "constrained JVP matches independent central solves")

    call sparse_direct_solve_constrained_vjp( &
        matrix, rhs, constrained, constrained_values, solution, solution_bar, &
        matrix_values_bar, rhs_bar, constrained_values_bar, status)
    call check(status == 0, "constrained solve VJP succeeds")
    lhs = dot_product(solution_bar, solution_dot)
    rhs_product = dot_product(matrix_values_bar, values_dot) + &
        dot_product(rhs_bar, rhs_dot) + &
        dot_product(constrained_values_bar, constrained_values_dot)
    call check(abs(lhs - rhs_product) < 2.0e-12_dp, &
        "constrained products satisfy the complete adjoint identity")

    if (failures > 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine make_matrix(matrix_values, output)
        real(dp), intent(in) :: matrix_values(10)
        type(csc_t), intent(out) :: output

        call csc_from_triplet( &
            n, n, rows, columns, matrix_values, output, sparse_status)
        call check(sparse_status%code == 0, "matrix fixture is valid CSC")
    end subroutine make_matrix

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_sparse_constrained_derivatives
