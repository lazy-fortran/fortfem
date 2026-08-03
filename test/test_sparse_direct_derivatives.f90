program test_sparse_direct_derivatives
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_factor_csc, &
        sparse_direct_factor_t, sparse_direct_factor_transpose_csc, &
        sparse_direct_free, sparse_direct_solve_csc, &
        sparse_direct_solve_factored, sparse_direct_solve_factored_jvp, &
        sparse_direct_solve_factored_vjp
    implicit none

    integer, parameter :: n = 3
    integer, parameter :: col_ptr(4) = [1, 3, 6, 8]
    integer, parameter :: row_ind(7) = [1, 2, 1, 2, 3, 2, 3]
    real(dp), parameter :: values(7) = [ &
        4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: values_dot(7) = [ &
        0.2_dp, 0.05_dp, -0.1_dp, 0.3_dp, 0.15_dp, -0.2_dp, 0.1_dp]
    real(dp), parameter :: step = 1.0e-6_dp
    real(dp) :: b(n), b_bar(n), b_dot(n), b_minus(n), b_plus(n)
    real(dp) :: values_bar(7), values_minus(7), values_plus(7)
    real(dp) :: x(n), x_bar(n), x_dot(n), x_minus(n), x_plus(n)
    real(dp) :: lhs, rhs
    type(sparse_direct_factor_t) :: factor, transpose_factor
    integer :: failures, status

    failures = 0
    b = [6.0_dp, 10.0_dp, 8.0_dp]
    b_dot = [0.3_dp, -0.4_dp, 0.2_dp]
    x_bar = [0.7_dp, -0.2_dp, 0.5_dp]

    call sparse_direct_factor_csc( &
        factor, n, col_ptr, row_ind, values, status)
    call check(status == 0, "primal factorization succeeds")
    call sparse_direct_solve_factored(factor, b, x, status)
    call check(status == 0, "primal solve succeeds")
    call sparse_direct_solve_factored_jvp( &
        factor, n, col_ptr, row_ind, values_dot, x, b_dot, x_dot, status)
    call check(status == 0, "solve JVP succeeds")

    values_plus = values + step*values_dot
    values_minus = values - step*values_dot
    b_plus = b + step*b_dot
    b_minus = b - step*b_dot
    call sparse_direct_solve_csc( &
        n, col_ptr, row_ind, values_plus, b_plus, x_plus, status)
    call sparse_direct_solve_csc( &
        n, col_ptr, row_ind, values_minus, b_minus, x_minus, status)
    call check(maxval(abs( &
        x_dot - (x_plus - x_minus)/(2.0_dp*step))) < 2.0e-9_dp, &
        "FortFEM solve JVP matches an independent central difference")

    call sparse_direct_factor_transpose_csc( &
        transpose_factor, n, col_ptr, row_ind, values, status)
    call check(status == 0, "transpose factorization succeeds")
    call sparse_direct_solve_factored_vjp( &
        transpose_factor, n, col_ptr, row_ind, x, x_bar, &
        b_bar, values_bar, status)
    call check(status == 0, "solve VJP succeeds")
    lhs = dot_product(x_bar, x_dot)
    rhs = dot_product(b_bar, b_dot) + dot_product(values_bar, values_dot)
    call check(abs(lhs - rhs) < 2.0e-12_dp, &
        "FortFEM solve JVP and VJP satisfy the adjoint identity")

    call sparse_direct_free(factor)
    call sparse_direct_free(transpose_factor)
    if (failures > 0) then
        write (error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write (error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_sparse_direct_derivatives
