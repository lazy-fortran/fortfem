program test_sparse_direct_complex_derivatives
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_sparse_direct, only: sparse_direct_factor_adjoint_csc, &
        sparse_direct_factor_csc, sparse_direct_factor_t, sparse_direct_free, &
        sparse_direct_solve_csc, sparse_direct_solve_factored, &
        sparse_direct_solve_factored_jvp, sparse_direct_solve_factored_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: n = 3
    integer, parameter :: col_ptr(4) = [1, 3, 6, 8]
    integer, parameter :: row_ind(7) = [1, 2, 1, 2, 3, 2, 3]
    complex(dp), parameter :: values(7) = [ &
        cmplx(4.0_dp, 0.2_dp, dp), cmplx(0.7_dp, 0.4_dp, dp), &
        cmplx(1.0_dp, -0.3_dp, dp), cmplx(3.0_dp, -0.1_dp, dp), &
        cmplx(0.8_dp, -0.5_dp, dp), cmplx(1.0_dp, 0.2_dp, dp), &
        cmplx(2.0_dp, 0.3_dp, dp)]
    complex(dp), parameter :: values_dot(7) = [ &
        cmplx(0.2_dp, -0.05_dp, dp), cmplx(0.05_dp, 0.03_dp, dp), &
        cmplx(-0.1_dp, 0.02_dp, dp), cmplx(0.3_dp, -0.04_dp, dp), &
        cmplx(0.15_dp, -0.02_dp, dp), cmplx(-0.2_dp, 0.06_dp, dp), &
        cmplx(0.1_dp, 0.05_dp, dp)]
    real(dp), parameter :: step = 1.0e-6_dp
    complex(dp) :: b(n), b_bar(n), b_dot(n), b_minus(n), b_plus(n)
    complex(dp) :: values_bar(7), values_minus(7), values_plus(7)
    complex(dp) :: x(n), x_bar(n), x_dot(n), x_minus(n), x_plus(n)
    real(dp) :: lhs, rhs
    type(sparse_direct_factor_t) :: adjoint_factor, factor
    integer :: failures, status

    failures = 0
    b = [ &
        cmplx(6.0_dp, 0.5_dp, dp), cmplx(10.0_dp, -0.7_dp, dp), &
        cmplx(8.0_dp, 0.2_dp, dp)]
    b_dot = [ &
        cmplx(0.3_dp, -0.1_dp, dp), cmplx(-0.4_dp, 0.2_dp, dp), &
        cmplx(0.2_dp, 0.05_dp, dp)]
    x_bar = [ &
        cmplx(0.7_dp, 0.1_dp, dp), cmplx(-0.2_dp, 0.3_dp, dp), &
        cmplx(0.5_dp, -0.4_dp, dp)]

    call sparse_direct_factor_csc( &
        factor, n, col_ptr, row_ind, values, status)
    call check(status == 0, "complex primal factorization succeeds")
    call sparse_direct_solve_factored(factor, b, x, status)
    call check(status == 0, "complex primal solve succeeds")
    call sparse_direct_solve_factored_jvp( &
        factor, n, col_ptr, row_ind, values_dot, x, b_dot, x_dot, status)
    call check(status == 0, "complex solve JVP succeeds")

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
        "FortFEM complex solve JVP matches central difference")

    call sparse_direct_factor_adjoint_csc( &
        adjoint_factor, n, col_ptr, row_ind, values, status)
    call check(status == 0, "complex adjoint factorization succeeds")
    call sparse_direct_solve_factored_vjp( &
        adjoint_factor, n, col_ptr, row_ind, x, x_bar, &
        b_bar, values_bar, status)
    call check(status == 0, "complex solve VJP succeeds")
    lhs = real(sum(conjg(x_bar)*x_dot), dp)
    rhs = real( &
        sum(conjg(b_bar)*b_dot) + sum(conjg(values_bar)*values_dot), dp)
    call check(abs(lhs - rhs) < 2.0e-12_dp, &
        "FortFEM complex solve products satisfy the real adjoint identity")

    call sparse_direct_free(factor)
    call sparse_direct_free(adjoint_factor)
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

end program test_sparse_direct_complex_derivatives
