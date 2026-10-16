program test_sparse_direct_backend
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_factor_t, &
        sparse_direct_factor_csc, sparse_direct_factor_transpose_csc, &
        sparse_direct_factor_adjoint_csc, sparse_direct_free, &
        sparse_direct_solve_constrained, sparse_direct_solve_csc, &
        sparse_direct_solve_factored
    use fortsparse, only: csc_t, FORTSPARSE_BACKEND_NATIVE
    use check, only: check_condition, check_summary
    implicit none

    logical :: all_passed

    all_passed = .true.
    call test_real_solve()
    call test_nonzero_constraint_elimination()
    call test_complex_solve()
    call test_reused_complex_factor()
    call test_pivoted_operations_and_lifetime()
    call test_failed_factor_recovery()
    call check_summary("Fortsparse direct backend")
    if (.not. all_passed) error stop 1

contains

    subroutine test_real_solve()
        integer, parameter :: col_ptr(3) = [1, 3, 5]
        integer, parameter :: row_ind(4) = [1, 2, 1, 2]
        real(dp), parameter :: values(4) = &
            [4.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]
        real(dp), parameter :: rhs(2) = [6.0_dp, 8.0_dp]
        real(dp), parameter :: expected(2) = [1.0_dp, 2.0_dp]
        real(dp) :: solution(2)
        integer :: status

        call sparse_direct_solve_csc( &
            2, col_ptr, row_ind, values, rhs, solution, status)

        call record_condition(status == 0, &
            "Real sparse factorization and solve succeeds")
        call record_condition(maxval(abs(solution - expected)) < 1.0e-13_dp, &
            "Real sparse solve matches the exact solution")
    end subroutine test_real_solve

    subroutine test_nonzero_constraint_elimination()
        logical, parameter :: constrained(2) = [.true., .false.]
        real(dp), parameter :: constrained_values(2) = [2.0_dp, 0.0_dp]
        real(dp), parameter :: expected(2) = [2.0_dp, 5.0_dp]
        real(dp), parameter :: rhs(2) = [13.0_dp, 19.0_dp]
        type(csc_t) :: matrix
        real(dp) :: solution(2)
        integer :: status

        matrix%nrow = 2
        matrix%ncol = 2
        matrix%col_ptr = [1, 3, 5]
        matrix%row_idx = [1, 2, 1, 2]
        matrix%val = [4.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]
        call sparse_direct_solve_constrained( &
            matrix, rhs, constrained, constrained_values, solution, status)

        call record_condition(status == 0, &
            "Nonzero sparse constraint elimination succeeds")
        call record_condition(maxval(abs(solution - expected)) < 1.0e-13_dp, &
            "Nonzero sparse constraints reproduce the exact linear solution")
    end subroutine test_nonzero_constraint_elimination

    subroutine test_complex_solve()
        integer, parameter :: col_ptr(3) = [1, 2, 3]
        integer, parameter :: row_ind(2) = [1, 2]
        complex(dp), parameter :: values(2) = &
            [cmplx(1.0_dp, 1.0_dp, dp), cmplx(2.0_dp, -1.0_dp, dp)]
        complex(dp), parameter :: rhs(2) = &
            [cmplx(3.0_dp, 1.0_dp, dp), cmplx(0.0_dp, 5.0_dp, dp)]
        complex(dp), parameter :: expected(2) = &
            [cmplx(2.0_dp, -1.0_dp, dp), cmplx(-1.0_dp, 2.0_dp, dp)]
        complex(dp) :: solution(2)
        integer :: status

        call sparse_direct_solve_csc( &
            2, col_ptr, row_ind, values, rhs, solution, status)

        call record_condition(status == 0, &
            "Complex sparse factorization and solve succeeds")
        call record_condition(maxval(abs(solution - expected)) < 1.0e-13_dp, &
            "Complex sparse solve matches the exact solution")
    end subroutine test_complex_solve

    subroutine test_reused_complex_factor()
        integer, parameter :: col_ptr(3) = [1, 3, 5]
        integer, parameter :: row_ind(4) = [1, 2, 1, 2]
        complex(dp), parameter :: values(4) = [ &
            cmplx(3.0_dp, 1.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
            cmplx(1.0_dp, 0.0_dp, dp), cmplx(2.0_dp, -1.0_dp, dp)]
        complex(dp), parameter :: expected_a(2) = [ &
            cmplx(1.0_dp, -1.0_dp, dp), cmplx(2.0_dp, 0.5_dp, dp)]
        complex(dp), parameter :: expected_b(2) = [ &
            cmplx(-2.0_dp, 0.25_dp, dp), cmplx(0.5_dp, 3.0_dp, dp)]
        complex(dp) :: rhs(2), solution(2)
        type(sparse_direct_factor_t) :: factor
        integer :: status

        call sparse_direct_factor_csc( &
            factor, 2, col_ptr, row_ind, values, status)
        call record_condition(status == 0, &
            "Complex sparse factorization for reuse succeeds")
        call record_condition(factor%backend_id == FORTSPARSE_BACKEND_NATIVE, &
            "Retained factor explicitly selects the native backend")

        rhs = matmul(reshape(values, [2, 2]), expected_a)
        call sparse_direct_solve_factored(factor, rhs, solution, status)
        call record_condition(status == 0 .and. &
            maxval(abs(solution - expected_a)) < 1.0e-13_dp, &
            "Retained factor solves the first exact right-hand side")

        rhs = matmul(reshape(values, [2, 2]), expected_b)
        call sparse_direct_solve_factored(factor, rhs, solution, status)
        call record_condition(status == 0 .and. &
            maxval(abs(solution - expected_b)) < 1.0e-13_dp, &
            "Retained factor solves a distinct exact right-hand side")
        call sparse_direct_free(factor)
    end subroutine test_reused_complex_factor

    subroutine test_pivoted_operations_and_lifetime()
        integer, parameter :: ptr(4) = [1, 4, 7, 10]
        integer, parameter :: rows(9) = [1, 2, 3, 1, 2, 3, 1, 2, 3]
        real(dp), parameter :: a(3, 3) = reshape([ &
            0.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 0.0_dp, 4.0_dp, &
            1.0_dp, 5.0_dp, 2.0_dp], [3, 3])
        real(dp), parameter :: expected(3) = [2.0_dp, -3.0_dp, 1.0_dp]
        complex(dp) :: z(3, 3), zx(3), zb(3), actual_z(3)
        real(dp) :: actual(3), rhs(3)
        real(dp), allocatable :: values(:)
        type(sparse_direct_factor_t) :: factor
        integer :: status

        values = reshape(a, [9])
        call sparse_direct_factor_csc(factor, 3, ptr, rows, values, status)
        call record_condition(status == 0, "Zero diagonal admits row pivoting")
        deallocate(values)
        rhs = matmul(a, expected)
        call sparse_direct_solve_factored(factor, rhs, actual, status)
        call record_condition(status == 0 .and. &
            maxval(abs(matmul(a, actual) - rhs)) < 1.0e-12_dp, &
            "Factors own their matrix after caller values are released")
        call sparse_direct_free(factor)
        call sparse_direct_free(factor)
        call sparse_direct_solve_factored(factor, rhs, actual, status)
        call record_condition(status /= 0, "Solve after repeated free fails safely")

        call sparse_direct_factor_transpose_csc( &
            factor, 3, ptr, rows, reshape(a, [9]), status)
        rhs = matmul(transpose(a), expected)
        call sparse_direct_solve_factored(factor, rhs, actual, status)
        call record_condition(status == 0 .and. &
            maxval(abs(actual - expected)) < 1.0e-12_dp, &
            "Pivoted transpose factor reproduces an independent exact vector")
        call sparse_direct_free(factor)

        z = cmplx(a, transpose(a), dp)
        zx = cmplx(expected, [1.0_dp, -2.0_dp, 0.5_dp], dp)
        call sparse_direct_factor_adjoint_csc( &
            factor, 3, ptr, rows, reshape(z, [9]), status)
        zb = matmul(conjg(transpose(z)), zx)
        call sparse_direct_solve_factored(factor, zb, actual_z, status)
        call record_condition(status == 0 .and. &
            maxval(abs(actual_z - zx)) < 1.0e-12_dp, &
            "Complex adjoint factor preserves conjugation and row pivots")
        call sparse_direct_free(factor)
    end subroutine test_pivoted_operations_and_lifetime

    subroutine test_failed_factor_recovery()
        integer, parameter :: ptr(3) = [1, 3, 5]
        integer, parameter :: rows(4) = [1, 2, 1, 2]
        type(sparse_direct_factor_t) :: factor
        real(dp) :: actual(2)
        integer :: status

        call sparse_direct_factor_csc(factor, 2, ptr, rows, &
            [1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp], status)
        call record_condition(status /= 0, "Singular matrix reports failure")
        call sparse_direct_solve_factored( &
            factor, [1.0_dp, 2.0_dp], actual, status)
        call record_condition(status /= 0, "Failed factors cannot be solved")
        call sparse_direct_factor_csc(factor, 2, ptr, rows, &
            [4.0_dp, 2.0_dp, 1.0_dp, 3.0_dp], status)
        call sparse_direct_solve_factored( &
            factor, [6.0_dp, 8.0_dp], actual, status)
        call record_condition(status == 0 .and. &
            maxval(abs(actual - [1.0_dp, 2.0_dp])) < 1.0e-13_dp, &
            "A singular factor attempt can be replaced by a valid matrix")
        call sparse_direct_free(factor)
    end subroutine test_failed_factor_recovery

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_sparse_direct_backend
