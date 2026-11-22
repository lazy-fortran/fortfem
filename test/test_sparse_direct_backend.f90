program test_sparse_direct_backend
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_factor_t, &
        sparse_direct_factor_csc, sparse_direct_factor_transpose_csc, &
        sparse_direct_factor_adjoint_csc, sparse_direct_free, &
        sparse_direct_solve_constrained, sparse_direct_solve_csc, &
        sparse_direct_solve_factored, sparse_direct_refactor_csc, &
        sparse_direct_factor_symmetric_csc, sparse_direct_inertia, &
        SPARSE_DIRECT_LOWER, SPARSE_DIRECT_UPPER
    use fortsparse, only: csc_t, FORTSPARSE_BACKEND_NATIVE, sparse_options_t, &
        native_ordering_available, FORTSPARSE_ORDER_AMD, FORTSPARSE_ORDER_ND, &
        FORTSPARSE_LU_SCALAR, FORTSPARSE_LU_SUPERNODAL
    use check, only: check_condition, check_summary
    implicit none

    logical :: all_passed

    all_passed = .true.
    call test_selection_and_limits()
    call test_real_solve()
    call test_nonzero_constraint_elimination()
    call test_complex_solve()
    call test_reused_complex_factor()
    call test_pivoted_operations_and_lifetime()
    call test_failed_factor_recovery()
    call test_real_refactor()
    call test_symmetric_one_triangle()
    call check_summary("Fortsparse direct backend")
    if (.not. all_passed) error stop 1

contains

    subroutine test_selection_and_limits()
        integer, parameter :: ptr(3) = [1, 3, 5]
        integer, parameter :: rows(4) = [1, 2, 1, 2]
        real(dp), parameter :: values(4) = [4.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]
        type(sparse_direct_factor_t) :: factor
        type(sparse_options_t) :: options
        real(dp) :: solution(2)
        integer :: status, expected_method

        expected_method = FORTSPARSE_LU_SCALAR
        if (native_ordering_available(FORTSPARSE_ORDER_AMD) .or. &
                native_ordering_available(FORTSPARSE_ORDER_ND)) then
            expected_method = FORTSPARSE_LU_SUPERNODAL
        end if
        factor%options%nthreads = 1
        call sparse_direct_factor_csc(factor, 2, ptr, rows, values, status)
        call record_condition(status == 0, "Default selected factor succeeds")
        call record_condition(factor%options%lu_method == expected_method, &
            "Default uses supernodal LU when graph ordering is available")
        call record_condition(factor%options%nthreads == 1, &
            "Default selection preserves the caller thread budget")
        call sparse_direct_free(factor)

        ! Explicit options retain the scalar, natural-order control.
        options%lu_method = FORTSPARSE_LU_SCALAR
        call sparse_direct_factor_csc( &
            factor, 2, ptr, rows, values, status, options)
        call sparse_direct_solve_factored( &
            factor, [6.0_dp, 8.0_dp], solution, status)
        call record_condition(status == 0 .and. &
            maxval(abs(solution - [1.0_dp, 2.0_dp])) < 1.0e-13_dp, &
            "Explicit scalar natural-order solve preserves the exact solution")
        call sparse_direct_free(factor)

        options%host_memory_limit_bytes = 1
        call sparse_direct_solve_csc(2, ptr, rows, values, &
            [6.0_dp, 8.0_dp], solution, status, options)
        call record_condition(status /= 0, &
            "One-shot options enforce an insufficient host payload cap")
        factor%options%host_memory_limit_bytes = 1
        call sparse_direct_factor_csc(factor, 2, ptr, rows, values, status)
        call record_condition(status /= 0, &
            "Default selection preserves an insufficient retained host cap")
        call sparse_direct_free(factor)
    end subroutine test_selection_and_limits

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

    subroutine test_real_refactor()
        !! Refactor with new values on the same pattern (a Newton loop):
        !! the solution equals that of a fresh factorisation, and the
        !! ordering and symbolic analysis are kept.
        integer, parameter :: m = 12, n = m*m
        integer, allocatable :: ptr(:), rows(:)
        real(dp), allocatable :: values(:), b(:), x1(:), x2(:)
        type(sparse_direct_factor_t) :: factor, fresh
        integer :: status, i, j, col, nz
        logical :: reused

        allocate (ptr(n + 1), rows(5*n), values(5*n), b(n), x1(n), x2(n))
        nz = 0
        ptr(1) = 1
        do j = 1, m
            do i = 1, m
                col = i + (j - 1)*m
                if (j > 1) call push(nz, rows, values, col - m, -1.0_dp)
                if (i > 1) call push(nz, rows, values, col - 1, -1.1_dp)
                call push(nz, rows, values, col, 4.5_dp)
                if (i < m) call push(nz, rows, values, col + 1, -1.1_dp)
                if (j < m) call push(nz, rows, values, col + m, -1.0_dp)
                ptr(col + 1) = nz + 1
            end do
        end do
        b = [(sin(0.1_dp*i), i=1, n)]
        call sparse_direct_factor_csc(factor, n, ptr, rows(:nz), values(:nz), status)
        call record_condition(status == 0, "Refactor test: first factor succeeds")
        values(:nz) = 2*values(:nz)
        do col = 1, n
            do i = ptr(col), ptr(col + 1) - 1
                if (rows(i) == col) values(i) = values(i) + 0.3_dp*col/n
            end do
        end do
        call sparse_direct_refactor_csc(factor, n, ptr, rows(:nz), values(:nz), status, reused)
        call record_condition(status == 0, "Refactor on the same pattern succeeds")
        call record_condition(reused, "Refactor keeps the ordering and symbolic analysis")
        call sparse_direct_solve_factored(factor, b, x1, status)
        call sparse_direct_factor_csc(fresh, n, ptr, rows(:nz), values(:nz), status)
        call sparse_direct_solve_factored(fresh, b, x2, status)
        call record_condition(maxval(abs(x1 - x2)) < 1.0e-12_dp*maxval(abs(x2)), &
            "Refactored solve equals a fresh factorisation")
        call sparse_direct_free(factor)
        call sparse_direct_free(fresh)
    end subroutine test_real_refactor

    subroutine test_symmetric_one_triangle()
        !! A symmetric saddle-point matrix [[H, B^T], [B, 0]] (zero block:
        !! 2x2 pivots) given as one triangle and factorised by LDL^T: the
        !! solution equals that of the full matrix by LU, the inertia is
        !! (mb, 0, nh), and a refactor with new values keeps the analysis.
        integer, parameter :: nh = 40, mb = 8, n = nh + mb
        real(dp) :: a(n, n), b(n), x_lu(n), x_up(n), x_low(n), x_new(n), x_ref(n)
        integer :: ptr(n + 1), rows(n*n), full_ptr(n + 1), full_rows(n*n)
        integer :: status, i, nz, full_nz, negative, zero, positive
        real(dp) :: values(n*n), full_values(n*n)
        type(sparse_direct_factor_t) :: lu, up, low, fresh
        logical :: reused, reliable

        call saddle_matrix(1.0_dp, a)
        b = [(cos(0.3_dp*i) + 0.5_dp, i=1, n)]
        call dense_to_csc(a, 0, full_ptr, full_rows, full_values, full_nz)
        call sparse_direct_factor_csc(lu, n, full_ptr, full_rows(:full_nz), &
            full_values(:full_nz), status)
        call sparse_direct_solve_factored(lu, b, x_lu, status)

        call dense_to_csc(a, 2, ptr, rows, values, nz)
        call sparse_direct_factor_symmetric_csc(up, n, ptr, rows(:nz), values(:nz), &
            status, SPARSE_DIRECT_UPPER)
        call record_condition(status == 0, "Symmetric upper-triangle factor succeeds")
        call sparse_direct_solve_factored(up, b, x_up, status)
        call record_condition(maxval(abs(x_up - x_lu)) < 1.0e-10_dp*maxval(abs(x_lu)), &
            "Upper triangle LDL^T solution equals the full-matrix LU solution")
        call sparse_direct_inertia(up, negative, zero, positive, status, reliable)
        call record_condition(status == 0 .and. negative == mb .and. zero == 0 &
            .and. positive == nh, "Inertia of the saddle-point matrix is (mb, 0, nh)")
        call record_condition(reliable, "Inertia is reliable")

        call dense_to_csc(a, 1, ptr, rows, values, nz)
        call sparse_direct_factor_symmetric_csc(low, n, ptr, rows(:nz), values(:nz), &
            status, SPARSE_DIRECT_LOWER)
        call sparse_direct_solve_factored(low, b, x_low, status)
        call record_condition(all(x_low == x_up), &
            "Lower and upper triangle give the same solution")

        call saddle_matrix(1.7_dp, a)
        call dense_to_csc(a, 2, ptr, rows, values, nz)
        call sparse_direct_refactor_csc(up, n, ptr, rows(:nz), values(:nz), status, reused)
        call record_condition(status == 0 .and. reused, &
            "Symmetric refactor keeps the ordering and symbolic analysis")
        call sparse_direct_solve_factored(up, b, x_new, status)
        call sparse_direct_inertia(up, negative, zero, positive, status)
        call record_condition(status == 0 .and. negative == mb .and. positive == nh, &
            "Inertia after the refactor")
        call dense_to_csc(a, 0, full_ptr, full_rows, full_values, full_nz)
        call sparse_direct_factor_csc(fresh, n, full_ptr, full_rows(:full_nz), &
            full_values(:full_nz), status)
        call sparse_direct_solve_factored(fresh, b, x_ref, status)
        call record_condition(maxval(abs(x_new - x_ref)) < 1.0e-10_dp*maxval(abs(x_ref)), &
            "Refactored symmetric solve equals the full-matrix LU solve")

        call sparse_direct_inertia(lu, negative, zero, positive, status)
        call record_condition(status /= 0, "LU factor reports no inertia")
        call sparse_direct_free(lu)
        call sparse_direct_free(up)
        call sparse_direct_free(low)
        call sparse_direct_free(fresh)
    end subroutine test_symmetric_one_triangle

    subroutine saddle_matrix(scale, a)
        integer, parameter :: nh = 40, mb = 8, n = nh + mb
        real(dp), intent(in) :: scale
        real(dp), intent(out) :: a(n, n)
        integer :: i, c
        a = 0
        do i = 1, nh
            a(i, i) = 2.5_dp*scale
            if (i < nh) then
                a(i, i + 1) = -1.0_dp
                a(i + 1, i) = -1.0_dp
            end if
        end do
        do c = 1, mb
            a(nh + c, 4*c) = 1.0_dp + 0.1_dp*scale*c
            a(4*c, nh + c) = a(nh + c, 4*c)
            a(nh + c, 4*c + 3) = -0.5_dp
            a(4*c + 3, nh + c) = -0.5_dp
        end do
    end subroutine saddle_matrix

    subroutine dense_to_csc(a, mode, ptr, rows, values, nz)
        !! mode 0: full matrix, 1: lower triangle, 2: upper triangle.
        real(dp), intent(in) :: a(:, :)
        integer, intent(in) :: mode
        integer, intent(out) :: ptr(:), rows(:), nz
        real(dp), intent(out) :: values(:)
        integer :: i, j
        nz = 0
        ptr(1) = 1
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                if (a(i, j) == 0.0_dp) cycle
                if (mode == 1 .and. i < j) cycle
                if (mode == 2 .and. i > j) cycle
                call push(nz, rows, values, i, a(i, j))
            end do
            ptr(j + 1) = nz + 1
        end do
    end subroutine dense_to_csc

    subroutine push(nz, rows, values, r, v)
        integer, intent(inout) :: nz
        integer, intent(inout) :: rows(:)
        real(dp), intent(inout) :: values(:)
        integer, intent(in) :: r
        real(dp), intent(in) :: v
        nz = nz + 1
        rows(nz) = r
        values(nz) = v
    end subroutine push

end program test_sparse_direct_backend
