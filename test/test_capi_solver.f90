program test_capi_solver
    use iso_c_binding, only: c_double_complex, c_int
    use check, only: check_condition, check_summary
    implicit none

    integer(c_int), parameter :: n = 2_c_int
    integer(c_int), parameter :: nnz = 4_c_int
    integer(c_int), parameter :: col_ptr(3) = [0_c_int, 2_c_int, 4_c_int]
    integer(c_int), parameter :: row_ind(4) = [ &
        0_c_int, 1_c_int, 0_c_int, 1_c_int]
    complex(c_double_complex), parameter :: values(4) = [ &
        cmplx(3.0, 1.0, c_double_complex), &
        cmplx(1.0, 0.0, c_double_complex), &
        cmplx(1.0, 0.0, c_double_complex), &
        cmplx(2.0, -1.0, c_double_complex)]
    complex(c_double_complex), parameter :: expected_a(2) = [ &
        cmplx(1.0, -1.0, c_double_complex), &
        cmplx(2.0, 0.5, c_double_complex)]
    complex(c_double_complex), parameter :: expected_b(2) = [ &
        cmplx(-2.0, 0.25, c_double_complex), &
        cmplx(0.5, 3.0, c_double_complex)]
    complex(c_double_complex) :: rhs(2), solution(2)
    integer(c_int) :: handle, status
    logical :: all_passed

    interface
        subroutine fortfem_complex_factor_csc( &
                n, nnz, col_ptr, row_ind, values, handle, status) &
                bind(c, name="fortfem_complex_factor_csc")
            import c_double_complex, c_int
            integer(c_int), value :: n, nnz
            integer(c_int), intent(in) :: col_ptr(*), row_ind(*)
            complex(c_double_complex), intent(in) :: values(*)
            integer(c_int), intent(out) :: handle, status
        end subroutine fortfem_complex_factor_csc

        subroutine fortfem_complex_solve(handle, n, rhs, solution, status) &
                bind(c, name="fortfem_complex_solve")
            import c_double_complex, c_int
            integer(c_int), value :: handle, n
            complex(c_double_complex), intent(in) :: rhs(*)
            complex(c_double_complex), intent(out) :: solution(*)
            integer(c_int), intent(out) :: status
        end subroutine fortfem_complex_solve

        subroutine fortfem_factor_free(handle, status) &
                bind(c, name="fortfem_factor_free")
            import c_int
            integer(c_int), value :: handle
            integer(c_int), intent(out) :: status
        end subroutine fortfem_factor_free
    end interface

    all_passed = .true.
    call fortfem_complex_factor_csc( &
        n, nnz, col_ptr, row_ind, values, handle, status)
    call record_condition(status == 0_c_int .and. handle > 0_c_int, &
        "C complex factorization returns a live handle")

    rhs = matmul(reshape(values, [2, 2]), expected_a)
    call fortfem_complex_solve(handle, n, rhs, solution, status)
    call record_condition(status == 0_c_int .and. &
        maxval(abs(solution - expected_a)) < 1.0e-13, &
        "C factor handle solves the first exact right-hand side")

    rhs = matmul(reshape(values, [2, 2]), expected_b)
    call fortfem_complex_solve(handle, n, rhs, solution, status)
    call record_condition(status == 0_c_int .and. &
        maxval(abs(solution - expected_b)) < 1.0e-13, &
        "C factor handle reuses factors for a second right-hand side")

    call fortfem_factor_free(handle, status)
    call record_condition(status == 0_c_int, &
        "C factor handle releases its resources")
    call fortfem_complex_solve(handle, n, rhs, solution, status)
    call record_condition(status /= 0_c_int, &
        "C solver rejects a released factor handle")

    call check_summary("C sparse solver API")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_capi_solver
