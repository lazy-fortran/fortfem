program test_sparse_direct_backend
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_solve_csc
    use check, only: check_condition, check_summary
    implicit none

    logical :: all_passed

    all_passed = .true.
    call test_real_solve()
    call test_complex_solve()
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

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_sparse_direct_backend
