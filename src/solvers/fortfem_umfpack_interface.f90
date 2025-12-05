module fortfem_umfpack_interface
    use fortfem_kinds, only: dp
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, &
                                           c_null_ptr, c_associated
    implicit none
    private

    integer(c_int), parameter :: umfpack_sys_A = 0_c_int

    public :: umfpack_available
    public :: umfpack_solve_csc

    interface
        function umfpack_di_symbolic(n_row, n_col, Ap, Ai, Ax, symbolic, &
                                     control, info) bind(C, &
                                                         name="umfpack_di_symbolic")
            import :: c_int, c_double, c_ptr
            integer(c_int), value :: n_row, n_col
            integer(c_int), intent(in) :: Ap(*)
            integer(c_int), intent(in) :: Ai(*)
            real(c_double), intent(in) :: Ax(*)
            type(c_ptr), intent(inout) :: symbolic
            type(c_ptr), value :: control
            type(c_ptr), value :: info
            integer(c_int) :: umfpack_di_symbolic
        end function umfpack_di_symbolic

        function umfpack_di_numeric(Ap, Ai, Ax, symbolic, numeric, control, &
                                    info) bind(C, &
                                               name="umfpack_di_numeric")
            import :: c_int, c_double, c_ptr
            integer(c_int), intent(in) :: Ap(*)
            integer(c_int), intent(in) :: Ai(*)
            real(c_double), intent(in) :: Ax(*)
            type(c_ptr), value :: symbolic
            type(c_ptr), intent(inout) :: numeric
            type(c_ptr), value :: control
            type(c_ptr), value :: info
            integer(c_int) :: umfpack_di_numeric
        end function umfpack_di_numeric

        function umfpack_di_solve(sys, Ap, Ai, Ax, x, b, numeric, control, &
                                  info) bind(C, name="umfpack_di_solve")
            import :: c_int, c_double, c_ptr
            integer(c_int), value :: sys
            integer(c_int), intent(in) :: Ap(*)
            integer(c_int), intent(in) :: Ai(*)
            real(c_double), intent(in) :: Ax(*)
            real(c_double), intent(inout) :: x(*)
            real(c_double), intent(in) :: b(*)
            type(c_ptr), value :: numeric
            type(c_ptr), value :: control
            type(c_ptr), value :: info
            integer(c_int) :: umfpack_di_solve
        end function umfpack_di_solve

        subroutine umfpack_di_free_symbolic(symbolic) bind(C, &
                                                           name="umfpack_di_free_symbolic")
            import :: c_ptr
            type(c_ptr), intent(inout) :: symbolic
        end subroutine umfpack_di_free_symbolic

        subroutine umfpack_di_free_numeric(numeric) bind(C, &
                                                         name="umfpack_di_free_numeric")
            import :: c_ptr
            type(c_ptr), intent(inout) :: numeric
        end subroutine umfpack_di_free_numeric
    end interface

contains

    function umfpack_available() result(available)
        logical :: available

        available = .true.
    end function umfpack_available

    subroutine check_umfpack_kinds()
        logical, parameter :: dp_matches_c_double = (dp == c_double)

        if (.not. dp_matches_c_double) then
            error stop "UMFPACK interface requires dp == c_double"
        end if
    end subroutine check_umfpack_kinds

    subroutine build_umfpack_arrays(n, col_ptr, row_ind, values_csc, b, &
                                    Ap_c, Ai_c, Ax_c, b_c, x_c, n_c)
        integer, intent(in) :: n
        integer, intent(in) :: col_ptr(:)
        integer, intent(in) :: row_ind(:)
        real(dp), intent(in) :: values_csc(:)
        real(dp), intent(in) :: b(:)
        integer(c_int), allocatable, intent(out) :: Ap_c(:)
        integer(c_int), allocatable, intent(out) :: Ai_c(:)
        real(c_double), allocatable, intent(out) :: Ax_c(:)
        real(c_double), allocatable, intent(out) :: b_c(:)
        real(c_double), allocatable, intent(out) :: x_c(:)
        integer(c_int), intent(out) :: n_c

        integer :: nnz, j

        nnz = size(values_csc)
        n_c = int(n, kind=c_int)

        allocate (Ap_c(size(col_ptr)), Ai_c(nnz))
        allocate (Ax_c(nnz), b_c(n), x_c(n))

        Ap_c(1) = 0_c_int
        do j = 2, size(col_ptr)
            Ap_c(j) = int(col_ptr(j) - 1, kind=c_int)
        end do

        do j = 1, nnz
            Ai_c(j) = int(row_ind(j) - 1, kind=c_int)
            Ax_c(j) = real(values_csc(j), kind=c_double)
        end do

        do j = 1, n
            b_c(j) = real(b(j), kind=c_double)
            x_c(j) = 0.0_c_double
        end do
    end subroutine build_umfpack_arrays

    subroutine umfpack_factor_and_solve(n_c, Ap_c, Ai_c, Ax_c, b_c, x_c, &
                                        status)
        integer(c_int), intent(in) :: n_c
        integer(c_int), intent(in) :: Ap_c(:)
        integer(c_int), intent(in) :: Ai_c(:)
        real(c_double), intent(in) :: Ax_c(:)
        real(c_double), intent(in) :: b_c(:)
        real(c_double), intent(inout) :: x_c(:)
        integer, intent(out) :: status

        type(c_ptr) :: symbolic, numeric
        integer(c_int) :: stat_sym, stat_num, stat_solve

        symbolic = c_null_ptr
        numeric = c_null_ptr

        stat_sym = umfpack_di_symbolic(n_c, n_c, Ap_c, Ai_c, Ax_c, symbolic, &
                                       c_null_ptr, c_null_ptr)

        if (stat_sym /= 0_c_int) then
            status = stat_sym
            if (c_associated(symbolic)) then
                call umfpack_di_free_symbolic(symbolic)
            end if
            return
        end if

        stat_num = umfpack_di_numeric(Ap_c, Ai_c, Ax_c, symbolic, numeric, &
                                      c_null_ptr, c_null_ptr)

        call umfpack_di_free_symbolic(symbolic)

        if (stat_num /= 0_c_int) then
            status = stat_num
            if (c_associated(numeric)) then
                call umfpack_di_free_numeric(numeric)
            end if
            return
        end if

        stat_solve = umfpack_di_solve(umfpack_sys_A, Ap_c, Ai_c, Ax_c, x_c, &
                                      b_c, numeric, c_null_ptr, c_null_ptr)

        call umfpack_di_free_numeric(numeric)

        status = stat_solve
    end subroutine umfpack_factor_and_solve

    subroutine umfpack_solve_csc(n, col_ptr, row_ind, values_csc, b, x, &
                                 status)
        integer, intent(in) :: n
        integer, intent(in) :: col_ptr(:)
        integer, intent(in) :: row_ind(:)
        real(dp), intent(in) :: values_csc(:)
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: status

        integer :: j
        integer(c_int) :: n_c
        integer(c_int), allocatable :: Ap_c(:), Ai_c(:)
        real(c_double), allocatable :: Ax_c(:), b_c(:), x_c(:)
        integer :: status_local

        call check_umfpack_kinds()

        call build_umfpack_arrays(n, col_ptr, row_ind, values_csc, b, Ap_c, &
                                  Ai_c, Ax_c, b_c, x_c, n_c)

        call umfpack_factor_and_solve(n_c, Ap_c, Ai_c, Ax_c, b_c, x_c, &
                                      status_local)

        do j = 1, n
            x(j) = real(x_c(j), kind=dp)
        end do

        status = status_local

        deallocate (Ap_c, Ai_c, Ax_c, b_c, x_c)
    end subroutine umfpack_solve_csc

end module fortfem_umfpack_interface
