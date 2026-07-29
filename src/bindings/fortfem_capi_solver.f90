module fortfem_capi_solver
    use iso_c_binding, only: c_double_complex, c_int
    use fortfem_kinds, only: dp
    use fortfem_sparse_direct, only: sparse_direct_factor_t, &
        sparse_direct_factor_csc, sparse_direct_free, &
        sparse_direct_solve_factored
    implicit none
    private

    integer, parameter :: max_factors = 32
    type(sparse_direct_factor_t), save :: factors(max_factors)
    integer, save :: factor_sizes(max_factors) = 0
    logical, save :: factor_occupied(max_factors) = .false.

    public :: fortfem_complex_factor_csc
    public :: fortfem_complex_solve
    public :: fortfem_factor_free

contains

    subroutine fortfem_complex_factor_csc( &
            n, nnz, col_ptr, row_ind, values, handle, status) &
            bind(c, name="fortfem_complex_factor_csc")
        integer(c_int), value :: n, nnz
        integer(c_int), intent(in) :: col_ptr(*), row_ind(*)
        complex(c_double_complex), intent(in) :: values(*)
        integer(c_int), intent(out) :: handle, status

        integer, allocatable :: fortran_col_ptr(:), fortran_row_ind(:)
        complex(dp), allocatable :: fortran_values(:)
        integer :: i, slot, solver_status

        handle = 0_c_int
        status = -1_c_int
        if (.not. valid_csc_input(n, nnz, col_ptr, row_ind)) return

        slot = 0
        do i = 1, max_factors
            if (.not. factor_occupied(i)) then
                slot = i
                exit
            end if
        end do
        if (slot == 0) then
            status = -3_c_int
            return
        end if

        allocate(fortran_col_ptr(int(n) + 1))
        allocate(fortran_row_ind(int(nnz)), fortran_values(int(nnz)))
        do i = 1, int(n) + 1
            fortran_col_ptr(i) = int(col_ptr(i)) + 1
        end do
        do i = 1, int(nnz)
            fortran_row_ind(i) = int(row_ind(i)) + 1
            fortran_values(i) = cmplx( &
                real(values(i), dp), aimag(values(i)), dp)
        end do

        call sparse_direct_factor_csc(factors(slot), int(n), &
            fortran_col_ptr, fortran_row_ind, fortran_values, solver_status)
        status = int(solver_status, c_int)
        if (solver_status /= 0) then
            call sparse_direct_free(factors(slot))
            return
        end if

        factor_occupied(slot) = .true.
        factor_sizes(slot) = int(n)
        handle = int(slot, c_int)
    end subroutine fortfem_complex_factor_csc

    subroutine fortfem_complex_solve(handle, n, rhs, solution, status) &
            bind(c, name="fortfem_complex_solve")
        integer(c_int), value :: handle, n
        complex(c_double_complex), intent(in) :: rhs(*)
        complex(c_double_complex), intent(out) :: solution(*)
        integer(c_int), intent(out) :: status

        complex(dp), allocatable :: fortran_rhs(:), fortran_solution(:)
        integer :: i, slot, solver_status

        status = -3_c_int
        slot = int(handle)
        if (slot < 1 .or. slot > max_factors) return
        if (.not. factor_occupied(slot)) return
        if (int(n) /= factor_sizes(slot)) then
            status = -1_c_int
            return
        end if

        allocate(fortran_rhs(int(n)), fortran_solution(int(n)))
        do i = 1, int(n)
            fortran_rhs(i) = cmplx(real(rhs(i), dp), aimag(rhs(i)), dp)
        end do
        call sparse_direct_solve_factored( &
            factors(slot), fortran_rhs, fortran_solution, solver_status)
        status = int(solver_status, c_int)
        if (solver_status /= 0) return

        do i = 1, int(n)
            solution(i) = cmplx( &
                real(fortran_solution(i), c_double_complex), &
                aimag(fortran_solution(i)), c_double_complex)
        end do
    end subroutine fortfem_complex_solve

    subroutine fortfem_factor_free(handle, status) &
            bind(c, name="fortfem_factor_free")
        integer(c_int), value :: handle
        integer(c_int), intent(out) :: status

        integer :: slot

        status = -3_c_int
        slot = int(handle)
        if (slot < 1 .or. slot > max_factors) return
        if (.not. factor_occupied(slot)) return

        call sparse_direct_free(factors(slot))
        factor_occupied(slot) = .false.
        factor_sizes(slot) = 0
        status = 0_c_int
    end subroutine fortfem_factor_free

    logical function valid_csc_input(n, nnz, col_ptr, row_ind) result(valid)
        integer(c_int), intent(in) :: n, nnz
        integer(c_int), intent(in) :: col_ptr(*), row_ind(*)

        integer :: i

        valid = .false.
        if (n < 1_c_int .or. nnz < 1_c_int) return
        if (col_ptr(1) /= 0_c_int .or. col_ptr(int(n) + 1) /= nnz) return
        do i = 1, int(n)
            if (col_ptr(i + 1) < col_ptr(i)) return
        end do
        do i = 1, int(nnz)
            if (row_ind(i) < 0_c_int .or. row_ind(i) >= n) return
        end do
        valid = .true.
    end function valid_csc_input

end module fortfem_capi_solver
