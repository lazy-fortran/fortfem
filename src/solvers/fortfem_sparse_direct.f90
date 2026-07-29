module fortfem_sparse_direct
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_t, csc_z_t, fortsparse_status_t, &
        sparse_solve_once
    implicit none
    private

    public :: sparse_direct_solve_csc

    interface sparse_direct_solve_csc
        module procedure sparse_direct_solve_csc_real
        module procedure sparse_direct_solve_csc_complex
    end interface sparse_direct_solve_csc

    interface validate_csc_dimensions
        module procedure validate_csc_dimensions_real
        module procedure validate_csc_dimensions_complex
    end interface validate_csc_dimensions

contains

    subroutine sparse_direct_solve_csc_real( &
            n, col_ptr, row_ind, values, rhs, solution, status)
        integer, intent(in) :: n
        integer, intent(in) :: col_ptr(:), row_ind(:)
        real(dp), intent(in) :: values(:), rhs(:)
        real(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        type(csc_t) :: matrix
        type(fortsparse_status_t) :: solver_status

        call validate_csc_dimensions( &
            n, col_ptr, row_ind, size(values), rhs, solution, status)
        if (status /= 0) return

        matrix%nrow = n
        matrix%ncol = n
        matrix%nnz = size(values)
        matrix%col_ptr = col_ptr
        matrix%row_idx = row_ind
        matrix%val = values

        call sparse_solve_once(matrix, rhs, solution, solver_status)
        status = solver_status%code
    end subroutine sparse_direct_solve_csc_real

    subroutine sparse_direct_solve_csc_complex( &
            n, col_ptr, row_ind, values, rhs, solution, status)
        integer, intent(in) :: n
        integer, intent(in) :: col_ptr(:), row_ind(:)
        complex(dp), intent(in) :: values(:), rhs(:)
        complex(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        type(csc_z_t) :: matrix
        type(fortsparse_status_t) :: solver_status

        call validate_csc_dimensions( &
            n, col_ptr, row_ind, size(values), rhs, solution, status)
        if (status /= 0) return

        matrix%nrow = n
        matrix%ncol = n
        matrix%nnz = size(values)
        matrix%col_ptr = col_ptr
        matrix%row_idx = row_ind
        matrix%val = values

        call sparse_solve_once(matrix, rhs, solution, solver_status)
        status = solver_status%code
    end subroutine sparse_direct_solve_csc_complex

    subroutine validate_csc_dimensions_real( &
            n, col_ptr, row_ind, nvalues, rhs, solution, status)
        integer, intent(in) :: n, col_ptr(:), row_ind(:), nvalues
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        status = 0
        if (n < 1) status = -1
        if (size(col_ptr) /= n + 1) status = -1
        if (size(row_ind) /= nvalues) status = -1
        if (size(rhs) /= n) status = -1
        if (size(solution) /= n) status = -1
    end subroutine validate_csc_dimensions_real

    subroutine validate_csc_dimensions_complex( &
            n, col_ptr, row_ind, nvalues, rhs, solution, status)
        integer, intent(in) :: n, col_ptr(:), row_ind(:), nvalues
        complex(dp), intent(in) :: rhs(:)
        complex(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        status = 0
        if (n < 1) status = -1
        if (size(col_ptr) /= n + 1) status = -1
        if (size(row_ind) /= nvalues) status = -1
        if (size(rhs) /= n) status = -1
        if (size(solution) /= n) status = -1
    end subroutine validate_csc_dimensions_complex

end module fortfem_sparse_direct
