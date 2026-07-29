module fortfem_sparse_direct
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_t, csc_z_t, fortsparse_status_t, &
        sparse_direct_factor_t => sparse_solver_t, sparse_destroy, &
        sparse_factor, sparse_solve, sparse_solve_once
    implicit none
    private

    public :: sparse_direct_solve_csc
    public :: sparse_direct_factor_t
    public :: sparse_direct_factor_csc
    public :: sparse_direct_solve_factored
    public :: sparse_direct_free

    interface sparse_direct_solve_csc
        module procedure sparse_direct_solve_csc_real
        module procedure sparse_direct_solve_csc_complex
    end interface sparse_direct_solve_csc

    interface sparse_direct_factor_csc
        module procedure sparse_direct_factor_csc_real
        module procedure sparse_direct_factor_csc_complex
    end interface sparse_direct_factor_csc

    interface sparse_direct_solve_factored
        module procedure sparse_direct_solve_factored_real
        module procedure sparse_direct_solve_factored_complex
    end interface sparse_direct_solve_factored

    interface validate_csc_dimensions
        module procedure validate_csc_dimensions_real
        module procedure validate_csc_dimensions_complex
    end interface validate_csc_dimensions

contains

    subroutine sparse_direct_factor_csc_real( &
            factor, n, col_ptr, row_ind, values, status)
        type(sparse_direct_factor_t), intent(inout) :: factor
        integer, intent(in) :: n
        integer, intent(in) :: col_ptr(:), row_ind(:)
        real(dp), intent(in) :: values(:)
        integer, intent(out) :: status

        type(csc_t) :: matrix
        type(fortsparse_status_t) :: solver_status

        call initialize_csc_real( &
            n, col_ptr, row_ind, values, matrix, status)
        if (status /= 0) return
        call sparse_factor(factor, matrix, solver_status)
        status = solver_status%code
    end subroutine sparse_direct_factor_csc_real

    subroutine sparse_direct_factor_csc_complex( &
            factor, n, col_ptr, row_ind, values, status)
        type(sparse_direct_factor_t), intent(inout) :: factor
        integer, intent(in) :: n
        integer, intent(in) :: col_ptr(:), row_ind(:)
        complex(dp), intent(in) :: values(:)
        integer, intent(out) :: status

        type(csc_z_t) :: matrix
        type(fortsparse_status_t) :: solver_status

        call initialize_csc_complex( &
            n, col_ptr, row_ind, values, matrix, status)
        if (status /= 0) return
        call sparse_factor(factor, matrix, solver_status)
        status = solver_status%code
    end subroutine sparse_direct_factor_csc_complex

    subroutine sparse_direct_solve_factored_real( &
            factor, rhs, solution, status)
        type(sparse_direct_factor_t), intent(inout) :: factor
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        type(fortsparse_status_t) :: solver_status

        call sparse_solve(factor, rhs, solution, solver_status)
        status = solver_status%code
    end subroutine sparse_direct_solve_factored_real

    subroutine sparse_direct_solve_factored_complex( &
            factor, rhs, solution, status)
        type(sparse_direct_factor_t), intent(inout) :: factor
        complex(dp), intent(in) :: rhs(:)
        complex(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        type(fortsparse_status_t) :: solver_status

        call sparse_solve(factor, rhs, solution, solver_status)
        status = solver_status%code
    end subroutine sparse_direct_solve_factored_complex

    subroutine sparse_direct_free(factor)
        type(sparse_direct_factor_t), intent(inout) :: factor

        call sparse_destroy(factor)
    end subroutine sparse_direct_free

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

    subroutine initialize_csc_real( &
            n, col_ptr, row_ind, values, matrix, status)
        integer, intent(in) :: n
        integer, intent(in) :: col_ptr(:), row_ind(:)
        real(dp), intent(in) :: values(:)
        type(csc_t), intent(out) :: matrix
        integer, intent(out) :: status

        call validate_csc_matrix_dimensions( &
            n, col_ptr, row_ind, size(values), status)
        if (status /= 0) return
        matrix%nrow = n
        matrix%ncol = n
        matrix%nnz = size(values)
        matrix%col_ptr = col_ptr
        matrix%row_idx = row_ind
        matrix%val = values
    end subroutine initialize_csc_real

    subroutine initialize_csc_complex( &
            n, col_ptr, row_ind, values, matrix, status)
        integer, intent(in) :: n
        integer, intent(in) :: col_ptr(:), row_ind(:)
        complex(dp), intent(in) :: values(:)
        type(csc_z_t), intent(out) :: matrix
        integer, intent(out) :: status

        call validate_csc_matrix_dimensions( &
            n, col_ptr, row_ind, size(values), status)
        if (status /= 0) return
        matrix%nrow = n
        matrix%ncol = n
        matrix%nnz = size(values)
        matrix%col_ptr = col_ptr
        matrix%row_idx = row_ind
        matrix%val = values
    end subroutine initialize_csc_complex

    subroutine validate_csc_matrix_dimensions( &
            n, col_ptr, row_ind, nvalues, status)
        integer, intent(in) :: n, col_ptr(:), row_ind(:), nvalues
        integer, intent(out) :: status

        status = 0
        if (n < 1) status = -1
        if (size(col_ptr) /= n + 1) status = -1
        if (size(row_ind) /= nvalues) status = -1
    end subroutine validate_csc_matrix_dimensions

end module fortfem_sparse_direct
