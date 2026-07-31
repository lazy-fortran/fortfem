module fortfem_sparse_direct
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_t, csc_z_t, fortsparse_status_t, &
        sparse_direct_factor_t => sparse_solver_t, sparse_destroy, &
        csc_conjugate_transpose, csc_transpose, sparse_factor, sparse_solve, &
        sparse_solve_jvp, sparse_solve_once, sparse_solve_vjp
    implicit none
    private

    public :: sparse_direct_solve_csc
    public :: sparse_direct_solve_constrained
    public :: sparse_direct_solve_constrained_jvp
    public :: sparse_direct_solve_constrained_vjp
    public :: sparse_direct_solve_zero_constrained
    public :: sparse_direct_factor_t
    public :: sparse_direct_factor_csc
    public :: sparse_direct_factor_adjoint_csc
    public :: sparse_direct_factor_transpose_csc
    public :: sparse_direct_solve_factored
    public :: sparse_direct_solve_factored_jvp
    public :: sparse_direct_solve_factored_vjp
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

    interface sparse_direct_solve_factored_jvp
        module procedure sparse_direct_solve_factored_jvp_real
        module procedure sparse_direct_solve_factored_jvp_complex
    end interface sparse_direct_solve_factored_jvp

    interface sparse_direct_solve_factored_vjp
        module procedure sparse_direct_solve_factored_vjp_real
        module procedure sparse_direct_solve_factored_vjp_complex
    end interface sparse_direct_solve_factored_vjp

    interface validate_csc_dimensions
        module procedure validate_csc_dimensions_real
        module procedure validate_csc_dimensions_complex
    end interface validate_csc_dimensions

contains

    subroutine sparse_direct_solve_zero_constrained( &
            matrix, rhs, constrained, solution, status)
        type(csc_t), intent(in) :: matrix
        real(dp), intent(in) :: rhs(:)
        logical, intent(in) :: constrained(:)
        real(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        real(dp), allocatable :: constrained_values(:)

        status = -1
        allocate (constrained_values(size(constrained)))
        constrained_values = 0.0_dp
        call sparse_direct_solve_constrained( &
            matrix, rhs, constrained, constrained_values, solution, status)
    end subroutine sparse_direct_solve_zero_constrained

    subroutine sparse_direct_solve_constrained( &
            matrix, rhs, constrained, constrained_values, solution, status)
        type(csc_t), intent(in) :: matrix
        real(dp), intent(in) :: rhs(:)
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: constrained_values(:)
        real(dp), intent(out) :: solution(:)
        integer, intent(out) :: status

        integer, allocatable :: column_pointers(:), free_dofs(:)
        integer, allocatable :: free_index(:), row_indices(:)
        real(dp), allocatable :: reduced_rhs(:), reduced_solution(:)
        real(dp), allocatable :: values(:)
        integer :: column, entry, free_count, reduced_entry, row

        status = -1
        solution = 0.0_dp
        if (matrix%nrow /= matrix%ncol) return
        if (size(rhs) /= matrix%nrow) return
        if (size(constrained) /= matrix%nrow) return
        if (size(constrained_values) /= matrix%nrow) return
        if (size(solution) /= matrix%nrow) return
        solution = constrained_values
        free_count = count(.not. constrained)
        if (free_count == 0) then
            status = 0
            return
        end if
        allocate(free_dofs(free_count), free_index(matrix%nrow))
        free_count = 0
        do column = 1, matrix%ncol
            if (constrained(column)) cycle
            free_count = free_count + 1
            free_dofs(free_count) = column
        end do
        free_index = 0
        do column = 1, free_count
            free_index(free_dofs(column)) = column
        end do
        allocate(column_pointers(free_count + 1))
        reduced_entry = 0
        do column = 1, free_count
            column_pointers(column) = reduced_entry + 1
            do entry = matrix%col_ptr(free_dofs(column)), &
                    matrix%col_ptr(free_dofs(column) + 1) - 1
                if (.not. constrained(matrix%row_idx(entry))) then
                    reduced_entry = reduced_entry + 1
                end if
            end do
        end do
        column_pointers(free_count + 1) = reduced_entry + 1
        allocate(row_indices(reduced_entry), values(reduced_entry))
        reduced_entry = 0
        do column = 1, free_count
            do entry = matrix%col_ptr(free_dofs(column)), &
                    matrix%col_ptr(free_dofs(column) + 1) - 1
                row = matrix%row_idx(entry)
                if (constrained(row)) cycle
                reduced_entry = reduced_entry + 1
                row_indices(reduced_entry) = free_index(row)
                values(reduced_entry) = matrix%val(entry)
            end do
        end do
        allocate(reduced_rhs(free_count), reduced_solution(free_count))
        reduced_rhs = rhs(free_dofs)
        do column = 1, matrix%ncol
            if (.not. constrained(column)) cycle
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                if (constrained(row)) cycle
                reduced_rhs(free_index(row)) = &
                    reduced_rhs(free_index(row)) - &
                    matrix%val(entry)*constrained_values(column)
            end do
        end do
        call sparse_direct_solve_csc( &
            free_count, column_pointers, row_indices, values, reduced_rhs, &
            reduced_solution, status)
        if (status /= 0) return
        solution(free_dofs) = reduced_solution
    end subroutine sparse_direct_solve_constrained

    subroutine sparse_direct_solve_constrained_jvp( &
            matrix, rhs, constrained, constrained_values, matrix_values_dot, &
            rhs_dot, constrained_values_dot, solution_dot, status)
        type(csc_t), intent(in) :: matrix
        real(dp), intent(in) :: rhs(:), constrained_values(:)
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: matrix_values_dot(:), rhs_dot(:)
        real(dp), intent(in) :: constrained_values_dot(:)
        real(dp), intent(out) :: solution_dot(:)
        integer, intent(out) :: status

        type(csc_t) :: reduced
        type(sparse_direct_factor_t) :: factor
        integer, allocatable :: free_dofs(:), free_index(:)
        integer, allocatable :: original_entries(:)
        real(dp), allocatable :: reduced_rhs_dot(:), reduced_solution(:)
        real(dp), allocatable :: reduced_solution_dot(:), reduced_values_dot(:)
        real(dp), allocatable :: solution(:)
        integer :: column, entry, free_count, row

        status = -1
        solution_dot = 0.0_dp
        if (.not. valid_constrained_shapes( &
            matrix, rhs, constrained, constrained_values, solution_dot)) &
            return
        if (size(matrix_values_dot) /= matrix%nnz) return
        if (size(rhs_dot) /= size(rhs)) return
        if (size(constrained_values_dot) /= size(constrained_values)) return
        allocate(solution(matrix%nrow))
        call sparse_direct_solve_constrained( &
            matrix, rhs, constrained, constrained_values, solution, status)
        if (status /= 0) return
        solution_dot = constrained_values_dot
        call build_constrained_reduction( &
            matrix, constrained, reduced, free_dofs, free_index, &
            original_entries, status)
        if (status /= 0) return
        free_count = size(free_dofs)
        if (free_count == 0) then
            status = 0
            return
        end if
        allocate(reduced_solution(free_count))
        allocate(reduced_solution_dot(free_count))
        allocate(reduced_rhs_dot(free_count))
        allocate(reduced_values_dot(reduced%nnz))
        reduced_solution = solution(free_dofs)
        reduced_rhs_dot = rhs_dot(free_dofs)
        reduced_values_dot = matrix_values_dot(original_entries)
        do column = 1, matrix%ncol
            if (.not. constrained(column)) cycle
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                if (constrained(row)) cycle
                reduced_rhs_dot(free_index(row)) = &
                    reduced_rhs_dot(free_index(row)) - &
                    matrix_values_dot(entry)*constrained_values(column) - &
                    matrix%val(entry)*constrained_values_dot(column)
            end do
        end do
        call sparse_direct_factor_csc( &
            factor, free_count, reduced%col_ptr, reduced%row_idx, &
            reduced%val, status)
        if (status /= 0) return
        call sparse_direct_solve_factored_jvp( &
            factor, free_count, reduced%col_ptr, reduced%row_idx, &
            reduced_values_dot, reduced_solution, reduced_rhs_dot, &
            reduced_solution_dot, status)
        call sparse_direct_free(factor)
        if (status /= 0) return
        solution_dot(free_dofs) = reduced_solution_dot
    end subroutine sparse_direct_solve_constrained_jvp

    subroutine sparse_direct_solve_constrained_vjp( &
            matrix, rhs, constrained, constrained_values, solution, &
            solution_bar, matrix_values_bar, rhs_bar, &
            constrained_values_bar, status)
        type(csc_t), intent(in) :: matrix
        real(dp), intent(in) :: rhs(:), constrained_values(:), solution(:)
        logical, intent(in) :: constrained(:)
        real(dp), intent(in) :: solution_bar(:)
        real(dp), intent(out) :: matrix_values_bar(:), rhs_bar(:)
        real(dp), intent(out) :: constrained_values_bar(:)
        integer, intent(out) :: status

        type(csc_t) :: reduced
        type(sparse_direct_factor_t) :: transpose_factor
        integer, allocatable :: free_dofs(:), free_index(:)
        integer, allocatable :: original_entries(:)
        real(dp), allocatable :: reduced_rhs_bar(:), reduced_values_bar(:)
        integer :: column, entry, free_count, row

        status = -1
        matrix_values_bar = 0.0_dp
        rhs_bar = 0.0_dp
        constrained_values_bar = 0.0_dp
        if (.not. valid_constrained_shapes( &
            matrix, rhs, constrained, constrained_values, solution)) return
        if (size(solution_bar) /= size(solution)) return
        if (size(matrix_values_bar) /= matrix%nnz) return
        if (size(rhs_bar) /= size(rhs)) return
        if (size(constrained_values_bar) /= size(constrained_values)) return
        call build_constrained_reduction( &
            matrix, constrained, reduced, free_dofs, free_index, &
            original_entries, status)
        if (status /= 0) return
        constrained_values_bar = merge(solution_bar, 0.0_dp, constrained)
        free_count = size(free_dofs)
        if (free_count == 0) then
            status = 0
            return
        end if
        allocate(reduced_rhs_bar(free_count))
        allocate(reduced_values_bar(reduced%nnz))
        call sparse_direct_factor_transpose_csc( &
            transpose_factor, free_count, reduced%col_ptr, reduced%row_idx, &
            reduced%val, status)
        if (status /= 0) return
        call sparse_direct_solve_factored_vjp( &
            transpose_factor, free_count, reduced%col_ptr, reduced%row_idx, &
            solution(free_dofs), solution_bar(free_dofs), reduced_rhs_bar, &
            reduced_values_bar, status)
        call sparse_direct_free(transpose_factor)
        if (status /= 0) return
        rhs_bar(free_dofs) = reduced_rhs_bar
        matrix_values_bar(original_entries) = reduced_values_bar
        do column = 1, matrix%ncol
            if (.not. constrained(column)) cycle
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                if (constrained(row)) cycle
                matrix_values_bar(entry) = &
                    -reduced_rhs_bar(free_index(row))* &
                    constrained_values(column)
                constrained_values_bar(column) = &
                    constrained_values_bar(column) - &
                    matrix%val(entry)*reduced_rhs_bar(free_index(row))
            end do
        end do
    end subroutine sparse_direct_solve_constrained_vjp

    subroutine build_constrained_reduction( &
            matrix, constrained, reduced, free_dofs, free_index, &
            original_entries, status)
        type(csc_t), intent(in) :: matrix
        logical, intent(in) :: constrained(:)
        type(csc_t), intent(out) :: reduced
        integer, allocatable, intent(out) :: free_dofs(:), free_index(:)
        integer, allocatable, intent(out) :: original_entries(:)
        integer, intent(out) :: status

        integer :: column, entry, free_count, reduced_entry, row

        status = -1
        free_count = count(.not. constrained)
        allocate(free_dofs(free_count), free_index(matrix%nrow))
        free_count = 0
        do column = 1, matrix%ncol
            if (constrained(column)) cycle
            free_count = free_count + 1
            free_dofs(free_count) = column
        end do
        free_index = 0
        do column = 1, free_count
            free_index(free_dofs(column)) = column
        end do
        reduced_entry = 0
        do column = 1, free_count
            do entry = matrix%col_ptr(free_dofs(column)), &
                    matrix%col_ptr(free_dofs(column) + 1) - 1
                if (.not. constrained(matrix%row_idx(entry))) then
                    reduced_entry = reduced_entry + 1
                end if
            end do
        end do
        reduced%nrow = free_count
        reduced%ncol = free_count
        reduced%nnz = reduced_entry
        allocate(reduced%col_ptr(free_count + 1))
        allocate(reduced%row_idx(reduced_entry), reduced%val(reduced_entry))
        allocate(original_entries(reduced_entry))
        reduced_entry = 0
        do column = 1, free_count
            reduced%col_ptr(column) = reduced_entry + 1
            do entry = matrix%col_ptr(free_dofs(column)), &
                    matrix%col_ptr(free_dofs(column) + 1) - 1
                row = matrix%row_idx(entry)
                if (constrained(row)) cycle
                reduced_entry = reduced_entry + 1
                reduced%row_idx(reduced_entry) = free_index(row)
                reduced%val(reduced_entry) = matrix%val(entry)
                original_entries(reduced_entry) = entry
            end do
        end do
        reduced%col_ptr(free_count + 1) = reduced_entry + 1
        status = 0
    end subroutine build_constrained_reduction

    pure logical function valid_constrained_shapes( &
            matrix, rhs, constrained, constrained_values, solution) &
            result(valid)
        type(csc_t), intent(in) :: matrix
        real(dp), intent(in) :: rhs(:), constrained_values(:), solution(:)
        logical, intent(in) :: constrained(:)

        valid = matrix%nrow == matrix%ncol
        if (.not. valid) return
        valid = size(rhs) == matrix%nrow
        if (.not. valid) return
        valid = size(constrained) == matrix%nrow
        if (.not. valid) return
        valid = size(constrained_values) == matrix%nrow
        if (.not. valid) return
        valid = size(solution) == matrix%nrow
    end function valid_constrained_shapes

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

    subroutine sparse_direct_factor_transpose_csc( &
            factor, n, col_ptr, row_ind, values, status)
        type(sparse_direct_factor_t), intent(inout) :: factor
        integer, intent(in) :: n
        integer, intent(in) :: col_ptr(:), row_ind(:)
        real(dp), intent(in) :: values(:)
        integer, intent(out) :: status

        type(csc_t) :: matrix, transpose_matrix
        type(fortsparse_status_t) :: solver_status

        call initialize_csc_real( &
            n, col_ptr, row_ind, values, matrix, status)
        if (status /= 0) return
        call csc_transpose(matrix, transpose_matrix, solver_status)
        if (solver_status%code /= 0) then
            status = solver_status%code
            return
        end if
        call sparse_factor(factor, transpose_matrix, solver_status)
        status = solver_status%code
    end subroutine sparse_direct_factor_transpose_csc

    subroutine sparse_direct_factor_adjoint_csc( &
            factor, n, col_ptr, row_ind, values, status)
        type(sparse_direct_factor_t), intent(inout) :: factor
        integer, intent(in) :: n
        integer, intent(in) :: col_ptr(:), row_ind(:)
        complex(dp), intent(in) :: values(:)
        integer, intent(out) :: status

        type(csc_z_t) :: adjoint_matrix, matrix
        type(fortsparse_status_t) :: solver_status

        call initialize_csc_complex( &
            n, col_ptr, row_ind, values, matrix, status)
        if (status /= 0) return
        call csc_conjugate_transpose(matrix, adjoint_matrix, solver_status)
        if (solver_status%code /= 0) then
            status = solver_status%code
            return
        end if
        call sparse_factor(factor, adjoint_matrix, solver_status)
        status = solver_status%code
    end subroutine sparse_direct_factor_adjoint_csc

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

    subroutine sparse_direct_solve_factored_jvp_real( &
            factor, n, col_ptr, row_ind, values_dot, solution, rhs_dot, &
            solution_dot, status)
        type(sparse_direct_factor_t), intent(inout) :: factor
        integer, intent(in) :: n
        integer, intent(in) :: col_ptr(:), row_ind(:)
        real(dp), intent(in) :: values_dot(:)
        real(dp), target, contiguous, intent(in) :: solution(:), rhs_dot(:)
        real(dp), target, contiguous, intent(out) :: solution_dot(:)
        integer, intent(out) :: status

        type(csc_t) :: matrix_dot
        type(fortsparse_status_t) :: solver_status

        call initialize_csc_real( &
            n, col_ptr, row_ind, values_dot, matrix_dot, status)
        if (status /= 0) return
        call sparse_solve_jvp( &
            factor, matrix_dot, solution, rhs_dot, solution_dot, solver_status)
        status = solver_status%code
    end subroutine sparse_direct_solve_factored_jvp_real

    subroutine sparse_direct_solve_factored_vjp_real( &
            transpose_factor, n, col_ptr, row_ind, solution, solution_bar, &
            rhs_bar, values_bar, status)
        type(sparse_direct_factor_t), intent(inout) :: transpose_factor
        integer, intent(in) :: n
        integer, intent(in) :: col_ptr(:), row_ind(:)
        real(dp), intent(in) :: solution(:)
        real(dp), target, contiguous, intent(in) :: solution_bar(:)
        real(dp), target, contiguous, intent(out) :: rhs_bar(:)
        real(dp), intent(out) :: values_bar(:)
        integer, intent(out) :: status

        type(csc_t) :: matrix
        type(fortsparse_status_t) :: solver_status
        real(dp), allocatable :: zero_values(:)

        allocate (zero_values(size(values_bar)), source=0.0_dp)
        call initialize_csc_real( &
            n, col_ptr, row_ind, zero_values, matrix, status)
        if (status /= 0) return
        call sparse_solve_vjp( &
            transpose_factor, matrix, solution, solution_bar, rhs_bar, &
            values_bar, solver_status)
        status = solver_status%code
    end subroutine sparse_direct_solve_factored_vjp_real

    subroutine sparse_direct_solve_factored_jvp_complex( &
            factor, n, col_ptr, row_ind, values_dot, solution, rhs_dot, &
            solution_dot, status)
        type(sparse_direct_factor_t), intent(inout) :: factor
        integer, intent(in) :: n
        integer, intent(in) :: col_ptr(:), row_ind(:)
        complex(dp), intent(in) :: values_dot(:), solution(:), rhs_dot(:)
        complex(dp), intent(out) :: solution_dot(:)
        integer, intent(out) :: status

        type(csc_z_t) :: matrix_dot
        type(fortsparse_status_t) :: solver_status

        call initialize_csc_complex( &
            n, col_ptr, row_ind, values_dot, matrix_dot, status)
        if (status /= 0) return
        call sparse_solve_jvp( &
            factor, matrix_dot, solution, rhs_dot, solution_dot, solver_status)
        status = solver_status%code
    end subroutine sparse_direct_solve_factored_jvp_complex

    subroutine sparse_direct_solve_factored_vjp_complex( &
            adjoint_factor, n, col_ptr, row_ind, solution, solution_bar, &
            rhs_bar, values_bar, status)
        type(sparse_direct_factor_t), intent(inout) :: adjoint_factor
        integer, intent(in) :: n
        integer, intent(in) :: col_ptr(:), row_ind(:)
        complex(dp), intent(in) :: solution(:), solution_bar(:)
        complex(dp), intent(out) :: rhs_bar(:), values_bar(:)
        integer, intent(out) :: status

        type(csc_z_t) :: matrix
        type(fortsparse_status_t) :: solver_status
        complex(dp), allocatable :: zero_values(:)

        allocate (zero_values(size(values_bar)), &
            source=cmplx(0.0_dp, 0.0_dp, dp))
        call initialize_csc_complex( &
            n, col_ptr, row_ind, zero_values, matrix, status)
        if (status /= 0) return
        call sparse_solve_vjp( &
            adjoint_factor, matrix, solution, solution_bar, rhs_bar, &
            values_bar, solver_status)
        status = solver_status%code
    end subroutine sparse_direct_solve_factored_vjp_complex

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
