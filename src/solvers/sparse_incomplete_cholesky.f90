module fortfem_sparse_incomplete_cholesky
    !! Sparse incomplete-Cholesky(0) factorization and fixed-factor apply.
    !!
    !! The lower CSC pattern is inherited from the lower triangle of a
    !! symmetric positive-definite input matrix; no fill is introduced.  The
    !! factor is a preconditioner primitive and is deliberately separate from
    !! the structure-preserving space and tree--cotree gauge.  Its fixed-factor
    !! JVP/VJP actions are useful when a solver differentiates only through the
    !! right-hand side while keeping preconditioner construction inactive.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_is_valid, csc_t, fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, FORTSPARSE_SINGULAR
    implicit none
    private

    type, public :: sparse_incomplete_cholesky_factor_t
        private
        integer :: dimension = 0
        integer, allocatable :: col_ptr(:)
        integer, allocatable :: row_idx(:)
        real(dp), allocatable :: lower(:)
    end type sparse_incomplete_cholesky_factor_t

    public :: build_sparse_incomplete_cholesky
    public :: apply_sparse_incomplete_cholesky
    public :: apply_sparse_incomplete_cholesky_jvp
    public :: apply_sparse_incomplete_cholesky_vjp

contains

    subroutine build_sparse_incomplete_cholesky(matrix, factor, status)
        !! Build an IC(0) factor with the input's lower-triangle sparsity.
        type(csc_t), intent(in) :: matrix
        type(sparse_incomplete_cholesky_factor_t), intent(inout) :: factor
        type(fortsparse_status_t), intent(out) :: status

        integer :: dimension, column, entry, lower_nnz, position, diagonal
        integer :: row, inner, left, right
        real(dp) :: value, tolerance, scale

        call clear_factor(factor)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "sparse IC(0) received an invalid CSC matrix")
        if (.not. csc_is_valid(matrix) .or. matrix%nrow /= matrix%ncol) return
        if (matrix%nnz > 0 .and. .not. all(ieee_is_finite(matrix%val))) return
        dimension = matrix%nrow
        scale = max(1.0_dp, maxval(abs(matrix%val)))
        tolerance = 128.0_dp*epsilon(1.0_dp)*scale*real(dimension, dp)
        do column = 1, dimension
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                if (abs(matrix%val(entry) - csc_entry(matrix, column, row)) > &
                    tolerance) return
            end do
        end do
        lower_nnz = 0
        do column = 1, dimension
            diagonal = 0
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                if (row >= column) then
                    lower_nnz = lower_nnz + 1
                    if (row == column) diagonal = 1
                end if
            end do
            if (diagonal == 0) return
        end do
        factor%dimension = dimension
        allocate(factor%col_ptr(dimension + 1), factor%row_idx(lower_nnz), &
            factor%lower(lower_nnz))
        factor%col_ptr(1) = 1
        position = 0
        do column = 1, dimension
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                if (row < column) cycle
                position = position + 1
                factor%row_idx(position) = row
                factor%lower(position) = 0.0_dp
            end do
            factor%col_ptr(column + 1) = position + 1
        end do
        do column = 1, dimension
            diagonal = lower_position(factor, column, column)
            value = matrix_value(matrix, column, column)
            do inner = 1, column - 1
                left = lower_position(factor, column, inner)
                if (left > 0) value = value - factor%lower(left)**2
            end do
            if (value <= tolerance .or. .not. ieee_is_finite(value)) then
                call clear_factor(factor)
                call status_set(status, FORTSPARSE_SINGULAR, &
                    "sparse IC(0) encountered a non-positive pivot")
                return
            end if
            factor%lower(diagonal) = sqrt(value)
            do entry = factor%col_ptr(column), factor%col_ptr(column + 1) - 1
                row = factor%row_idx(entry)
                if (row <= column) cycle
                value = matrix_value(matrix, row, column)
                do inner = 1, column - 1
                    left = lower_position(factor, row, inner)
                    right = lower_position(factor, column, inner)
                    if (left > 0 .and. right > 0) value = value - &
                        factor%lower(left)*factor%lower(right)
                end do
                factor%lower(entry) = value/factor%lower(diagonal)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_sparse_incomplete_cholesky

    subroutine apply_sparse_incomplete_cholesky( &
            factor, right_hand_side, solution, status)
        !! Apply `(L L^T)^(-1)` with the fixed sparse factor.
        type(sparse_incomplete_cholesky_factor_t), intent(in) :: factor
        real(dp), intent(in) :: right_hand_side(:)
        real(dp), allocatable, intent(out) :: solution(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: dimension, column, entry, diagonal
        real(dp), allocatable :: intermediate(:)

        if (allocated(solution)) deallocate(solution)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "sparse IC(0) apply received an invalid factor")
        if (.not. valid_factor(factor) .or. size(right_hand_side) /= &
            factor%dimension .or. (size(right_hand_side) > 0 .and. &
            .not. all(ieee_is_finite(right_hand_side)))) then
            allocate(solution(0))
            return
        end if
        dimension = factor%dimension
        allocate(intermediate(dimension), solution(dimension))
        intermediate = right_hand_side
        do column = 1, dimension
            diagonal = lower_position(factor, column, column)
            if (factor%lower(diagonal) <= 0.0_dp) then
                deallocate(intermediate, solution)
                allocate(solution(0))
                call status_set(status, FORTSPARSE_SINGULAR, &
                    "sparse IC(0) apply encountered a zero diagonal")
                return
            end if
            intermediate(column) = intermediate(column)/factor%lower(diagonal)
            do entry = factor%col_ptr(column), factor%col_ptr(column + 1) - 1
                if (factor%row_idx(entry) > column) intermediate(factor%row_idx(entry)) = &
                    intermediate(factor%row_idx(entry)) - &
                    factor%lower(entry)*intermediate(column)
            end do
        end do
        solution = intermediate
        do column = dimension, 1, -1
            do entry = factor%col_ptr(column), factor%col_ptr(column + 1) - 1
                if (factor%row_idx(entry) > column) solution(column) = &
                    solution(column) - factor%lower(entry)*solution(factor%row_idx(entry))
            end do
            diagonal = lower_position(factor, column, column)
            solution(column) = solution(column)/factor%lower(diagonal)
        end do
        deallocate(intermediate)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_sparse_incomplete_cholesky

    subroutine apply_sparse_incomplete_cholesky_jvp( &
            factor, right_hand_side_dot, solution_dot, status)
        !! Fixed-factor JVP; factor construction is an inactive solver branch.
        type(sparse_incomplete_cholesky_factor_t), intent(in) :: factor
        real(dp), intent(in) :: right_hand_side_dot(:)
        real(dp), allocatable, intent(out) :: solution_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        call apply_sparse_incomplete_cholesky( &
            factor, right_hand_side_dot, solution_dot, status)
    end subroutine apply_sparse_incomplete_cholesky_jvp

    subroutine apply_sparse_incomplete_cholesky_vjp( &
            factor, solution_bar, right_hand_side_bar, status)
        !! Fixed-factor real VJP; IC(0) is symmetric by construction.
        type(sparse_incomplete_cholesky_factor_t), intent(in) :: factor
        real(dp), intent(in) :: solution_bar(:)
        real(dp), allocatable, intent(out) :: right_hand_side_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        call apply_sparse_incomplete_cholesky( &
            factor, solution_bar, right_hand_side_bar, status)
    end subroutine apply_sparse_incomplete_cholesky_vjp

    pure logical function valid_factor(factor) result(valid)
        type(sparse_incomplete_cholesky_factor_t), intent(in) :: factor
        integer :: column, entry

        valid = .false.
        if (factor%dimension < 1) return
        if (.not. allocated(factor%col_ptr) .or. &
            .not. allocated(factor%row_idx) .or. .not. allocated(factor%lower)) return
        if (size(factor%col_ptr) /= factor%dimension + 1 .or. &
            size(factor%row_idx) /= size(factor%lower)) return
        if (factor%col_ptr(1) /= 1 .or. factor%col_ptr(factor%dimension + 1) /= &
            size(factor%lower) + 1) return
        do column = 1, factor%dimension
            if (factor%col_ptr(column + 1) < factor%col_ptr(column)) return
            if (lower_position(factor, column, column) == 0) return
            do entry = factor%col_ptr(column), factor%col_ptr(column + 1) - 1
                if (factor%row_idx(entry) < column .or. &
                    factor%row_idx(entry) > factor%dimension) return
                if (.not. ieee_is_finite(factor%lower(entry))) return
            end do
        end do
        valid = .true.
    end function valid_factor

    pure integer function lower_position(factor, row, column) result(position)
        type(sparse_incomplete_cholesky_factor_t), intent(in) :: factor
        integer, intent(in) :: row, column
        integer :: entry

        position = 0
        if (.not. allocated(factor%col_ptr)) return
        if (column < 1 .or. column > factor%dimension) return
        do entry = factor%col_ptr(column), factor%col_ptr(column + 1) - 1
            if (factor%row_idx(entry) == row) then
                position = entry
                return
            end if
        end do
    end function lower_position

    pure real(dp) function matrix_value(matrix, row, column) result(value)
        type(csc_t), intent(in) :: matrix
        integer, intent(in) :: row, column
        integer :: entry

        value = 0.0_dp
        do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
            if (matrix%row_idx(entry) == row) then
                value = matrix%val(entry)
                return
            end if
        end do
    end function matrix_value

    pure real(dp) function csc_entry(matrix, row, column) result(value)
        type(csc_t), intent(in) :: matrix
        integer, intent(in) :: row, column
        integer :: entry

        value = 0.0_dp
        do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
            if (matrix%row_idx(entry) == row) then
                value = matrix%val(entry)
                return
            end if
        end do
    end function csc_entry

    subroutine clear_factor(factor)
        type(sparse_incomplete_cholesky_factor_t), intent(inout) :: factor

        if (allocated(factor%col_ptr)) deallocate(factor%col_ptr)
        if (allocated(factor%row_idx)) deallocate(factor%row_idx)
        if (allocated(factor%lower)) deallocate(factor%lower)
        factor%dimension = 0
    end subroutine clear_factor

end module fortfem_sparse_incomplete_cholesky
