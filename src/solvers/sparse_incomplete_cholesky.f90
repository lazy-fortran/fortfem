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
        integer :: max_fill_per_column = 0
        real(dp) :: drop_tolerance = 0.0_dp
        integer, allocatable :: col_ptr(:)
        integer, allocatable :: row_idx(:)
        real(dp), allocatable :: lower(:)
    end type sparse_incomplete_cholesky_factor_t

    public :: build_sparse_incomplete_cholesky
    public :: build_sparse_ichol
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

    subroutine build_sparse_ichol( &
            matrix, drop_tolerance, max_fill_per_column, factor, status)
        !! Build a drop- and fill-controlled incomplete Cholesky factor.
        !!
        !! The dense numeric phase is a deterministic reference construction;
        !! the retained lower factor and all apply paths remain sparse CSC.
        type(csc_t), intent(in) :: matrix
        real(dp), intent(in) :: drop_tolerance
        integer, intent(in) :: max_fill_per_column
        type(sparse_incomplete_cholesky_factor_t), intent(inout) :: factor
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: work(:, :)
        logical, allocatable :: keep(:, :)
        real(dp) :: scale, tolerance, value, diagonal
        integer :: dimension, column, entry, row, inner, lower_nnz, position

        call clear_factor(factor)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "controlled sparse ICHOL received an invalid matrix or policy")
        if (.not. csc_is_valid(matrix) .or. matrix%nrow /= matrix%ncol) return
        if (.not. ieee_is_finite(drop_tolerance) .or. &
            drop_tolerance < 0.0_dp .or. max_fill_per_column < 0) return
        if (matrix%nnz > 0) then
            if (.not. all(ieee_is_finite(matrix%val))) return
        end if
        dimension = matrix%nrow
        if (dimension < 1) return
        scale = 1.0_dp
        if (matrix%nnz > 0) scale = max(scale, maxval(abs(matrix%val)))
        tolerance = 128.0_dp*epsilon(1.0_dp)*scale*real(dimension, dp)
        allocate(work(dimension, dimension))
        work = 0.0_dp
        do column = 1, dimension
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                work(row, column) = matrix%val(entry)
            end do
        end do
        if (maxval(abs(work - transpose(work))) > tolerance) then
            deallocate(work)
            return
        end if

        ! Complete lower Cholesky numeric phase.
        do column = 1, dimension
            diagonal = work(column, column)
            do inner = 1, column - 1
                diagonal = diagonal - work(column, inner)**2
            end do
            if (.not. ieee_is_finite(diagonal) .or. diagonal <= tolerance) then
                deallocate(work)
                call status_set(status, FORTSPARSE_SINGULAR, &
                    "controlled sparse ICHOL encountered a non-positive pivot")
                return
            end if
            work(column, column) = sqrt(diagonal)
            do row = column + 1, dimension
                value = work(row, column)
                do inner = 1, column - 1
                    value = value - work(row, inner)*work(column, inner)
                end do
                work(row, column) = value/work(column, column)
                if (.not. ieee_is_finite(work(row, column))) then
                    deallocate(work)
                    call status_set(status, FORTSPARSE_SINGULAR, &
                        "controlled sparse ICHOL produced a non-finite entry")
                    return
                end if
            end do
        end do

        allocate(keep(dimension, dimension))
        keep = .false.
        do column = 1, dimension
            call select_lower_entries(work, column, column + 1, dimension, &
                drop_tolerance*scale, max_fill_per_column, keep)
            keep(column, column) = .true.
        end do
        lower_nnz = count(keep)
        factor%dimension = dimension
        factor%max_fill_per_column = max_fill_per_column
        factor%drop_tolerance = drop_tolerance
        allocate(factor%col_ptr(dimension + 1), factor%row_idx(lower_nnz), &
            factor%lower(lower_nnz))
        factor%col_ptr(1) = 1
        position = 0
        do column = 1, dimension
            do row = column, dimension
                if (.not. keep(row, column)) cycle
                position = position + 1
                factor%row_idx(position) = row
                factor%lower(position) = work(row, column)
            end do
            factor%col_ptr(column + 1) = position + 1
        end do
        deallocate(work, keep)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_sparse_ichol

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

    subroutine select_lower_entries( &
            work, column, first_row, last_row, threshold, max_fill, keep)
        real(dp), intent(in) :: work(:, :)
        integer, intent(in) :: column, first_row, last_row, max_fill
        real(dp), intent(in) :: threshold
        logical, intent(inout) :: keep(:, :)

        logical, allocatable :: selected(:)
        integer :: row, candidate_count, target_count, selected_count, best_row
        real(dp) :: best_value, value

        if (first_row > last_row .or. max_fill == 0) return
        candidate_count = 0
        do row = first_row, last_row
            if (abs(work(row, column)) > threshold) candidate_count = &
                candidate_count + 1
        end do
        if (candidate_count == 0) return
        target_count = min(max_fill, candidate_count)
        allocate(selected(size(work, 1)))
        selected = .false.
        if (target_count == candidate_count) then
            do row = first_row, last_row
                if (abs(work(row, column)) > threshold) keep(row, column) = .true.
            end do
            deallocate(selected)
            return
        end if
        do selected_count = 1, target_count
            best_row = 0
            best_value = -1.0_dp
            do row = first_row, last_row
                value = abs(work(row, column))
                if (.not. selected(row) .and. value > threshold .and. &
                    value > best_value) then
                    best_row = row
                    best_value = value
                end if
            end do
            if (best_row == 0) exit
            selected(best_row) = .true.
            keep(best_row, column) = .true.
        end do
        deallocate(selected)
    end subroutine select_lower_entries

    pure logical function valid_factor(factor) result(valid)
        type(sparse_incomplete_cholesky_factor_t), intent(in) :: factor
        integer :: column, entry

        valid = .false.
        if (factor%dimension < 1) return
        if (factor%max_fill_per_column < 0 .or. &
            .not. ieee_is_finite(factor%drop_tolerance) .or. &
            factor%drop_tolerance < 0.0_dp) return
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
        factor%max_fill_per_column = 0
        factor%drop_tolerance = 0.0_dp
    end subroutine clear_factor

end module fortfem_sparse_incomplete_cholesky
