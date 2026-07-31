module fortfem_sparse_incomplete_lu
    !! Sparse incomplete-LU(0) factorization and fixed-factor apply.
    !!
    !! The strict lower and upper CSC patterns are inherited from the input;
    !! no fill is introduced.  This is the explicit nonsymmetric counterpart
    !! to sparse IC(0), intended for nonsymmetric or response blocks for which
    !! an SPD Cholesky preconditioner is not valid.  Factor construction is an
    !! inactive solver branch; the right-hand-side JVP/VJP actions retain a
    !! fixed pattern and fixed factors.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_is_valid, csc_t, fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, FORTSPARSE_SINGULAR
    implicit none
    private

    type, public :: sparse_incomplete_lu_factor_t
        private
        integer :: dimension = 0
        integer, allocatable :: lower_col_ptr(:), lower_row_idx(:)
        real(dp), allocatable :: lower(:)
        integer, allocatable :: upper_col_ptr(:), upper_row_idx(:)
        real(dp), allocatable :: upper(:)
    end type sparse_incomplete_lu_factor_t

    public :: build_sparse_incomplete_lu
    public :: apply_sparse_incomplete_lu
    public :: apply_sparse_incomplete_lu_jvp
    public :: apply_sparse_incomplete_lu_vjp

contains

    subroutine build_sparse_incomplete_lu(matrix, factor, status)
        !! Build an ILU(0) factor with the input's fixed lower/upper pattern.
        type(csc_t), intent(in) :: matrix
        type(sparse_incomplete_lu_factor_t), intent(inout) :: factor
        type(fortsparse_status_t), intent(out) :: status

        integer :: dimension, column, entry, lower_nnz, upper_nnz
        integer :: position, row, diagonal, inner, left, right
        real(dp) :: value, tolerance, scale

        call clear_factor(factor)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "sparse ILU(0) received an invalid CSC matrix")
        if (.not. csc_is_valid(matrix) .or. matrix%nrow /= matrix%ncol) return
        if (matrix%nnz > 0 .and. .not. all(ieee_is_finite(matrix%val))) return
        dimension = matrix%nrow
        scale = max(1.0_dp, maxval(abs(matrix%val)))
        tolerance = 128.0_dp*epsilon(1.0_dp)*scale*real(dimension, dp)
        lower_nnz = 0
        upper_nnz = 0
        do column = 1, dimension
            diagonal = 0
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                if (row > column) then
                    lower_nnz = lower_nnz + 1
                else
                    upper_nnz = upper_nnz + 1
                    if (row == column) diagonal = 1
                end if
            end do
            if (diagonal == 0) return
        end do
        factor%dimension = dimension
        allocate(factor%lower_col_ptr(dimension + 1), &
            factor%lower_row_idx(lower_nnz), factor%lower(lower_nnz), &
            factor%upper_col_ptr(dimension + 1), &
            factor%upper_row_idx(upper_nnz), factor%upper(upper_nnz))
        factor%lower_col_ptr(1) = 1
        factor%upper_col_ptr(1) = 1
        lower_nnz = 0
        upper_nnz = 0
        do column = 1, dimension
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                if (row > column) then
                    lower_nnz = lower_nnz + 1
                    factor%lower_row_idx(lower_nnz) = row
                    factor%lower(lower_nnz) = 0.0_dp
                else
                    upper_nnz = upper_nnz + 1
                    factor%upper_row_idx(upper_nnz) = row
                    factor%upper(upper_nnz) = 0.0_dp
                end if
            end do
            factor%lower_col_ptr(column + 1) = lower_nnz + 1
            factor%upper_col_ptr(column + 1) = upper_nnz + 1
        end do

        do column = 1, dimension
            ! U(i,column), i <= column.  Rows are sorted by CSC input, so
            ! U(k,column) is available when row i is processed.
            do entry = factor%upper_col_ptr(column), &
                    factor%upper_col_ptr(column + 1) - 1
                row = factor%upper_row_idx(entry)
                if (row == column) cycle
                value = matrix_value(matrix, row, column)
                do inner = 1, row - 1
                    left = lower_position(factor, row, inner)
                    right = upper_position(factor, inner, column)
                    if (left > 0 .and. right > 0) value = value - &
                        factor%lower(left)*factor%upper(right)
                end do
                factor%upper(entry) = value
            end do
            right = upper_position(factor, column, column)
            value = matrix_value(matrix, column, column)
            do inner = 1, column - 1
                left = lower_position(factor, column, inner)
                right = upper_position(factor, inner, column)
                if (left > 0 .and. right > 0) value = value - &
                    factor%lower(left)*factor%upper(right)
            end do
            right = upper_position(factor, column, column)
            factor%upper(right) = value
            if (abs(value) <= tolerance .or. .not. ieee_is_finite(value)) then
                call clear_factor(factor)
                call status_set(status, FORTSPARSE_SINGULAR, &
                    "sparse ILU(0) encountered a zero pivot")
                return
            end if
            ! L(i,column), i > column, using the completed U(*,column).
            do entry = factor%lower_col_ptr(column), &
                    factor%lower_col_ptr(column + 1) - 1
                row = factor%lower_row_idx(entry)
                if (row == column) cycle
                value = matrix_value(matrix, row, column)
                do inner = 1, column - 1
                    left = lower_position(factor, row, inner)
                    right = upper_position(factor, inner, column)
                    if (left > 0 .and. right > 0) value = value - &
                        factor%lower(left)*factor%upper(right)
                end do
                factor%lower(entry) = value/factor%upper( &
                    upper_position(factor, column, column))
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_sparse_incomplete_lu

    subroutine apply_sparse_incomplete_lu( &
            factor, right_hand_side, solution, status)
        !! Apply `(L U)^(-1)` with the fixed sparse factors.
        type(sparse_incomplete_lu_factor_t), intent(in) :: factor
        real(dp), intent(in) :: right_hand_side(:)
        real(dp), allocatable, intent(out) :: solution(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: dimension, column, entry, row, diagonal, future_column
        real(dp), allocatable :: intermediate(:)

        if (allocated(solution)) deallocate(solution)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "sparse ILU(0) apply received an invalid factor")
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
            do entry = factor%lower_col_ptr(column), &
                    factor%lower_col_ptr(column + 1) - 1
                row = factor%lower_row_idx(entry)
                if (row > column) intermediate(row) = intermediate(row) - &
                    factor%lower(entry)*intermediate(column)
            end do
        end do
        solution = intermediate
        do column = dimension, 1, -1
            do future_column = column + 1, dimension
                entry = upper_position(factor, column, future_column)
                if (entry > 0) solution(column) = solution(column) - &
                    factor%upper(entry)*solution(future_column)
            end do
            diagonal = upper_position(factor, column, column)
            if (abs(factor%upper(diagonal)) <= tiny(1.0_dp)) then
                deallocate(intermediate, solution)
                allocate(solution(0))
                call status_set(status, FORTSPARSE_SINGULAR, &
                    "sparse ILU(0) apply encountered a zero pivot")
                return
            end if
            solution(column) = solution(column)/factor%upper(diagonal)
        end do
        deallocate(intermediate)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_sparse_incomplete_lu

    subroutine apply_sparse_incomplete_lu_jvp( &
            factor, right_hand_side_dot, solution_dot, status)
        !! Fixed-factor JVP; factor construction is an inactive solver branch.
        type(sparse_incomplete_lu_factor_t), intent(in) :: factor
        real(dp), intent(in) :: right_hand_side_dot(:)
        real(dp), allocatable, intent(out) :: solution_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        call apply_sparse_incomplete_lu( &
            factor, right_hand_side_dot, solution_dot, status)
    end subroutine apply_sparse_incomplete_lu_jvp

    subroutine apply_sparse_incomplete_lu_vjp( &
            factor, solution_bar, right_hand_side_bar, status)
        !! Apply the transpose solve `(L U)^(-T)` to a cotangent.
        type(sparse_incomplete_lu_factor_t), intent(in) :: factor
        real(dp), intent(in) :: solution_bar(:)
        real(dp), allocatable, intent(out) :: right_hand_side_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: dimension, column, entry, row, diagonal
        real(dp), allocatable :: intermediate(:)

        if (allocated(right_hand_side_bar)) deallocate(right_hand_side_bar)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "sparse ILU(0) VJP received an invalid factor")
        if (.not. valid_factor(factor) .or. size(solution_bar) /= &
            factor%dimension .or. (size(solution_bar) > 0 .and. &
            .not. all(ieee_is_finite(solution_bar)))) then
            allocate(right_hand_side_bar(0))
            return
        end if
        dimension = factor%dimension
        allocate(intermediate(dimension), right_hand_side_bar(dimension))
        intermediate = solution_bar
        ! U^T z = solution_bar, a forward solve with the upper columns.
        do column = 1, dimension
            do entry = factor%upper_col_ptr(column), &
                    factor%upper_col_ptr(column + 1) - 1
                row = factor%upper_row_idx(entry)
                if (row < column) intermediate(column) = intermediate(column) - &
                    factor%upper(entry)*intermediate(row)
            end do
            diagonal = upper_position(factor, column, column)
            if (abs(factor%upper(diagonal)) <= tiny(1.0_dp)) then
                deallocate(intermediate, right_hand_side_bar)
                allocate(right_hand_side_bar(0))
                call status_set(status, FORTSPARSE_SINGULAR, &
                    "sparse ILU(0) VJP encountered a zero pivot")
                return
            end if
            intermediate(column) = intermediate(column)/factor%upper(diagonal)
        end do
        right_hand_side_bar = intermediate
        ! L^T rhs_bar = z; L has an implicit unit diagonal.
        do column = dimension, 1, -1
            do entry = factor%lower_col_ptr(column), &
                    factor%lower_col_ptr(column + 1) - 1
                row = factor%lower_row_idx(entry)
                if (row > column) right_hand_side_bar(column) = &
                    right_hand_side_bar(column) - factor%lower(entry)* &
                    right_hand_side_bar(row)
            end do
        end do
        deallocate(intermediate)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_sparse_incomplete_lu_vjp

    pure logical function valid_factor(factor) result(valid)
        type(sparse_incomplete_lu_factor_t), intent(in) :: factor
        integer :: column, entry, row

        valid = .false.
        if (factor%dimension < 1) return
        if (.not. allocated(factor%lower_col_ptr) .or. &
            .not. allocated(factor%lower_row_idx) .or. &
            .not. allocated(factor%lower) .or. &
            .not. allocated(factor%upper_col_ptr) .or. &
            .not. allocated(factor%upper_row_idx) .or. &
            .not. allocated(factor%upper)) return
        if (size(factor%lower_col_ptr) /= factor%dimension + 1 .or. &
            size(factor%lower_row_idx) /= size(factor%lower) .or. &
            size(factor%upper_col_ptr) /= factor%dimension + 1 .or. &
            size(factor%upper_row_idx) /= size(factor%upper)) return
        if (factor%lower_col_ptr(1) /= 1 .or. &
            factor%lower_col_ptr(factor%dimension + 1) /= size(factor%lower) + 1 .or. &
            factor%upper_col_ptr(1) /= 1 .or. &
            factor%upper_col_ptr(factor%dimension + 1) /= size(factor%upper) + 1) return
        do column = 1, factor%dimension
            if (factor%lower_col_ptr(column + 1) < factor%lower_col_ptr(column) .or. &
                factor%upper_col_ptr(column + 1) < factor%upper_col_ptr(column)) return
            if (upper_position(factor, column, column) == 0) return
            do entry = factor%lower_col_ptr(column), &
                    factor%lower_col_ptr(column + 1) - 1
                row = factor%lower_row_idx(entry)
                if (row <= column .or. row > factor%dimension .or. &
                    .not. ieee_is_finite(factor%lower(entry))) return
            end do
            do entry = factor%upper_col_ptr(column), &
                    factor%upper_col_ptr(column + 1) - 1
                row = factor%upper_row_idx(entry)
                if (row > column .or. row < 1 .or. &
                    .not. ieee_is_finite(factor%upper(entry))) return
            end do
        end do
        valid = .true.
    end function valid_factor

    pure integer function lower_position(factor, row, column) result(position)
        type(sparse_incomplete_lu_factor_t), intent(in) :: factor
        integer, intent(in) :: row, column
        integer :: entry

        position = 0
        if (.not. allocated(factor%lower_col_ptr)) return
        if (column < 1 .or. column > factor%dimension) return
        do entry = factor%lower_col_ptr(column), &
                factor%lower_col_ptr(column + 1) - 1
            if (factor%lower_row_idx(entry) == row) then
                position = entry
                return
            end if
        end do
    end function lower_position

    pure integer function upper_position(factor, row, column) result(position)
        type(sparse_incomplete_lu_factor_t), intent(in) :: factor
        integer, intent(in) :: row, column
        integer :: entry

        position = 0
        if (.not. allocated(factor%upper_col_ptr)) return
        if (column < 1 .or. column > factor%dimension) return
        do entry = factor%upper_col_ptr(column), &
                factor%upper_col_ptr(column + 1) - 1
            if (factor%upper_row_idx(entry) == row) then
                position = entry
                return
            end if
        end do
    end function upper_position

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

    subroutine clear_factor(factor)
        type(sparse_incomplete_lu_factor_t), intent(inout) :: factor

        if (allocated(factor%lower_col_ptr)) deallocate(factor%lower_col_ptr)
        if (allocated(factor%lower_row_idx)) deallocate(factor%lower_row_idx)
        if (allocated(factor%lower)) deallocate(factor%lower)
        if (allocated(factor%upper_col_ptr)) deallocate(factor%upper_col_ptr)
        if (allocated(factor%upper_row_idx)) deallocate(factor%upper_row_idx)
        if (allocated(factor%upper)) deallocate(factor%upper)
        factor%dimension = 0
    end subroutine clear_factor

end module fortfem_sparse_incomplete_lu
