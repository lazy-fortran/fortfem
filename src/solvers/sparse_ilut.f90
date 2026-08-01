module fortfem_sparse_ilut
    !! Drop- and fill-controlled sparse ILU preconditioner.
    !!
    !! The deterministic dense constructor is retained as a reference path for
    !! fixed-factor preconditioning.  The row-oriented constructor below uses
    !! O(n + nnz) work storage, while exporting the same sparse CSC factors and
    !! apply/JVP/VJP contract.  The two construction paths provide an
    !! independent small-matrix oracle and a memory-scalable production path.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_is_valid, csc_t, fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, FORTSPARSE_SINGULAR
    implicit none
    private

    type, public :: sparse_ilut_factor_t
        private
        integer :: dimension = 0
        integer :: max_fill_per_column = 0
        real(dp) :: drop_tolerance = 0.0_dp
        integer, allocatable :: lower_col_ptr(:), lower_row_idx(:)
        real(dp), allocatable :: lower(:)
        integer, allocatable :: upper_col_ptr(:), upper_row_idx(:)
        real(dp), allocatable :: upper(:)
    end type sparse_ilut_factor_t

    type :: ilut_row_entries_t
        integer :: count = 0
        integer, allocatable :: index(:)
        real(dp), allocatable :: value(:)
    end type ilut_row_entries_t

    public :: build_sparse_ilut
    public :: build_sparse_ilut_row
    public :: apply_sparse_ilut
    public :: apply_sparse_ilut_jvp
    public :: apply_sparse_ilut_vjp

contains

    subroutine build_sparse_ilut( &
            matrix, drop_tolerance, max_fill_per_column, factor, status)
        !! Build an ILUT factor with deterministic drop and fill policies.
        type(csc_t), intent(in) :: matrix
        real(dp), intent(in) :: drop_tolerance
        integer, intent(in) :: max_fill_per_column
        type(sparse_ilut_factor_t), intent(inout) :: factor
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: work(:, :)
        logical, allocatable :: lower_keep(:, :), upper_keep(:, :)
        real(dp) :: scale, threshold, pivot_tolerance, pivot
        integer :: dimension, column, entry, row, inner, lower_nnz
        integer :: upper_nnz, position, diagonal

        call clear_factor(factor)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "sparse ILUT received an invalid CSC matrix or policy")
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
        pivot_tolerance = 128.0_dp*epsilon(1.0_dp)*scale*real(dimension, dp)
        threshold = drop_tolerance*scale
        allocate(work(dimension, dimension))
        work = 0.0_dp
        do column = 1, dimension
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                work(row, column) = matrix%val(entry)
            end do
        end do

        ! Complete dense LU numeric phase.  Dropping is applied only after
        ! this phase, so the result is reproducible for equivalent CSC maps.
        do column = 1, dimension
            pivot = work(column, column)
            if (.not. ieee_is_finite(pivot) .or. &
                abs(pivot) <= pivot_tolerance) then
                deallocate(work)
                call status_set(status, FORTSPARSE_SINGULAR, &
                    "sparse ILUT encountered a zero pivot")
                return
            end if
            do row = column + 1, dimension
                work(row, column) = work(row, column)/pivot
                if (.not. ieee_is_finite(work(row, column))) then
                    deallocate(work)
                    call status_set(status, FORTSPARSE_SINGULAR, &
                        "sparse ILUT produced a non-finite multiplier")
                    return
                end if
                do inner = column + 1, dimension
                    work(row, inner) = work(row, inner) - &
                        work(row, column)*work(column, inner)
                end do
            end do
        end do
        if (.not. all(ieee_is_finite(work))) then
            deallocate(work)
            call status_set(status, FORTSPARSE_SINGULAR, &
                "sparse ILUT produced a non-finite factor")
            return
        end if

        allocate(lower_keep(dimension, dimension), upper_keep(dimension, dimension))
        lower_keep = .false.
        upper_keep = .false.
        do column = 1, dimension
            call select_column_entries(work, column, column + 1, dimension, &
                threshold, max_fill_per_column, lower_keep)
            call select_column_entries(work, column, 1, column - 1, threshold, &
                max_fill_per_column, upper_keep)
            upper_keep(column, column) = .true.
        end do

        lower_nnz = count(lower_keep)
        upper_nnz = count(upper_keep)
        factor%dimension = dimension
        factor%max_fill_per_column = max_fill_per_column
        factor%drop_tolerance = drop_tolerance
        allocate(factor%lower_col_ptr(dimension + 1), &
            factor%lower_row_idx(lower_nnz), factor%lower(lower_nnz), &
            factor%upper_col_ptr(dimension + 1), &
            factor%upper_row_idx(upper_nnz), factor%upper(upper_nnz))
        factor%lower_col_ptr(1) = 1
        factor%upper_col_ptr(1) = 1
        position = 0
        do column = 1, dimension
            do row = column + 1, dimension
                if (.not. lower_keep(row, column)) cycle
                position = position + 1
                factor%lower_row_idx(position) = row
                factor%lower(position) = work(row, column)
            end do
            factor%lower_col_ptr(column + 1) = position + 1
        end do
        position = 0
        do column = 1, dimension
            do row = 1, column
                if (.not. upper_keep(row, column)) cycle
                position = position + 1
                factor%upper_row_idx(position) = row
                factor%upper(position) = work(row, column)
            end do
            factor%upper_col_ptr(column + 1) = position + 1
        end do

        deallocate(work, lower_keep, upper_keep)
        do column = 1, dimension
            diagonal = upper_position(factor, column, column)
            if (diagonal == 0 .or. abs(factor%upper(diagonal)) <= &
                pivot_tolerance) then
                call clear_factor(factor)
                call status_set(status, FORTSPARSE_SINGULAR, &
                    "sparse ILUT retained a zero diagonal")
                return
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_sparse_ilut

    subroutine build_sparse_ilut_row( &
            matrix, drop_tolerance, max_fill_per_row, factor, status)
        !! Build ILUT with row-oriented O(n + nnz) work storage.
        !!
        !! The factorization is the usual no-pivot ILUT sweep: each row is
        !! eliminated against the already retained U rows, then the largest
        !! `max_fill_per_row` lower and upper off-diagonal entries are kept.
        !! The exported factor is converted to the same CSC representation as
        !! `build_sparse_ilut`, so apply and derivative clients are unchanged.
        type(csc_t), intent(in) :: matrix
        real(dp), intent(in) :: drop_tolerance
        integer, intent(in) :: max_fill_per_row
        type(sparse_ilut_factor_t), intent(inout) :: factor
        type(fortsparse_status_t), intent(out) :: status

        type(ilut_row_entries_t), allocatable :: lower_rows(:), upper_rows(:)
        integer, allocatable :: row_ptr(:), row_cursor(:), row_columns(:)
        integer, allocatable :: lower_counts(:), upper_counts(:)
        integer, allocatable :: lower_next(:), upper_next(:)
        integer :: dimension, row, column, entry, pivot_row, position
        integer :: lower_nnz, upper_nnz
        integer :: diagonal
        integer, allocatable :: marker(:)
        real(dp), allocatable :: row_values(:), work(:), diagonal_values(:)
        real(dp) :: scale, threshold, pivot_tolerance, pivot, multiplier

        call clear_factor(factor)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "row-oriented ILUT received an invalid matrix or policy")
        if (.not. csc_is_valid(matrix) .or. matrix%nrow /= matrix%ncol) return
        if (.not. ieee_is_finite(drop_tolerance) .or. &
            drop_tolerance < 0.0_dp .or. max_fill_per_row < 0) return
        if (matrix%nnz > 0) then
            if (.not. all(ieee_is_finite(matrix%val))) return
        end if

        dimension = matrix%nrow
        if (dimension < 1) return
        scale = 1.0_dp
        if (matrix%nnz > 0) scale = max(scale, maxval(abs(matrix%val)))
        threshold = drop_tolerance*scale
        pivot_tolerance = 128.0_dp*epsilon(1.0_dp)*scale* &
            real(dimension, dp)

        allocate(row_ptr(dimension + 1), row_cursor(dimension), &
            row_columns(matrix%nnz), row_values(matrix%nnz))
        row_ptr = 1
        do column = 1, dimension
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                row_ptr(row + 1) = row_ptr(row + 1) + 1
            end do
        end do
        row_ptr(1) = 1
        do row = 1, dimension
            row_ptr(row + 1) = row_ptr(row + 1) + row_ptr(row) - 1
        end do
        row_cursor = row_ptr(1:dimension)
        do column = 1, dimension
            do entry = matrix%col_ptr(column), matrix%col_ptr(column + 1) - 1
                row = matrix%row_idx(entry)
                position = row_cursor(row)
                row_columns(position) = column
                row_values(position) = matrix%val(entry)
                row_cursor(row) = position + 1
            end do
        end do

        allocate(lower_rows(dimension), upper_rows(dimension), &
            marker(dimension), work(dimension), diagonal_values(dimension))
        allocate(lower_counts(dimension), upper_counts(dimension), &
            lower_next(dimension), upper_next(dimension))
        marker = 0
        diagonal_values = 0.0_dp

        do row = 1, dimension
            work = 0.0_dp
            do entry = row_ptr(row), row_ptr(row + 1) - 1
                column = row_columns(entry)
                if (marker(column) /= row) then
                    marker(column) = row
                    work(column) = row_values(entry)
                else
                    work(column) = work(column) + row_values(entry)
                end if
            end do

            do pivot_row = 1, row - 1
                if (marker(pivot_row) /= row) cycle
                if (.not. ieee_is_finite(work(pivot_row))) then
                    call clear_factor(factor)
                    call status_set(status, FORTSPARSE_SINGULAR, &
                        "row-oriented ILUT produced a non-finite multiplier")
                    return
                end if
                multiplier = work(pivot_row)/diagonal_values(pivot_row)
                work(pivot_row) = multiplier
                if (.not. ieee_is_finite(multiplier)) then
                    call clear_factor(factor)
                    call status_set(status, FORTSPARSE_SINGULAR, &
                        "row-oriented ILUT produced a non-finite multiplier")
                    return
                end if
                do entry = 1, upper_rows(pivot_row)%count
                    column = upper_rows(pivot_row)%index(entry)
                    if (column <= pivot_row) cycle
                    if (marker(column) /= row) then
                        marker(column) = row
                        work(column) = 0.0_dp
                    end if
                    work(column) = work(column) - multiplier* &
                        upper_rows(pivot_row)%value(entry)
                end do
            end do

            pivot = 0.0_dp
            if (marker(row) == row) pivot = work(row)
            if (.not. ieee_is_finite(pivot) .or. abs(pivot) <= pivot_tolerance) then
                call clear_factor(factor)
                call status_set(status, FORTSPARSE_SINGULAR, &
                    "row-oriented ILUT encountered a zero pivot")
                return
            end if
            diagonal_values(row) = pivot
            call retain_row_side(work, 1, row - 1, threshold, &
                max_fill_per_row, lower_rows(row))
            call retain_row_side(work, row + 1, dimension, threshold, &
                max_fill_per_row, upper_rows(row))
            call prepend_row_diagonal(upper_rows(row), row, pivot)
        end do

        lower_counts = 0
        upper_counts = 0
        do row = 1, dimension
            do entry = 1, lower_rows(row)%count
                column = lower_rows(row)%index(entry)
                lower_counts(column) = lower_counts(column) + 1
            end do
            do entry = 1, upper_rows(row)%count
                column = upper_rows(row)%index(entry)
                upper_counts(column) = upper_counts(column) + 1
            end do
        end do
        lower_nnz = sum(lower_counts)
        upper_nnz = sum(upper_counts)

        factor%dimension = dimension
        factor%max_fill_per_column = max_fill_per_row
        factor%drop_tolerance = drop_tolerance
        allocate(factor%lower_col_ptr(dimension + 1), &
            factor%lower_row_idx(lower_nnz), factor%lower(lower_nnz), &
            factor%upper_col_ptr(dimension + 1), &
            factor%upper_row_idx(upper_nnz), factor%upper(upper_nnz))
        factor%lower_col_ptr(1) = 1
        factor%upper_col_ptr(1) = 1
        do column = 1, dimension
            factor%lower_col_ptr(column + 1) = &
                factor%lower_col_ptr(column) + lower_counts(column)
            factor%upper_col_ptr(column + 1) = &
                factor%upper_col_ptr(column) + upper_counts(column)
        end do
        lower_next = factor%lower_col_ptr(1:dimension)
        upper_next = factor%upper_col_ptr(1:dimension)
        do row = 1, dimension
            do entry = 1, lower_rows(row)%count
                column = lower_rows(row)%index(entry)
                position = lower_next(column)
                factor%lower_row_idx(position) = row
                factor%lower(position) = lower_rows(row)%value(entry)
                lower_next(column) = position + 1
            end do
            do entry = 1, upper_rows(row)%count
                column = upper_rows(row)%index(entry)
                position = upper_next(column)
                factor%upper_row_idx(position) = row
                factor%upper(position) = upper_rows(row)%value(entry)
                upper_next(column) = position + 1
            end do
        end do

        do column = 1, dimension
            diagonal = upper_position(factor, column, column)
            if (diagonal == 0) then
                call clear_factor(factor)
                call status_set(status, FORTSPARSE_SINGULAR, &
                    "row-oriented ILUT retained a zero diagonal")
                return
            end if
            if (abs(factor%upper(diagonal)) <= pivot_tolerance) then
                call clear_factor(factor)
                call status_set(status, FORTSPARSE_SINGULAR, &
                    "row-oriented ILUT retained a zero diagonal")
                return
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_sparse_ilut_row

    subroutine apply_sparse_ilut( &
            factor, right_hand_side, solution, status)
        !! Apply `(L U)^(-1)` with the fixed sparse ILUT factors.
        type(sparse_ilut_factor_t), intent(in) :: factor
        real(dp), intent(in) :: right_hand_side(:)
        real(dp), allocatable, intent(out) :: solution(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: intermediate(:)
        integer :: dimension, column, entry, row, diagonal

        if (allocated(solution)) deallocate(solution)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "sparse ILUT apply received an invalid factor")
        if (.not. valid_factor(factor)) then
            allocate(solution(0))
            return
        end if
        dimension = factor%dimension
        if (size(right_hand_side) /= dimension) then
            allocate(solution(0))
            return
        end if
        if (dimension > 0) then
            if (.not. all(ieee_is_finite(right_hand_side))) then
                allocate(solution(0))
                return
            end if
        end if

        allocate(intermediate(dimension), solution(dimension))
        intermediate = right_hand_side
        do column = 1, dimension
            do entry = factor%lower_col_ptr(column), &
                    factor%lower_col_ptr(column + 1) - 1
                row = factor%lower_row_idx(entry)
                intermediate(row) = intermediate(row) - &
                    factor%lower(entry)*intermediate(column)
            end do
        end do
        solution = intermediate
        do column = dimension, 1, -1
            diagonal = upper_position(factor, column, column)
            solution(column) = solution(column)/factor%upper(diagonal)
            do entry = factor%upper_col_ptr(column), &
                    factor%upper_col_ptr(column + 1) - 1
                row = factor%upper_row_idx(entry)
                if (row < column) solution(row) = solution(row) - &
                    factor%upper(entry)*solution(column)
            end do
        end do
        deallocate(intermediate)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_sparse_ilut

    subroutine apply_sparse_ilut_jvp( &
            factor, right_hand_side_dot, solution_dot, status)
        !! Fixed-factor JVP; factor construction is an inactive branch.
        type(sparse_ilut_factor_t), intent(in) :: factor
        real(dp), intent(in) :: right_hand_side_dot(:)
        real(dp), allocatable, intent(out) :: solution_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        call apply_sparse_ilut(factor, right_hand_side_dot, solution_dot, status)
    end subroutine apply_sparse_ilut_jvp

    subroutine apply_sparse_ilut_vjp( &
            factor, solution_bar, right_hand_side_bar, status)
        !! Fixed-factor real VJP, applying the transpose triangular solves.
        type(sparse_ilut_factor_t), intent(in) :: factor
        real(dp), intent(in) :: solution_bar(:)
        real(dp), allocatable, intent(out) :: right_hand_side_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: intermediate(:)
        integer :: dimension, column, entry, row, diagonal

        if (allocated(right_hand_side_bar)) deallocate(right_hand_side_bar)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "sparse ILUT VJP received an invalid factor")
        if (.not. valid_factor(factor)) then
            allocate(right_hand_side_bar(0))
            return
        end if
        dimension = factor%dimension
        if (size(solution_bar) /= dimension) then
            allocate(right_hand_side_bar(0))
            return
        end if
        if (dimension > 0) then
            if (.not. all(ieee_is_finite(solution_bar))) then
                allocate(right_hand_side_bar(0))
                return
            end if
        end if

        allocate(intermediate(dimension), right_hand_side_bar(dimension))
        intermediate = solution_bar
        do column = 1, dimension
            do entry = factor%upper_col_ptr(column), &
                    factor%upper_col_ptr(column + 1) - 1
                row = factor%upper_row_idx(entry)
                if (row < column) intermediate(column) = intermediate(column) - &
                    factor%upper(entry)*intermediate(row)
            end do
            diagonal = upper_position(factor, column, column)
            intermediate(column) = intermediate(column)/factor%upper(diagonal)
        end do
        right_hand_side_bar = intermediate
        do column = dimension, 1, -1
            do entry = factor%lower_col_ptr(column), &
                    factor%lower_col_ptr(column + 1) - 1
                row = factor%lower_row_idx(entry)
                right_hand_side_bar(column) = right_hand_side_bar(column) - &
                    factor%lower(entry)*right_hand_side_bar(row)
            end do
        end do
        deallocate(intermediate)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_sparse_ilut_vjp

    subroutine retain_row_side( &
            work, first_column, last_column, threshold, max_fill, entries)
        real(dp), intent(in) :: work(:)
        integer, intent(in) :: first_column, last_column, max_fill
        real(dp), intent(in) :: threshold
        type(ilut_row_entries_t), intent(inout) :: entries

        logical, allocatable :: selected(:)
        integer :: candidate_count, target_count, selected_count
        integer :: column, best_column
        real(dp) :: best_value, value

        call clear_row_entries(entries)
        if (first_column > last_column .or. max_fill == 0) return
        candidate_count = 0
        do column = first_column, last_column
            if (abs(work(column)) > threshold) candidate_count = &
                candidate_count + 1
        end do
        if (candidate_count == 0) return
        target_count = min(max_fill, candidate_count)
        allocate(entries%index(target_count), entries%value(target_count))
        entries%count = target_count
        if (target_count == candidate_count) then
            selected_count = 0
            do column = first_column, last_column
                if (abs(work(column)) <= threshold) cycle
                selected_count = selected_count + 1
                entries%index(selected_count) = column
                entries%value(selected_count) = work(column)
            end do
            return
        end if

        allocate(selected(size(work)))
        selected = .false.
        do selected_count = 1, target_count
            best_column = 0
            best_value = -1.0_dp
            do column = first_column, last_column
                value = abs(work(column))
                if (selected(column)) cycle
                if (value <= threshold .or. value <= best_value) cycle
                best_column = column
                best_value = value
            end do
            if (best_column == 0) exit
            selected(best_column) = .true.
            entries%index(selected_count) = best_column
            entries%value(selected_count) = work(best_column)
        end do
        deallocate(selected)
    end subroutine retain_row_side

    subroutine prepend_row_diagonal(entries, row, value)
        type(ilut_row_entries_t), intent(inout) :: entries
        integer, intent(in) :: row
        real(dp), intent(in) :: value

        integer, allocatable :: index(:)
        real(dp), allocatable :: values(:)
        integer :: entry, old_count

        old_count = entries%count
        allocate(index(old_count + 1), values(old_count + 1))
        index(1) = row
        values(1) = value
        do entry = 1, old_count
            index(entry + 1) = entries%index(entry)
            values(entry + 1) = entries%value(entry)
        end do
        if (allocated(entries%index)) deallocate(entries%index)
        if (allocated(entries%value)) deallocate(entries%value)
        call move_alloc(index, entries%index)
        call move_alloc(values, entries%value)
        entries%count = old_count + 1
    end subroutine prepend_row_diagonal

    subroutine clear_row_entries(entries)
        type(ilut_row_entries_t), intent(inout) :: entries

        if (allocated(entries%index)) deallocate(entries%index)
        if (allocated(entries%value)) deallocate(entries%value)
        entries%count = 0
    end subroutine clear_row_entries

    subroutine select_column_entries( &
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
    end subroutine select_column_entries

    pure logical function valid_factor(factor) result(valid)
        type(sparse_ilut_factor_t), intent(in) :: factor
        integer :: column, entry, row

        valid = .false.
        if (factor%dimension < 1 .or. factor%max_fill_per_column < 0 .or. &
            .not. ieee_is_finite(factor%drop_tolerance) .or. &
            factor%drop_tolerance < 0.0_dp) return
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
                if (row < 1 .or. row > column .or. &
                    .not. ieee_is_finite(factor%upper(entry))) return
            end do
        end do
        valid = .true.
    end function valid_factor

    pure integer function lower_position(factor, row, column) result(position)
        type(sparse_ilut_factor_t), intent(in) :: factor
        integer, intent(in) :: row, column
        integer :: entry

        position = 0
        if (.not. allocated(factor%lower_col_ptr)) return
        if (column < 1 .or. column > factor%dimension) return
        do entry = factor%lower_col_ptr(column), factor%lower_col_ptr(column + 1) - 1
            if (factor%lower_row_idx(entry) == row) then
                position = entry
                return
            end if
        end do
    end function lower_position

    pure integer function upper_position(factor, row, column) result(position)
        type(sparse_ilut_factor_t), intent(in) :: factor
        integer, intent(in) :: row, column
        integer :: entry

        position = 0
        if (.not. allocated(factor%upper_col_ptr)) return
        if (column < 1 .or. column > factor%dimension) return
        do entry = factor%upper_col_ptr(column), factor%upper_col_ptr(column + 1) - 1
            if (factor%upper_row_idx(entry) == row) then
                position = entry
                return
            end if
        end do
    end function upper_position

    subroutine clear_factor(factor)
        type(sparse_ilut_factor_t), intent(inout) :: factor

        if (allocated(factor%lower_col_ptr)) deallocate(factor%lower_col_ptr)
        if (allocated(factor%lower_row_idx)) deallocate(factor%lower_row_idx)
        if (allocated(factor%lower)) deallocate(factor%lower)
        if (allocated(factor%upper_col_ptr)) deallocate(factor%upper_col_ptr)
        if (allocated(factor%upper_row_idx)) deallocate(factor%upper_row_idx)
        if (allocated(factor%upper)) deallocate(factor%upper)
        factor%dimension = 0
        factor%max_fill_per_column = 0
        factor%drop_tolerance = 0.0_dp
    end subroutine clear_factor

end module fortfem_sparse_ilut
