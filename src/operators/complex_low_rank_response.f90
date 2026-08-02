module fortfem_complex_low_rank_response
    !! Neutral complex low-rank response storage and matrix-free actions.
    !!
    !! Cross approximation is deterministic for a fixed dense matrix and
    !! fixed tolerance/rank cap.  The retained factors use A ~= U V^T (not a
    !! Hermitian factorization), which is the natural convention for complex
    !! frequency response blocks.  Factor perturbations have complete JVP and
    !! real-complex VJP products; pivot selection is a fixed-topology build
    !! step and is not differentiated.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_kinds, only: dp
    implicit none
    private

    character(len=*), parameter :: schema_version = "fortfem-complex-low-rank-1"
    integer(int64), parameter :: max_factor_values = 8_int64*1000_int64*1000_int64

    type, public :: complex_low_rank_matrix_t
        character(len=32) :: schema = schema_version
        integer :: row_count = 0
        integer :: column_count = 0
        integer :: rank = 0
        real(dp) :: relative_error_bound = 0.0_dp
        complex(dp), allocatable :: left(:, :)
        complex(dp), allocatable :: right(:, :)
    contains
        procedure, private :: assign_complex_low_rank_matrix
        generic :: assignment(=) => assign_complex_low_rank_matrix
    end type complex_low_rank_matrix_t

    public :: apply_complex_low_rank_matrix
    public :: apply_complex_low_rank_matrix_jvp
    public :: apply_complex_low_rank_matrix_vjp
    public :: compress_complex_matrix_cross
    public :: initialize_complex_low_rank_matrix
    public :: materialize_complex_low_rank_matrix
    public :: validate_complex_low_rank_matrix

    interface finite_complex
        module procedure finite_complex_matrix
        module procedure finite_complex_vector
    end interface finite_complex

contains

    subroutine initialize_complex_low_rank_matrix( &
            left, right, relative_error_bound, low_rank, status)
        complex(dp), intent(in) :: left(:, :), right(:, :)
        real(dp), intent(in) :: relative_error_bound
        type(complex_low_rank_matrix_t), intent(out) :: low_rank
        integer, intent(out) :: status

        call clear_complex_low_rank_matrix(low_rank)
        status = 1
        if (.not. valid_factor_inputs(left, right, relative_error_bound)) return
        low_rank%row_count = size(left, 1)
        low_rank%column_count = size(right, 1)
        low_rank%rank = size(left, 2)
        low_rank%relative_error_bound = relative_error_bound
        allocate(low_rank%left, source=left)
        allocate(low_rank%right, source=right)
        status = 0
    end subroutine initialize_complex_low_rank_matrix

    subroutine compress_complex_matrix_cross(matrix, tolerance, max_rank, low_rank, status)
        complex(dp), intent(in) :: matrix(:, :)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_rank
        type(complex_low_rank_matrix_t), intent(out) :: low_rank
        integer, intent(out) :: status

        complex(dp), allocatable :: residual(:, :), work_left(:, :), work_right(:, :)
        real(dp) :: scale, pivot_abs
        complex(dp) :: pivot
        integer :: row, column, pivot_row, pivot_column, retained_rank

        call clear_complex_low_rank_matrix(low_rank)
        status = 1
        if (size(matrix, 1) < 1 .or. size(matrix, 2) < 1 .or. &
            tolerance < 0.0_dp .or. .not. ieee_is_finite(tolerance) .or. &
            max_rank < 0 .or. max_rank > min(size(matrix, 1), size(matrix, 2)) .or. &
            .not. finite_complex(matrix)) return
        if (int(size(matrix, 1), int64)*int(max_rank, int64) > max_factor_values .or. &
            int(size(matrix, 2), int64)*int(max_rank, int64) > max_factor_values) return

        allocate(residual, source=matrix)
        allocate(work_left(size(matrix, 1), max_rank), work_right(size(matrix, 2), max_rank))
        work_left = cmplx(0.0_dp, 0.0_dp, dp)
        work_right = cmplx(0.0_dp, 0.0_dp, dp)
        scale = max(1.0_dp, maxval(abs(matrix)))
        retained_rank = 0
        do while (retained_rank < max_rank)
            pivot_abs = 0.0_dp
            pivot_row = 1
            pivot_column = 1
            do column = 1, size(matrix, 2)
                do row = 1, size(matrix, 1)
                    if (abs(residual(row, column)) > pivot_abs) then
                        pivot_abs = abs(residual(row, column))
                        pivot_row = row
                        pivot_column = column
                    end if
                end do
            end do
            if (pivot_abs <= tolerance*scale) exit
            retained_rank = retained_rank + 1
            pivot = residual(pivot_row, pivot_column)
            work_left(:, retained_rank) = residual(:, pivot_column)
            work_right(:, retained_rank) = residual(pivot_row, :)/pivot
            do column = 1, size(matrix, 2)
                do row = 1, size(matrix, 1)
                    residual(row, column) = residual(row, column) - &
                        work_left(row, retained_rank)*work_right(column, retained_rank)
                end do
            end do
        end do

        low_rank%row_count = size(matrix, 1)
        low_rank%column_count = size(matrix, 2)
        low_rank%rank = retained_rank
        low_rank%relative_error_bound = maxval(abs(residual))/scale
        allocate(low_rank%left(size(matrix, 1), retained_rank), &
            low_rank%right(size(matrix, 2), retained_rank))
        if (retained_rank > 0) then
            low_rank%left = work_left(:, 1:retained_rank)
            low_rank%right = work_right(:, 1:retained_rank)
        end if
        status = 0
    end subroutine compress_complex_matrix_cross

    subroutine materialize_complex_low_rank_matrix(low_rank, matrix, status)
        type(complex_low_rank_matrix_t), intent(in) :: low_rank
        complex(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. validate_complex_low_rank_matrix(low_rank, status) .or. &
            any(shape(matrix) /= [low_rank%row_count, low_rank%column_count])) return
        if (low_rank%rank > 0) matrix = matmul(low_rank%left, transpose(low_rank%right))
        status = 0
    end subroutine materialize_complex_low_rank_matrix

    subroutine apply_complex_low_rank_matrix(low_rank, x, y, status)
        type(complex_low_rank_matrix_t), intent(in) :: low_rank
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: status

        integer :: factor
        complex(dp) :: coefficient

        y = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. validate_complex_low_rank_matrix(low_rank, status) .or. &
            size(x) /= low_rank%column_count .or. size(y) /= low_rank%row_count .or. &
            .not. finite_complex(x)) return
        do factor = 1, low_rank%rank
            coefficient = sum(low_rank%right(:, factor)*x)
            y = y + low_rank%left(:, factor)*coefficient
        end do
        status = 0
    end subroutine apply_complex_low_rank_matrix

    subroutine apply_complex_low_rank_matrix_jvp( &
            low_rank, left_dot, right_dot, x, x_dot, y_dot, status)
        type(complex_low_rank_matrix_t), intent(in) :: low_rank
        complex(dp), intent(in) :: left_dot(:, :), right_dot(:, :)
        complex(dp), intent(in) :: x(:), x_dot(:)
        complex(dp), intent(out) :: y_dot(:)
        integer, intent(out) :: status

        integer :: factor
        complex(dp) :: coefficient, coefficient_dot

        y_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. validate_complex_low_rank_matrix(low_rank, status) .or. &
            any(shape(left_dot) /= shape(low_rank%left)) .or. &
            any(shape(right_dot) /= shape(low_rank%right)) .or. &
            size(x) /= low_rank%column_count .or. size(x_dot) /= size(x) .or. &
            size(y_dot) /= low_rank%row_count .or. .not. finite_complex(left_dot) .or. &
            .not. finite_complex(right_dot) .or. .not. finite_complex(x) .or. &
            .not. finite_complex(x_dot)) return
        do factor = 1, low_rank%rank
            coefficient = sum(low_rank%right(:, factor)*x)
            coefficient_dot = sum(right_dot(:, factor)*x) + &
                sum(low_rank%right(:, factor)*x_dot)
            y_dot = y_dot + left_dot(:, factor)*coefficient + &
                low_rank%left(:, factor)*coefficient_dot
        end do
        status = 0
    end subroutine apply_complex_low_rank_matrix_jvp

    subroutine apply_complex_low_rank_matrix_vjp( &
            low_rank, x, y_bar, left_bar, right_bar, x_bar, status)
        type(complex_low_rank_matrix_t), intent(in) :: low_rank
        complex(dp), intent(in) :: x(:), y_bar(:)
        complex(dp), allocatable, intent(out) :: left_bar(:, :), right_bar(:, :)
        complex(dp), intent(out) :: x_bar(:)
        integer, intent(out) :: status

        integer :: factor
        complex(dp) :: coefficient, coefficient_bar

        if (allocated(left_bar)) deallocate(left_bar)
        if (allocated(right_bar)) deallocate(right_bar)
        status = 1
        if (.not. validate_complex_low_rank_matrix(low_rank, status)) then
            allocate(left_bar(0, 0), right_bar(0, 0))
            x_bar = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if
        allocate(left_bar(low_rank%row_count, low_rank%rank), &
            right_bar(low_rank%column_count, low_rank%rank))
        left_bar = cmplx(0.0_dp, 0.0_dp, dp)
        right_bar = cmplx(0.0_dp, 0.0_dp, dp)
        x_bar = cmplx(0.0_dp, 0.0_dp, dp)
        if (size(x) /= low_rank%column_count .or. size(y_bar) /= low_rank%row_count .or. &
            size(x_bar) /= low_rank%column_count .or. .not. finite_complex(x) .or. &
            .not. finite_complex(y_bar)) return
        do factor = 1, low_rank%rank
            coefficient = sum(low_rank%right(:, factor)*x)
            coefficient_bar = sum(conjg(low_rank%left(:, factor))*y_bar)
            left_bar(:, factor) = y_bar*conjg(coefficient)
            right_bar(:, factor) = coefficient_bar*conjg(x)
            x_bar = x_bar + coefficient_bar*conjg(low_rank%right(:, factor))
        end do
        status = 0
    end subroutine apply_complex_low_rank_matrix_vjp

    logical function validate_complex_low_rank_matrix(low_rank, status) result(valid)
        type(complex_low_rank_matrix_t), intent(in) :: low_rank
        integer, intent(out) :: status

        valid = .false.
        status = 1
        if (low_rank%schema /= schema_version .or. low_rank%row_count < 1 .or. &
            low_rank%column_count < 1 .or. low_rank%rank < 0 .or. &
            low_rank%rank > min(low_rank%row_count, low_rank%column_count) .or. &
            .not. ieee_is_finite(low_rank%relative_error_bound) .or. &
            low_rank%relative_error_bound < 0.0_dp) return
        if (.not. allocated(low_rank%left) .or. .not. allocated(low_rank%right)) return
        if (any(shape(low_rank%left) /= [low_rank%row_count, low_rank%rank]) .or. &
            any(shape(low_rank%right) /= [low_rank%column_count, low_rank%rank]) .or. &
            .not. finite_complex(low_rank%left) .or. .not. finite_complex(low_rank%right)) return
        valid = .true.
        status = 0
    end function validate_complex_low_rank_matrix

    subroutine assign_complex_low_rank_matrix(lhs, rhs)
        class(complex_low_rank_matrix_t), intent(out) :: lhs
        type(complex_low_rank_matrix_t), intent(in) :: rhs

        call clear_complex_low_rank_matrix(lhs)
        lhs%schema = rhs%schema
        lhs%row_count = rhs%row_count
        lhs%column_count = rhs%column_count
        lhs%rank = rhs%rank
        lhs%relative_error_bound = rhs%relative_error_bound
        if (allocated(rhs%left)) allocate(lhs%left, source=rhs%left)
        if (allocated(rhs%right)) allocate(lhs%right, source=rhs%right)
    end subroutine assign_complex_low_rank_matrix

    subroutine clear_complex_low_rank_matrix(low_rank)
        type(complex_low_rank_matrix_t), intent(inout) :: low_rank

        if (allocated(low_rank%left)) deallocate(low_rank%left)
        if (allocated(low_rank%right)) deallocate(low_rank%right)
        low_rank%schema = schema_version
        low_rank%row_count = 0
        low_rank%column_count = 0
        low_rank%rank = 0
        low_rank%relative_error_bound = 0.0_dp
    end subroutine clear_complex_low_rank_matrix

    logical function valid_factor_inputs(left, right, relative_error_bound) result(valid)
        complex(dp), intent(in) :: left(:, :), right(:, :)
        real(dp), intent(in) :: relative_error_bound

        valid = size(left, 1) > 0 .and. size(right, 1) > 0 .and. &
            size(left, 2) == size(right, 2) .and. &
            int(size(left), int64) <= max_factor_values .and. &
            int(size(right), int64) <= max_factor_values .and. &
            relative_error_bound >= 0.0_dp .and. ieee_is_finite(relative_error_bound) .and. &
            finite_complex(left) .and. finite_complex(right)
    end function valid_factor_inputs

    logical function finite_complex_matrix(value) result(valid)
        complex(dp), intent(in) :: value(:, :)

        valid = all(ieee_is_finite(real(value, dp))) .and. &
            all(ieee_is_finite(aimag(value)))
    end function finite_complex_matrix

    logical function finite_complex_vector(value) result(valid)
        complex(dp), intent(in) :: value(:)

        valid = all(ieee_is_finite(real(value, dp))) .and. &
            all(ieee_is_finite(aimag(value)))
    end function finite_complex_vector

end module fortfem_complex_low_rank_response
