module fortfem_incomplete_cholesky
    !! Dense incomplete-Cholesky(0) factor and fixed-factor apply.
    !!
    !! The lower-triangular sparsity pattern is inherited from the lower
    !! triangle of the input matrix; no fill is introduced. The factorization
    !! is a solver/preconditioner primitive, not a replacement for a
    !! structure-preserving space or a tree-cotree gauge.
    use fortfem_kinds, only: dp
    implicit none
    private

    type, public :: incomplete_cholesky_factor_t
        private
        real(dp), allocatable :: lower(:, :)
    end type incomplete_cholesky_factor_t

    public :: build_incomplete_cholesky
    public :: apply_incomplete_cholesky

contains

    subroutine build_incomplete_cholesky(matrix, factor, status)
        real(dp), intent(in) :: matrix(:, :)
        type(incomplete_cholesky_factor_t), intent(inout) :: factor
        integer, intent(out) :: status
        real(dp) :: scale, tolerance, value
        integer :: dimension, row, column, inner

        call clear_factor(factor)
        status = 1
        dimension = size(matrix, 1)
        if (dimension < 1 .or. size(matrix, 2) /= dimension) return
        if (.not. finite_values(matrix)) return
        scale = max(1.0_dp, maxval(abs(matrix)))
        tolerance = 128.0_dp*epsilon(1.0_dp)*scale*real(dimension, dp)
        if (maxval(abs(matrix - transpose(matrix))) > tolerance) return

        allocate(factor%lower(dimension, dimension))
        factor%lower = 0.0_dp
        do row = 1, dimension
            do column = 1, row
                if (row /= column .and. &
                    abs(matrix(row, column)) <= tolerance) cycle
                value = matrix(row, column)
                do inner = 1, column - 1
                    value = value - factor%lower(row, inner)* &
                        factor%lower(column, inner)
                end do
                if (row == column) then
                    if (value <= tolerance) then
                        call clear_factor(factor)
                        status = 2
                        return
                    end if
                    factor%lower(row, column) = sqrt(value)
                else
                    if (factor%lower(column, column) <= 0.0_dp) then
                        call clear_factor(factor)
                        status = 3
                        return
                    end if
                    factor%lower(row, column) = value/ &
                        factor%lower(column, column)
                end if
            end do
        end do
        status = 0
    end subroutine build_incomplete_cholesky

    subroutine apply_incomplete_cholesky( &
            factor, right_hand_side, solution, status)
        type(incomplete_cholesky_factor_t), intent(in) :: factor
        real(dp), intent(in) :: right_hand_side(:)
        real(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: status
        real(dp), allocatable :: intermediate(:)
        integer :: dimension, row, column

        if (allocated(solution)) deallocate(solution)
        status = 1
        if (.not. allocated(factor%lower)) then
            allocate(solution(0))
            return
        end if
        dimension = size(factor%lower, 1)
        if (size(factor%lower, 2) /= dimension .or. &
            size(right_hand_side) /= dimension) then
            allocate(solution(0))
            return
        end if
        allocate(intermediate(dimension), solution(dimension))
        do row = 1, dimension
            intermediate(row) = right_hand_side(row)
            do column = 1, row - 1
                intermediate(row) = intermediate(row) - &
                    factor%lower(row, column)*intermediate(column)
            end do
            if (factor%lower(row, row) <= 0.0_dp) then
                deallocate(intermediate, solution)
                allocate(solution(0))
                return
            end if
            intermediate(row) = intermediate(row)/factor%lower(row, row)
        end do
        do row = dimension, 1, -1
            solution(row) = intermediate(row)
            do column = row + 1, dimension
                solution(row) = solution(row) - factor%lower(column, row)* &
                    solution(column)
            end do
            solution(row) = solution(row)/factor%lower(row, row)
        end do
        status = 0
    end subroutine apply_incomplete_cholesky

    subroutine clear_factor(factor)
        type(incomplete_cholesky_factor_t), intent(inout) :: factor

        if (allocated(factor%lower)) deallocate(factor%lower)
    end subroutine clear_factor

    pure logical function finite_values(values)
        real(dp), intent(in) :: values(:, :)
        integer :: row, column

        finite_values = .true.
        do column = 1, size(values, 2)
            do row = 1, size(values, 1)
                if (.not. (values(row, column) == values(row, column))) then
                    finite_values = .false.
                    return
                end if
                if (abs(values(row, column)) > huge(1.0_dp)) then
                    finite_values = .false.
                    return
                end if
            end do
        end do
    end function finite_values

end module fortfem_incomplete_cholesky
