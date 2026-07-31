module fortfem_interface_jump_penalty
    !! Symmetric trace-jump penalty block for fitted interfaces and DG.
    !!
    !! With plus-minus orientation, the block is
    !! A = eta * integral [T_plus, -T_minus]^T [T_plus, -T_minus] dS.
    !! It is therefore symmetric positive semidefinite for eta >= 0 and can
    !! be used as the penalty portion of Nitsche or interior-penalty forms.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_interface_jump_penalty

contains

    subroutine assemble_interface_jump_penalty( &
            plus_trace, minus_trace, surface_weights, penalty, matrix, status)
        !! Assemble a scalar trace-jump penalty block.
        real(dp), intent(in) :: plus_trace(:, :), minus_trace(:, :)
        real(dp), intent(in) :: surface_weights(:)
        real(dp), intent(in) :: penalty
        real(dp), intent(out) :: matrix(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, plus_count, minus_count, total_count
        integer :: quadrature, row, column
        real(dp), allocatable :: jump_trace(:)
        real(dp) :: scale

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "interface jump penalty received incompatible arrays")
        matrix = 0.0_dp
        quadrature_count = size(plus_trace, 1)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        total_count = plus_count + minus_count
        if (quadrature_count < 1 .or. plus_count < 1 .or. minus_count < 1) return
        if (size(minus_trace, 1) /= quadrature_count .or. &
            size(surface_weights) /= quadrature_count .or. &
            size(matrix, 1) /= total_count .or. &
            size(matrix, 2) /= total_count) return
        if (.not. ieee_is_finite(penalty) .or. penalty < 0.0_dp) return
        if (any(.not. ieee_is_finite(plus_trace)) .or. &
            any(.not. ieee_is_finite(minus_trace)) .or. &
            any(.not. ieee_is_finite(surface_weights))) return
        if (any(surface_weights <= 0.0_dp)) return

        allocate(jump_trace(total_count))
        do quadrature = 1, quadrature_count
            jump_trace(1:plus_count) = plus_trace(quadrature, :)
            jump_trace(plus_count + 1:total_count) = &
                -minus_trace(quadrature, :)
            scale = penalty*surface_weights(quadrature)
            do row = 1, total_count
                do column = 1, total_count
                    matrix(row, column) = matrix(row, column) + &
                        scale*jump_trace(row)*jump_trace(column)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_interface_jump_penalty

end module fortfem_interface_jump_penalty
