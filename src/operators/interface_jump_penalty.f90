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
    public :: assemble_interface_jump_penalty_jvp
    public :: assemble_interface_jump_penalty_vjp

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

        matrix = 0.0_dp
        call validate_jump_penalty_inputs( &
            plus_trace, minus_trace, surface_weights, penalty, matrix, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(plus_trace, 1)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        total_count = plus_count + minus_count

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

    subroutine assemble_interface_jump_penalty_jvp( &
            plus_trace, minus_trace, surface_weights, penalty, plus_trace_dot, &
            minus_trace_dot, surface_weights_dot, penalty_dot, matrix_dot, status)
        !! Apply the product-rule JVP of the symmetric jump penalty block.
        real(dp), intent(in) :: plus_trace(:, :), minus_trace(:, :)
        real(dp), intent(in) :: surface_weights(:), penalty
        real(dp), intent(in) :: plus_trace_dot(:, :), minus_trace_dot(:, :)
        real(dp), intent(in) :: surface_weights_dot(:), penalty_dot
        real(dp), intent(out) :: matrix_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, plus_count, minus_count, total_count
        integer :: quadrature, row, column
        real(dp), allocatable :: jump_trace(:), jump_trace_dot(:)
        real(dp) :: scale, scale_dot

        matrix_dot = 0.0_dp
        call validate_jump_penalty_inputs( &
            plus_trace, minus_trace, surface_weights, penalty, matrix_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(plus_trace, 1)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        total_count = plus_count + minus_count
        if (size(plus_trace_dot, 1) /= quadrature_count .or. &
            size(plus_trace_dot, 2) /= plus_count .or. &
            size(minus_trace_dot, 1) /= quadrature_count .or. &
            size(minus_trace_dot, 2) /= minus_count .or. &
            size(surface_weights_dot) /= quadrature_count .or. &
            .not. ieee_is_finite(penalty_dot) .or. &
            any(.not. ieee_is_finite(plus_trace_dot)) .or. &
            any(.not. ieee_is_finite(minus_trace_dot)) .or. &
            any(.not. ieee_is_finite(surface_weights_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "interface jump penalty JVP has incompatible increments")
            return
        end if
        allocate(jump_trace(total_count), jump_trace_dot(total_count))
        do quadrature = 1, quadrature_count
            jump_trace(1:plus_count) = plus_trace(quadrature, :)
            jump_trace(plus_count + 1:total_count) = -minus_trace(quadrature, :)
            jump_trace_dot(1:plus_count) = plus_trace_dot(quadrature, :)
            jump_trace_dot(plus_count + 1:total_count) = &
                -minus_trace_dot(quadrature, :)
            scale = penalty*surface_weights(quadrature)
            scale_dot = penalty_dot*surface_weights(quadrature) + &
                penalty*surface_weights_dot(quadrature)
            do row = 1, total_count
                do column = 1, total_count
                    matrix_dot(row, column) = matrix_dot(row, column) + &
                        scale_dot*jump_trace(row)*jump_trace(column) + &
                        scale*(jump_trace_dot(row)*jump_trace(column) + &
                        jump_trace(row)*jump_trace_dot(column))
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_interface_jump_penalty_jvp

    subroutine assemble_interface_jump_penalty_vjp( &
            plus_trace, minus_trace, surface_weights, penalty, matrix_bar, &
            plus_trace_bar, minus_trace_bar, surface_weights_bar, penalty_bar, &
            status)
        !! Apply the reverse product of the symmetric jump penalty block.
        real(dp), intent(in) :: plus_trace(:, :), minus_trace(:, :)
        real(dp), intent(in) :: surface_weights(:), penalty, matrix_bar(:, :)
        real(dp), intent(out) :: plus_trace_bar(:, :), minus_trace_bar(:, :)
        real(dp), intent(out) :: surface_weights_bar(:), penalty_bar
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, plus_count, minus_count, total_count
        integer :: quadrature, row, column
        real(dp), allocatable :: jump_trace(:), jump_bar(:)
        real(dp) :: matrix_bar_quadratic

        plus_trace_bar = 0.0_dp
        minus_trace_bar = 0.0_dp
        surface_weights_bar = 0.0_dp
        penalty_bar = 0.0_dp
        call validate_jump_penalty_inputs( &
            plus_trace, minus_trace, surface_weights, penalty, matrix_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(plus_trace, 1)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        total_count = plus_count + minus_count
        if (size(plus_trace_bar, 1) /= quadrature_count .or. &
            size(plus_trace_bar, 2) /= plus_count .or. &
            size(minus_trace_bar, 1) /= quadrature_count .or. &
            size(minus_trace_bar, 2) /= minus_count .or. &
            size(surface_weights_bar) /= quadrature_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "interface jump penalty VJP has incompatible cotangents")
            return
        end if
        allocate(jump_trace(total_count), jump_bar(total_count))
        do quadrature = 1, quadrature_count
            jump_trace(1:plus_count) = plus_trace(quadrature, :)
            jump_trace(plus_count + 1:total_count) = -minus_trace(quadrature, :)
            jump_bar = 0.0_dp
            do row = 1, total_count
                do column = 1, total_count
                    jump_bar(row) = jump_bar(row) + &
                        (matrix_bar(row, column) + matrix_bar(column, row))* &
                        jump_trace(column)
                end do
            end do
            jump_bar = penalty*surface_weights(quadrature)*jump_bar
            plus_trace_bar(quadrature, :) = jump_bar(1:plus_count)
            minus_trace_bar(quadrature, :) = -jump_bar(plus_count + 1:total_count)
            matrix_bar_quadratic = 0.0_dp
            do row = 1, total_count
                do column = 1, total_count
                    matrix_bar_quadratic = matrix_bar_quadratic + &
                        matrix_bar(row, column)*jump_trace(row)*jump_trace(column)
                end do
            end do
            surface_weights_bar(quadrature) = penalty*matrix_bar_quadratic
            penalty_bar = penalty_bar + &
                surface_weights(quadrature)*matrix_bar_quadratic
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_interface_jump_penalty_vjp

    subroutine validate_jump_penalty_inputs( &
            plus_trace, minus_trace, surface_weights, penalty, matrix, status)
        real(dp), intent(in) :: plus_trace(:, :), minus_trace(:, :)
        real(dp), intent(in) :: surface_weights(:), penalty, matrix(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: quadrature_count, plus_count, minus_count, total_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "interface jump penalty received incompatible arrays")
        quadrature_count = size(plus_trace, 1)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        total_count = plus_count + minus_count
        if (quadrature_count < 1 .or. plus_count < 1 .or. minus_count < 1) return
        if (size(minus_trace, 1) /= quadrature_count .or. &
            size(surface_weights) /= quadrature_count .or. &
            size(matrix, 1) /= total_count .or. size(matrix, 2) /= total_count) return
        if (.not. ieee_is_finite(penalty) .or. penalty < 0.0_dp .or. &
            any(.not. ieee_is_finite(plus_trace)) .or. &
            any(.not. ieee_is_finite(minus_trace)) .or. &
            any(.not. ieee_is_finite(surface_weights)) .or. &
            any(surface_weights <= 0.0_dp) .or. &
            any(.not. ieee_is_finite(matrix))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_jump_penalty_inputs

end module fortfem_interface_jump_penalty
