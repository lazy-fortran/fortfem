module fortfem_vector_jump_penalty
    !! Tensor-weighted vector trace-jump penalty for broken FEEC fields.
    !!
    !! The component metric is supplied by the caller.  It can therefore be a
    !! tangential or normal projector, an anisotropic constitutive tensor, or a
    !! generated H(curl)/H(div) trace metric without this neutral operator
    !! choosing a physical law.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_vector_jump_penalty
    public :: assemble_vector_jump_penalty_jvp
    public :: assemble_vector_jump_penalty_vjp

contains

    subroutine assemble_vector_jump_penalty( &
            plus_trace, minus_trace, surface_weights, penalty, component_metric, &
            matrix, status)
        !! Assemble penalty * integral J M J^T dS for vector traces.
        real(dp), intent(in) :: plus_trace(:, :, :), minus_trace(:, :, :)
        real(dp), intent(in) :: surface_weights(:), penalty
        real(dp), intent(in) :: component_metric(:, :, :)
        real(dp), intent(out) :: matrix(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, plus_count, minus_count, component_count
        integer :: total_count, quadrature, row, column, a, b
        real(dp), allocatable :: jump(:, :)
        real(dp) :: scale

        matrix = 0.0_dp
        call validate_vector_penalty_inputs( &
            plus_trace, minus_trace, surface_weights, penalty, component_metric, &
            matrix, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(plus_trace, 1)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        component_count = size(plus_trace, 3)
        total_count = plus_count + minus_count
        allocate(jump(total_count, component_count))
        do quadrature = 1, quadrature_count
            call fill_jump(plus_trace(quadrature, :, :), &
                minus_trace(quadrature, :, :), jump)
            scale = penalty*surface_weights(quadrature)
            do row = 1, total_count
                do column = 1, total_count
                    do a = 1, component_count
                        do b = 1, component_count
                            matrix(row, column) = matrix(row, column) + &
                                scale*jump(row, a)*component_metric(quadrature, a, b)* &
                                jump(column, b)
                        end do
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_vector_jump_penalty

    subroutine assemble_vector_jump_penalty_jvp( &
            plus_trace, minus_trace, surface_weights, penalty, component_metric, &
            plus_trace_dot, minus_trace_dot, surface_weights_dot, penalty_dot, &
            component_metric_dot, matrix_dot, status)
        !! Apply the product-rule JVP of a tensor-weighted vector penalty.
        real(dp), intent(in) :: plus_trace(:, :, :), minus_trace(:, :, :)
        real(dp), intent(in) :: surface_weights(:), penalty
        real(dp), intent(in) :: component_metric(:, :, :)
        real(dp), intent(in) :: plus_trace_dot(:, :, :), minus_trace_dot(:, :, :)
        real(dp), intent(in) :: surface_weights_dot(:), penalty_dot
        real(dp), intent(in) :: component_metric_dot(:, :, :)
        real(dp), intent(out) :: matrix_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, plus_count, minus_count, component_count
        integer :: total_count, quadrature, row, column, a, b
        real(dp), allocatable :: jump(:, :), jump_dot(:, :)
        real(dp) :: scale, scale_dot, kernel, kernel_dot

        matrix_dot = 0.0_dp
        call validate_vector_penalty_inputs( &
            plus_trace, minus_trace, surface_weights, penalty, component_metric, &
            matrix_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(plus_trace, 1)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        component_count = size(plus_trace, 3)
        total_count = plus_count + minus_count
        if (.not. valid_vector_penalty_direction( &
            plus_trace_dot, minus_trace_dot, surface_weights_dot, penalty_dot, &
            component_metric_dot, quadrature_count, plus_count, minus_count, &
            component_count)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "vector jump penalty JVP has incompatible increments")
            return
        end if
        allocate(jump(total_count, component_count), jump_dot(total_count, component_count))
        do quadrature = 1, quadrature_count
            call fill_jump(plus_trace(quadrature, :, :), &
                minus_trace(quadrature, :, :), jump)
            call fill_jump(plus_trace_dot(quadrature, :, :), &
                minus_trace_dot(quadrature, :, :), jump_dot)
            scale = penalty*surface_weights(quadrature)
            scale_dot = penalty_dot*surface_weights(quadrature) + &
                penalty*surface_weights_dot(quadrature)
            do row = 1, total_count
                do column = 1, total_count
                    kernel = 0.0_dp
                    kernel_dot = 0.0_dp
                    do a = 1, component_count
                        do b = 1, component_count
                            kernel = kernel + jump(row, a)* &
                                component_metric(quadrature, a, b)*jump(column, b)
                            kernel_dot = kernel_dot + &
                                jump_dot(row, a)*component_metric(quadrature, a, b)* &
                                jump(column, b) + &
                                jump(row, a)*component_metric_dot(quadrature, a, b)* &
                                jump(column, b) + &
                                jump(row, a)*component_metric(quadrature, a, b)* &
                                jump_dot(column, b)
                        end do
                    end do
                    matrix_dot(row, column) = matrix_dot(row, column) + &
                        scale_dot*kernel + scale*kernel_dot
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_vector_jump_penalty_jvp

    subroutine assemble_vector_jump_penalty_vjp( &
            plus_trace, minus_trace, surface_weights, penalty, component_metric, &
            matrix_bar, plus_trace_bar, minus_trace_bar, surface_weights_bar, &
            penalty_bar, component_metric_bar, status)
        !! Apply the reverse product of a tensor-weighted vector penalty.
        real(dp), intent(in) :: plus_trace(:, :, :), minus_trace(:, :, :)
        real(dp), intent(in) :: surface_weights(:), penalty
        real(dp), intent(in) :: component_metric(:, :, :), matrix_bar(:, :)
        real(dp), intent(out) :: plus_trace_bar(:, :, :), minus_trace_bar(:, :, :)
        real(dp), intent(out) :: surface_weights_bar(:), penalty_bar
        real(dp), intent(out) :: component_metric_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, plus_count, minus_count, component_count
        integer :: total_count, quadrature, row, column, a, b
        real(dp), allocatable :: jump(:, :), jump_bar(:, :)
        real(dp) :: scale, kernel_contraction

        plus_trace_bar = 0.0_dp
        minus_trace_bar = 0.0_dp
        surface_weights_bar = 0.0_dp
        penalty_bar = 0.0_dp
        component_metric_bar = 0.0_dp
        call validate_vector_penalty_inputs( &
            plus_trace, minus_trace, surface_weights, penalty, component_metric, &
            matrix_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(plus_trace, 1)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        component_count = size(plus_trace, 3)
        total_count = plus_count + minus_count
        if (size(plus_trace_bar, 1) /= quadrature_count .or. &
            size(plus_trace_bar, 2) /= plus_count .or. &
            size(plus_trace_bar, 3) /= component_count .or. &
            size(minus_trace_bar, 1) /= quadrature_count .or. &
            size(minus_trace_bar, 2) /= minus_count .or. &
            size(minus_trace_bar, 3) /= component_count .or. &
            size(surface_weights_bar) /= quadrature_count .or. &
            size(component_metric_bar, 1) /= quadrature_count .or. &
            size(component_metric_bar, 2) /= component_count .or. &
            size(component_metric_bar, 3) /= component_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "vector jump penalty VJP has incompatible cotangents")
            return
        end if

        allocate(jump(total_count, component_count), jump_bar(total_count, component_count))
        do quadrature = 1, quadrature_count
            call fill_jump(plus_trace(quadrature, :, :), &
                minus_trace(quadrature, :, :), jump)
            jump_bar = 0.0_dp
            scale = penalty*surface_weights(quadrature)
            do row = 1, total_count
                do column = 1, total_count
                    do a = 1, component_count
                        do b = 1, component_count
                            jump_bar(row, a) = jump_bar(row, a) + scale*( &
                                matrix_bar(row, column)*component_metric(quadrature, a, b)* &
                                jump(column, b) + matrix_bar(column, row)* &
                                jump(column, b)*component_metric(quadrature, b, a))
                            component_metric_bar(quadrature, a, b) = &
                                component_metric_bar(quadrature, a, b) + scale* &
                                matrix_bar(row, column)*jump(row, a)*jump(column, b)
                        end do
                    end do
                end do
            end do
            plus_trace_bar(quadrature, :, :) = jump_bar(1:plus_count, :)
            minus_trace_bar(quadrature, :, :) = -jump_bar(plus_count + 1:total_count, :)
            kernel_contraction = 0.0_dp
            do row = 1, total_count
                do column = 1, total_count
                    do a = 1, component_count
                        do b = 1, component_count
                            kernel_contraction = kernel_contraction + matrix_bar(row, column)* &
                                jump(row, a)*component_metric(quadrature, a, b)* &
                                jump(column, b)
                        end do
                    end do
                end do
            end do
            surface_weights_bar(quadrature) = penalty*kernel_contraction
            penalty_bar = penalty_bar + surface_weights(quadrature)*kernel_contraction
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_vector_jump_penalty_vjp

    subroutine fill_jump(plus_trace, minus_trace, jump)
        real(dp), intent(in) :: plus_trace(:, :), minus_trace(:, :)
        real(dp), intent(out) :: jump(:, :)
        integer :: plus_count

        plus_count = size(plus_trace, 1)
        jump(1:plus_count, :) = plus_trace
        jump(plus_count + 1:, :) = -minus_trace
    end subroutine fill_jump

    logical function valid_vector_penalty_direction( &
            plus_trace_dot, minus_trace_dot, surface_weights_dot, penalty_dot, &
            component_metric_dot, quadrature_count, plus_count, minus_count, &
            component_count)
        real(dp), intent(in) :: plus_trace_dot(:, :, :), minus_trace_dot(:, :, :)
        real(dp), intent(in) :: surface_weights_dot(:), penalty_dot
        real(dp), intent(in) :: component_metric_dot(:, :, :)
        integer, intent(in) :: quadrature_count, plus_count, minus_count, component_count

        valid_vector_penalty_direction = &
            size(plus_trace_dot, 1) == quadrature_count .and. &
            size(plus_trace_dot, 2) == plus_count .and. &
            size(plus_trace_dot, 3) == component_count .and. &
            size(minus_trace_dot, 1) == quadrature_count .and. &
            size(minus_trace_dot, 2) == minus_count .and. &
            size(minus_trace_dot, 3) == component_count .and. &
            size(surface_weights_dot) == quadrature_count .and. &
            size(component_metric_dot, 1) == quadrature_count .and. &
            size(component_metric_dot, 2) == component_count .and. &
            size(component_metric_dot, 3) == component_count .and. &
            ieee_is_finite(penalty_dot) .and. &
            all(ieee_is_finite(plus_trace_dot)) .and. &
            all(ieee_is_finite(minus_trace_dot)) .and. &
            all(ieee_is_finite(surface_weights_dot)) .and. &
            all(ieee_is_finite(component_metric_dot))
    end function valid_vector_penalty_direction

    subroutine validate_vector_penalty_inputs( &
            plus_trace, minus_trace, surface_weights, penalty, component_metric, &
            matrix, status)
        real(dp), intent(in) :: plus_trace(:, :, :), minus_trace(:, :, :)
        real(dp), intent(in) :: surface_weights(:), penalty
        real(dp), intent(in) :: component_metric(:, :, :), matrix(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, plus_count, minus_count, component_count
        integer :: total_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "vector jump penalty received incompatible arrays")
        quadrature_count = size(plus_trace, 1)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        component_count = size(plus_trace, 3)
        total_count = plus_count + minus_count
        if (quadrature_count < 1 .or. plus_count < 1 .or. minus_count < 1 .or. &
            component_count < 1) return
        if (size(minus_trace, 1) /= quadrature_count .or. &
            size(minus_trace, 3) /= component_count .or. &
            size(surface_weights) /= quadrature_count .or. &
            size(component_metric, 1) /= quadrature_count .or. &
            size(component_metric, 2) /= component_count .or. &
            size(component_metric, 3) /= component_count .or. &
            size(matrix, 1) /= total_count .or. size(matrix, 2) /= total_count) return
        if (.not. ieee_is_finite(penalty) .or. penalty < 0.0_dp .or. &
            any(surface_weights <= 0.0_dp) .or. &
            any(.not. ieee_is_finite(plus_trace)) .or. &
            any(.not. ieee_is_finite(minus_trace)) .or. &
            any(.not. ieee_is_finite(surface_weights)) .or. &
            any(.not. ieee_is_finite(component_metric)) .or. &
            any(.not. ieee_is_finite(matrix))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_vector_penalty_inputs

end module fortfem_vector_jump_penalty
