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
    public :: assemble_vector_sipg_interface
    public :: assemble_vector_sipg_interface_jvp
    public :: assemble_vector_sipg_interface_vjp

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

    subroutine assemble_vector_sipg_interface( &
            test_plus_trace, test_minus_trace, trial_plus_trace, &
            trial_minus_trace, test_plus_flux, test_minus_flux, &
            trial_plus_flux, trial_minus_flux, surface_weights, penalty, &
            consistency_sign, component_metric, matrix, status)
        !! Assemble a tensor-weighted vector SIPG interface block.
        real(dp), intent(in) :: test_plus_trace(:, :, :), test_minus_trace(:, :, :)
        real(dp), intent(in) :: trial_plus_trace(:, :, :), trial_minus_trace(:, :, :)
        real(dp), intent(in) :: test_plus_flux(:, :, :), test_minus_flux(:, :, :)
        real(dp), intent(in) :: trial_plus_flux(:, :, :), trial_minus_flux(:, :, :)
        real(dp), intent(in) :: surface_weights(:), penalty
        integer, intent(in) :: consistency_sign
        real(dp), intent(in) :: component_metric(:, :, :)
        real(dp), intent(out) :: matrix(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_plus_count, test_minus_count
        integer :: trial_plus_count, trial_minus_count, component_count
        integer :: test_count, trial_count, quadrature
        real(dp), allocatable :: test_jump(:, :), trial_jump(:, :)
        real(dp), allocatable :: test_average(:, :), trial_average(:, :)

        matrix = 0.0_dp
        call validate_vector_sipg_inputs( &
            test_plus_trace, test_minus_trace, trial_plus_trace, trial_minus_trace, &
            test_plus_flux, test_minus_flux, trial_plus_flux, trial_minus_flux, &
            surface_weights, penalty, consistency_sign, component_metric, matrix, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(test_plus_trace, 1)
        test_plus_count = size(test_plus_trace, 2)
        test_minus_count = size(test_minus_trace, 2)
        trial_plus_count = size(trial_plus_trace, 2)
        trial_minus_count = size(trial_minus_trace, 2)
        component_count = size(test_plus_trace, 3)
        test_count = test_plus_count + test_minus_count
        trial_count = trial_plus_count + trial_minus_count
        allocate(test_jump(test_count, component_count), &
            trial_jump(trial_count, component_count), &
            test_average(test_count, component_count), &
            trial_average(trial_count, component_count))
        do quadrature = 1, quadrature_count
            call fill_vector_jump_average( &
                test_plus_trace(quadrature, :, :), test_minus_trace(quadrature, :, :), &
                test_plus_count, test_jump, test_average, &
                test_plus_flux(quadrature, :, :), test_minus_flux(quadrature, :, :))
            call fill_vector_jump_average( &
                trial_plus_trace(quadrature, :, :), trial_minus_trace(quadrature, :, :), &
                trial_plus_count, trial_jump, trial_average, &
                trial_plus_flux(quadrature, :, :), trial_minus_flux(quadrature, :, :))
            call add_vector_sipg_kernel(test_jump, trial_jump, test_average, &
                trial_average, surface_weights(quadrature), penalty, consistency_sign, &
                component_metric(quadrature, :, :), matrix)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_vector_sipg_interface

    subroutine assemble_vector_sipg_interface_jvp( &
            test_plus_trace, test_minus_trace, trial_plus_trace, &
            trial_minus_trace, test_plus_flux, test_minus_flux, &
            trial_plus_flux, trial_minus_flux, surface_weights, penalty, &
            consistency_sign, component_metric, test_plus_trace_dot, &
            test_minus_trace_dot, trial_plus_trace_dot, trial_minus_trace_dot, &
            test_plus_flux_dot, test_minus_flux_dot, trial_plus_flux_dot, &
            trial_minus_flux_dot, surface_weights_dot, penalty_dot, &
            component_metric_dot, matrix_dot, status)
        !! Apply the product-rule JVP of a tensor-weighted vector SIPG block.
        real(dp), intent(in) :: test_plus_trace(:, :, :), test_minus_trace(:, :, :)
        real(dp), intent(in) :: trial_plus_trace(:, :, :), trial_minus_trace(:, :, :)
        real(dp), intent(in) :: test_plus_flux(:, :, :), test_minus_flux(:, :, :)
        real(dp), intent(in) :: trial_plus_flux(:, :, :), trial_minus_flux(:, :, :)
        real(dp), intent(in) :: surface_weights(:), penalty
        integer, intent(in) :: consistency_sign
        real(dp), intent(in) :: component_metric(:, :, :)
        real(dp), intent(in) :: test_plus_trace_dot(:, :, :), test_minus_trace_dot(:, :, :)
        real(dp), intent(in) :: trial_plus_trace_dot(:, :, :), trial_minus_trace_dot(:, :, :)
        real(dp), intent(in) :: test_plus_flux_dot(:, :, :), test_minus_flux_dot(:, :, :)
        real(dp), intent(in) :: trial_plus_flux_dot(:, :, :), trial_minus_flux_dot(:, :, :)
        real(dp), intent(in) :: surface_weights_dot(:), penalty_dot
        real(dp), intent(in) :: component_metric_dot(:, :, :)
        real(dp), intent(out) :: matrix_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_plus_count, test_minus_count
        integer :: trial_plus_count, trial_minus_count, component_count
        integer :: test_count, trial_count, quadrature
        real(dp), allocatable :: test_jump(:, :), trial_jump(:, :)
        real(dp), allocatable :: test_average(:, :), trial_average(:, :)
        real(dp), allocatable :: test_jump_dot(:, :), trial_jump_dot(:, :)
        real(dp), allocatable :: test_average_dot(:, :), trial_average_dot(:, :)

        matrix_dot = 0.0_dp
        call validate_vector_sipg_inputs( &
            test_plus_trace, test_minus_trace, trial_plus_trace, trial_minus_trace, &
            test_plus_flux, test_minus_flux, trial_plus_flux, trial_minus_flux, &
            surface_weights, penalty, consistency_sign, component_metric, matrix_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(test_plus_trace, 1)
        test_plus_count = size(test_plus_trace, 2)
        test_minus_count = size(test_minus_trace, 2)
        trial_plus_count = size(trial_plus_trace, 2)
        trial_minus_count = size(trial_minus_trace, 2)
        component_count = size(test_plus_trace, 3)
        test_count = test_plus_count + test_minus_count
        trial_count = trial_plus_count + trial_minus_count
        if (.not. valid_vector_sipg_direction( &
            test_plus_trace_dot, test_minus_trace_dot, trial_plus_trace_dot, &
            trial_minus_trace_dot, test_plus_flux_dot, test_minus_flux_dot, &
            trial_plus_flux_dot, trial_minus_flux_dot, surface_weights_dot, &
            penalty_dot, component_metric_dot, quadrature_count, test_plus_count, &
            test_minus_count, trial_plus_count, trial_minus_count, component_count)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "vector SIPG interface JVP has incompatible increments")
            return
        end if
        allocate(test_jump(test_count, component_count), &
            trial_jump(trial_count, component_count), &
            test_average(test_count, component_count), &
            trial_average(trial_count, component_count), &
            test_jump_dot(test_count, component_count), &
            trial_jump_dot(trial_count, component_count), &
            test_average_dot(test_count, component_count), &
            trial_average_dot(trial_count, component_count))
        do quadrature = 1, quadrature_count
            call fill_vector_jump_average( &
                test_plus_trace(quadrature, :, :), test_minus_trace(quadrature, :, :), &
                test_plus_count, test_jump, test_average, &
                test_plus_flux(quadrature, :, :), test_minus_flux(quadrature, :, :))
            call fill_vector_jump_average( &
                trial_plus_trace(quadrature, :, :), trial_minus_trace(quadrature, :, :), &
                trial_plus_count, trial_jump, trial_average, &
                trial_plus_flux(quadrature, :, :), trial_minus_flux(quadrature, :, :))
            call fill_vector_jump_average( &
                test_plus_trace_dot(quadrature, :, :), test_minus_trace_dot(quadrature, :, :), &
                test_plus_count, test_jump_dot, test_average_dot, &
                test_plus_flux_dot(quadrature, :, :), test_minus_flux_dot(quadrature, :, :))
            call fill_vector_jump_average( &
                trial_plus_trace_dot(quadrature, :, :), trial_minus_trace_dot(quadrature, :, :), &
                trial_plus_count, trial_jump_dot, trial_average_dot, &
                trial_plus_flux_dot(quadrature, :, :), trial_minus_flux_dot(quadrature, :, :))
            call add_vector_sipg_kernel_jvp( &
                test_jump, trial_jump, test_average, trial_average, test_jump_dot, &
                trial_jump_dot, test_average_dot, trial_average_dot, &
                surface_weights(quadrature), surface_weights_dot(quadrature), penalty, &
                penalty_dot, consistency_sign, component_metric(quadrature, :, :), &
                component_metric_dot(quadrature, :, :), matrix_dot)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_vector_sipg_interface_jvp

    subroutine assemble_vector_sipg_interface_vjp( &
            test_plus_trace, test_minus_trace, trial_plus_trace, &
            trial_minus_trace, test_plus_flux, test_minus_flux, &
            trial_plus_flux, trial_minus_flux, surface_weights, penalty, &
            consistency_sign, component_metric, matrix_bar, test_plus_trace_bar, &
            test_minus_trace_bar, trial_plus_trace_bar, trial_minus_trace_bar, &
            test_plus_flux_bar, test_minus_flux_bar, trial_plus_flux_bar, &
            trial_minus_flux_bar, surface_weights_bar, penalty_bar, &
            component_metric_bar, status)
        !! Apply the reverse product of a tensor-weighted vector SIPG block.
        real(dp), intent(in) :: test_plus_trace(:, :, :), test_minus_trace(:, :, :)
        real(dp), intent(in) :: trial_plus_trace(:, :, :), trial_minus_trace(:, :, :)
        real(dp), intent(in) :: test_plus_flux(:, :, :), test_minus_flux(:, :, :)
        real(dp), intent(in) :: trial_plus_flux(:, :, :), trial_minus_flux(:, :, :)
        real(dp), intent(in) :: surface_weights(:), penalty
        integer, intent(in) :: consistency_sign
        real(dp), intent(in) :: component_metric(:, :, :), matrix_bar(:, :)
        real(dp), intent(out) :: test_plus_trace_bar(:, :, :), test_minus_trace_bar(:, :, :)
        real(dp), intent(out) :: trial_plus_trace_bar(:, :, :), trial_minus_trace_bar(:, :, :)
        real(dp), intent(out) :: test_plus_flux_bar(:, :, :), test_minus_flux_bar(:, :, :)
        real(dp), intent(out) :: trial_plus_flux_bar(:, :, :), trial_minus_flux_bar(:, :, :)
        real(dp), intent(out) :: surface_weights_bar(:), penalty_bar
        real(dp), intent(out) :: component_metric_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_plus_count, test_minus_count
        integer :: trial_plus_count, trial_minus_count, component_count
        integer :: test_count, trial_count, quadrature, row, column, a, b
        real(dp), allocatable :: test_jump(:, :), trial_jump(:, :)
        real(dp), allocatable :: test_average(:, :), trial_average(:, :)
        real(dp), allocatable :: test_jump_bar(:, :), trial_jump_bar(:, :)
        real(dp), allocatable :: test_average_bar(:, :), trial_average_bar(:, :)
        real(dp) :: scale, kernel_contraction, penalty_contraction

        test_plus_trace_bar = 0.0_dp
        test_minus_trace_bar = 0.0_dp
        trial_plus_trace_bar = 0.0_dp
        trial_minus_trace_bar = 0.0_dp
        test_plus_flux_bar = 0.0_dp
        test_minus_flux_bar = 0.0_dp
        trial_plus_flux_bar = 0.0_dp
        trial_minus_flux_bar = 0.0_dp
        surface_weights_bar = 0.0_dp
        penalty_bar = 0.0_dp
        component_metric_bar = 0.0_dp
        call validate_vector_sipg_inputs( &
            test_plus_trace, test_minus_trace, trial_plus_trace, trial_minus_trace, &
            test_plus_flux, test_minus_flux, trial_plus_flux, trial_minus_flux, &
            surface_weights, penalty, consistency_sign, component_metric, matrix_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(test_plus_trace, 1)
        test_plus_count = size(test_plus_trace, 2)
        test_minus_count = size(test_minus_trace, 2)
        trial_plus_count = size(trial_plus_trace, 2)
        trial_minus_count = size(trial_minus_trace, 2)
        component_count = size(test_plus_trace, 3)
        test_count = test_plus_count + test_minus_count
        trial_count = trial_plus_count + trial_minus_count
        if (.not. valid_vector_sipg_cotangent( &
            test_plus_trace_bar, test_minus_trace_bar, trial_plus_trace_bar, &
            trial_minus_trace_bar, test_plus_flux_bar, test_minus_flux_bar, &
            trial_plus_flux_bar, trial_minus_flux_bar, surface_weights_bar, &
            component_metric_bar, quadrature_count, test_plus_count, test_minus_count, &
            trial_plus_count, trial_minus_count, component_count)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "vector SIPG interface VJP has incompatible cotangents")
            return
        end if
        allocate(test_jump(test_count, component_count), &
            trial_jump(trial_count, component_count), &
            test_average(test_count, component_count), &
            trial_average(trial_count, component_count), &
            test_jump_bar(test_count, component_count), &
            trial_jump_bar(trial_count, component_count), &
            test_average_bar(test_count, component_count), &
            trial_average_bar(trial_count, component_count))
        do quadrature = 1, quadrature_count
            call fill_vector_jump_average( &
                test_plus_trace(quadrature, :, :), test_minus_trace(quadrature, :, :), &
                test_plus_count, test_jump, test_average, &
                test_plus_flux(quadrature, :, :), test_minus_flux(quadrature, :, :))
            call fill_vector_jump_average( &
                trial_plus_trace(quadrature, :, :), trial_minus_trace(quadrature, :, :), &
                trial_plus_count, trial_jump, trial_average, &
                trial_plus_flux(quadrature, :, :), trial_minus_flux(quadrature, :, :))
            test_jump_bar = 0.0_dp
            trial_jump_bar = 0.0_dp
            test_average_bar = 0.0_dp
            trial_average_bar = 0.0_dp
            scale = surface_weights(quadrature)
            do row = 1, test_count
                do column = 1, trial_count
                    do a = 1, component_count
                        do b = 1, component_count
                            test_jump_bar(row, a) = test_jump_bar(row, a) + scale*( &
                                -consistency_sign*matrix_bar(row, column)* &
                                component_metric(quadrature, a, b)*trial_average(column, b) + &
                                penalty*matrix_bar(row, column)*component_metric( &
                                quadrature, a, b)*trial_jump(column, b))
                            test_average_bar(row, a) = test_average_bar(row, a) - scale* &
                                matrix_bar(row, column)*component_metric(quadrature, a, b)* &
                                trial_jump(column, b)
                            trial_jump_bar(column, b) = trial_jump_bar(column, b) + scale*( &
                                -matrix_bar(row, column)*test_average(row, a)* &
                                component_metric(quadrature, a, b) + penalty*matrix_bar( &
                                row, column)*test_jump(row, a)*component_metric( &
                                quadrature, a, b))
                            trial_average_bar(column, b) = trial_average_bar(column, b) - &
                                scale*consistency_sign*matrix_bar(row, column)*test_jump(row, a)* &
                                component_metric(quadrature, a, b)
                            component_metric_bar(quadrature, a, b) = &
                                component_metric_bar(quadrature, a, b) + scale*matrix_bar( &
                                row, column)*(-test_average(row, a)*trial_jump(column, b) - &
                                consistency_sign*test_jump(row, a)*trial_average(column, b) + &
                                penalty*test_jump(row, a)*trial_jump(column, b))
                        end do
                    end do
                end do
            end do
            test_plus_trace_bar(quadrature, :, :) = test_jump_bar(1:test_plus_count, :)
            test_minus_trace_bar(quadrature, :, :) = -test_jump_bar( &
                test_plus_count + 1:test_count, :)
            trial_plus_trace_bar(quadrature, :, :) = trial_jump_bar(1:trial_plus_count, :)
            trial_minus_trace_bar(quadrature, :, :) = -trial_jump_bar( &
                trial_plus_count + 1:trial_count, :)
            test_plus_flux_bar(quadrature, :, :) = 0.5_dp*test_average_bar(1:test_plus_count, :)
            test_minus_flux_bar(quadrature, :, :) = 0.5_dp*test_average_bar( &
                test_plus_count + 1:test_count, :)
            trial_plus_flux_bar(quadrature, :, :) = 0.5_dp*trial_average_bar(1:trial_plus_count, :)
            trial_minus_flux_bar(quadrature, :, :) = 0.5_dp*trial_average_bar( &
                trial_plus_count + 1:trial_count, :)
            kernel_contraction = 0.0_dp
            penalty_contraction = 0.0_dp
            do row = 1, test_count
                do column = 1, trial_count
                    do a = 1, component_count
                        do b = 1, component_count
                            kernel_contraction = kernel_contraction + matrix_bar(row, column)*( &
                                -test_average(row, a)*component_metric(quadrature, a, b)* &
                                trial_jump(column, b) - consistency_sign*test_jump(row, a)* &
                                component_metric(quadrature, a, b)*trial_average(column, b))
                            penalty_contraction = penalty_contraction + matrix_bar(row, column)* &
                                test_jump(row, a)*component_metric(quadrature, a, b)* &
                                trial_jump(column, b)
                        end do
                    end do
                end do
            end do
            surface_weights_bar(quadrature) = kernel_contraction + penalty*penalty_contraction
            penalty_bar = penalty_bar + surface_weights(quadrature)*penalty_contraction
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_vector_sipg_interface_vjp

    subroutine fill_vector_jump_average( &
            plus_trace, minus_trace, plus_count, jump, average, plus_flux, minus_flux)
        real(dp), intent(in) :: plus_trace(:, :), minus_trace(:, :)
        integer, intent(in) :: plus_count
        real(dp), intent(out) :: jump(:, :), average(:, :)
        real(dp), intent(in) :: plus_flux(:, :), minus_flux(:, :)
        integer :: minus_count

        minus_count = size(minus_trace, 1)
        jump(1:plus_count, :) = plus_trace
        jump(plus_count + 1:plus_count + minus_count, :) = -minus_trace
        average(1:plus_count, :) = 0.5_dp*plus_flux
        average(plus_count + 1:plus_count + minus_count, :) = 0.5_dp*minus_flux
    end subroutine fill_vector_jump_average

    subroutine add_vector_sipg_kernel( &
            test_jump, trial_jump, test_average, trial_average, scale, penalty, &
            consistency_sign, component_metric, matrix)
        real(dp), intent(in) :: test_jump(:, :), trial_jump(:, :), test_average(:, :)
        real(dp), intent(in) :: trial_average(:, :), scale, penalty
        integer, intent(in) :: consistency_sign
        real(dp), intent(in) :: component_metric(:, :)
        real(dp), intent(inout) :: matrix(:, :)
        integer :: row, column, a, b

        do row = 1, size(test_jump, 1)
            do column = 1, size(trial_jump, 1)
                do a = 1, size(test_jump, 2)
                    do b = 1, size(test_jump, 2)
                        matrix(row, column) = matrix(row, column) + scale*( &
                            -test_average(row, a)*component_metric(a, b)*trial_jump(column, b) - &
                            consistency_sign*test_jump(row, a)*component_metric(a, b)* &
                            trial_average(column, b) + penalty*test_jump(row, a)* &
                            component_metric(a, b)*trial_jump(column, b))
                    end do
                end do
            end do
        end do
    end subroutine add_vector_sipg_kernel

    subroutine add_vector_sipg_kernel_jvp( &
            test_jump, trial_jump, test_average, trial_average, test_jump_dot, &
            trial_jump_dot, test_average_dot, trial_average_dot, scale, scale_dot, &
            penalty, penalty_dot, consistency_sign, component_metric, &
            component_metric_dot, matrix_dot)
        real(dp), intent(in) :: test_jump(:, :), trial_jump(:, :), test_average(:, :)
        real(dp), intent(in) :: trial_average(:, :), test_jump_dot(:, :)
        real(dp), intent(in) :: trial_jump_dot(:, :), test_average_dot(:, :)
        real(dp), intent(in) :: trial_average_dot(:, :), scale, scale_dot, penalty
        real(dp), intent(in) :: penalty_dot, component_metric(:, :)
        real(dp), intent(in) :: component_metric_dot(:, :)
        integer, intent(in) :: consistency_sign
        real(dp), intent(inout) :: matrix_dot(:, :)
        integer :: row, column, a, b
        real(dp) :: kernel, kernel_dot

        do row = 1, size(test_jump, 1)
            do column = 1, size(trial_jump, 1)
                kernel = 0.0_dp
                kernel_dot = 0.0_dp
                do a = 1, size(test_jump, 2)
                    do b = 1, size(test_jump, 2)
                        kernel = kernel - test_average(row, a)*component_metric(a, b)* &
                            trial_jump(column, b) - consistency_sign*test_jump(row, a)* &
                            component_metric(a, b)*trial_average(column, b) + penalty* &
                            test_jump(row, a)*component_metric(a, b)*trial_jump(column, b)
                        kernel_dot = kernel_dot - test_average_dot(row, a)*component_metric(a, b)* &
                            trial_jump(column, b) - test_average(row, a)*component_metric_dot(a, b)* &
                            trial_jump(column, b) - test_average(row, a)*component_metric(a, b)* &
                            trial_jump_dot(column, b) - consistency_sign*(test_jump_dot(row, a)* &
                            component_metric(a, b)*trial_average(column, b) + test_jump(row, a)* &
                            component_metric_dot(a, b)*trial_average(column, b) + test_jump(row, a)* &
                            component_metric(a, b)*trial_average_dot(column, b)) + penalty_dot* &
                            test_jump(row, a)*component_metric(a, b)*trial_jump(column, b) + &
                            penalty*test_jump_dot(row, a)*component_metric(a, b)*trial_jump(column, b) + &
                            penalty*test_jump(row, a)*component_metric_dot(a, b)*trial_jump(column, b) + &
                            penalty*test_jump(row, a)*component_metric(a, b)*trial_jump_dot(column, b)
                    end do
                end do
                matrix_dot(row, column) = matrix_dot(row, column) + &
                    scale_dot*kernel + scale*kernel_dot
            end do
        end do
    end subroutine add_vector_sipg_kernel_jvp

    logical function valid_vector_sipg_direction( &
            test_plus_trace_dot, test_minus_trace_dot, trial_plus_trace_dot, &
            trial_minus_trace_dot, test_plus_flux_dot, test_minus_flux_dot, &
            trial_plus_flux_dot, trial_minus_flux_dot, surface_weights_dot, &
            penalty_dot, component_metric_dot, quadrature_count, test_plus_count, &
            test_minus_count, trial_plus_count, trial_minus_count, component_count)
        real(dp), intent(in) :: test_plus_trace_dot(:, :, :), test_minus_trace_dot(:, :, :)
        real(dp), intent(in) :: trial_plus_trace_dot(:, :, :), trial_minus_trace_dot(:, :, :)
        real(dp), intent(in) :: test_plus_flux_dot(:, :, :), test_minus_flux_dot(:, :, :)
        real(dp), intent(in) :: trial_plus_flux_dot(:, :, :), trial_minus_flux_dot(:, :, :)
        real(dp), intent(in) :: surface_weights_dot(:), penalty_dot
        real(dp), intent(in) :: component_metric_dot(:, :, :)
        integer, intent(in) :: quadrature_count, test_plus_count, test_minus_count
        integer, intent(in) :: trial_plus_count, trial_minus_count, component_count

        valid_vector_sipg_direction = &
            size(test_plus_trace_dot, 1) == quadrature_count .and. &
            size(test_plus_trace_dot, 2) == test_plus_count .and. &
            size(test_plus_trace_dot, 3) == component_count .and. &
            size(test_minus_trace_dot, 1) == quadrature_count .and. &
            size(test_minus_trace_dot, 2) == test_minus_count .and. &
            size(test_minus_trace_dot, 3) == component_count .and. &
            size(trial_plus_trace_dot, 1) == quadrature_count .and. &
            size(trial_plus_trace_dot, 2) == trial_plus_count .and. &
            size(trial_plus_trace_dot, 3) == component_count .and. &
            size(trial_minus_trace_dot, 1) == quadrature_count .and. &
            size(trial_minus_trace_dot, 2) == trial_minus_count .and. &
            size(trial_minus_trace_dot, 3) == component_count .and. &
            size(test_plus_flux_dot, 1) == quadrature_count .and. &
            size(test_plus_flux_dot, 2) == test_plus_count .and. &
            size(test_plus_flux_dot, 3) == component_count .and. &
            size(test_minus_flux_dot, 1) == quadrature_count .and. &
            size(test_minus_flux_dot, 2) == test_minus_count .and. &
            size(test_minus_flux_dot, 3) == component_count .and. &
            size(trial_plus_flux_dot, 1) == quadrature_count .and. &
            size(trial_plus_flux_dot, 2) == trial_plus_count .and. &
            size(trial_plus_flux_dot, 3) == component_count .and. &
            size(trial_minus_flux_dot, 1) == quadrature_count .and. &
            size(trial_minus_flux_dot, 2) == trial_minus_count .and. &
            size(trial_minus_flux_dot, 3) == component_count .and. &
            size(surface_weights_dot) == quadrature_count .and. &
            size(component_metric_dot, 1) == quadrature_count .and. &
            size(component_metric_dot, 2) == component_count .and. &
            size(component_metric_dot, 3) == component_count .and. &
            ieee_is_finite(penalty_dot) .and. all(ieee_is_finite(test_plus_trace_dot)) .and. &
            all(ieee_is_finite(test_minus_trace_dot)) .and. all(ieee_is_finite(trial_plus_trace_dot)) .and. &
            all(ieee_is_finite(trial_minus_trace_dot)) .and. all(ieee_is_finite(test_plus_flux_dot)) .and. &
            all(ieee_is_finite(test_minus_flux_dot)) .and. all(ieee_is_finite(trial_plus_flux_dot)) .and. &
            all(ieee_is_finite(trial_minus_flux_dot)) .and. all(ieee_is_finite(surface_weights_dot)) .and. &
            all(ieee_is_finite(component_metric_dot))
    end function valid_vector_sipg_direction

    logical function valid_vector_sipg_cotangent( &
            test_plus_trace_bar, test_minus_trace_bar, trial_plus_trace_bar, &
            trial_minus_trace_bar, test_plus_flux_bar, test_minus_flux_bar, &
            trial_plus_flux_bar, trial_minus_flux_bar, surface_weights_bar, &
            component_metric_bar, quadrature_count, test_plus_count, test_minus_count, &
            trial_plus_count, trial_minus_count, component_count)
        real(dp), intent(in) :: test_plus_trace_bar(:, :, :), test_minus_trace_bar(:, :, :)
        real(dp), intent(in) :: trial_plus_trace_bar(:, :, :), trial_minus_trace_bar(:, :, :)
        real(dp), intent(in) :: test_plus_flux_bar(:, :, :), test_minus_flux_bar(:, :, :)
        real(dp), intent(in) :: trial_plus_flux_bar(:, :, :), trial_minus_flux_bar(:, :, :)
        real(dp), intent(in) :: surface_weights_bar(:), component_metric_bar(:, :, :)
        integer, intent(in) :: quadrature_count, test_plus_count, test_minus_count
        integer, intent(in) :: trial_plus_count, trial_minus_count, component_count

        valid_vector_sipg_cotangent = &
            size(test_plus_trace_bar, 1) == quadrature_count .and. &
            size(test_plus_trace_bar, 2) == test_plus_count .and. &
            size(test_plus_trace_bar, 3) == component_count .and. &
            size(test_minus_trace_bar, 1) == quadrature_count .and. &
            size(test_minus_trace_bar, 2) == test_minus_count .and. &
            size(test_minus_trace_bar, 3) == component_count .and. &
            size(trial_plus_trace_bar, 1) == quadrature_count .and. &
            size(trial_plus_trace_bar, 2) == trial_plus_count .and. &
            size(trial_plus_trace_bar, 3) == component_count .and. &
            size(trial_minus_trace_bar, 1) == quadrature_count .and. &
            size(trial_minus_trace_bar, 2) == trial_minus_count .and. &
            size(trial_minus_trace_bar, 3) == component_count .and. &
            size(test_plus_flux_bar, 1) == quadrature_count .and. &
            size(test_plus_flux_bar, 2) == test_plus_count .and. &
            size(test_plus_flux_bar, 3) == component_count .and. &
            size(test_minus_flux_bar, 1) == quadrature_count .and. &
            size(test_minus_flux_bar, 2) == test_minus_count .and. &
            size(test_minus_flux_bar, 3) == component_count .and. &
            size(trial_plus_flux_bar, 1) == quadrature_count .and. &
            size(trial_plus_flux_bar, 2) == trial_plus_count .and. &
            size(trial_plus_flux_bar, 3) == component_count .and. &
            size(trial_minus_flux_bar, 1) == quadrature_count .and. &
            size(trial_minus_flux_bar, 2) == trial_minus_count .and. &
            size(trial_minus_flux_bar, 3) == component_count .and. &
            size(surface_weights_bar) == quadrature_count .and. &
            size(component_metric_bar, 1) == quadrature_count .and. &
            size(component_metric_bar, 2) == component_count .and. &
            size(component_metric_bar, 3) == component_count
    end function valid_vector_sipg_cotangent

    subroutine validate_vector_sipg_inputs( &
            test_plus_trace, test_minus_trace, trial_plus_trace, trial_minus_trace, &
            test_plus_flux, test_minus_flux, trial_plus_flux, trial_minus_flux, &
            surface_weights, penalty, consistency_sign, component_metric, matrix, status)
        real(dp), intent(in) :: test_plus_trace(:, :, :), test_minus_trace(:, :, :)
        real(dp), intent(in) :: trial_plus_trace(:, :, :), trial_minus_trace(:, :, :)
        real(dp), intent(in) :: test_plus_flux(:, :, :), test_minus_flux(:, :, :)
        real(dp), intent(in) :: trial_plus_flux(:, :, :), trial_minus_flux(:, :, :)
        real(dp), intent(in) :: surface_weights(:), penalty
        integer, intent(in) :: consistency_sign
        real(dp), intent(in) :: component_metric(:, :, :), matrix(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_plus_count, test_minus_count
        integer :: trial_plus_count, trial_minus_count, component_count
        integer :: test_count, trial_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "vector SIPG interface received incompatible arrays")
        quadrature_count = size(test_plus_trace, 1)
        test_plus_count = size(test_plus_trace, 2)
        test_minus_count = size(test_minus_trace, 2)
        trial_plus_count = size(trial_plus_trace, 2)
        trial_minus_count = size(trial_minus_trace, 2)
        component_count = size(test_plus_trace, 3)
        test_count = test_plus_count + test_minus_count
        trial_count = trial_plus_count + trial_minus_count
        if (quadrature_count < 1 .or. test_plus_count < 1 .or. test_minus_count < 1 .or. &
            trial_plus_count < 1 .or. trial_minus_count < 1 .or. component_count < 1) return
        if (size(test_minus_trace, 1) /= quadrature_count .or. &
            size(test_minus_trace, 3) /= component_count .or. &
            size(trial_plus_trace, 1) /= quadrature_count .or. &
            size(trial_plus_trace, 3) /= component_count .or. &
            size(trial_minus_trace, 1) /= quadrature_count .or. &
            size(trial_minus_trace, 3) /= component_count .or. &
            size(test_plus_flux, 1) /= quadrature_count .or. &
            size(test_plus_flux, 2) /= test_plus_count .or. &
            size(test_plus_flux, 3) /= component_count .or. &
            size(test_minus_flux, 1) /= quadrature_count .or. &
            size(test_minus_flux, 2) /= test_minus_count .or. &
            size(test_minus_flux, 3) /= component_count .or. &
            size(trial_plus_flux, 1) /= quadrature_count .or. &
            size(trial_plus_flux, 2) /= trial_plus_count .or. &
            size(trial_plus_flux, 3) /= component_count .or. &
            size(trial_minus_flux, 1) /= quadrature_count .or. &
            size(trial_minus_flux, 2) /= trial_minus_count .or. &
            size(trial_minus_flux, 3) /= component_count .or. &
            size(surface_weights) /= quadrature_count .or. &
            size(component_metric, 1) /= quadrature_count .or. &
            size(component_metric, 2) /= component_count .or. &
            size(component_metric, 3) /= component_count .or. &
            size(matrix, 1) /= test_count .or. size(matrix, 2) /= trial_count) return
        if (consistency_sign < -1 .or. consistency_sign > 1 .or. &
            .not. ieee_is_finite(penalty) .or. penalty < 0.0_dp .or. &
            any(surface_weights <= 0.0_dp) .or. &
            any(.not. ieee_is_finite(test_plus_trace)) .or. &
            any(.not. ieee_is_finite(test_minus_trace)) .or. &
            any(.not. ieee_is_finite(trial_plus_trace)) .or. &
            any(.not. ieee_is_finite(trial_minus_trace)) .or. &
            any(.not. ieee_is_finite(test_plus_flux)) .or. &
            any(.not. ieee_is_finite(test_minus_flux)) .or. &
            any(.not. ieee_is_finite(trial_plus_flux)) .or. &
            any(.not. ieee_is_finite(trial_minus_flux)) .or. &
            any(.not. ieee_is_finite(surface_weights)) .or. &
            any(.not. ieee_is_finite(component_metric)) .or. &
            any(.not. ieee_is_finite(matrix))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_vector_sipg_inputs

end module fortfem_vector_jump_penalty
