module fortfem_nitsche_interface
    !! Symmetric scalar Nitsche interface consistency and penalty block.
    !!
    !! The supplied plus/minus flux traces use one common normal convention.
    !! For jump J=[T_plus,-T_minus] and average flux A=[(F_plus+F_minus)/2],
    !! the assembled block is
    !! -A^T W J - J^T W A + eta J^T W J.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_symmetric_nitsche_interface
    public :: assemble_symmetric_nitsche_interface_jvp
    public :: assemble_symmetric_nitsche_interface_vjp

contains

    subroutine assemble_symmetric_nitsche_interface( &
            plus_trace, minus_trace, plus_flux, minus_flux, surface_weights, &
            penalty, matrix, status)
        !! Assemble the symmetric scalar Nitsche interface block.
        real(dp), intent(in) :: plus_trace(:, :), minus_trace(:, :)
        real(dp), intent(in) :: plus_flux(:, :), minus_flux(:, :)
        real(dp), intent(in) :: surface_weights(:), penalty
        real(dp), intent(out) :: matrix(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, plus_count, minus_count, total_count
        integer :: quadrature, row, column
        real(dp), allocatable :: jump_trace(:), average_flux(:)
        real(dp) :: scale

        matrix = 0.0_dp
        call validate_nitsche_inputs( &
            plus_trace, minus_trace, plus_flux, minus_flux, surface_weights, &
            penalty, matrix, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(plus_trace, 1)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        total_count = plus_count + minus_count

        allocate(jump_trace(total_count), average_flux(total_count))
        do quadrature = 1, quadrature_count
            jump_trace(1:plus_count) = plus_trace(quadrature, :)
            jump_trace(plus_count + 1:total_count) = &
                -minus_trace(quadrature, :)
            average_flux(1:plus_count) = 0.5_dp*plus_flux(quadrature, :)
            average_flux(plus_count + 1:total_count) = &
                0.5_dp*minus_flux(quadrature, :)
            scale = surface_weights(quadrature)
            do row = 1, total_count
                do column = 1, total_count
                    matrix(row, column) = matrix(row, column) + scale*( &
                        -average_flux(row)*jump_trace(column) - &
                        jump_trace(row)*average_flux(column) + &
                        penalty*jump_trace(row)*jump_trace(column))
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_symmetric_nitsche_interface

    subroutine assemble_symmetric_nitsche_interface_jvp( &
            plus_trace, minus_trace, plus_flux, minus_flux, surface_weights, &
            penalty, plus_trace_dot, minus_trace_dot, plus_flux_dot, &
            minus_flux_dot, surface_weights_dot, penalty_dot, matrix_dot, status)
        !! Apply the product-rule JVP of the symmetric Nitsche block.
        real(dp), intent(in) :: plus_trace(:, :), minus_trace(:, :)
        real(dp), intent(in) :: plus_flux(:, :), minus_flux(:, :)
        real(dp), intent(in) :: surface_weights(:), penalty
        real(dp), intent(in) :: plus_trace_dot(:, :), minus_trace_dot(:, :)
        real(dp), intent(in) :: plus_flux_dot(:, :), minus_flux_dot(:, :)
        real(dp), intent(in) :: surface_weights_dot(:), penalty_dot
        real(dp), intent(out) :: matrix_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, plus_count, minus_count, total_count
        integer :: quadrature, row, column
        real(dp), allocatable :: jump_trace(:), jump_trace_dot(:)
        real(dp), allocatable :: average_flux(:), average_flux_dot(:)
        real(dp) :: scale, scale_dot

        matrix_dot = 0.0_dp
        call validate_nitsche_inputs( &
            plus_trace, minus_trace, plus_flux, minus_flux, surface_weights, &
            penalty, matrix_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(plus_trace, 1)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        total_count = plus_count + minus_count
        if (size(plus_trace_dot, 1) /= quadrature_count .or. &
            size(plus_trace_dot, 2) /= plus_count .or. &
            size(minus_trace_dot, 1) /= quadrature_count .or. &
            size(minus_trace_dot, 2) /= minus_count .or. &
            size(plus_flux_dot, 1) /= quadrature_count .or. &
            size(plus_flux_dot, 2) /= plus_count .or. &
            size(minus_flux_dot, 1) /= quadrature_count .or. &
            size(minus_flux_dot, 2) /= minus_count .or. &
            size(surface_weights_dot) /= quadrature_count .or. &
            .not. ieee_is_finite(penalty_dot) .or. &
            any(.not. ieee_is_finite(plus_trace_dot)) .or. &
            any(.not. ieee_is_finite(minus_trace_dot)) .or. &
            any(.not. ieee_is_finite(plus_flux_dot)) .or. &
            any(.not. ieee_is_finite(minus_flux_dot)) .or. &
            any(.not. ieee_is_finite(surface_weights_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "symmetric Nitsche interface JVP has incompatible increments")
            return
        end if

        allocate(jump_trace(total_count), jump_trace_dot(total_count), &
            average_flux(total_count), average_flux_dot(total_count))
        do quadrature = 1, quadrature_count
            jump_trace(1:plus_count) = plus_trace(quadrature, :)
            jump_trace(plus_count + 1:total_count) = &
                -minus_trace(quadrature, :)
            jump_trace_dot(1:plus_count) = plus_trace_dot(quadrature, :)
            jump_trace_dot(plus_count + 1:total_count) = &
                -minus_trace_dot(quadrature, :)
            average_flux(1:plus_count) = 0.5_dp*plus_flux(quadrature, :)
            average_flux(plus_count + 1:total_count) = &
                0.5_dp*minus_flux(quadrature, :)
            average_flux_dot(1:plus_count) = &
                0.5_dp*plus_flux_dot(quadrature, :)
            average_flux_dot(plus_count + 1:total_count) = &
                0.5_dp*minus_flux_dot(quadrature, :)
            scale = surface_weights(quadrature)
            scale_dot = surface_weights_dot(quadrature)
            do row = 1, total_count
                do column = 1, total_count
                    matrix_dot(row, column) = matrix_dot(row, column) + &
                        scale_dot*( &
                        -average_flux(row)*jump_trace(column) - &
                        jump_trace(row)*average_flux(column) + &
                        penalty*jump_trace(row)*jump_trace(column)) + &
                        scale*( &
                        -average_flux_dot(row)*jump_trace(column) - &
                        average_flux(row)*jump_trace_dot(column) - &
                        jump_trace_dot(row)*average_flux(column) - &
                        jump_trace(row)*average_flux_dot(column) + &
                        penalty_dot*jump_trace(row)*jump_trace(column) + &
                        penalty*jump_trace_dot(row)*jump_trace(column) + &
                        penalty*jump_trace(row)*jump_trace_dot(column))
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_symmetric_nitsche_interface_jvp

    subroutine assemble_symmetric_nitsche_interface_vjp( &
            plus_trace, minus_trace, plus_flux, minus_flux, surface_weights, &
            penalty, matrix_bar, plus_trace_bar, minus_trace_bar, &
            plus_flux_bar, minus_flux_bar, surface_weights_bar, penalty_bar, &
            status)
        !! Apply the reverse product of the symmetric Nitsche block.
        real(dp), intent(in) :: plus_trace(:, :), minus_trace(:, :)
        real(dp), intent(in) :: plus_flux(:, :), minus_flux(:, :)
        real(dp), intent(in) :: surface_weights(:), penalty, matrix_bar(:, :)
        real(dp), intent(out) :: plus_trace_bar(:, :), minus_trace_bar(:, :)
        real(dp), intent(out) :: plus_flux_bar(:, :), minus_flux_bar(:, :)
        real(dp), intent(out) :: surface_weights_bar(:), penalty_bar
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, plus_count, minus_count, total_count
        integer :: quadrature, row, column
        real(dp), allocatable :: jump_trace(:), jump_bar(:)
        real(dp), allocatable :: average_flux(:), average_flux_bar(:)
        real(dp) :: matrix_bar_quadratic, penalty_quadratic

        plus_trace_bar = 0.0_dp
        minus_trace_bar = 0.0_dp
        plus_flux_bar = 0.0_dp
        minus_flux_bar = 0.0_dp
        surface_weights_bar = 0.0_dp
        penalty_bar = 0.0_dp
        call validate_nitsche_inputs( &
            plus_trace, minus_trace, plus_flux, minus_flux, surface_weights, &
            penalty, matrix_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(plus_trace, 1)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        total_count = plus_count + minus_count
        if (size(plus_trace_bar, 1) /= quadrature_count .or. &
            size(plus_trace_bar, 2) /= plus_count .or. &
            size(minus_trace_bar, 1) /= quadrature_count .or. &
            size(minus_trace_bar, 2) /= minus_count .or. &
            size(plus_flux_bar, 1) /= quadrature_count .or. &
            size(plus_flux_bar, 2) /= plus_count .or. &
            size(minus_flux_bar, 1) /= quadrature_count .or. &
            size(minus_flux_bar, 2) /= minus_count .or. &
            size(surface_weights_bar) /= quadrature_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "symmetric Nitsche interface VJP has incompatible cotangents")
            return
        end if

        allocate(jump_trace(total_count), jump_bar(total_count), &
            average_flux(total_count), average_flux_bar(total_count))
        do quadrature = 1, quadrature_count
            jump_trace(1:plus_count) = plus_trace(quadrature, :)
            jump_trace(plus_count + 1:total_count) = -minus_trace(quadrature, :)
            average_flux(1:plus_count) = 0.5_dp*plus_flux(quadrature, :)
            average_flux(plus_count + 1:total_count) = &
                0.5_dp*minus_flux(quadrature, :)
            jump_bar = 0.0_dp
            average_flux_bar = 0.0_dp
            do row = 1, total_count
                do column = 1, total_count
                    jump_bar(row) = jump_bar(row) + &
                        surface_weights(quadrature)*( &
                        -matrix_bar(row, column)*average_flux(column) - &
                        matrix_bar(column, row)*average_flux(column) + &
                        penalty*(matrix_bar(row, column) + &
                        matrix_bar(column, row))*jump_trace(column))
                    average_flux_bar(row) = average_flux_bar(row) + &
                        surface_weights(quadrature)*( &
                        -matrix_bar(row, column)*jump_trace(column) - &
                        matrix_bar(column, row)*jump_trace(column))
                end do
            end do
            plus_trace_bar(quadrature, :) = jump_bar(1:plus_count)
            minus_trace_bar(quadrature, :) = - &
                jump_bar(plus_count + 1:total_count)
            plus_flux_bar(quadrature, :) = &
                0.5_dp*average_flux_bar(1:plus_count)
            minus_flux_bar(quadrature, :) = &
                0.5_dp*average_flux_bar(plus_count + 1:total_count)
            matrix_bar_quadratic = 0.0_dp
            penalty_quadratic = 0.0_dp
            do row = 1, total_count
                do column = 1, total_count
                    matrix_bar_quadratic = matrix_bar_quadratic + &
                        matrix_bar(row, column)*( &
                        -average_flux(row)*jump_trace(column) - &
                        jump_trace(row)*average_flux(column))
                    penalty_quadratic = penalty_quadratic + &
                        matrix_bar(row, column)*jump_trace(row)*jump_trace(column)
                end do
            end do
            surface_weights_bar(quadrature) = matrix_bar_quadratic + &
                penalty*penalty_quadratic
            penalty_bar = penalty_bar + &
                surface_weights(quadrature)*penalty_quadratic
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_symmetric_nitsche_interface_vjp

    subroutine validate_nitsche_inputs( &
            plus_trace, minus_trace, plus_flux, minus_flux, surface_weights, &
            penalty, matrix, status)
        real(dp), intent(in) :: plus_trace(:, :), minus_trace(:, :)
        real(dp), intent(in) :: plus_flux(:, :), minus_flux(:, :)
        real(dp), intent(in) :: surface_weights(:), penalty, matrix(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, plus_count, minus_count, total_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "symmetric Nitsche interface received incompatible arrays")
        quadrature_count = size(plus_trace, 1)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        total_count = plus_count + minus_count
        if (quadrature_count < 1 .or. plus_count < 1 .or. minus_count < 1) return
        if (size(minus_trace, 1) /= quadrature_count .or. &
            size(plus_flux, 1) /= quadrature_count .or. &
            size(plus_flux, 2) /= plus_count .or. &
            size(minus_flux, 1) /= quadrature_count .or. &
            size(minus_flux, 2) /= minus_count .or. &
            size(surface_weights) /= quadrature_count .or. &
            size(matrix, 1) /= total_count .or. &
            size(matrix, 2) /= total_count) return
        if (.not. ieee_is_finite(penalty) .or. penalty < 0.0_dp .or. &
            any(.not. ieee_is_finite(plus_trace)) .or. &
            any(.not. ieee_is_finite(minus_trace)) .or. &
            any(.not. ieee_is_finite(plus_flux)) .or. &
            any(.not. ieee_is_finite(minus_flux)) .or. &
            any(.not. ieee_is_finite(surface_weights)) .or. &
            any(.not. ieee_is_finite(matrix)) .or. &
            any(surface_weights <= 0.0_dp)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_nitsche_inputs

end module fortfem_nitsche_interface
