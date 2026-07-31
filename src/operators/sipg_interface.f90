module fortfem_scalar_sipg_interface
    !! Scalar SIPG/NIPG/IIPG interface blocks with fixed normal orientation.
    !!
    !! For test traces V and trial traces U, the returned test-by-trial block
    !! is the weighted form
    !! -{q}[u] - theta[v]{p} + eta[v][u], where theta=1,0,-1 selects
    !! symmetric, incomplete, and nonsymmetric interior penalty coupling.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_scalar_sipg_interface
    public :: assemble_scalar_sipg_interface_jvp
    public :: assemble_scalar_sipg_interface_vjp

contains

    subroutine assemble_scalar_sipg_interface( &
            test_plus_trace, test_minus_trace, trial_plus_trace, &
            trial_minus_trace, test_plus_flux, test_minus_flux, &
            trial_plus_flux, trial_minus_flux, surface_weights, penalty, &
            consistency_sign, matrix, status)
        !! Assemble one scalar SIPG interface block.
        real(dp), intent(in) :: test_plus_trace(:, :), test_minus_trace(:, :)
        real(dp), intent(in) :: trial_plus_trace(:, :), trial_minus_trace(:, :)
        real(dp), intent(in) :: test_plus_flux(:, :), test_minus_flux(:, :)
        real(dp), intent(in) :: trial_plus_flux(:, :), trial_minus_flux(:, :)
        real(dp), intent(in) :: surface_weights(:), penalty
        integer, intent(in) :: consistency_sign
        real(dp), intent(out) :: matrix(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_plus_count, test_minus_count
        integer :: trial_plus_count, trial_minus_count, test_count, trial_count
        integer :: quadrature, row, column
        real(dp), allocatable :: test_jump(:), trial_jump(:)
        real(dp), allocatable :: test_average(:), trial_average(:)
        real(dp) :: scale

        matrix = 0.0_dp
        call validate_sipg_inputs( &
            test_plus_trace, test_minus_trace, trial_plus_trace, &
            trial_minus_trace, test_plus_flux, test_minus_flux, &
            trial_plus_flux, trial_minus_flux, surface_weights, penalty, &
            consistency_sign, matrix, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(test_plus_trace, 1)
        test_plus_count = size(test_plus_trace, 2)
        test_minus_count = size(test_minus_trace, 2)
        trial_plus_count = size(trial_plus_trace, 2)
        trial_minus_count = size(trial_minus_trace, 2)
        test_count = test_plus_count + test_minus_count
        trial_count = trial_plus_count + trial_minus_count

        allocate(test_jump(test_count), trial_jump(trial_count), &
            test_average(test_count), trial_average(trial_count))
        do quadrature = 1, quadrature_count
            call fill_jump_average( &
                test_plus_trace(quadrature, :), test_minus_trace(quadrature, :), &
                test_plus_count, test_jump, test_average, &
                test_plus_flux(quadrature, :), test_minus_flux(quadrature, :))
            call fill_jump_average( &
                trial_plus_trace(quadrature, :), trial_minus_trace(quadrature, :), &
                trial_plus_count, trial_jump, trial_average, &
                trial_plus_flux(quadrature, :), trial_minus_flux(quadrature, :))
            scale = surface_weights(quadrature)
            call add_sipg_kernel( &
                test_jump, trial_jump, test_average, trial_average, scale, &
                penalty, consistency_sign, matrix)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_scalar_sipg_interface

    subroutine assemble_scalar_sipg_interface_jvp( &
            test_plus_trace, test_minus_trace, trial_plus_trace, &
            trial_minus_trace, test_plus_flux, test_minus_flux, &
            trial_plus_flux, trial_minus_flux, surface_weights, penalty, &
            consistency_sign, test_plus_trace_dot, test_minus_trace_dot, &
            trial_plus_trace_dot, trial_minus_trace_dot, test_plus_flux_dot, &
            test_minus_flux_dot, trial_plus_flux_dot, trial_minus_flux_dot, &
            surface_weights_dot, penalty_dot, matrix_dot, status)
        !! Apply the product-rule JVP of a scalar SIPG interface block.
        real(dp), intent(in) :: test_plus_trace(:, :), test_minus_trace(:, :)
        real(dp), intent(in) :: trial_plus_trace(:, :), trial_minus_trace(:, :)
        real(dp), intent(in) :: test_plus_flux(:, :), test_minus_flux(:, :)
        real(dp), intent(in) :: trial_plus_flux(:, :), trial_minus_flux(:, :)
        real(dp), intent(in) :: surface_weights(:), penalty
        integer, intent(in) :: consistency_sign
        real(dp), intent(in) :: test_plus_trace_dot(:, :), test_minus_trace_dot(:, :)
        real(dp), intent(in) :: trial_plus_trace_dot(:, :), trial_minus_trace_dot(:, :)
        real(dp), intent(in) :: test_plus_flux_dot(:, :), test_minus_flux_dot(:, :)
        real(dp), intent(in) :: trial_plus_flux_dot(:, :), trial_minus_flux_dot(:, :)
        real(dp), intent(in) :: surface_weights_dot(:), penalty_dot
        real(dp), intent(out) :: matrix_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_plus_count, test_minus_count
        integer :: trial_plus_count, trial_minus_count, test_count, trial_count
        integer :: quadrature
        real(dp), allocatable :: test_jump(:), trial_jump(:), test_jump_dot(:)
        real(dp), allocatable :: trial_jump_dot(:), test_average(:), trial_average(:)
        real(dp), allocatable :: test_average_dot(:), trial_average_dot(:)
        real(dp) :: scale, scale_dot

        matrix_dot = 0.0_dp
        call validate_sipg_inputs( &
            test_plus_trace, test_minus_trace, trial_plus_trace, &
            trial_minus_trace, test_plus_flux, test_minus_flux, &
            trial_plus_flux, trial_minus_flux, surface_weights, penalty, &
            consistency_sign, matrix_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(test_plus_trace, 1)
        test_plus_count = size(test_plus_trace, 2)
        test_minus_count = size(test_minus_trace, 2)
        trial_plus_count = size(trial_plus_trace, 2)
        trial_minus_count = size(trial_minus_trace, 2)
        test_count = test_plus_count + test_minus_count
        trial_count = trial_plus_count + trial_minus_count
        if (.not. valid_sipg_direction( &
            test_plus_trace_dot, test_minus_trace_dot, trial_plus_trace_dot, &
            trial_minus_trace_dot, test_plus_flux_dot, test_minus_flux_dot, &
            trial_plus_flux_dot, trial_minus_flux_dot, surface_weights_dot, &
            penalty_dot, quadrature_count, test_plus_count, test_minus_count, &
            trial_plus_count, trial_minus_count)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "scalar SIPG interface JVP has incompatible increments")
            return
        end if

        allocate(test_jump(test_count), trial_jump(trial_count), &
            test_jump_dot(test_count), trial_jump_dot(trial_count), &
            test_average(test_count), trial_average(trial_count), &
            test_average_dot(test_count), trial_average_dot(trial_count))
        do quadrature = 1, quadrature_count
            call fill_jump_average( &
                test_plus_trace(quadrature, :), test_minus_trace(quadrature, :), &
                test_plus_count, test_jump, test_average, &
                test_plus_flux(quadrature, :), test_minus_flux(quadrature, :))
            call fill_jump_average( &
                trial_plus_trace(quadrature, :), trial_minus_trace(quadrature, :), &
                trial_plus_count, trial_jump, trial_average, &
                trial_plus_flux(quadrature, :), trial_minus_flux(quadrature, :))
            call fill_jump_average( &
                test_plus_trace_dot(quadrature, :), test_minus_trace_dot(quadrature, :), &
                test_plus_count, test_jump_dot, test_average_dot, &
                test_plus_flux_dot(quadrature, :), test_minus_flux_dot(quadrature, :))
            call fill_jump_average( &
                trial_plus_trace_dot(quadrature, :), trial_minus_trace_dot(quadrature, :), &
                trial_plus_count, trial_jump_dot, trial_average_dot, &
                trial_plus_flux_dot(quadrature, :), trial_minus_flux_dot(quadrature, :))
            scale = surface_weights(quadrature)
            scale_dot = surface_weights_dot(quadrature)
            call add_sipg_kernel_jvp( &
                test_jump, trial_jump, test_average, trial_average, &
                test_jump_dot, trial_jump_dot, test_average_dot, trial_average_dot, &
                scale, scale_dot, penalty, penalty_dot, consistency_sign, matrix_dot)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_scalar_sipg_interface_jvp

    subroutine assemble_scalar_sipg_interface_vjp( &
            test_plus_trace, test_minus_trace, trial_plus_trace, &
            trial_minus_trace, test_plus_flux, test_minus_flux, &
            trial_plus_flux, trial_minus_flux, surface_weights, penalty, &
            consistency_sign, matrix_bar, test_plus_trace_bar, &
            test_minus_trace_bar, trial_plus_trace_bar, trial_minus_trace_bar, &
            test_plus_flux_bar, test_minus_flux_bar, trial_plus_flux_bar, &
            trial_minus_flux_bar, surface_weights_bar, penalty_bar, status)
        !! Apply the reverse product of a scalar SIPG interface block.
        real(dp), intent(in) :: test_plus_trace(:, :), test_minus_trace(:, :)
        real(dp), intent(in) :: trial_plus_trace(:, :), trial_minus_trace(:, :)
        real(dp), intent(in) :: test_plus_flux(:, :), test_minus_flux(:, :)
        real(dp), intent(in) :: trial_plus_flux(:, :), trial_minus_flux(:, :)
        real(dp), intent(in) :: surface_weights(:), penalty
        integer, intent(in) :: consistency_sign
        real(dp), intent(in) :: matrix_bar(:, :)
        real(dp), intent(out) :: test_plus_trace_bar(:, :), test_minus_trace_bar(:, :)
        real(dp), intent(out) :: trial_plus_trace_bar(:, :), trial_minus_trace_bar(:, :)
        real(dp), intent(out) :: test_plus_flux_bar(:, :), test_minus_flux_bar(:, :)
        real(dp), intent(out) :: trial_plus_flux_bar(:, :), trial_minus_flux_bar(:, :)
        real(dp), intent(out) :: surface_weights_bar(:), penalty_bar
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_plus_count, test_minus_count
        integer :: trial_plus_count, trial_minus_count, test_count, trial_count
        integer :: quadrature, row, column
        real(dp), allocatable :: test_jump(:), trial_jump(:), test_jump_bar(:)
        real(dp), allocatable :: trial_jump_bar(:), test_average(:), trial_average(:)
        real(dp), allocatable :: test_average_bar(:), trial_average_bar(:)
        real(dp) :: kernel_contraction, penalty_contraction

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
        call validate_sipg_inputs( &
            test_plus_trace, test_minus_trace, trial_plus_trace, &
            trial_minus_trace, test_plus_flux, test_minus_flux, &
            trial_plus_flux, trial_minus_flux, surface_weights, penalty, &
            consistency_sign, matrix_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(test_plus_trace, 1)
        test_plus_count = size(test_plus_trace, 2)
        test_minus_count = size(test_minus_trace, 2)
        trial_plus_count = size(trial_plus_trace, 2)
        trial_minus_count = size(trial_minus_trace, 2)
        test_count = test_plus_count + test_minus_count
        trial_count = trial_plus_count + trial_minus_count
        if (.not. valid_sipg_cotangent( &
            test_plus_trace_bar, test_minus_trace_bar, trial_plus_trace_bar, &
            trial_minus_trace_bar, test_plus_flux_bar, test_minus_flux_bar, &
            trial_plus_flux_bar, trial_minus_flux_bar, surface_weights_bar, &
            quadrature_count, test_plus_count, test_minus_count, &
            trial_plus_count, trial_minus_count)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "scalar SIPG interface VJP has incompatible cotangents")
            return
        end if

        allocate(test_jump(test_count), trial_jump(trial_count), &
            test_jump_bar(test_count), trial_jump_bar(trial_count), &
            test_average(test_count), trial_average(trial_count), &
            test_average_bar(test_count), trial_average_bar(trial_count))
        do quadrature = 1, quadrature_count
            call fill_jump_average( &
                test_plus_trace(quadrature, :), test_minus_trace(quadrature, :), &
                test_plus_count, test_jump, test_average, &
                test_plus_flux(quadrature, :), test_minus_flux(quadrature, :))
            call fill_jump_average( &
                trial_plus_trace(quadrature, :), trial_minus_trace(quadrature, :), &
                trial_plus_count, trial_jump, trial_average, &
                trial_plus_flux(quadrature, :), trial_minus_flux(quadrature, :))
            test_jump_bar = 0.0_dp
            trial_jump_bar = 0.0_dp
            test_average_bar = 0.0_dp
            trial_average_bar = 0.0_dp
            do row = 1, test_count
                do column = 1, trial_count
                    test_jump_bar(row) = test_jump_bar(row) + &
                        surface_weights(quadrature)*( &
                        -consistency_sign*matrix_bar(row, column)* &
                        trial_average(column) + penalty*matrix_bar(row, column)* &
                        trial_jump(column))
                    test_average_bar(row) = test_average_bar(row) - &
                        surface_weights(quadrature)*matrix_bar(row, column)* &
                        trial_jump(column)
                    trial_jump_bar(column) = trial_jump_bar(column) + &
                        surface_weights(quadrature)*( &
                        -matrix_bar(row, column)*test_average(row) + &
                        penalty*matrix_bar(row, column)*test_jump(row))
                    trial_average_bar(column) = trial_average_bar(column) - &
                        surface_weights(quadrature)*consistency_sign* &
                        matrix_bar(row, column)*test_jump(row)
                end do
            end do
            test_plus_trace_bar(quadrature, :) = test_jump_bar(1:test_plus_count)
            test_minus_trace_bar(quadrature, :) = - &
                test_jump_bar(test_plus_count + 1:test_count)
            trial_plus_trace_bar(quadrature, :) = trial_jump_bar(1:trial_plus_count)
            trial_minus_trace_bar(quadrature, :) = - &
                trial_jump_bar(trial_plus_count + 1:trial_count)
            test_plus_flux_bar(quadrature, :) = &
                0.5_dp*test_average_bar(1:test_plus_count)
            test_minus_flux_bar(quadrature, :) = &
                0.5_dp*test_average_bar(test_plus_count + 1:test_count)
            trial_plus_flux_bar(quadrature, :) = &
                0.5_dp*trial_average_bar(1:trial_plus_count)
            trial_minus_flux_bar(quadrature, :) = &
                0.5_dp*trial_average_bar(trial_plus_count + 1:trial_count)
            kernel_contraction = 0.0_dp
            penalty_contraction = 0.0_dp
            do row = 1, test_count
                do column = 1, trial_count
                    kernel_contraction = kernel_contraction + matrix_bar(row, column)*( &
                        -test_average(row)*trial_jump(column) - &
                        consistency_sign*test_jump(row)*trial_average(column))
                    penalty_contraction = penalty_contraction + &
                        matrix_bar(row, column)*test_jump(row)*trial_jump(column)
                end do
            end do
            surface_weights_bar(quadrature) = kernel_contraction + &
                penalty*penalty_contraction
            penalty_bar = penalty_bar + &
                surface_weights(quadrature)*penalty_contraction
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_scalar_sipg_interface_vjp

    subroutine fill_jump_average( &
            plus_trace, minus_trace, plus_count, jump, average, plus_flux, minus_flux)
        real(dp), intent(in) :: plus_trace(:), minus_trace(:)
        integer, intent(in) :: plus_count
        real(dp), intent(out) :: jump(:), average(:)
        real(dp), intent(in) :: plus_flux(:), minus_flux(:)
        integer :: minus_count

        minus_count = size(minus_trace)
        jump(1:plus_count) = plus_trace
        jump(plus_count + 1:plus_count + minus_count) = -minus_trace
        average(1:plus_count) = 0.5_dp*plus_flux
        average(plus_count + 1:plus_count + minus_count) = 0.5_dp*minus_flux
    end subroutine fill_jump_average

    subroutine add_sipg_kernel( &
            test_jump, trial_jump, test_average, trial_average, scale, penalty, &
            consistency_sign, matrix)
        real(dp), intent(in) :: test_jump(:), trial_jump(:), test_average(:)
        real(dp), intent(in) :: trial_average(:), scale, penalty
        integer, intent(in) :: consistency_sign
        real(dp), intent(inout) :: matrix(:, :)
        integer :: row, column

        do row = 1, size(test_jump)
            do column = 1, size(trial_jump)
                matrix(row, column) = matrix(row, column) + scale*( &
                    -test_average(row)*trial_jump(column) - &
                    consistency_sign*test_jump(row)*trial_average(column) + &
                    penalty*test_jump(row)*trial_jump(column))
            end do
        end do
    end subroutine add_sipg_kernel

    subroutine add_sipg_kernel_jvp( &
            test_jump, trial_jump, test_average, trial_average, test_jump_dot, &
            trial_jump_dot, test_average_dot, trial_average_dot, scale, scale_dot, &
            penalty, penalty_dot, consistency_sign, matrix_dot)
        real(dp), intent(in) :: test_jump(:), trial_jump(:), test_average(:)
        real(dp), intent(in) :: trial_average(:), test_jump_dot(:), trial_jump_dot(:)
        real(dp), intent(in) :: test_average_dot(:), trial_average_dot(:)
        real(dp), intent(in) :: scale, scale_dot, penalty, penalty_dot
        integer, intent(in) :: consistency_sign
        real(dp), intent(inout) :: matrix_dot(:, :)
        integer :: row, column
        real(dp) :: kernel, kernel_dot

        do row = 1, size(test_jump)
            do column = 1, size(trial_jump)
                kernel = -test_average(row)*trial_jump(column) - &
                    consistency_sign*test_jump(row)*trial_average(column) + &
                    penalty*test_jump(row)*trial_jump(column)
                kernel_dot = -test_average_dot(row)*trial_jump(column) - &
                    test_average(row)*trial_jump_dot(column) - &
                    consistency_sign*(test_jump_dot(row)*trial_average(column) + &
                    test_jump(row)*trial_average_dot(column)) + &
                    penalty_dot*test_jump(row)*trial_jump(column) + &
                    penalty*test_jump_dot(row)*trial_jump(column) + &
                    penalty*test_jump(row)*trial_jump_dot(column)
                matrix_dot(row, column) = matrix_dot(row, column) + &
                    scale_dot*kernel + scale*kernel_dot
            end do
        end do
    end subroutine add_sipg_kernel_jvp

    logical function valid_sipg_direction( &
            test_plus_trace_dot, test_minus_trace_dot, trial_plus_trace_dot, &
            trial_minus_trace_dot, test_plus_flux_dot, test_minus_flux_dot, &
            trial_plus_flux_dot, trial_minus_flux_dot, surface_weights_dot, &
            penalty_dot, quadrature_count, test_plus_count, test_minus_count, &
            trial_plus_count, trial_minus_count)
        real(dp), intent(in) :: test_plus_trace_dot(:, :), test_minus_trace_dot(:, :)
        real(dp), intent(in) :: trial_plus_trace_dot(:, :), trial_minus_trace_dot(:, :)
        real(dp), intent(in) :: test_plus_flux_dot(:, :), test_minus_flux_dot(:, :)
        real(dp), intent(in) :: trial_plus_flux_dot(:, :), trial_minus_flux_dot(:, :)
        real(dp), intent(in) :: surface_weights_dot(:), penalty_dot
        integer, intent(in) :: quadrature_count, test_plus_count, test_minus_count
        integer, intent(in) :: trial_plus_count, trial_minus_count

        valid_sipg_direction = size(test_plus_trace_dot, 1) == quadrature_count .and. &
            size(test_plus_trace_dot, 2) == test_plus_count .and. &
            size(test_minus_trace_dot, 1) == quadrature_count .and. &
            size(test_minus_trace_dot, 2) == test_minus_count .and. &
            size(trial_plus_trace_dot, 1) == quadrature_count .and. &
            size(trial_plus_trace_dot, 2) == trial_plus_count .and. &
            size(trial_minus_trace_dot, 1) == quadrature_count .and. &
            size(trial_minus_trace_dot, 2) == trial_minus_count .and. &
            size(test_plus_flux_dot, 1) == quadrature_count .and. &
            size(test_plus_flux_dot, 2) == test_plus_count .and. &
            size(test_minus_flux_dot, 1) == quadrature_count .and. &
            size(test_minus_flux_dot, 2) == test_minus_count .and. &
            size(trial_plus_flux_dot, 1) == quadrature_count .and. &
            size(trial_plus_flux_dot, 2) == trial_plus_count .and. &
            size(trial_minus_flux_dot, 1) == quadrature_count .and. &
            size(trial_minus_flux_dot, 2) == trial_minus_count .and. &
            size(surface_weights_dot) == quadrature_count .and. &
            ieee_is_finite(penalty_dot) .and. &
            all(ieee_is_finite(test_plus_trace_dot)) .and. &
            all(ieee_is_finite(test_minus_trace_dot)) .and. &
            all(ieee_is_finite(trial_plus_trace_dot)) .and. &
            all(ieee_is_finite(trial_minus_trace_dot)) .and. &
            all(ieee_is_finite(test_plus_flux_dot)) .and. &
            all(ieee_is_finite(test_minus_flux_dot)) .and. &
            all(ieee_is_finite(trial_plus_flux_dot)) .and. &
            all(ieee_is_finite(trial_minus_flux_dot)) .and. &
            all(ieee_is_finite(surface_weights_dot))
    end function valid_sipg_direction

    logical function valid_sipg_cotangent( &
            test_plus_trace_bar, test_minus_trace_bar, trial_plus_trace_bar, &
            trial_minus_trace_bar, test_plus_flux_bar, test_minus_flux_bar, &
            trial_plus_flux_bar, trial_minus_flux_bar, surface_weights_bar, &
            quadrature_count, test_plus_count, test_minus_count, trial_plus_count, &
            trial_minus_count)
        real(dp), intent(in) :: test_plus_trace_bar(:, :), test_minus_trace_bar(:, :)
        real(dp), intent(in) :: trial_plus_trace_bar(:, :), trial_minus_trace_bar(:, :)
        real(dp), intent(in) :: test_plus_flux_bar(:, :), test_minus_flux_bar(:, :)
        real(dp), intent(in) :: trial_plus_flux_bar(:, :), trial_minus_flux_bar(:, :)
        real(dp), intent(in) :: surface_weights_bar(:)
        integer, intent(in) :: quadrature_count, test_plus_count, test_minus_count
        integer, intent(in) :: trial_plus_count, trial_minus_count

        valid_sipg_cotangent = size(test_plus_trace_bar, 1) == quadrature_count .and. &
            size(test_plus_trace_bar, 2) == test_plus_count .and. &
            size(test_minus_trace_bar, 1) == quadrature_count .and. &
            size(test_minus_trace_bar, 2) == test_minus_count .and. &
            size(trial_plus_trace_bar, 1) == quadrature_count .and. &
            size(trial_plus_trace_bar, 2) == trial_plus_count .and. &
            size(trial_minus_trace_bar, 1) == quadrature_count .and. &
            size(trial_minus_trace_bar, 2) == trial_minus_count .and. &
            size(test_plus_flux_bar, 1) == quadrature_count .and. &
            size(test_plus_flux_bar, 2) == test_plus_count .and. &
            size(test_minus_flux_bar, 1) == quadrature_count .and. &
            size(test_minus_flux_bar, 2) == test_minus_count .and. &
            size(trial_plus_flux_bar, 1) == quadrature_count .and. &
            size(trial_plus_flux_bar, 2) == trial_plus_count .and. &
            size(trial_minus_flux_bar, 1) == quadrature_count .and. &
            size(trial_minus_flux_bar, 2) == trial_minus_count .and. &
            size(surface_weights_bar) == quadrature_count
    end function valid_sipg_cotangent

    subroutine validate_sipg_inputs( &
            test_plus_trace, test_minus_trace, trial_plus_trace, &
            trial_minus_trace, test_plus_flux, test_minus_flux, &
            trial_plus_flux, trial_minus_flux, surface_weights, penalty, &
            consistency_sign, matrix, status)
        real(dp), intent(in) :: test_plus_trace(:, :), test_minus_trace(:, :)
        real(dp), intent(in) :: trial_plus_trace(:, :), trial_minus_trace(:, :)
        real(dp), intent(in) :: test_plus_flux(:, :), test_minus_flux(:, :)
        real(dp), intent(in) :: trial_plus_flux(:, :), trial_minus_flux(:, :)
        real(dp), intent(in) :: surface_weights(:), penalty, matrix(:, :)
        integer, intent(in) :: consistency_sign
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, test_plus_count, test_minus_count
        integer :: trial_plus_count, trial_minus_count, test_count, trial_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "scalar SIPG interface received incompatible arrays")
        quadrature_count = size(test_plus_trace, 1)
        test_plus_count = size(test_plus_trace, 2)
        test_minus_count = size(test_minus_trace, 2)
        trial_plus_count = size(trial_plus_trace, 2)
        trial_minus_count = size(trial_minus_trace, 2)
        test_count = test_plus_count + test_minus_count
        trial_count = trial_plus_count + trial_minus_count
        if (quadrature_count < 1 .or. test_plus_count < 1 .or. &
            test_minus_count < 1 .or. trial_plus_count < 1 .or. &
            trial_minus_count < 1) return
        if (size(test_minus_trace, 1) /= quadrature_count .or. &
            size(trial_plus_trace, 1) /= quadrature_count .or. &
            size(trial_minus_trace, 1) /= quadrature_count .or. &
            size(test_plus_flux, 1) /= quadrature_count .or. &
            size(test_plus_flux, 2) /= test_plus_count .or. &
            size(test_minus_flux, 1) /= quadrature_count .or. &
            size(test_minus_flux, 2) /= test_minus_count .or. &
            size(trial_plus_flux, 1) /= quadrature_count .or. &
            size(trial_plus_flux, 2) /= trial_plus_count .or. &
            size(trial_minus_flux, 1) /= quadrature_count .or. &
            size(trial_minus_flux, 2) /= trial_minus_count .or. &
            size(surface_weights) /= quadrature_count .or. &
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
            any(.not. ieee_is_finite(matrix))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_sipg_inputs

end module fortfem_scalar_sipg_interface
