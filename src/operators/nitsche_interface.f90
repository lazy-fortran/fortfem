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

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "symmetric Nitsche interface received incompatible arrays")
        matrix = 0.0_dp
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
        if (.not. ieee_is_finite(penalty) .or. penalty < 0.0_dp) return
        if (any(.not. ieee_is_finite(plus_trace)) .or. &
            any(.not. ieee_is_finite(minus_trace)) .or. &
            any(.not. ieee_is_finite(plus_flux)) .or. &
            any(.not. ieee_is_finite(minus_flux)) .or. &
            any(.not. ieee_is_finite(surface_weights))) return
        if (any(surface_weights <= 0.0_dp)) return

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

end module fortfem_nitsche_interface
