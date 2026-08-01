module fortfem_fitted_trace_constraint
    !! Signed multiplier coupling for fitted duplicated interface traces.
    !!
    !! The block is the weak constraint
    !!
    !!   integral lambda * (trace_plus - trace_minus) dS.
    !!
    !! Multiplier, plus, and minus traces may all have different numbers of
    !! degrees of freedom.  Their quadrature rows are the only shared
    !! topology.  Constitutive jump laws, metric projections, and global
    !! numbering remain caller-owned.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_fitted_trace_constraint
    public :: assemble_fitted_trace_constraint_jvp
    public :: assemble_fitted_trace_constraint_vjp

contains

    subroutine assemble_fitted_trace_constraint( &
            multiplier_trace, plus_trace, minus_trace, surface_weights, matrix, &
            status)
        !! Assemble the signed plus/minus multiplier cross-mass block.
        real(dp), intent(in) :: multiplier_trace(:, :), plus_trace(:, :)
        real(dp), intent(in) :: minus_trace(:, :), surface_weights(:)
        real(dp), intent(out) :: matrix(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature, multiplier, plus, minus
        integer :: multiplier_count, plus_count, minus_count
        real(dp) :: scale

        matrix = 0.0_dp
        call validate_inputs( &
            multiplier_trace, plus_trace, minus_trace, surface_weights, matrix, &
            status)
        if (status%code /= FORTSPARSE_OK) return
        multiplier_count = size(multiplier_trace, 2)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        do quadrature = 1, size(surface_weights)
            scale = surface_weights(quadrature)
            do multiplier = 1, multiplier_count
                do plus = 1, plus_count
                    matrix(multiplier, plus) = matrix(multiplier, plus) + &
                        scale*multiplier_trace(quadrature, multiplier)* &
                        plus_trace(quadrature, plus)
                end do
                do minus = 1, minus_count
                    matrix(multiplier, plus_count + minus) = &
                        matrix(multiplier, plus_count + minus) - &
                        scale*multiplier_trace(quadrature, multiplier)* &
                        minus_trace(quadrature, minus)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fitted_trace_constraint

    subroutine assemble_fitted_trace_constraint_jvp( &
            multiplier_trace, plus_trace, minus_trace, surface_weights, &
            multiplier_trace_dot, plus_trace_dot, minus_trace_dot, &
            surface_weights_dot, matrix_dot, status)
        !! Apply the product-rule JVP of the signed trace constraint.
        real(dp), intent(in) :: multiplier_trace(:, :), plus_trace(:, :)
        real(dp), intent(in) :: minus_trace(:, :), surface_weights(:)
        real(dp), intent(in) :: multiplier_trace_dot(:, :), plus_trace_dot(:, :)
        real(dp), intent(in) :: minus_trace_dot(:, :), surface_weights_dot(:)
        real(dp), intent(out) :: matrix_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature, multiplier, plus, minus
        integer :: multiplier_count, plus_count, minus_count
        real(dp) :: scale, scale_dot, multiplier_value, multiplier_dot

        matrix_dot = 0.0_dp
        call validate_inputs( &
            multiplier_trace, plus_trace, minus_trace, surface_weights, matrix_dot, &
            status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_direction( &
            multiplier_trace_dot, plus_trace_dot, minus_trace_dot, &
            surface_weights_dot, multiplier_trace, plus_trace, minus_trace, &
            surface_weights)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "fitted trace constraint JVP has incompatible increments")
            return
        end if
        multiplier_count = size(multiplier_trace, 2)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        do quadrature = 1, size(surface_weights)
            scale = surface_weights(quadrature)
            scale_dot = surface_weights_dot(quadrature)
            do multiplier = 1, multiplier_count
                multiplier_value = multiplier_trace(quadrature, multiplier)
                multiplier_dot = multiplier_trace_dot(quadrature, multiplier)
                do plus = 1, plus_count
                    matrix_dot(multiplier, plus) = &
                        matrix_dot(multiplier, plus) + scale_dot*multiplier_value* &
                        plus_trace(quadrature, plus) + scale*multiplier_dot* &
                        plus_trace(quadrature, plus) + scale*multiplier_value* &
                        plus_trace_dot(quadrature, plus)
                end do
                do minus = 1, minus_count
                    matrix_dot(multiplier, plus_count + minus) = &
                        matrix_dot(multiplier, plus_count + minus) - &
                        scale_dot*multiplier_value*minus_trace(quadrature, minus) - &
                        scale*multiplier_dot*minus_trace(quadrature, minus) - &
                        scale*multiplier_value*minus_trace_dot(quadrature, minus)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fitted_trace_constraint_jvp

    subroutine assemble_fitted_trace_constraint_vjp( &
            multiplier_trace, plus_trace, minus_trace, surface_weights, matrix_bar, &
            multiplier_trace_bar, plus_trace_bar, minus_trace_bar, &
            surface_weights_bar, status)
        !! Apply the real reverse product of the signed trace constraint.
        real(dp), intent(in) :: multiplier_trace(:, :), plus_trace(:, :)
        real(dp), intent(in) :: minus_trace(:, :), surface_weights(:)
        real(dp), intent(in) :: matrix_bar(:, :)
        real(dp), intent(out) :: multiplier_trace_bar(:, :), plus_trace_bar(:, :)
        real(dp), intent(out) :: minus_trace_bar(:, :), surface_weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature, multiplier, plus, minus
        integer :: multiplier_count, plus_count, minus_count
        real(dp) :: scale, multiplier_value, plus_value, minus_value

        multiplier_trace_bar = 0.0_dp
        plus_trace_bar = 0.0_dp
        minus_trace_bar = 0.0_dp
        surface_weights_bar = 0.0_dp
        call validate_inputs( &
            multiplier_trace, plus_trace, minus_trace, surface_weights, matrix_bar, &
            status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(multiplier_trace_bar, 1) /= size(multiplier_trace, 1) .or. &
            size(multiplier_trace_bar, 2) /= size(multiplier_trace, 2) .or. &
            size(plus_trace_bar, 1) /= size(plus_trace, 1) .or. &
            size(plus_trace_bar, 2) /= size(plus_trace, 2) .or. &
            size(minus_trace_bar, 1) /= size(minus_trace, 1) .or. &
            size(minus_trace_bar, 2) /= size(minus_trace, 2) .or. &
            size(surface_weights_bar) /= size(surface_weights)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "fitted trace constraint VJP has incompatible cotangents")
            return
        end if
        multiplier_count = size(multiplier_trace, 2)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        do quadrature = 1, size(surface_weights)
            scale = surface_weights(quadrature)
            do multiplier = 1, multiplier_count
                multiplier_value = multiplier_trace(quadrature, multiplier)
                do plus = 1, plus_count
                    plus_value = plus_trace(quadrature, plus)
                    multiplier_trace_bar(quadrature, multiplier) = &
                        multiplier_trace_bar(quadrature, multiplier) + &
                        scale*matrix_bar(multiplier, plus)*plus_value
                    plus_trace_bar(quadrature, plus) = &
                        plus_trace_bar(quadrature, plus) + &
                        scale*matrix_bar(multiplier, plus)*multiplier_value
                    surface_weights_bar(quadrature) = &
                        surface_weights_bar(quadrature) + &
                        matrix_bar(multiplier, plus)*multiplier_value*plus_value
                end do
                do minus = 1, minus_count
                    minus_value = minus_trace(quadrature, minus)
                    multiplier_trace_bar(quadrature, multiplier) = &
                        multiplier_trace_bar(quadrature, multiplier) - &
                        scale*matrix_bar(multiplier, plus_count + minus)*minus_value
                    minus_trace_bar(quadrature, minus) = &
                        minus_trace_bar(quadrature, minus) - &
                        scale*matrix_bar(multiplier, plus_count + minus)* &
                        multiplier_value
                    surface_weights_bar(quadrature) = &
                        surface_weights_bar(quadrature) - &
                        matrix_bar(multiplier, plus_count + minus)*multiplier_value* &
                        minus_value
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fitted_trace_constraint_vjp

    subroutine validate_inputs( &
            multiplier_trace, plus_trace, minus_trace, surface_weights, matrix, &
            status)
        real(dp), intent(in) :: multiplier_trace(:, :), plus_trace(:, :)
        real(dp), intent(in) :: minus_trace(:, :), surface_weights(:), matrix(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: quadrature_count, multiplier_count, plus_count, minus_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "fitted trace constraint received incompatible arrays")
        quadrature_count = size(multiplier_trace, 1)
        multiplier_count = size(multiplier_trace, 2)
        plus_count = size(plus_trace, 2)
        minus_count = size(minus_trace, 2)
        if (quadrature_count < 1 .or. multiplier_count < 1 .or. &
            plus_count < 1 .or. minus_count < 1) return
        if (size(plus_trace, 1) /= quadrature_count .or. &
            size(minus_trace, 1) /= quadrature_count .or. &
            size(surface_weights) /= quadrature_count .or. &
            size(matrix, 1) /= multiplier_count .or. &
            size(matrix, 2) /= plus_count + minus_count) return
        if (any(.not. ieee_is_finite(multiplier_trace)) .or. &
            any(.not. ieee_is_finite(plus_trace)) .or. &
            any(.not. ieee_is_finite(minus_trace)) .or. &
            any(.not. ieee_is_finite(surface_weights)) .or. &
            any(.not. ieee_is_finite(matrix)) .or. &
            any(surface_weights <= 0.0_dp)) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_inputs

    pure logical function valid_direction( &
            multiplier_trace_dot, plus_trace_dot, minus_trace_dot, &
            surface_weights_dot, multiplier_trace, plus_trace, minus_trace, &
            surface_weights) result(valid)
        real(dp), intent(in) :: multiplier_trace_dot(:, :), plus_trace_dot(:, :)
        real(dp), intent(in) :: minus_trace_dot(:, :), surface_weights_dot(:)
        real(dp), intent(in) :: multiplier_trace(:, :), plus_trace(:, :)
        real(dp), intent(in) :: minus_trace(:, :), surface_weights(:)

        valid = all(shape(multiplier_trace_dot) == shape(multiplier_trace)) .and. &
            all(shape(plus_trace_dot) == shape(plus_trace)) .and. &
            all(shape(minus_trace_dot) == shape(minus_trace)) .and. &
            size(surface_weights_dot) == size(surface_weights) .and. &
            all(ieee_is_finite(multiplier_trace_dot)) .and. &
            all(ieee_is_finite(plus_trace_dot)) .and. &
            all(ieee_is_finite(minus_trace_dot)) .and. &
            all(ieee_is_finite(surface_weights_dot))
    end function valid_direction

end module fortfem_fitted_trace_constraint
