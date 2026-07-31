module fortfem_feec_exact_sequence
    !! Differentiable algebraic diagnostics for a discrete de Rham sequence.
    !!
    !! The incidence maps are metric-independent and may come from simplicial
    !! FEEC, tensor-product IGA, multipatch quotients, or periodic complexes.
    !! The composition defects expose curl(grad) and div(curl) without hiding
    !! them in a solver or a material/Hodge matrix.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_feec_exact_sequence
    public :: assemble_feec_exact_sequence_jvp
    public :: assemble_feec_exact_sequence_vjp

contains

    subroutine assemble_feec_exact_sequence( &
            gradient, curl, divergence, curl_gradient, divergence_curl, status)
        !! Return the two algebraic exact-sequence composition defects.
        real(dp), intent(in) :: gradient(:, :), curl(:, :), divergence(:, :)
        real(dp), intent(out) :: curl_gradient(:, :), divergence_curl(:, :)
        type(fortsparse_status_t), intent(out) :: status

        curl_gradient = 0.0_dp
        divergence_curl = 0.0_dp
        call validate_sequence_inputs( &
            gradient, curl, divergence, curl_gradient, divergence_curl, status)
        if (status%code /= FORTSPARSE_OK) return
        curl_gradient = matmul(curl, gradient)
        divergence_curl = matmul(divergence, curl)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_feec_exact_sequence

    subroutine assemble_feec_exact_sequence_jvp( &
            gradient, curl, divergence, gradient_dot, curl_dot, divergence_dot, &
            curl_gradient_dot, divergence_curl_dot, status)
        !! Apply the product-rule JVP of both composition defects.
        real(dp), intent(in) :: gradient(:, :), curl(:, :), divergence(:, :)
        real(dp), intent(in) :: gradient_dot(:, :), curl_dot(:, :), divergence_dot(:, :)
        real(dp), intent(out) :: curl_gradient_dot(:, :), divergence_curl_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        curl_gradient_dot = 0.0_dp
        divergence_curl_dot = 0.0_dp
        call validate_sequence_inputs( &
            gradient, curl, divergence, curl_gradient_dot, divergence_curl_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_sequence_direction( &
            gradient_dot, curl_dot, divergence_dot, gradient, curl, divergence)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FEEC exact-sequence JVP has incompatible increments")
            return
        end if
        curl_gradient_dot = matmul(curl_dot, gradient) + matmul(curl, gradient_dot)
        divergence_curl_dot = matmul(divergence_dot, curl) + &
            matmul(divergence, curl_dot)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_feec_exact_sequence_jvp

    subroutine assemble_feec_exact_sequence_vjp( &
            gradient, curl, divergence, curl_gradient_bar, divergence_curl_bar, &
            gradient_bar, curl_bar, divergence_bar, status)
        !! Apply the real reverse product of both composition defects.
        real(dp), intent(in) :: gradient(:, :), curl(:, :), divergence(:, :)
        real(dp), intent(in) :: curl_gradient_bar(:, :), divergence_curl_bar(:, :)
        real(dp), intent(out) :: gradient_bar(:, :), curl_bar(:, :), divergence_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        gradient_bar = 0.0_dp
        curl_bar = 0.0_dp
        divergence_bar = 0.0_dp
        call validate_sequence_inputs( &
            gradient, curl, divergence, curl_gradient_bar, divergence_curl_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(gradient_bar, 1) /= size(gradient, 1) .or. &
            size(gradient_bar, 2) /= size(gradient, 2) .or. &
            size(curl_bar, 1) /= size(curl, 1) .or. &
            size(curl_bar, 2) /= size(curl, 2) .or. &
            size(divergence_bar, 1) /= size(divergence, 1) .or. &
            size(divergence_bar, 2) /= size(divergence, 2)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FEEC exact-sequence VJP has incompatible cotangents")
            return
        end if
        gradient_bar = matmul(transpose(curl), curl_gradient_bar)
        curl_bar = matmul(curl_gradient_bar, transpose(gradient)) + &
            matmul(transpose(divergence), divergence_curl_bar)
        divergence_bar = matmul(divergence_curl_bar, transpose(curl))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_feec_exact_sequence_vjp

    pure logical function valid_sequence_direction( &
            gradient_dot, curl_dot, divergence_dot, gradient, curl, divergence) result(valid)
        real(dp), intent(in) :: gradient_dot(:, :), curl_dot(:, :), divergence_dot(:, :)
        real(dp), intent(in) :: gradient(:, :), curl(:, :), divergence(:, :)

        valid = all(shape(gradient_dot) == shape(gradient)) .and. &
            all(shape(curl_dot) == shape(curl)) .and. &
            all(shape(divergence_dot) == shape(divergence)) .and. &
            all(ieee_is_finite(gradient_dot)) .and. &
            all(ieee_is_finite(curl_dot)) .and. &
            all(ieee_is_finite(divergence_dot))
    end function valid_sequence_direction

    subroutine validate_sequence_inputs( &
            gradient, curl, divergence, curl_gradient, divergence_curl, status)
        real(dp), intent(in) :: gradient(:, :), curl(:, :), divergence(:, :)
        real(dp), intent(in) :: curl_gradient(:, :), divergence_curl(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: scalar_count, hcurl_count, hdiv_count, l2_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FEEC exact-sequence maps have incompatible dimensions")
        hcurl_count = size(gradient, 1)
        scalar_count = size(gradient, 2)
        hdiv_count = size(curl, 1)
        l2_count = size(divergence, 1)
        if (scalar_count < 1 .or. hcurl_count < 1 .or. hdiv_count < 1 .or. &
            l2_count < 1 .or. size(curl, 2) /= hcurl_count .or. &
            size(divergence, 2) /= hdiv_count .or. &
            size(curl_gradient, 1) /= hdiv_count .or. &
            size(curl_gradient, 2) /= scalar_count .or. &
            size(divergence_curl, 1) /= l2_count .or. &
            size(divergence_curl, 2) /= hcurl_count) return
        if (any(.not. ieee_is_finite(gradient)) .or. &
            any(.not. ieee_is_finite(curl)) .or. &
            any(.not. ieee_is_finite(divergence)) .or. &
            any(.not. ieee_is_finite(curl_gradient)) .or. &
            any(.not. ieee_is_finite(divergence_curl))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_sequence_inputs

end module fortfem_feec_exact_sequence
