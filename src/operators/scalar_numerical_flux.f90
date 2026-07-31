module fortfem_scalar_numerical_flux
    !! Conservative scalar central, upwind, and Lax--Friedrichs fluxes.
    !!
    !! The normal speed and Lax--Friedrichs dissipation are caller-owned
    !! coefficients.  The interface returns equal-and-opposite weighted
    !! residuals on the two sides and the quadratic-entropy production
    !! diagnostic w*alpha*(u_plus-u_minus)^2.  Upwind absolute values are
    !! differentiated only on a fixed nonzero-speed topology.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    integer, parameter, public :: NUMERICAL_FLUX_CENTRAL = 0
    integer, parameter, public :: NUMERICAL_FLUX_UPWIND = 1
    integer, parameter, public :: NUMERICAL_FLUX_LAX_FRIEDRICHS = 2

    public :: assemble_scalar_numerical_flux
    public :: assemble_scalar_numerical_flux_jvp
    public :: assemble_scalar_numerical_flux_vjp

contains

    subroutine assemble_scalar_numerical_flux( &
            plus_state, minus_state, normal_speed, dissipation, surface_weights, &
            flux_kind, plus_residual, minus_residual, entropy_production, status)
        !! Assemble weighted conservative scalar interface fluxes.
        real(dp), intent(in) :: plus_state(:), minus_state(:), normal_speed(:)
        real(dp), intent(in) :: dissipation(:), surface_weights(:)
        integer, intent(in) :: flux_kind
        real(dp), intent(out) :: plus_residual(:), minus_residual(:)
        real(dp), intent(out) :: entropy_production(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, q
        real(dp) :: jump, alpha, flux

        plus_residual = 0.0_dp
        minus_residual = 0.0_dp
        entropy_production = 0.0_dp
        call validate_flux_inputs( &
            plus_state, minus_state, normal_speed, dissipation, surface_weights, &
            flux_kind, plus_residual, minus_residual, entropy_production, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(plus_state)
        do q = 1, quadrature_count
            jump = plus_state(q) - minus_state(q)
            call select_dissipation( &
                flux_kind, normal_speed(q), dissipation(q), alpha, status)
            if (status%code /= FORTSPARSE_OK) return
            flux = 0.5_dp*normal_speed(q)*(plus_state(q) + minus_state(q)) - &
                0.5_dp*alpha*jump
            plus_residual(q) = surface_weights(q)*flux
            minus_residual(q) = -plus_residual(q)
            entropy_production(q) = surface_weights(q)*alpha*jump**2
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_scalar_numerical_flux

    subroutine assemble_scalar_numerical_flux_jvp( &
            plus_state, minus_state, normal_speed, dissipation, surface_weights, &
            flux_kind, plus_state_dot, minus_state_dot, normal_speed_dot, &
            dissipation_dot, surface_weights_dot, plus_residual_dot, &
            minus_residual_dot, entropy_production_dot, status)
        !! Apply the fixed-topology product-rule JVP of the numerical flux.
        real(dp), intent(in) :: plus_state(:), minus_state(:), normal_speed(:)
        real(dp), intent(in) :: dissipation(:), surface_weights(:)
        integer, intent(in) :: flux_kind
        real(dp), intent(in) :: plus_state_dot(:), minus_state_dot(:)
        real(dp), intent(in) :: normal_speed_dot(:), dissipation_dot(:)
        real(dp), intent(in) :: surface_weights_dot(:)
        real(dp), intent(out) :: plus_residual_dot(:), minus_residual_dot(:)
        real(dp), intent(out) :: entropy_production_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, q
        real(dp) :: jump, jump_dot, alpha, alpha_dot, flux, flux_dot

        plus_residual_dot = 0.0_dp
        minus_residual_dot = 0.0_dp
        entropy_production_dot = 0.0_dp
        call validate_flux_inputs( &
            plus_state, minus_state, normal_speed, dissipation, surface_weights, &
            flux_kind, plus_residual_dot, minus_residual_dot, &
            entropy_production_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(plus_state)
        if (.not. valid_flux_direction( &
            plus_state_dot, minus_state_dot, normal_speed_dot, dissipation_dot, &
            surface_weights_dot, quadrature_count)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "scalar numerical flux JVP has incompatible increments")
            return
        end if
        do q = 1, quadrature_count
            jump = plus_state(q) - minus_state(q)
            jump_dot = plus_state_dot(q) - minus_state_dot(q)
            call select_dissipation( &
                flux_kind, normal_speed(q), dissipation(q), alpha, status)
            if (status%code /= FORTSPARSE_OK) return
            call select_dissipation_dot( &
                flux_kind, normal_speed(q), normal_speed_dot(q), dissipation_dot(q), &
                alpha_dot)
            flux = 0.5_dp*normal_speed(q)*(plus_state(q) + minus_state(q)) - &
                0.5_dp*alpha*jump
            flux_dot = 0.5_dp*normal_speed_dot(q)*(plus_state(q) + minus_state(q)) + &
                0.5_dp*normal_speed(q)*(plus_state_dot(q) + minus_state_dot(q)) - &
                0.5_dp*alpha_dot*jump - 0.5_dp*alpha*jump_dot
            plus_residual_dot(q) = surface_weights_dot(q)*flux + &
                surface_weights(q)*flux_dot
            minus_residual_dot(q) = -plus_residual_dot(q)
            entropy_production_dot(q) = surface_weights_dot(q)*alpha*jump**2 + &
                surface_weights(q)*alpha_dot*jump**2 + &
                2.0_dp*surface_weights(q)*alpha*jump*jump_dot
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_scalar_numerical_flux_jvp

    subroutine assemble_scalar_numerical_flux_vjp( &
            plus_state, minus_state, normal_speed, dissipation, surface_weights, &
            flux_kind, plus_residual_bar, minus_residual_bar, entropy_bar, &
            plus_state_bar, minus_state_bar, normal_speed_bar, dissipation_bar, &
            surface_weights_bar, status)
        !! Apply the real reverse product of the numerical flux and entropy output.
        real(dp), intent(in) :: plus_state(:), minus_state(:), normal_speed(:)
        real(dp), intent(in) :: dissipation(:), surface_weights(:)
        integer, intent(in) :: flux_kind
        real(dp), intent(in) :: plus_residual_bar(:), minus_residual_bar(:)
        real(dp), intent(in) :: entropy_bar(:)
        real(dp), intent(out) :: plus_state_bar(:), minus_state_bar(:)
        real(dp), intent(out) :: normal_speed_bar(:), dissipation_bar(:)
        real(dp), intent(out) :: surface_weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, q
        real(dp) :: jump, alpha, flux, residual_bar, alpha_speed_factor

        plus_state_bar = 0.0_dp
        minus_state_bar = 0.0_dp
        normal_speed_bar = 0.0_dp
        dissipation_bar = 0.0_dp
        surface_weights_bar = 0.0_dp
        call validate_flux_inputs( &
            plus_state, minus_state, normal_speed, dissipation, surface_weights, &
            flux_kind, plus_residual_bar, minus_residual_bar, entropy_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(plus_state)
        if (size(plus_state_bar) /= quadrature_count .or. &
            size(minus_state_bar) /= quadrature_count .or. &
            size(normal_speed_bar) /= quadrature_count .or. &
            size(dissipation_bar) /= quadrature_count .or. &
            size(surface_weights_bar) /= quadrature_count .or. &
            any(.not. ieee_is_finite(plus_residual_bar)) .or. &
            any(.not. ieee_is_finite(minus_residual_bar)) .or. &
            any(.not. ieee_is_finite(entropy_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "scalar numerical flux VJP has incompatible cotangents")
            return
        end if
        do q = 1, quadrature_count
            jump = plus_state(q) - minus_state(q)
            call select_dissipation( &
                flux_kind, normal_speed(q), dissipation(q), alpha, status)
            if (status%code /= FORTSPARSE_OK) return
            call select_dissipation_speed_factor( &
                flux_kind, normal_speed(q), alpha_speed_factor)
            flux = 0.5_dp*normal_speed(q)*(plus_state(q) + minus_state(q)) - &
                0.5_dp*alpha*jump
            residual_bar = surface_weights(q)*( &
                plus_residual_bar(q) - minus_residual_bar(q))
            plus_state_bar(q) = residual_bar*( &
                0.5_dp*normal_speed(q) - 0.5_dp*alpha) + &
                entropy_bar(q)*2.0_dp*surface_weights(q)*alpha*jump
            minus_state_bar(q) = residual_bar*( &
                0.5_dp*normal_speed(q) + 0.5_dp*alpha) - &
                entropy_bar(q)*2.0_dp*surface_weights(q)*alpha*jump
            normal_speed_bar(q) = residual_bar*( &
                0.5_dp*(plus_state(q) + minus_state(q)) - &
                0.5_dp*alpha_speed_factor*jump) + &
                entropy_bar(q)*surface_weights(q)*alpha_speed_factor*jump**2
            dissipation_bar(q) = -0.5_dp*residual_bar*jump + &
                entropy_bar(q)*surface_weights(q)*jump**2
            if (flux_kind /= NUMERICAL_FLUX_LAX_FRIEDRICHS) then
                dissipation_bar(q) = 0.0_dp
            end if
            surface_weights_bar(q) = ( &
                plus_residual_bar(q) - minus_residual_bar(q))*flux + &
                entropy_bar(q)*alpha*jump**2
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_scalar_numerical_flux_vjp

    subroutine select_dissipation(flux_kind, speed, supplied, alpha, status)
        integer, intent(in) :: flux_kind
        real(dp), intent(in) :: speed, supplied
        real(dp), intent(out) :: alpha
        type(fortsparse_status_t), intent(inout) :: status

        alpha = 0.0_dp
        select case (flux_kind)
        case (NUMERICAL_FLUX_CENTRAL)
            return
        case (NUMERICAL_FLUX_UPWIND)
            if (speed == 0.0_dp) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "upwind numerical flux has a zero-speed topology event")
                return
            end if
            alpha = abs(speed)
        case (NUMERICAL_FLUX_LAX_FRIEDRICHS)
            alpha = supplied
        case default
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "scalar numerical flux received an unknown flux kind")
        end select
    end subroutine select_dissipation

    pure subroutine select_dissipation_dot( &
            flux_kind, speed, speed_dot, supplied_dot, alpha_dot)
        integer, intent(in) :: flux_kind
        real(dp), intent(in) :: speed, speed_dot, supplied_dot
        real(dp), intent(out) :: alpha_dot

        select case (flux_kind)
        case (NUMERICAL_FLUX_UPWIND)
            alpha_dot = sign(1.0_dp, speed)*speed_dot
        case (NUMERICAL_FLUX_LAX_FRIEDRICHS)
            alpha_dot = supplied_dot
        case default
            alpha_dot = 0.0_dp
        end select
    end subroutine select_dissipation_dot

    pure subroutine select_dissipation_speed_factor(flux_kind, speed, factor)
        integer, intent(in) :: flux_kind
        real(dp), intent(in) :: speed
        real(dp), intent(out) :: factor

        if (flux_kind == NUMERICAL_FLUX_UPWIND) then
            factor = sign(1.0_dp, speed)
        else
            factor = 0.0_dp
        end if
    end subroutine select_dissipation_speed_factor

    pure logical function valid_flux_direction( &
            plus_state_dot, minus_state_dot, normal_speed_dot, dissipation_dot, &
            surface_weights_dot, quadrature_count) result(valid)
        real(dp), intent(in) :: plus_state_dot(:), minus_state_dot(:)
        real(dp), intent(in) :: normal_speed_dot(:), dissipation_dot(:)
        real(dp), intent(in) :: surface_weights_dot(:)
        integer, intent(in) :: quadrature_count

        valid = size(plus_state_dot) == quadrature_count .and. &
            size(minus_state_dot) == quadrature_count .and. &
            size(normal_speed_dot) == quadrature_count .and. &
            size(dissipation_dot) == quadrature_count .and. &
            size(surface_weights_dot) == quadrature_count .and. &
            all(ieee_is_finite(plus_state_dot)) .and. &
            all(ieee_is_finite(minus_state_dot)) .and. &
            all(ieee_is_finite(normal_speed_dot)) .and. &
            all(ieee_is_finite(dissipation_dot)) .and. &
            all(ieee_is_finite(surface_weights_dot))
    end function valid_flux_direction

    subroutine validate_flux_inputs( &
            plus_state, minus_state, normal_speed, dissipation, surface_weights, &
            flux_kind, plus_residual, minus_residual, entropy_production, status)
        real(dp), intent(in) :: plus_state(:), minus_state(:), normal_speed(:)
        real(dp), intent(in) :: dissipation(:), surface_weights(:)
        integer, intent(in) :: flux_kind
        real(dp), intent(in) :: plus_residual(:), minus_residual(:), entropy_production(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "scalar numerical flux received incompatible arrays")
        quadrature_count = size(plus_state)
        if (quadrature_count < 1 .or. size(minus_state) /= quadrature_count .or. &
            size(normal_speed) /= quadrature_count .or. &
            size(dissipation) /= quadrature_count .or. &
            size(surface_weights) /= quadrature_count .or. &
            size(plus_residual) /= quadrature_count .or. &
            size(minus_residual) /= quadrature_count .or. &
            size(entropy_production) /= quadrature_count) return
        if (flux_kind < NUMERICAL_FLUX_CENTRAL .or. &
            flux_kind > NUMERICAL_FLUX_LAX_FRIEDRICHS .or. &
            any(surface_weights <= 0.0_dp) .or. &
            any(dissipation < 0.0_dp) .or. &
            any(.not. ieee_is_finite(plus_state)) .or. &
            any(.not. ieee_is_finite(minus_state)) .or. &
            any(.not. ieee_is_finite(normal_speed)) .or. &
            any(.not. ieee_is_finite(dissipation)) .or. &
            any(.not. ieee_is_finite(surface_weights)) .or. &
            any(.not. ieee_is_finite(plus_residual)) .or. &
            any(.not. ieee_is_finite(minus_residual)) .or. &
            any(.not. ieee_is_finite(entropy_production))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_flux_inputs

end module fortfem_scalar_numerical_flux
