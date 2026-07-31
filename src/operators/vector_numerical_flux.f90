module fortfem_vector_numerical_flux
    !! Conservative numerical fluxes for vector-valued states.
    !!
    !! A caller-supplied component metric maps the normal average and jump
    !! terms.  The residuals on the two sides are equal and opposite, while
    !! the returned quadratic entropy production is
    !! w*alpha*jump^T*M*jump.  Upwind absolute values are differentiated only
    !! on a fixed nonzero-speed topology.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_scalar_numerical_flux, only: &
        NUMERICAL_FLUX_CENTRAL, NUMERICAL_FLUX_UPWIND, &
        NUMERICAL_FLUX_LAX_FRIEDRICHS
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_vector_numerical_flux
    public :: assemble_vector_numerical_flux_jvp
    public :: assemble_vector_numerical_flux_vjp
    public :: assemble_vector_entropy_stable_flux
    public :: assemble_vector_entropy_stable_flux_jvp
    public :: assemble_vector_entropy_stable_flux_vjp

contains

    subroutine assemble_vector_entropy_stable_flux( &
            plus_state, minus_state, normal_speed, dissipation, component_metric, &
            surface_weights, flux_kind, plus_residual, minus_residual, &
            entropy_production, status)
        !! Assemble a vector flux with an SPD entropy metric.
        real(dp), intent(in) :: plus_state(:, :), minus_state(:, :)
        real(dp), intent(in) :: normal_speed(:), dissipation(:)
        real(dp), intent(in) :: component_metric(:, :, :), surface_weights(:)
        integer, intent(in) :: flux_kind
        real(dp), intent(out) :: plus_residual(:, :), minus_residual(:, :)
        real(dp), intent(out) :: entropy_production(:)
        type(fortsparse_status_t), intent(out) :: status

        call validate_entropy_metric(component_metric, status)
        if (status%code /= FORTSPARSE_OK) then
            plus_residual = 0.0_dp
            minus_residual = 0.0_dp
            entropy_production = 0.0_dp
            return
        end if
        call assemble_vector_numerical_flux( &
            plus_state, minus_state, normal_speed, dissipation, component_metric, &
            surface_weights, flux_kind, plus_residual, minus_residual, &
            entropy_production, status)
    end subroutine assemble_vector_entropy_stable_flux

    subroutine assemble_vector_entropy_stable_flux_jvp( &
            plus_state, minus_state, normal_speed, dissipation, component_metric, &
            surface_weights, flux_kind, plus_state_dot, minus_state_dot, &
            normal_speed_dot, dissipation_dot, component_metric_dot, &
            surface_weights_dot, plus_residual_dot, minus_residual_dot, &
            entropy_production_dot, status)
        !! Apply the fixed-SPD-topology JVP of the entropy-stable flux.
        real(dp), intent(in) :: plus_state(:, :), minus_state(:, :)
        real(dp), intent(in) :: normal_speed(:), dissipation(:)
        real(dp), intent(in) :: component_metric(:, :, :), surface_weights(:)
        integer, intent(in) :: flux_kind
        real(dp), intent(in) :: plus_state_dot(:, :), minus_state_dot(:, :)
        real(dp), intent(in) :: normal_speed_dot(:), dissipation_dot(:)
        real(dp), intent(in) :: component_metric_dot(:, :, :)
        real(dp), intent(in) :: surface_weights_dot(:)
        real(dp), intent(out) :: plus_residual_dot(:, :), minus_residual_dot(:, :)
        real(dp), intent(out) :: entropy_production_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        call validate_entropy_metric(component_metric, status)
        if (status%code /= FORTSPARSE_OK) then
            plus_residual_dot = 0.0_dp
            minus_residual_dot = 0.0_dp
            entropy_production_dot = 0.0_dp
            return
        end if
        call assemble_vector_numerical_flux_jvp( &
            plus_state, minus_state, normal_speed, dissipation, component_metric, &
            surface_weights, flux_kind, plus_state_dot, minus_state_dot, &
            normal_speed_dot, dissipation_dot, component_metric_dot, &
            surface_weights_dot, plus_residual_dot, minus_residual_dot, &
            entropy_production_dot, status)
    end subroutine assemble_vector_entropy_stable_flux_jvp

    subroutine assemble_vector_entropy_stable_flux_vjp( &
            plus_state, minus_state, normal_speed, dissipation, component_metric, &
            surface_weights, flux_kind, plus_residual_bar, minus_residual_bar, &
            entropy_bar, plus_state_bar, minus_state_bar, normal_speed_bar, &
            dissipation_bar, component_metric_bar, surface_weights_bar, status)
        !! Apply the fixed-SPD-topology VJP of the entropy-stable flux.
        real(dp), intent(in) :: plus_state(:, :), minus_state(:, :)
        real(dp), intent(in) :: normal_speed(:), dissipation(:)
        real(dp), intent(in) :: component_metric(:, :, :), surface_weights(:)
        integer, intent(in) :: flux_kind
        real(dp), intent(in) :: plus_residual_bar(:, :), minus_residual_bar(:, :)
        real(dp), intent(in) :: entropy_bar(:)
        real(dp), intent(out) :: plus_state_bar(:, :), minus_state_bar(:, :)
        real(dp), intent(out) :: normal_speed_bar(:), dissipation_bar(:)
        real(dp), intent(out) :: component_metric_bar(:, :, :)
        real(dp), intent(out) :: surface_weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        call validate_entropy_metric(component_metric, status)
        if (status%code /= FORTSPARSE_OK) then
            plus_state_bar = 0.0_dp
            minus_state_bar = 0.0_dp
            normal_speed_bar = 0.0_dp
            dissipation_bar = 0.0_dp
            component_metric_bar = 0.0_dp
            surface_weights_bar = 0.0_dp
            return
        end if
        call assemble_vector_numerical_flux_vjp( &
            plus_state, minus_state, normal_speed, dissipation, component_metric, &
            surface_weights, flux_kind, plus_residual_bar, minus_residual_bar, &
            entropy_bar, plus_state_bar, minus_state_bar, normal_speed_bar, &
            dissipation_bar, component_metric_bar, surface_weights_bar, status)
    end subroutine assemble_vector_entropy_stable_flux_vjp

    subroutine assemble_vector_numerical_flux( &
            plus_state, minus_state, normal_speed, dissipation, component_metric, &
            surface_weights, flux_kind, plus_residual, minus_residual, &
            entropy_production, status)
        !! Assemble weighted conservative vector interface fluxes.
        real(dp), intent(in) :: plus_state(:, :), minus_state(:, :)
        real(dp), intent(in) :: normal_speed(:), dissipation(:)
        real(dp), intent(in) :: component_metric(:, :, :), surface_weights(:)
        integer, intent(in) :: flux_kind
        real(dp), intent(out) :: plus_residual(:, :), minus_residual(:, :)
        real(dp), intent(out) :: entropy_production(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, component_count, q, a, b
        real(dp) :: jump(size(plus_state, 2)), average(size(plus_state, 2))
        real(dp) :: metric_jump(size(plus_state, 2))
        real(dp) :: inner(size(plus_state, 2)), flux(size(plus_state, 2))
        real(dp) :: alpha, quadratic

        plus_residual = 0.0_dp
        minus_residual = 0.0_dp
        entropy_production = 0.0_dp
        call validate_flux_inputs( &
            plus_state, minus_state, normal_speed, dissipation, component_metric, &
            surface_weights, flux_kind, plus_residual, minus_residual, &
            entropy_production, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(plus_state, 1)
        component_count = size(plus_state, 2)
        do q = 1, quadrature_count
            jump = plus_state(q, :) - minus_state(q, :)
            average = 0.5_dp*(plus_state(q, :) + minus_state(q, :))
            call select_dissipation( &
                flux_kind, normal_speed(q), dissipation(q), alpha, status)
            if (status%code /= FORTSPARSE_OK) return
            inner = normal_speed(q)*average - 0.5_dp*alpha*jump
            metric_jump = 0.0_dp
            flux = 0.0_dp
            do a = 1, component_count
                do b = 1, component_count
                    metric_jump(a) = metric_jump(a) + &
                        component_metric(q, a, b)*jump(b)
                    flux(a) = flux(a) + component_metric(q, a, b)*inner(b)
                end do
            end do
            quadratic = dot_product(jump, metric_jump)
            plus_residual(q, :) = surface_weights(q)*flux
            minus_residual(q, :) = -plus_residual(q, :)
            entropy_production(q) = surface_weights(q)*alpha*quadratic
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_vector_numerical_flux

    subroutine assemble_vector_numerical_flux_jvp( &
            plus_state, minus_state, normal_speed, dissipation, component_metric, &
            surface_weights, flux_kind, plus_state_dot, minus_state_dot, &
            normal_speed_dot, dissipation_dot, component_metric_dot, &
            surface_weights_dot, plus_residual_dot, minus_residual_dot, &
            entropy_production_dot, status)
        !! Apply the fixed-topology product-rule JVP of the vector flux.
        real(dp), intent(in) :: plus_state(:, :), minus_state(:, :)
        real(dp), intent(in) :: normal_speed(:), dissipation(:)
        real(dp), intent(in) :: component_metric(:, :, :), surface_weights(:)
        integer, intent(in) :: flux_kind
        real(dp), intent(in) :: plus_state_dot(:, :), minus_state_dot(:, :)
        real(dp), intent(in) :: normal_speed_dot(:), dissipation_dot(:)
        real(dp), intent(in) :: component_metric_dot(:, :, :)
        real(dp), intent(in) :: surface_weights_dot(:)
        real(dp), intent(out) :: plus_residual_dot(:, :), minus_residual_dot(:, :)
        real(dp), intent(out) :: entropy_production_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, component_count, q, a, b
        real(dp) :: jump(size(plus_state, 2)), average(size(plus_state, 2))
        real(dp) :: jump_dot(size(plus_state, 2)), average_dot(size(plus_state, 2))
        real(dp) :: metric_jump(size(plus_state, 2))
        real(dp) :: metric_jump_dot(size(plus_state, 2))
        real(dp) :: inner(size(plus_state, 2)), inner_dot(size(plus_state, 2))
        real(dp) :: flux(size(plus_state, 2)), flux_dot(size(plus_state, 2))
        real(dp) :: alpha, alpha_dot, quadratic, quadratic_dot

        plus_residual_dot = 0.0_dp
        minus_residual_dot = 0.0_dp
        entropy_production_dot = 0.0_dp
        call validate_flux_inputs( &
            plus_state, minus_state, normal_speed, dissipation, component_metric, &
            surface_weights, flux_kind, plus_residual_dot, minus_residual_dot, &
            entropy_production_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(plus_state, 1)
        component_count = size(plus_state, 2)
        if (.not. valid_flux_direction( &
            plus_state_dot, minus_state_dot, normal_speed_dot, dissipation_dot, &
            component_metric_dot, surface_weights_dot, plus_residual_dot, &
            minus_residual_dot, entropy_production_dot, quadrature_count, &
            component_count)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "vector numerical flux JVP has incompatible increments")
            return
        end if
        do q = 1, quadrature_count
            jump = plus_state(q, :) - minus_state(q, :)
            average = 0.5_dp*(plus_state(q, :) + minus_state(q, :))
            jump_dot = plus_state_dot(q, :) - minus_state_dot(q, :)
            average_dot = 0.5_dp*(plus_state_dot(q, :) + minus_state_dot(q, :))
            call select_dissipation( &
                flux_kind, normal_speed(q), dissipation(q), alpha, status)
            if (status%code /= FORTSPARSE_OK) return
            call select_dissipation_dot( &
                flux_kind, normal_speed(q), normal_speed_dot(q), &
                dissipation_dot(q), alpha_dot)
            inner = normal_speed(q)*average - 0.5_dp*alpha*jump
            inner_dot = normal_speed_dot(q)*average + &
                normal_speed(q)*average_dot - 0.5_dp*alpha_dot*jump - &
                0.5_dp*alpha*jump_dot
            metric_jump = 0.0_dp
            metric_jump_dot = 0.0_dp
            flux = 0.0_dp
            flux_dot = 0.0_dp
            do a = 1, component_count
                do b = 1, component_count
                    metric_jump(a) = metric_jump(a) + &
                        component_metric(q, a, b)*jump(b)
                    metric_jump_dot(a) = metric_jump_dot(a) + &
                        component_metric_dot(q, a, b)*jump(b) + &
                        component_metric(q, a, b)*jump_dot(b)
                    flux(a) = flux(a) + component_metric(q, a, b)*inner(b)
                    flux_dot(a) = flux_dot(a) + &
                        component_metric_dot(q, a, b)*inner(b) + &
                        component_metric(q, a, b)*inner_dot(b)
                end do
            end do
            quadratic = dot_product(jump, metric_jump)
            quadratic_dot = dot_product(jump_dot, metric_jump) + &
                dot_product(jump, metric_jump_dot)
            plus_residual_dot(q, :) = surface_weights_dot(q)*flux + &
                surface_weights(q)*flux_dot
            minus_residual_dot(q, :) = -plus_residual_dot(q, :)
            entropy_production_dot(q) = surface_weights_dot(q)*alpha*quadratic + &
                surface_weights(q)*alpha_dot*quadratic + &
                surface_weights(q)*alpha*quadratic_dot
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_vector_numerical_flux_jvp

    subroutine assemble_vector_numerical_flux_vjp( &
            plus_state, minus_state, normal_speed, dissipation, component_metric, &
            surface_weights, flux_kind, plus_residual_bar, minus_residual_bar, &
            entropy_bar, plus_state_bar, minus_state_bar, normal_speed_bar, &
            dissipation_bar, component_metric_bar, surface_weights_bar, status)
        !! Apply the real reverse product of the vector flux and entropy output.
        real(dp), intent(in) :: plus_state(:, :), minus_state(:, :)
        real(dp), intent(in) :: normal_speed(:), dissipation(:)
        real(dp), intent(in) :: component_metric(:, :, :), surface_weights(:)
        integer, intent(in) :: flux_kind
        real(dp), intent(in) :: plus_residual_bar(:, :), minus_residual_bar(:, :)
        real(dp), intent(in) :: entropy_bar(:)
        real(dp), intent(out) :: plus_state_bar(:, :), minus_state_bar(:, :)
        real(dp), intent(out) :: normal_speed_bar(:), dissipation_bar(:)
        real(dp), intent(out) :: component_metric_bar(:, :, :)
        real(dp), intent(out) :: surface_weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, component_count, q, a, b
        real(dp) :: jump(size(plus_state, 2)), average(size(plus_state, 2))
        real(dp) :: metric_jump(size(plus_state, 2))
        real(dp) :: inner(size(plus_state, 2)), flux(size(plus_state, 2))
        real(dp) :: average_bar(size(plus_state, 2)), jump_bar(size(plus_state, 2))
        real(dp) :: flux_bar(size(plus_state, 2)), residual_bar(size(plus_state, 2))
        real(dp) :: alpha, alpha_bar, quadratic, alpha_speed_factor

        plus_state_bar = 0.0_dp
        minus_state_bar = 0.0_dp
        normal_speed_bar = 0.0_dp
        dissipation_bar = 0.0_dp
        component_metric_bar = 0.0_dp
        surface_weights_bar = 0.0_dp
        call validate_flux_inputs( &
            plus_state, minus_state, normal_speed, dissipation, component_metric, &
            surface_weights, flux_kind, plus_residual_bar, minus_residual_bar, &
            entropy_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        quadrature_count = size(plus_state, 1)
        component_count = size(plus_state, 2)
        if (size(plus_state_bar, 1) /= quadrature_count .or. &
            size(plus_state_bar, 2) /= component_count .or. &
            size(minus_state_bar, 1) /= quadrature_count .or. &
            size(minus_state_bar, 2) /= component_count .or. &
            size(normal_speed_bar) /= quadrature_count .or. &
            size(dissipation_bar) /= quadrature_count .or. &
            size(component_metric_bar, 1) /= quadrature_count .or. &
            size(component_metric_bar, 2) /= component_count .or. &
            size(component_metric_bar, 3) /= component_count .or. &
            size(surface_weights_bar) /= quadrature_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "vector numerical flux VJP has incompatible cotangents")
            return
        end if
        do q = 1, quadrature_count
            jump = plus_state(q, :) - minus_state(q, :)
            average = 0.5_dp*(plus_state(q, :) + minus_state(q, :))
            call select_dissipation( &
                flux_kind, normal_speed(q), dissipation(q), alpha, status)
            if (status%code /= FORTSPARSE_OK) return
            call select_dissipation_speed_factor( &
                flux_kind, normal_speed(q), alpha_speed_factor)
            inner = normal_speed(q)*average - 0.5_dp*alpha*jump
            metric_jump = 0.0_dp
            flux = 0.0_dp
            do a = 1, component_count
                do b = 1, component_count
                    metric_jump(a) = metric_jump(a) + &
                        component_metric(q, a, b)*jump(b)
                    flux(a) = flux(a) + component_metric(q, a, b)*inner(b)
                end do
            end do
            quadratic = dot_product(jump, metric_jump)
            residual_bar = surface_weights(q)*( &
                plus_residual_bar(q, :) - minus_residual_bar(q, :))
            flux_bar = residual_bar
            average_bar = 0.0_dp
            jump_bar = 0.0_dp
            alpha_bar = -0.5_dp*dot_product(flux_bar, metric_jump) + &
                entropy_bar(q)*surface_weights(q)*quadratic
            do a = 1, component_count
                do b = 1, component_count
                    average_bar(b) = average_bar(b) + &
                        normal_speed(q)*component_metric(q, a, b)*flux_bar(a)
                    jump_bar(b) = jump_bar(b) - 0.5_dp*alpha* &
                        component_metric(q, a, b)*flux_bar(a)
                    component_metric_bar(q, a, b) = flux_bar(a)*inner(b) + &
                        entropy_bar(q)*surface_weights(q)*alpha*jump(a)*jump(b)
                    jump_bar(b) = jump_bar(b) + entropy_bar(q)*surface_weights(q)* &
                        alpha*(component_metric(q, a, b) + &
                        component_metric(q, b, a))*jump(a)
                end do
            end do
            plus_state_bar(q, :) = 0.5_dp*average_bar + jump_bar
            minus_state_bar(q, :) = 0.5_dp*average_bar - jump_bar
            normal_speed_bar(q) = dot_product(flux_bar, matvec( &
                component_metric(q, :, :), average)) + &
                alpha_bar*alpha_speed_factor
            if (flux_kind == NUMERICAL_FLUX_LAX_FRIEDRICHS) then
                dissipation_bar(q) = alpha_bar
            end if
            surface_weights_bar(q) = dot_product( &
                plus_residual_bar(q, :) - minus_residual_bar(q, :), flux) + &
                entropy_bar(q)*alpha*quadratic
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_vector_numerical_flux_vjp

    pure function matvec(matrix, vector) result(product)
        real(dp), intent(in) :: matrix(:, :), vector(:)
        real(dp) :: product(size(vector))
        integer :: a, b

        product = 0.0_dp
        do a = 1, size(vector)
            do b = 1, size(vector)
                product(a) = product(a) + matrix(a, b)*vector(b)
            end do
        end do
    end function matvec

    subroutine validate_entropy_metric(component_metric, status)
        !! Validate the fixed SPD topology required by the entropy contract.
        real(dp), intent(in) :: component_metric(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, component_count, q, a, b, k
        real(dp), allocatable :: factor(:, :)
        real(dp) :: value, accumulation, scale

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "entropy-stable vector flux requires an SPD component metric")
        quadrature_count = size(component_metric, 1)
        component_count = size(component_metric, 2)
        if (quadrature_count < 1 .or. component_count < 1 .or. &
            size(component_metric, 3) /= component_count .or. &
            any(.not. ieee_is_finite(component_metric))) return
        allocate(factor(component_count, component_count))
        do q = 1, quadrature_count
            do a = 1, component_count
                do b = a + 1, component_count
                    scale = max(1.0_dp, abs(component_metric(q, a, b)), &
                        abs(component_metric(q, b, a)))
                    if (abs(component_metric(q, a, b) - &
                        component_metric(q, b, a)) > 1.0e-12_dp*scale) return
                end do
            end do
            factor = 0.0_dp
            do a = 1, component_count
                do b = 1, a
                    accumulation = 0.0_dp
                    do k = 1, b - 1
                        accumulation = accumulation + factor(a, k)*factor(b, k)
                    end do
                    value = component_metric(q, a, b) - accumulation
                    if (a == b) then
                        if (value <= 0.0_dp) return
                        factor(a, b) = sqrt(value)
                    else
                        if (factor(b, b) <= 0.0_dp) return
                        factor(a, b) = value/factor(b, b)
                    end if
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_entropy_metric

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
                "vector numerical flux received an unknown flux kind")
        end select
    end subroutine select_dissipation

    pure subroutine select_dissipation_dot( &
            flux_kind, speed, speed_dot, supplied_dot, alpha_dot)
        integer, intent(in) :: flux_kind
        real(dp), intent(in) :: speed, speed_dot, supplied_dot
        real(dp), intent(out) :: alpha_dot

        select case (flux_kind)
        case (NUMERICAL_FLUX_UPWIND)
            ! The base speed is validated nonzero before this derivative is used.
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
            component_metric_dot, surface_weights_dot, plus_residual_dot, &
            minus_residual_dot, entropy_production_dot, quadrature_count, &
            component_count) result(valid)
        real(dp), intent(in) :: plus_state_dot(:, :), minus_state_dot(:, :)
        real(dp), intent(in) :: normal_speed_dot(:), dissipation_dot(:)
        real(dp), intent(in) :: component_metric_dot(:, :, :)
        real(dp), intent(in) :: surface_weights_dot(:)
        real(dp), intent(in) :: plus_residual_dot(:, :), minus_residual_dot(:, :)
        real(dp), intent(in) :: entropy_production_dot(:)
        integer, intent(in) :: quadrature_count, component_count

        valid = size(plus_state_dot, 1) == quadrature_count .and. &
            size(plus_state_dot, 2) == component_count .and. &
            size(minus_state_dot, 1) == quadrature_count .and. &
            size(minus_state_dot, 2) == component_count .and. &
            size(normal_speed_dot) == quadrature_count .and. &
            size(dissipation_dot) == quadrature_count .and. &
            size(component_metric_dot, 1) == quadrature_count .and. &
            size(component_metric_dot, 2) == component_count .and. &
            size(component_metric_dot, 3) == component_count .and. &
            size(surface_weights_dot) == quadrature_count .and. &
            size(plus_residual_dot, 1) == quadrature_count .and. &
            size(plus_residual_dot, 2) == component_count .and. &
            size(minus_residual_dot, 1) == quadrature_count .and. &
            size(minus_residual_dot, 2) == component_count .and. &
            size(entropy_production_dot) == quadrature_count .and. &
            all(ieee_is_finite(plus_state_dot)) .and. &
            all(ieee_is_finite(minus_state_dot)) .and. &
            all(ieee_is_finite(normal_speed_dot)) .and. &
            all(ieee_is_finite(dissipation_dot)) .and. &
            all(ieee_is_finite(component_metric_dot)) .and. &
            all(ieee_is_finite(surface_weights_dot))
    end function valid_flux_direction

    subroutine validate_flux_inputs( &
            plus_state, minus_state, normal_speed, dissipation, component_metric, &
            surface_weights, flux_kind, plus_residual, minus_residual, &
            entropy_production, status)
        real(dp), intent(in) :: plus_state(:, :), minus_state(:, :)
        real(dp), intent(in) :: normal_speed(:), dissipation(:)
        real(dp), intent(in) :: component_metric(:, :, :), surface_weights(:)
        integer, intent(in) :: flux_kind
        real(dp), intent(in) :: plus_residual(:, :), minus_residual(:, :)
        real(dp), intent(in) :: entropy_production(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature_count, component_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "vector numerical flux received incompatible arrays")
        quadrature_count = size(plus_state, 1)
        component_count = size(plus_state, 2)
        if (quadrature_count < 1 .or. component_count < 1 .or. &
            size(minus_state, 1) /= quadrature_count .or. &
            size(minus_state, 2) /= component_count .or. &
            size(normal_speed) /= quadrature_count .or. &
            size(dissipation) /= quadrature_count .or. &
            size(surface_weights) /= quadrature_count .or. &
            size(component_metric, 1) /= quadrature_count .or. &
            size(component_metric, 2) /= component_count .or. &
            size(component_metric, 3) /= component_count .or. &
            size(plus_residual, 1) /= quadrature_count .or. &
            size(plus_residual, 2) /= component_count .or. &
            size(minus_residual, 1) /= quadrature_count .or. &
            size(minus_residual, 2) /= component_count .or. &
            size(entropy_production) /= quadrature_count) return
        if (flux_kind < NUMERICAL_FLUX_CENTRAL .or. &
            flux_kind > NUMERICAL_FLUX_LAX_FRIEDRICHS .or. &
            any(surface_weights <= 0.0_dp) .or. any(dissipation < 0.0_dp) .or. &
            any(.not. ieee_is_finite(plus_state)) .or. &
            any(.not. ieee_is_finite(minus_state)) .or. &
            any(.not. ieee_is_finite(normal_speed)) .or. &
            any(.not. ieee_is_finite(dissipation)) .or. &
            any(.not. ieee_is_finite(component_metric)) .or. &
            any(.not. ieee_is_finite(surface_weights)) .or. &
            any(.not. ieee_is_finite(plus_residual)) .or. &
            any(.not. ieee_is_finite(minus_residual)) .or. &
            any(.not. ieee_is_finite(entropy_production))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_flux_inputs

end module fortfem_vector_numerical_flux
