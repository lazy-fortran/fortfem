program test_scalar_numerical_flux_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_scalar_numerical_flux, &
        assemble_scalar_numerical_flux_jvp, assemble_scalar_numerical_flux_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: plus_state(3) = [1.2_dp, -0.3_dp, 0.7_dp]
    real(dp), parameter :: minus_state(3) = [-0.4_dp, 0.8_dp, 0.1_dp]
    real(dp), parameter :: normal_speed(3) = [2.0_dp, -1.5_dp, 0.8_dp]
    real(dp), parameter :: dissipation(3) = [2.4_dp, 2.0_dp, 1.1_dp]
    real(dp), parameter :: surface_weights(3) = [1.0_dp, 2.0_dp, 0.5_dp]
    real(dp), parameter :: plus_state_dot(3) = [0.1_dp, -0.2_dp, 0.3_dp]
    real(dp), parameter :: minus_state_dot(3) = [-0.05_dp, 0.4_dp, -0.1_dp]
    real(dp), parameter :: normal_speed_dot(3) = [0.2_dp, -0.1_dp, 0.05_dp]
    real(dp), parameter :: dissipation_dot(3) = [0.1_dp, 0.3_dp, -0.2_dp]
    real(dp), parameter :: surface_weights_dot(3) = [0.05_dp, -0.1_dp, 0.2_dp]
    real(dp), parameter :: plus_residual_bar(3) = [0.4_dp, -0.3_dp, 0.7_dp]
    real(dp), parameter :: minus_residual_bar(3) = [-0.2_dp, 0.6_dp, -0.5_dp]
    real(dp), parameter :: entropy_bar(3) = [0.3_dp, -0.4_dp, 0.2_dp]
    real(dp) :: plus_residual(3), minus_residual(3), entropy(3)
    real(dp) :: plus_residual_dot(3), minus_residual_dot(3), entropy_dot(3)
    real(dp) :: plus_state_bar(3), minus_state_bar(3), normal_speed_bar(3)
    real(dp) :: dissipation_bar(3), surface_weights_bar(3), lhs, rhs
    real(dp) :: plus_expected(3), minus_expected(3), entropy_expected(3)
    type(fortsparse_status_t) :: status

    call assemble_scalar_numerical_flux( &
        plus_state, minus_state, normal_speed, dissipation, surface_weights, 0, &
        plus_residual, minus_residual, entropy, status)
    call oracle_flux(0, plus_expected, minus_expected, entropy_expected)
    call check_condition(status%code == 0 .and. &
        maxval(abs(plus_residual - plus_expected)) < 1.0e-14_dp .and. &
        maxval(abs(minus_residual - minus_expected)) < 1.0e-14_dp .and. &
        maxval(abs(entropy - entropy_expected)) < 1.0e-14_dp, &
        "central numerical flux matches the conservative oracle")

    call assemble_scalar_numerical_flux( &
        plus_state, minus_state, normal_speed, dissipation, surface_weights, 1, &
        plus_residual, minus_residual, entropy, status)
    call oracle_flux(1, plus_expected, minus_expected, entropy_expected)
    call check_condition(status%code == 0 .and. &
        maxval(abs(plus_residual - plus_expected)) < 1.0e-14_dp .and. &
        maxval(abs(minus_residual - minus_expected)) < 1.0e-14_dp .and. &
        maxval(abs(entropy - entropy_expected)) < 1.0e-14_dp .and. &
        maxval(abs(plus_residual + minus_residual)) < 1.0e-14_dp, &
        "upwind numerical flux is conservative and entropy producing")

    call assemble_scalar_numerical_flux_jvp( &
        plus_state, minus_state, normal_speed, dissipation, surface_weights, 1, &
        plus_state_dot, minus_state_dot, normal_speed_dot, dissipation_dot, &
        surface_weights_dot, plus_residual_dot, minus_residual_dot, entropy_dot, status)
    call check_condition(status%code == 0, &
        "upwind numerical flux JVP accepts state, speed, and weight directions")

    call assemble_scalar_numerical_flux_vjp( &
        plus_state, minus_state, normal_speed, dissipation, surface_weights, 1, &
        plus_residual_bar, minus_residual_bar, entropy_bar, plus_state_bar, &
        minus_state_bar, normal_speed_bar, dissipation_bar, surface_weights_bar, status)
    lhs = sum(plus_residual_bar*plus_residual_dot) + &
        sum(minus_residual_bar*minus_residual_dot) + sum(entropy_bar*entropy_dot)
    rhs = sum(plus_state_bar*plus_state_dot) + sum(minus_state_bar*minus_state_dot) + &
        sum(normal_speed_bar*normal_speed_dot) + sum(dissipation_bar*dissipation_dot) + &
        sum(surface_weights_bar*surface_weights_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "upwind numerical flux VJP satisfies the real dot-product identity")

    call assemble_scalar_numerical_flux( &
        plus_state, minus_state, normal_speed, dissipation, surface_weights, 2, &
        plus_residual, minus_residual, entropy, status)
    call oracle_flux(2, plus_expected, minus_expected, entropy_expected)
    call check_condition(status%code == 0 .and. &
        maxval(abs(plus_residual - plus_expected)) < 1.0e-14_dp, &
        "Lax-Friedrichs numerical flux accepts caller dissipation")

    call assemble_scalar_numerical_flux( &
        plus_state, minus_state, [0.0_dp, -1.5_dp, 0.8_dp], dissipation, &
        surface_weights, 1, plus_residual, minus_residual, entropy, status)
    call check_condition(status%code /= 0, &
        "upwind numerical flux rejects a zero-speed nondifferentiable event")
    call assemble_scalar_numerical_flux( &
        plus_state, minus_state, normal_speed, [-1.0_dp, 2.0_dp, 1.1_dp], &
        surface_weights, 2, plus_residual, minus_residual, entropy, status)
    call check_condition(status%code /= 0, &
        "Lax-Friedrichs numerical flux rejects negative dissipation")
    call check_summary("scalar numerical flux AD")

contains

    subroutine oracle_flux(flux_kind, plus_result, minus_result, entropy_result)
        integer, intent(in) :: flux_kind
        real(dp), intent(out) :: plus_result(:), minus_result(:), entropy_result(:)
        integer :: q
        real(dp) :: jump, alpha, flux

        do q = 1, 3
            jump = plus_state(q) - minus_state(q)
            if (flux_kind == 0) then
                alpha = 0.0_dp
            else if (flux_kind == 1) then
                alpha = abs(normal_speed(q))
            else
                alpha = dissipation(q)
            end if
            flux = 0.5_dp*normal_speed(q)*(plus_state(q) + minus_state(q)) - &
                0.5_dp*alpha*jump
            plus_result(q) = surface_weights(q)*flux
            minus_result(q) = -plus_result(q)
            entropy_result(q) = surface_weights(q)*alpha*jump**2
        end do
    end subroutine oracle_flux

end program test_scalar_numerical_flux_ad
