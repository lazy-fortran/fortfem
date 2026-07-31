program test_vector_entropy_flux_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_vector_entropy_stable_flux, &
        assemble_vector_entropy_stable_flux_jvp, &
        assemble_vector_entropy_stable_flux_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: plus_state(2, 2) = reshape([1.2_dp, -0.3_dp, 0.7_dp, 0.4_dp], [2, 2])
    real(dp), parameter :: minus_state(2, 2) = reshape([-0.4_dp, 0.8_dp, 0.1_dp, -0.2_dp], [2, 2])
    real(dp), parameter :: normal_speed(2) = [2.0_dp, -1.5_dp]
    real(dp), parameter :: dissipation(2) = [2.4_dp, 2.0_dp]
    real(dp), parameter :: surface_weights(2) = [1.0_dp, 2.0_dp]
    real(dp), parameter :: component_metric(2, 2, 2) = reshape([ &
        2.0_dp, 1.2_dp, 0.3_dp, -0.2_dp, 0.3_dp, -0.2_dp, &
        1.5_dp, 0.8_dp], [2, 2, 2])
    real(dp), parameter :: non_spd_metric(2, 2, 2) = reshape([ &
        1.0_dp, 1.0_dp, 2.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, &
        1.0_dp, 1.0_dp], [2, 2, 2])
    real(dp), parameter :: plus_state_dot(2, 2) = reshape([0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp], [2, 2])
    real(dp), parameter :: minus_state_dot(2, 2) = reshape([-0.05_dp, 0.4_dp, -0.1_dp, 0.2_dp], [2, 2])
    real(dp), parameter :: normal_speed_dot(2) = [0.2_dp, -0.1_dp]
    real(dp), parameter :: dissipation_dot(2) = [0.1_dp, 0.3_dp]
    real(dp), parameter :: surface_weights_dot(2) = [0.05_dp, -0.1_dp]
    real(dp), parameter :: component_metric_dot(2, 2, 2) = reshape([ &
        -0.1_dp, 0.2_dp, 0.05_dp, -0.4_dp, 0.3_dp, 0.1_dp, &
        -0.2_dp, 0.15_dp], [2, 2, 2])
    real(dp), parameter :: plus_residual_bar(2, 2) = reshape([0.4_dp, -0.3_dp, 0.7_dp, -0.2_dp], [2, 2])
    real(dp), parameter :: minus_residual_bar(2, 2) = reshape([-0.2_dp, 0.6_dp, -0.5_dp, 0.3_dp], [2, 2])
    real(dp), parameter :: entropy_bar(2) = [0.3_dp, -0.4_dp]
    real(dp) :: plus_residual(2, 2), minus_residual(2, 2), entropy(2)
    real(dp) :: plus_residual_dot(2, 2), minus_residual_dot(2, 2), entropy_dot(2)
    real(dp) :: plus_state_bar(2, 2), minus_state_bar(2, 2), normal_speed_bar(2)
    real(dp) :: dissipation_bar(2), surface_weights_bar(2), component_metric_bar(2, 2, 2)
    real(dp) :: plus_expected(2, 2), minus_expected(2, 2), entropy_expected(2)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_vector_entropy_stable_flux( &
        plus_state, minus_state, normal_speed, dissipation, component_metric, &
        surface_weights, 1, plus_residual, minus_residual, entropy, status)
    call oracle_flux(plus_expected, minus_expected, entropy_expected)
    call check_condition(status%code == 0 .and. &
        maxval(abs(plus_residual - plus_expected)) < 1.0e-14_dp .and. &
        maxval(abs(minus_residual - minus_expected)) < 1.0e-14_dp .and. &
        maxval(abs(entropy - entropy_expected)) < 1.0e-14_dp .and. &
        minval(entropy) >= 0.0_dp, &
        "SPD-metric entropy-stable flux matches the independent oracle")

    call assemble_vector_entropy_stable_flux_jvp( &
        plus_state, minus_state, normal_speed, dissipation, component_metric, &
        surface_weights, 1, plus_state_dot, minus_state_dot, normal_speed_dot, &
        dissipation_dot, component_metric_dot, surface_weights_dot, &
        plus_residual_dot, minus_residual_dot, entropy_dot, status)
    call check_condition(status%code == 0, &
        "entropy-stable vector flux JVP accepts metric directions")

    call assemble_vector_entropy_stable_flux_vjp( &
        plus_state, minus_state, normal_speed, dissipation, component_metric, &
        surface_weights, 1, plus_residual_bar, minus_residual_bar, entropy_bar, &
        plus_state_bar, minus_state_bar, normal_speed_bar, dissipation_bar, &
        component_metric_bar, surface_weights_bar, status)
    lhs = sum(plus_residual_bar*plus_residual_dot) + &
        sum(minus_residual_bar*minus_residual_dot) + sum(entropy_bar*entropy_dot)
    rhs = sum(plus_state_bar*plus_state_dot) + sum(minus_state_bar*minus_state_dot) + &
        sum(normal_speed_bar*normal_speed_dot) + sum(dissipation_bar*dissipation_dot) + &
        sum(component_metric_bar*component_metric_dot) + &
        sum(surface_weights_bar*surface_weights_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "entropy-stable vector flux VJP satisfies the real dot-product identity")

    call assemble_vector_entropy_stable_flux( &
        plus_state, minus_state, normal_speed, dissipation, non_spd_metric, &
        surface_weights, 1, plus_residual, minus_residual, entropy, status)
    call check_condition(status%code /= 0, &
        "entropy-stable vector flux rejects a non-SPD entropy metric")
    call check_summary("vector entropy flux AD")

contains

    subroutine oracle_flux(plus_result, minus_result, entropy_result)
        real(dp), intent(out) :: plus_result(:, :), minus_result(:, :), entropy_result(:)
        integer :: q, a, b
        real(dp) :: jump(2), average(2), metric_jump(2), flux(2), entropy_value

        do q = 1, 2
            jump = plus_state(q, :) - minus_state(q, :)
            average = 0.5_dp*(plus_state(q, :) + minus_state(q, :))
            metric_jump = 0.0_dp
            flux = 0.0_dp
            do a = 1, 2
                do b = 1, 2
                    metric_jump(a) = metric_jump(a) + &
                        component_metric(q, a, b)*jump(b)
                    flux(a) = flux(a) + component_metric(q, a, b)* &
                        (normal_speed(q)*average(b) - 0.5_dp*abs(normal_speed(q))*jump(b))
                end do
            end do
            entropy_value = dot_product(jump, metric_jump)
            plus_result(q, :) = surface_weights(q)*flux
            minus_result(q, :) = -plus_result(q, :)
            entropy_result(q) = surface_weights(q)*abs(normal_speed(q))*entropy_value
        end do
    end subroutine oracle_flux

end program test_vector_entropy_flux_ad
