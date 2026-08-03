program test_singular_layer_matching
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        assemble_singular_layer_matching, &
        assemble_singular_layer_matching_jvp, &
        assemble_singular_layer_matching_vjp
    implicit none

    complex(dp) :: outer_trace(2, 2), inner_trace(2, 1)
    complex(dp) :: outer_state(2), inner_state(1), jump(2)
    complex(dp) :: residual(2), residual_plus(2), residual_minus(2)
    complex(dp) :: outer_trace_dot(2, 2), inner_trace_dot(2, 1)
    complex(dp) :: outer_state_dot(2), inner_state_dot(1), jump_dot(2)
    complex(dp) :: residual_dot(2), residual_bar(2)
    complex(dp) :: outer_trace_bar(2, 2), inner_trace_bar(2, 1)
    complex(dp) :: outer_state_bar(2), inner_state_bar(1), jump_bar(2)
    real(dp) :: weights(2), weights_dot(2), weights_bar(2)
    real(dp) :: lhs, rhs, step
    integer :: status, status_plus, status_minus
    logical :: all_passed

    all_passed = .true.
    outer_trace = reshape([ &
        cmplx(1.0_dp, 0.2_dp, dp), cmplx(-0.4_dp, 0.1_dp, dp), &
        cmplx(0.7_dp, -0.3_dp, dp), cmplx(1.3_dp, 0.4_dp, dp)], [2, 2])
    inner_trace = reshape([ &
        cmplx(0.5_dp, -0.2_dp, dp), cmplx(-0.8_dp, 0.3_dp, dp)], [2, 1])
    outer_state = [cmplx(0.3_dp, -0.6_dp, dp), cmplx(-0.2_dp, 0.5_dp, dp)]
    inner_state = [cmplx(0.4_dp, 0.1_dp, dp)]
    jump = [cmplx(0.1_dp, -0.2_dp, dp), cmplx(-0.3_dp, 0.2_dp, dp)]
    weights = [1.5_dp, 2.0_dp]

    call assemble_singular_layer_matching( &
        outer_trace, inner_trace, weights, outer_state, inner_state, jump, &
        residual, status)
    call record_condition(status == 0, "singular-layer residual assembles")
    call record_condition(maxval(abs(residual - weights*( &
        matmul(outer_trace, outer_state) - matmul(inner_trace, inner_state) - &
        jump))) < 1.0e-14_dp, "residual matches the independent trace oracle")

    outer_trace_dot = reshape([ &
        cmplx(-0.2_dp, 0.1_dp, dp), cmplx(0.3_dp, -0.4_dp, dp), &
        cmplx(0.5_dp, 0.2_dp, dp), cmplx(-0.1_dp, 0.6_dp, dp)], [2, 2])
    inner_trace_dot = reshape([ &
        cmplx(0.2_dp, 0.4_dp, dp), cmplx(-0.3_dp, 0.1_dp, dp)], [2, 1])
    outer_state_dot = [cmplx(-0.7_dp, 0.2_dp, dp), cmplx(0.1_dp, -0.5_dp, dp)]
    inner_state_dot = [cmplx(0.6_dp, -0.3_dp, dp)]
    jump_dot = [cmplx(0.2_dp, 0.1_dp, dp), cmplx(-0.4_dp, -0.2_dp, dp)]
    weights_dot = [0.2_dp, -0.3_dp]
    call assemble_singular_layer_matching_jvp( &
        outer_trace, inner_trace, weights, outer_state, inner_state, jump, &
        outer_trace_dot, inner_trace_dot, weights_dot, outer_state_dot, &
        inner_state_dot, jump_dot, residual_dot, status)
    call record_condition(status == 0, "singular-layer JVP assembles")

    step = 1.0e-6_dp
    call assemble_singular_layer_matching( &
        outer_trace + step*outer_trace_dot, inner_trace + step*inner_trace_dot, &
        weights + step*weights_dot, outer_state + step*outer_state_dot, &
        inner_state + step*inner_state_dot, jump + step*jump_dot, residual_plus, &
        status_plus)
    call assemble_singular_layer_matching( &
        outer_trace - step*outer_trace_dot, inner_trace - step*inner_trace_dot, &
        weights - step*weights_dot, outer_state - step*outer_state_dot, &
        inner_state - step*inner_state_dot, jump - step*jump_dot, residual_minus, &
        status_minus)
    call record_condition(status_plus == 0 .and. status_minus == 0, &
        "perturbed singular-layer residuals assemble")
    call record_condition(maxval(abs(residual_dot - &
        (residual_plus - residual_minus)/(2.0_dp*step))) < 2.0e-9_dp, &
        "singular-layer JVP matches central differences")

    residual_bar = [cmplx(0.8_dp, -0.4_dp, dp), cmplx(-0.2_dp, 0.7_dp, dp)]
    call assemble_singular_layer_matching_vjp( &
        outer_trace, inner_trace, weights, outer_state, inner_state, jump, &
        residual_bar, outer_trace_bar, inner_trace_bar, weights_bar, &
        outer_state_bar, inner_state_bar, jump_bar, status)
    call record_condition(status == 0, "singular-layer VJP assembles")
    lhs = real(sum(conjg(residual_bar)*residual_dot), dp)
    rhs = real(sum(conjg(outer_trace_bar)*outer_trace_dot), dp) + &
        real(sum(conjg(inner_trace_bar)*inner_trace_dot), dp) + &
        dot_product(weights_bar, weights_dot) + &
        real(sum(conjg(outer_state_bar)*outer_state_dot), dp) + &
        real(sum(conjg(inner_state_bar)*inner_state_dot), dp) + &
        real(sum(conjg(jump_bar)*jump_dot), dp)
    call record_condition(abs(lhs - rhs) < 2.0e-12_dp, &
        "singular-layer JVP and VJP satisfy the real complex adjoint identity")

    weights(1) = 0.0_dp
    call assemble_singular_layer_matching( &
        outer_trace, inner_trace, weights, outer_state, inner_state, jump, &
        residual, status)
    call record_condition(status /= 0, "non-positive singular-layer weight rejects")

    call check_summary("singular-layer trace matching")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_singular_layer_matching
