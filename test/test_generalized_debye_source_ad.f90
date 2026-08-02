program test_generalized_debye_source_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_generalized_debye_source_residual, &
        assemble_generalized_debye_source_residual_jvp, &
        assemble_generalized_debye_source_residual_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: tangent_count = 3, scalar_count = 2
    integer, parameter :: harmonic_count = 1, source_count = 2, period_count = 1
    real(dp), parameter :: step = 1.0e-6_dp
    complex(dp) :: gradient_lift(tangent_count, scalar_count)
    complex(dp) :: cogradient_lift(tangent_count, scalar_count)
    complex(dp) :: harmonic_basis(tangent_count, harmonic_count)
    complex(dp) :: source_operator(source_count, tangent_count)
    complex(dp) :: period_operator(period_count, tangent_count)
    complex(dp) :: gradient_coefficients(scalar_count)
    complex(dp) :: cogradient_coefficients(scalar_count)
    complex(dp) :: harmonic_coefficients(harmonic_count)
    complex(dp) :: source_target(source_count), period_target(period_count)
    complex(dp) :: gradient_lift_dot(tangent_count, scalar_count)
    complex(dp) :: cogradient_lift_dot(tangent_count, scalar_count)
    complex(dp) :: harmonic_basis_dot(tangent_count, harmonic_count)
    complex(dp) :: source_operator_dot(source_count, tangent_count)
    complex(dp) :: period_operator_dot(period_count, tangent_count)
    complex(dp) :: gradient_coefficients_dot(scalar_count)
    complex(dp) :: cogradient_coefficients_dot(scalar_count)
    complex(dp) :: harmonic_coefficients_dot(harmonic_count)
    complex(dp) :: source_target_dot(source_count), period_target_dot(period_count)
    complex(dp) :: surface_current(tangent_count), source_residual(source_count)
    complex(dp) :: period_residual(period_count)
    complex(dp) :: surface_current_dot(tangent_count), source_residual_dot(source_count)
    complex(dp) :: period_residual_dot(period_count)
    complex(dp) :: surface_current_plus(tangent_count), source_residual_plus(source_count)
    complex(dp) :: period_residual_plus(period_count)
    complex(dp) :: surface_current_minus(tangent_count), source_residual_minus(source_count)
    complex(dp) :: period_residual_minus(period_count)
    complex(dp) :: source_residual_bar(source_count), period_residual_bar(period_count)
    complex(dp) :: surface_current_bar(tangent_count)
    complex(dp) :: gradient_lift_bar(tangent_count, scalar_count)
    complex(dp) :: cogradient_lift_bar(tangent_count, scalar_count)
    complex(dp) :: harmonic_basis_bar(tangent_count, harmonic_count)
    complex(dp) :: source_operator_bar(source_count, tangent_count)
    complex(dp) :: period_operator_bar(period_count, tangent_count)
    complex(dp) :: gradient_coefficients_bar(scalar_count)
    complex(dp) :: cogradient_coefficients_bar(scalar_count)
    complex(dp) :: harmonic_coefficients_bar(harmonic_count)
    complex(dp) :: source_target_bar(source_count), period_target_bar(period_count)
    real(dp) :: derivative_error, adjoint_error, lhs, rhs
    integer :: status
    logical :: all_passed

    all_passed = .true.
    gradient_lift = reshape([ &
        cmplx(1.0_dp, 0.1_dp, dp), cmplx(-0.2_dp, 0.3_dp, dp), &
        cmplx(0.4_dp, -0.1_dp, dp), cmplx(0.7_dp, 0.2_dp, dp), &
        cmplx(-0.3_dp, 0.2_dp, dp), cmplx(0.5_dp, -0.4_dp, dp)], [3, 2])
    cogradient_lift = reshape([ &
        cmplx(0.2_dp, -0.3_dp, dp), cmplx(0.6_dp, 0.1_dp, dp), &
        cmplx(-0.5_dp, 0.2_dp, dp), cmplx(0.3_dp, -0.2_dp, dp), &
        cmplx(0.4_dp, 0.5_dp, dp), cmplx(-0.1_dp, 0.6_dp, dp)], [3, 2])
    harmonic_basis(:, 1) = [ &
        cmplx(0.3_dp, -0.2_dp, dp), cmplx(-0.4_dp, 0.1_dp, dp), &
        cmplx(0.2_dp, 0.5_dp, dp)]
    source_operator = reshape([ &
        cmplx(0.7_dp, 0.1_dp, dp), cmplx(-0.2_dp, 0.4_dp, dp), &
        cmplx(0.3_dp, -0.5_dp, dp), cmplx(-0.6_dp, 0.2_dp, dp), &
        cmplx(0.4_dp, 0.3_dp, dp), cmplx(0.1_dp, -0.7_dp, dp)], [2, 3])
    period_operator(1, :) = [ &
        cmplx(0.2_dp, 0.3_dp, dp), cmplx(-0.4_dp, 0.2_dp, dp), &
        cmplx(0.6_dp, -0.1_dp, dp)]
    gradient_coefficients = [cmplx(0.5_dp, -0.1_dp, dp), cmplx(-0.2_dp, 0.4_dp, dp)]
    cogradient_coefficients = [cmplx(-0.3_dp, 0.2_dp, dp), cmplx(0.6_dp, 0.1_dp, dp)]
    harmonic_coefficients = [cmplx(0.4_dp, -0.5_dp, dp)]
    source_target = [cmplx(0.2_dp, -0.1_dp, dp), cmplx(-0.4_dp, 0.3_dp, dp)]
    period_target = [cmplx(0.1_dp, 0.2_dp, dp)]

    gradient_lift_dot = cmplx(reshape([0.03_dp, -0.02_dp, 0.01_dp, 0.04_dp, &
        -0.01_dp, 0.02_dp], [3, 2]), 0.0_dp, dp)
    cogradient_lift_dot = cmplx(reshape([-0.02_dp, 0.01_dp, 0.03_dp, -0.01_dp, &
        0.04_dp, -0.03_dp], [3, 2]), 0.0_dp, dp)
    harmonic_basis_dot(:, 1) = cmplx([0.02_dp, -0.01_dp, 0.03_dp], 0.0_dp, dp)
    source_operator_dot = cmplx(reshape([0.01_dp, -0.02_dp, 0.03_dp, -0.01_dp, &
        0.02_dp, 0.01_dp], [2, 3]), 0.0_dp, dp)
    period_operator_dot(1, :) = cmplx([0.03_dp, -0.02_dp, 0.01_dp], 0.0_dp, dp)
    gradient_coefficients_dot = cmplx([0.02_dp, -0.03_dp], 0.0_dp, dp)
    cogradient_coefficients_dot = cmplx([-0.01_dp, 0.04_dp], 0.0_dp, dp)
    harmonic_coefficients_dot = cmplx([0.03_dp], 0.0_dp, dp)
    source_target_dot = cmplx([0.01_dp, -0.02_dp], 0.0_dp, dp)
    period_target_dot = cmplx([0.02_dp], 0.0_dp, dp)

    call assemble_generalized_debye_source_residual( &
        gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
        period_operator, gradient_coefficients, cogradient_coefficients, &
        harmonic_coefficients, source_target, period_target, surface_current, &
        source_residual, period_residual, status)
    call record_condition(status == 0, &
        "generalized Debye-source lift accepts explicit harmonic coordinates")
    call assemble_generalized_debye_source_residual_jvp( &
        gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
        period_operator, gradient_coefficients, cogradient_coefficients, &
        harmonic_coefficients, source_target, period_target, gradient_lift_dot, &
        cogradient_lift_dot, harmonic_basis_dot, source_operator_dot, &
        period_operator_dot, gradient_coefficients_dot, cogradient_coefficients_dot, &
        harmonic_coefficients_dot, source_target_dot, period_target_dot, &
        surface_current_dot, source_residual_dot, period_residual_dot, status)
    call assemble_generalized_debye_source_residual( &
        gradient_lift + step*gradient_lift_dot, cogradient_lift + step*cogradient_lift_dot, &
        harmonic_basis + step*harmonic_basis_dot, source_operator + step*source_operator_dot, &
        period_operator + step*period_operator_dot, gradient_coefficients + &
        step*gradient_coefficients_dot, cogradient_coefficients + step*cogradient_coefficients_dot, &
        harmonic_coefficients + step*harmonic_coefficients_dot, source_target + &
        step*source_target_dot, period_target + step*period_target_dot, surface_current_plus, &
        source_residual_plus, period_residual_plus, status)
    call assemble_generalized_debye_source_residual( &
        gradient_lift - step*gradient_lift_dot, cogradient_lift - step*cogradient_lift_dot, &
        harmonic_basis - step*harmonic_basis_dot, source_operator - step*source_operator_dot, &
        period_operator - step*period_operator_dot, gradient_coefficients - &
        step*gradient_coefficients_dot, cogradient_coefficients - step*cogradient_coefficients_dot, &
        harmonic_coefficients - step*harmonic_coefficients_dot, source_target - &
        step*source_target_dot, period_target - step*period_target_dot, surface_current_minus, &
        source_residual_minus, period_residual_minus, status)
    derivative_error = max( &
        maxval(abs(surface_current_dot - &
        (surface_current_plus - surface_current_minus)/(2.0_dp*step))), &
        max(maxval(abs(source_residual_dot - &
        (source_residual_plus - source_residual_minus)/(2.0_dp*step))), &
        maxval(abs(period_residual_dot - &
        (period_residual_plus - period_residual_minus)/(2.0_dp*step)))))
    call record_condition(status == 0 .and. derivative_error < 2.0e-7_dp, &
        "generalized Debye-source JVP matches independent central reassembly")

    source_residual_bar = [cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.1_dp, 0.5_dp, dp)]
    period_residual_bar = [cmplx(-0.4_dp, 0.3_dp, dp)]
    surface_current_bar = [cmplx(0.2_dp, 0.6_dp, dp), cmplx(-0.5_dp, 0.1_dp, dp), &
        cmplx(0.3_dp, -0.2_dp, dp)]
    call assemble_generalized_debye_source_residual_vjp( &
        gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
        period_operator, gradient_coefficients, cogradient_coefficients, &
        harmonic_coefficients, source_target, period_target, source_residual_bar, &
        period_residual_bar, surface_current_bar, gradient_lift_bar, cogradient_lift_bar, &
        harmonic_basis_bar, source_operator_bar, period_operator_bar, &
        gradient_coefficients_bar, cogradient_coefficients_bar, harmonic_coefficients_bar, &
        source_target_bar, period_target_bar, status)
    lhs = real(sum(conjg(surface_current_bar)*surface_current_dot) + &
        sum(conjg(source_residual_bar)*source_residual_dot) + &
        sum(conjg(period_residual_bar)*period_residual_dot), dp)
    rhs = real(sum(conjg(gradient_lift_bar)*gradient_lift_dot) + &
        sum(conjg(cogradient_lift_bar)*cogradient_lift_dot) + &
        sum(conjg(harmonic_basis_bar)*harmonic_basis_dot) + &
        sum(conjg(source_operator_bar)*source_operator_dot) + &
        sum(conjg(period_operator_bar)*period_operator_dot) + &
        sum(conjg(gradient_coefficients_bar)*gradient_coefficients_dot) + &
        sum(conjg(cogradient_coefficients_bar)*cogradient_coefficients_dot) + &
        sum(conjg(harmonic_coefficients_bar)*harmonic_coefficients_dot) + &
        sum(conjg(source_target_bar)*source_target_dot) + &
        sum(conjg(period_target_bar)*period_target_dot), dp)
    adjoint_error = abs(lhs - rhs)
    call record_condition(status == 0 .and. adjoint_error < &
        5.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "generalized Debye-source VJP satisfies the real complex adjoint")

    call check_summary("generalized Debye-source layer")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_generalized_debye_source_ad
