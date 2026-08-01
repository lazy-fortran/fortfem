program test_wave_reflection_diagnostics_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_weighted_complex_error, &
        evaluate_weighted_complex_error_jvp, &
        evaluate_weighted_complex_error_vjp, &
        evaluate_weighted_reflection_coefficient, &
        evaluate_weighted_reflection_coefficient_jvp, &
        evaluate_weighted_reflection_coefficient_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: sample_count = 7
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: coordinates(sample_count), weights(sample_count)
    real(dp) :: weights_dot(sample_count), weights_bar(sample_count)
    complex(dp) :: incident(1, sample_count), total(1, sample_count)
    complex(dp) :: incident_dot(1, sample_count), total_dot(1, sample_count)
    complex(dp) :: incident_bar(1, sample_count), total_bar(1, sample_count)
    complex(dp) :: reference_bar(1, sample_count), candidate_bar(1, sample_count)
    complex(dp) :: reflected(1, sample_count)
    real(dp) :: absolute_error, relative_error, absolute_error_dot
    real(dp) :: relative_error_dot, coefficient, coefficient_dot
    real(dp) :: coefficient_plus, coefficient_minus, expected, lhs, rhs
    real(dp) :: absolute_plus, absolute_minus
    real(dp) :: reflection_bar, step_wave, wave_number
    integer :: index, status, status_plus, status_minus
    logical :: all_passed

    all_passed = .true.
    wave_number = 1.3_dp
    coordinates = [0.0_dp, 0.2_dp, 0.5_dp, 0.9_dp, 1.4_dp, 2.0_dp, 2.7_dp]
    weights = [0.4_dp, 0.7_dp, 0.5_dp, 0.9_dp, 0.6_dp, 0.8_dp, 0.3_dp]
    weights_dot = [0.01_dp, -0.02_dp, 0.015_dp, 0.01_dp, &
        -0.01_dp, 0.02_dp, -0.015_dp]
    do index = 1, sample_count
        incident(1, index) = exp(cmplx(0.0_dp, wave_number*coordinates(index), dp))
        reflected(1, index) = cmplx(0.2_dp, 0.1_dp, dp)*exp( &
            cmplx(0.0_dp, -wave_number*coordinates(index), dp))
        total(1, index) = incident(1, index) + reflected(1, index)
        incident_dot(1, index) = cmplx(0.01_dp*index, -0.004_dp*index, dp)
        total_dot(1, index) = cmplx(-0.003_dp*index, 0.008_dp*index, dp)
    end do

    call evaluate_weighted_complex_error( &
        incident, total, weights, absolute_error, relative_error, status)
    expected = sqrt(0.2_dp**2 + 0.1_dp**2)
    call record_condition(status == 0, "weighted complex error accepts wave samples")
    call record_condition(abs(absolute_error - expected*sqrt(sum(weights))) < &
        2.0e-13_dp, "weighted error matches the plane-wave reflection oracle")
    call record_condition(abs(relative_error - expected) < 2.0e-13_dp, &
        "weighted relative error is the reflection coefficient")

    call evaluate_weighted_complex_error_jvp( &
        incident, total, weights, incident_dot, total_dot, weights_dot, &
        absolute_error_dot, relative_error_dot, status)
    step_wave = step
    call evaluate_weighted_reflection_coefficient( &
        incident - step_wave*incident_dot, total - step_wave*total_dot, &
        weights - step_wave*weights_dot, coefficient_minus, status_minus)
    call evaluate_weighted_reflection_coefficient( &
        incident + step_wave*incident_dot, total + step_wave*total_dot, &
        weights + step_wave*weights_dot, coefficient_plus, status_plus)
    call evaluate_weighted_complex_error( &
        incident - step_wave*incident_dot, total - step_wave*total_dot, &
        weights - step_wave*weights_dot, absolute_minus, coefficient_minus, &
        status_minus)
    call evaluate_weighted_complex_error( &
        incident + step_wave*incident_dot, total + step_wave*total_dot, &
        weights + step_wave*weights_dot, absolute_plus, coefficient_plus, &
        status_plus)
    call record_condition(status == 0 .and. status_minus == 0 .and. &
        status_plus == 0, "weighted reflection JVP states assemble")
    call record_condition(abs(relative_error_dot - &
        (coefficient_plus - coefficient_minus)/(2.0_dp*step_wave)) < 2.0e-8_dp, &
        "weighted reflection JVP matches central differences")
    call record_condition(abs(absolute_error_dot - &
        (absolute_plus - absolute_minus)/(2.0_dp*step_wave)) < 2.0e-8_dp, &
        "weighted absolute-error JVP matches central differences")

    reflection_bar = 1.0_dp
    call evaluate_weighted_reflection_coefficient_vjp( &
        incident, total, weights, reflection_bar, incident_bar, total_bar, &
        weights_bar, status)
    lhs = relative_error_dot
    rhs = real(sum(conjg(incident_bar)*incident_dot) + &
        sum(conjg(total_bar)*total_dot), dp) + sum(weights_bar*weights_dot)
    call record_condition(status == 0, "weighted reflection VJP succeeds")
    call record_condition(abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "weighted reflection JVP and VJP satisfy the real complex adjoint identity")

    call evaluate_weighted_complex_error_vjp( &
        incident, total, weights, 0.7_dp, -0.2_dp, reference_bar, &
        candidate_bar, weights_bar, status)
    call record_condition(status == 0, "weighted error VJP accepts both output bars")

    weights(3) = 0.0_dp
    call evaluate_weighted_reflection_coefficient( &
        incident, total, weights, coefficient, status)
    call record_condition(status /= 0, &
        "weighted reflection rejects nonpositive weights")

    call check_summary("Weighted wave reflection diagnostics")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_wave_reflection_diagnostics_ad
