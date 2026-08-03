program test_flux_surface_average
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_flux_surface_average, &
        evaluate_flux_surface_average_jvp, &
        evaluate_flux_surface_average_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t
    implicit none

    integer, parameter :: component_count = 2, sample_count = 4
    real(dp), parameter :: epsilon_fd = 1.0e-7_dp
    real(dp), parameter :: tolerance = 2.0e-12_dp
    real(dp), parameter :: derivative_tolerance = 3.0e-8_dp
    real(dp) :: samples(component_count, sample_count)
    real(dp) :: samples_dot(component_count, sample_count)
    real(dp) :: samples_plus(component_count, sample_count)
    real(dp) :: samples_minus(component_count, sample_count)
    real(dp) :: weights(sample_count), weights_dot(sample_count)
    real(dp) :: weights_plus(sample_count), weights_minus(sample_count)
    real(dp) :: factors(sample_count), factors_dot(sample_count)
    real(dp) :: factors_plus(sample_count), factors_minus(sample_count)
    real(dp) :: average(component_count), average_dot(component_count)
    real(dp) :: average_plus(component_count), average_minus(component_count)
    real(dp) :: samples_bar(component_count, sample_count)
    real(dp) :: weights_bar(sample_count), factors_bar(sample_count)
    real(dp) :: average_bar(component_count)
    real(dp) :: denominator, denominator_dot, denominator_plus, denominator_minus
    real(dp) :: denominator_bar, lhs, rhs, effective_weights(sample_count)
    real(dp) :: expected_average(component_count), expected_dot(component_count)
    real(dp) :: expected_denominator, effective_weights_dot(sample_count)
    real(dp) :: scalar_sample(1, sample_count), scalar_average(1)
    real(dp) :: scalar_denominator
    real(dp) :: empty_samples(1, 0), empty_weights(0), empty_average(1)
    type(fortsparse_status_t) :: status
    logical :: all_passed
    integer :: component, sample

    all_passed = .true.
    samples = reshape([ &
        0.4_dp, -0.7_dp, 1.1_dp, 0.2_dp, -0.3_dp, 0.9_dp, 0.8_dp, -1.2_dp], &
        shape(samples))
    samples_dot = reshape([ &
        -0.2_dp, 0.3_dp, 0.4_dp, -0.1_dp, 0.5_dp, -0.6_dp, 0.7_dp, 0.8_dp], &
        shape(samples_dot))
    weights = [1.2_dp, 0.7_dp, 1.6_dp, 0.9_dp]
    weights_dot = [0.1_dp, -0.2_dp, 0.3_dp, -0.15_dp]
    factors = [0.8_dp, 1.4_dp, 0.6_dp, 1.1_dp]
    factors_dot = [0.05_dp, -0.08_dp, 0.04_dp, 0.06_dp]
    effective_weights = weights*factors
    effective_weights_dot = weights_dot*factors + weights*factors_dot
    expected_denominator = sum(effective_weights)
    expected_average = 0.0_dp
    expected_dot = 0.0_dp
    do component = 1, component_count
        do sample = 1, sample_count
            expected_average(component) = expected_average(component) + &
                effective_weights(sample)*samples(component, sample)
            expected_dot(component) = expected_dot(component) + &
                effective_weights(sample)*samples_dot(component, sample) + &
                effective_weights_dot(sample)*samples(component, sample)
        end do
        expected_average(component) = expected_average(component)/expected_denominator
    end do
    expected_dot = (expected_dot - expected_average*sum(effective_weights_dot))/ &
        expected_denominator

    call evaluate_flux_surface_average( &
        samples, weights, average, denominator, status, factors)
    call record_condition(status%code == FORTSPARSE_OK .and. &
        maxval(abs(average - expected_average)) < tolerance .and. &
        abs(denominator - expected_denominator) < tolerance, &
        "flux-surface average matches independent weighted oracle")

    call evaluate_flux_surface_average_jvp( &
        samples, weights, samples_dot, weights_dot, average_dot, denominator_dot, &
        status, factors, factors_dot)
    samples_plus = samples + epsilon_fd*samples_dot
    samples_minus = samples - epsilon_fd*samples_dot
    weights_plus = weights + epsilon_fd*weights_dot
    weights_minus = weights - epsilon_fd*weights_dot
    factors_plus = factors + epsilon_fd*factors_dot
    factors_minus = factors - epsilon_fd*factors_dot
    call evaluate_flux_surface_average( &
        samples_plus, weights_plus, average_plus, denominator_plus, status, &
        factors_plus)
    call evaluate_flux_surface_average( &
        samples_minus, weights_minus, average_minus, denominator_minus, status, &
        factors_minus)
    call record_condition(status%code == FORTSPARSE_OK .and. &
        maxval(abs(average_dot - expected_dot)) < tolerance .and. &
        maxval(abs(average_dot - (average_plus - average_minus)/ &
        (2.0_dp*epsilon_fd))) < derivative_tolerance .and. &
        abs(denominator_dot - (denominator_plus - denominator_minus)/ &
        (2.0_dp*epsilon_fd)) < derivative_tolerance, &
        "flux-surface average JVP matches analytic and central-difference oracles")

    average_bar = [0.6_dp, -0.4_dp]
    denominator_bar = 0.35_dp
    call evaluate_flux_surface_average_vjp( &
        samples, weights, average_bar, denominator_bar, samples_bar, weights_bar, &
        status, factors, factors_bar)
    lhs = dot_product(average_bar, average_dot) + denominator_bar*denominator_dot
    rhs = sum(samples_bar*samples_dot) + dot_product(weights_bar, weights_dot) + &
        dot_product(factors_bar, factors_dot)
    call record_condition(status%code == FORTSPARSE_OK .and. &
        abs(lhs - rhs) < derivative_tolerance, &
        "flux-surface average VJP satisfies the real dot-product identity")

    scalar_sample = reshape([0.2_dp, -0.5_dp, 0.9_dp, 1.1_dp], shape(scalar_sample))
    call evaluate_flux_surface_average( &
        scalar_sample, weights, scalar_average, scalar_denominator, status)
    call record_condition(status%code == FORTSPARSE_OK .and. &
        abs(scalar_average(1) - sum(scalar_sample(1, :)*weights)/ &
        sum(weights)) < tolerance, &
        "scalar samples use the same component-major contract without factors")

    call evaluate_flux_surface_average( &
        empty_samples, empty_weights, empty_average, scalar_denominator, status)
    call record_condition(status%code == FORTSPARSE_INVALID_MATRIX, &
        "empty surface samples are rejected")

    factors(2) = 0.0_dp
    call evaluate_flux_surface_average( &
        samples, weights, average, denominator, status, factors)
    call record_condition(status%code == FORTSPARSE_INVALID_MATRIX, &
        "non-positive geometric factors are rejected")

    call check_summary("Flux-surface average")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_flux_surface_average
