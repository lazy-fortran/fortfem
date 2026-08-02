program test_regularized_surface_current_layer
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use check, only: check_condition, check_summary
    use fortfem_interop, only: evaluate_regularized_surface_current_layer, &
        evaluate_regularized_surface_current_layer_jvp, &
        evaluate_regularized_surface_current_layer_vjp, &
        evaluate_regularized_surface_current_integral
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: sample_count = 5, component_count = 2
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: distance(sample_count), current(sample_count, component_count)
    real(dp) :: distance_dot(sample_count), current_dot(sample_count, component_count)
    real(dp) :: volume_current(sample_count, component_count)
    real(dp) :: volume_current_dot(sample_count, component_count)
    real(dp) :: volume_current_plus(sample_count, component_count)
    real(dp) :: volume_current_minus(sample_count, component_count)
    real(dp) :: volume_current_bar(sample_count, component_count)
    real(dp) :: distance_bar(sample_count), current_bar(sample_count, component_count)
    real(dp) :: reference(sample_count, component_count)
    real(dp) :: reference_dot(sample_count, component_count)
    real(dp) :: weights(sample_count), integrated_current(component_count)
    real(dp) :: expected_integrated(component_count)
    real(dp) :: epsilon, epsilon_dot, epsilon_bar, normalization
    real(dp) :: expected_normalization
    real(dp) :: profile, profile_dot, pi, lhs, rhs, nan_value
    type(fortsparse_status_t) :: status
    integer :: sample, component

    epsilon = 0.4_dp
    epsilon_dot = -0.07_dp
    distance = [-0.8_dp, -0.4_dp, 0.0_dp, 0.4_dp, 0.8_dp]
    current = reshape([ &
        0.7_dp, -0.2_dp, 0.9_dp, 1.1_dp, -0.4_dp, &
        -0.3_dp, 0.8_dp, 0.5_dp, -0.6_dp, 1.2_dp], shape(current))
    distance_dot = [0.2_dp, -0.1_dp, 0.3_dp, -0.4_dp, 0.15_dp]
    current_dot = reshape([ &
        -0.2_dp, 0.4_dp, 0.1_dp, -0.5_dp, 0.6_dp, &
        0.3_dp, -0.1_dp, 0.7_dp, 0.2_dp, -0.4_dp], shape(current_dot))
    volume_current_bar = reshape([ &
        0.3_dp, -0.6_dp, 0.8_dp, 0.4_dp, -0.2_dp, &
        0.5_dp, 0.7_dp, -0.1_dp, 0.9_dp, -0.3_dp], &
        shape(volume_current_bar))
    pi = acos(-1.0_dp)

    do sample = 1, sample_count
        profile = exp(-(distance(sample)/epsilon)**2)/(sqrt(pi)*epsilon)
        profile_dot = profile*( &
            -2.0_dp*distance(sample)*distance_dot(sample)/epsilon**2 + &
            (-1.0_dp/epsilon + 2.0_dp*distance(sample)**2/epsilon**3)* &
            epsilon_dot)
        do component = 1, component_count
            reference(sample, component) = profile*current(sample, component)
            reference_dot(sample, component) = &
                profile*current_dot(sample, component) + &
                profile_dot*current(sample, component)
        end do
    end do

    call evaluate_regularized_surface_current_layer( &
        distance, current, epsilon, volume_current, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(volume_current - reference)) < 2.0e-15_dp, &
        "regularized current matches an independent Gaussian oracle")

    call evaluate_regularized_surface_current_layer_jvp( &
        distance, current, epsilon, distance_dot, current_dot, epsilon_dot, &
        volume_current_dot, status)
    call evaluate_regularized_surface_current_layer( &
        distance + step*distance_dot, current + step*current_dot, &
        epsilon + step*epsilon_dot, volume_current_plus, status)
    call evaluate_regularized_surface_current_layer( &
        distance - step*distance_dot, current - step*current_dot, &
        epsilon - step*epsilon_dot, volume_current_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(volume_current_dot - reference_dot)) < 3.0e-14_dp .and. &
        maxval(abs(volume_current_dot - &
        (volume_current_plus - volume_current_minus)/(2.0_dp*step))) < &
        2.0e-9_dp, &
        "regularized current JVP matches oracle and central difference")

    call evaluate_regularized_surface_current_layer_vjp( &
        distance, current, epsilon, volume_current_bar, distance_bar, &
        current_bar, epsilon_bar, status)
    lhs = sum(volume_current_bar*volume_current_dot)
    rhs = dot_product(distance_bar, distance_dot) + &
        sum(current_bar*current_dot) + epsilon_bar*epsilon_dot
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 2.0e-13_dp, &
        "regularized current VJP satisfies the real adjoint identity")

    weights = [0.2_dp, 0.4_dp, 0.4_dp, 0.4_dp, 0.2_dp]
    expected_integrated = 0.0_dp
    expected_normalization = 0.0_dp
    do sample = 1, sample_count
        expected_normalization = expected_normalization + weights(sample)* &
            exp(-(distance(sample)/epsilon)**2)/(sqrt(pi)*epsilon)
    end do
    do component = 1, component_count
        expected_integrated(component) = &
            dot_product(weights, reference(:, component))
    end do
    call evaluate_regularized_surface_current_integral( &
        distance, weights, current, epsilon, normalization, &
        integrated_current, status)
    call check_condition(status%code == 0 .and. &
        abs(normalization - expected_normalization) < 2.0e-15_dp .and. &
        maxval(abs(integrated_current - expected_integrated)) < 2.0e-15_dp, &
        "regularized current diagnostic reports profile and current integrals")

    call evaluate_regularized_surface_current_layer( &
        distance, current, 0.0_dp, volume_current, status)
    call check_condition(status%code /= 0, &
        "regularized current rejects zero thickness")

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    distance(2) = nan_value
    call evaluate_regularized_surface_current_layer( &
        distance, current, epsilon, volume_current, status)
    call check_condition(status%code /= 0, &
        "regularized current rejects non-finite distances")

    call check_summary("regularized surface-current layer")
end program test_regularized_surface_current_layer
