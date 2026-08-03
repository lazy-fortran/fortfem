program test_complex_interchange_sample_derivatives
    use check, only: check_condition, check_summary
    use fortfem_interop, only: &
        compare_complex_interchange_samples, &
        compare_complex_interchange_samples_jvp, &
        compare_complex_interchange_samples_vjp, &
        complex_interchange_sample_set_t, initialize_complex_interchange_samples
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: component_count = 2, sample_count = 3
    real(dp), parameter :: step = 1.0e-6_dp
    real(dp), parameter :: coordinates(2, sample_count) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, sample_count])
    complex(dp), parameter :: reference_values(component_count, sample_count) = reshape([ &
        cmplx(1.0_dp, 0.2_dp, dp), cmplx(2.0_dp, -0.3_dp, dp), &
        cmplx(0.5_dp, 0.7_dp, dp), cmplx(1.5_dp, 0.1_dp, dp), &
        cmplx(-0.4_dp, 0.6_dp, dp), cmplx(2.4_dp, -0.8_dp, dp)], &
        [component_count, sample_count])
    complex(dp), parameter :: candidate_values(component_count, sample_count) = reshape([ &
        cmplx(1.2_dp, 0.0_dp, dp), cmplx(1.7_dp, -0.1_dp, dp), &
        cmplx(0.2_dp, 0.9_dp, dp), cmplx(1.8_dp, 0.4_dp, dp), &
        cmplx(-0.1_dp, 0.4_dp, dp), cmplx(2.8_dp, -0.5_dp, dp)], &
        [component_count, sample_count])
    real(dp), parameter :: weights(sample_count) = [0.5_dp, 1.0_dp, 1.5_dp]
    complex(dp), parameter :: reference_values_dot(component_count, sample_count) = reshape([ &
        cmplx(-0.1_dp, 0.3_dp, dp), cmplx(0.2_dp, -0.2_dp, dp), &
        cmplx(0.3_dp, 0.1_dp, dp), cmplx(-0.4_dp, 0.2_dp, dp), &
        cmplx(0.5_dp, -0.5_dp, dp), cmplx(0.6_dp, 0.4_dp, dp)], &
        [component_count, sample_count])
    complex(dp), parameter :: candidate_values_dot(component_count, sample_count) = reshape([ &
        cmplx(0.4_dp, -0.1_dp, dp), cmplx(-0.3_dp, 0.5_dp, dp), &
        cmplx(0.2_dp, 0.6_dp, dp), cmplx(0.1_dp, -0.4_dp, dp), &
        cmplx(-0.5_dp, 0.3_dp, dp), cmplx(0.7_dp, -0.2_dp, dp)], &
        [component_count, sample_count])
    real(dp), parameter :: weights_dot(sample_count) = [0.1_dp, -0.2_dp, 0.3_dp]
    real(dp), parameter :: absolute_error_bar = 0.8_dp
    real(dp), parameter :: relative_error_bar = -0.6_dp

    type(complex_interchange_sample_set_t) :: reference, candidate
    type(complex_interchange_sample_set_t) :: reference_plus, reference_minus
    type(complex_interchange_sample_set_t) :: candidate_plus, candidate_minus
    real(dp) :: absolute_error, relative_error, maximum_error
    real(dp) :: absolute_error_dot, relative_error_dot
    real(dp) :: absolute_error_plus, relative_error_plus, maximum_error_plus
    real(dp) :: absolute_error_minus, relative_error_minus, maximum_error_minus
    complex(dp) :: reference_values_bar(component_count, sample_count)
    complex(dp) :: candidate_values_bar(component_count, sample_count)
    real(dp) :: weights_bar(sample_count)
    real(dp) :: lhs, rhs
    integer :: status, status_plus, status_minus

    call initialize_complex_interchange_samples( &
        reference, coordinates, reference_values, weights, "reference", &
        "analytic-v1", status)
    call initialize_complex_interchange_samples( &
        candidate, coordinates, candidate_values, weights, "candidate", &
        "adapter-v1", status)
    call check_condition(status == 0, &
        "complex derivative comparison accepts matching coordinates and weights")

    call compare_complex_interchange_samples( &
        reference, candidate, 1.0e-14_dp, absolute_error, relative_error, &
        maximum_error, status)
    call compare_complex_interchange_samples_jvp( &
        reference, candidate, 1.0e-14_dp, reference_values_dot, &
        candidate_values_dot, weights_dot, absolute_error_dot, relative_error_dot, &
        status)
    call initialize_complex_interchange_samples( &
        reference_plus, coordinates, reference_values + step*reference_values_dot, &
        weights + step*weights_dot, "reference", "analytic-v1", status_plus)
    call initialize_complex_interchange_samples( &
        reference_minus, coordinates, reference_values - step*reference_values_dot, &
        weights - step*weights_dot, "reference", "analytic-v1", status_minus)
    call initialize_complex_interchange_samples( &
        candidate_plus, coordinates, candidate_values + step*candidate_values_dot, &
        weights + step*weights_dot, "candidate", "adapter-v1", status_plus)
    call initialize_complex_interchange_samples( &
        candidate_minus, coordinates, candidate_values - step*candidate_values_dot, &
        weights - step*weights_dot, "candidate", "adapter-v1", status_minus)
    call compare_complex_interchange_samples( &
        reference_plus, candidate_plus, 1.0e-14_dp, absolute_error_plus, &
        relative_error_plus, maximum_error_plus, status_plus)
    call compare_complex_interchange_samples( &
        reference_minus, candidate_minus, 1.0e-14_dp, absolute_error_minus, &
        relative_error_minus, maximum_error_minus, status_minus)
    call check_condition(status == 0 .and. status_plus == 0 .and. status_minus == 0 .and. &
        abs(absolute_error_dot - (absolute_error_plus - absolute_error_minus)/ &
        (2.0_dp*step)) < 5.0e-8_dp .and. &
        abs(relative_error_dot - (relative_error_plus - relative_error_minus)/ &
        (2.0_dp*step)) < 5.0e-8_dp, &
        "complex sample comparison JVP matches central differences")

    call compare_complex_interchange_samples_vjp( &
        reference, candidate, 1.0e-14_dp, absolute_error_bar, relative_error_bar, &
        reference_values_bar, candidate_values_bar, weights_bar, status)
    lhs = absolute_error_bar*absolute_error_dot + relative_error_bar*relative_error_dot
    rhs = real(sum(conjg(reference_values_bar)*reference_values_dot) + &
        sum(conjg(candidate_values_bar)*candidate_values_dot), dp) + &
        sum(weights_bar*weights_dot)
    call check_condition(status == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "complex sample comparison VJP satisfies the real-part adjoint identity")

    call initialize_complex_interchange_samples( &
        candidate, coordinates, reference_values, weights, "candidate", &
        "zero-error", status)
    call compare_complex_interchange_samples_jvp( &
        reference, candidate, 1.0e-14_dp, reference_values_dot, &
        candidate_values_dot, weights_dot, absolute_error_dot, relative_error_dot, &
        status)
    call check_condition(status /= 0, &
        "complex sample comparison derivative rejects the zero-error point")
    call check_summary("differentiable complex interchange sample comparison")
end program test_complex_interchange_sample_derivatives
