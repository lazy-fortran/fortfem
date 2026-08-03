program test_complex_interchange_samples
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_interop, only: &
        compare_complex_interchange_samples, &
        complex_interchange_sample_set_t, &
        initialize_complex_interchange_samples, &
        validate_complex_interchange_samples
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: spatial_dimension = 3, component_count = 2, sample_count = 3
    type(complex_interchange_sample_set_t) :: reference, candidate, copy, invalid
    real(dp) :: coordinates(spatial_dimension, sample_count), weights(sample_count)
    complex(dp) :: values(component_count, sample_count), candidate_values(component_count, sample_count)
    real(dp) :: absolute_error, relative_error, maximum_error
    integer :: status

    coordinates = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp], shape(coordinates))
    weights = [0.5_dp, 1.0_dp, 1.5_dp]
    values = reshape([ &
        cmplx(1.0_dp, 0.2_dp, dp), cmplx(-0.4_dp, 0.1_dp, dp), &
        cmplx(0.3_dp, -0.5_dp, dp), cmplx(0.2_dp, 0.8_dp, dp), &
        cmplx(-0.7_dp, 0.4_dp, dp), cmplx(0.6_dp, -0.3_dp, dp)], shape(values))
    candidate_values = values
    candidate_values(1, 2) = candidate_values(1, 2) + cmplx(0.1_dp, -0.2_dp, dp)

    call initialize_complex_interchange_samples( &
        reference, coordinates, values, weights, "analytic", "complex-oracle", status)
    call check_condition(status == 0, "complex samples initialize")
    call initialize_complex_interchange_samples( &
        candidate, coordinates, candidate_values, weights, "fem", "same-grid", status)
    call check_condition(status == 0, "complex candidate initializes")
    call check_condition(validate_complex_interchange_samples(reference, status) .and. status == 0, &
        "complex samples validate")
    copy = reference
    call check_condition(validate_complex_interchange_samples(copy, status) .and. &
        maxval(abs(copy%values - reference%values)) < 1.0e-14_dp, &
        "complex sample assignment deep-copies values")

    call compare_complex_interchange_samples( &
        reference, candidate, 1.0e-14_dp, absolute_error, relative_error, maximum_error, status)
    call check_condition(status == 0 .and. absolute_error > 0.0_dp .and. &
        relative_error > 0.0_dp .and. maximum_error > 0.0_dp, &
        "weighted complex sample error is reported")
    invalid = candidate
    invalid%coordinates(1, 1) = invalid%coordinates(1, 1) + 1.0e-3_dp
    call compare_complex_interchange_samples( &
        reference, invalid, 1.0e-6_dp, absolute_error, relative_error, maximum_error, status)
    call check_condition(status /= 0, "complex sample comparison rejects mismatched points")
    call check_summary("complex interchange samples")
end program test_complex_interchange_samples
