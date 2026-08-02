program test_interop_facade
    !! Smoke-test the canonical interoperability facade without the umbrella API.
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use check, only: check_condition, check_summary
    use fortfem_interop, only: &
        BOUNDARY_OPERATOR_BACKEND_BEM, BOUNDARY_OPERATOR_BACKEND_FEM, &
        BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL, &
        boundary_operator_contract_t, boundary_operator_parity_t, &
        compare_boundary_operator_parity, &
        compare_boundary_operator_parity_jvp, &
        compare_complex_interchange_samples, &
        complex_interchange_sample_set_t, &
        compare_interchange_samples_facade => compare_interchange_samples, &
        initialize_boundary_operator_contract, &
        initialize_boundary_operator_trace_metadata, &
        initialize_complex_interchange_samples, initialize_interchange_samples, &
        initialize_oracle_manifest, interchange_sample_set_t, &
        oracle_manifest_t, oracle_normalization_t, oracle_timing_t, &
        oracle_tolerance_t, validate_boundary_operator_contract, &
        validate_boundary_operator_parity, validate_complex_interchange_samples, &
        validate_interchange_samples, validate_oracle_manifest
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: sample_count = 3, backend_count = 2
    real(dp) :: coordinates(2, sample_count), values(1, sample_count)
    real(dp) :: weights(sample_count), absolute_error, relative_error
    real(dp) :: maximum_error
    complex(dp) :: complex_values(1, sample_count)
    complex(dp) :: reference(1, sample_count), candidates(1, sample_count, backend_count)
    complex(dp) :: reference_dot(1, sample_count)
    complex(dp) :: candidates_dot(1, sample_count, backend_count)
    real(dp) :: weights_dot(sample_count), reference_norm_dot
    real(dp) :: absolute_error_dot(backend_count), relative_error_dot(backend_count)
    integer :: backend_kinds(backend_count)
    type(interchange_sample_set_t) :: samples, candidate_samples
    type(complex_interchange_sample_set_t) :: complex_samples, complex_candidate
    type(oracle_manifest_t) :: manifest
    type(oracle_normalization_t) :: normalization
    type(oracle_tolerance_t) :: tolerances
    type(oracle_timing_t) :: timing
    type(boundary_operator_contract_t) :: contracts(backend_count)
    type(boundary_operator_parity_t) :: parity
    integer :: backend, status

    coordinates = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], &
        shape(coordinates))
    values = reshape([1.0_dp, 2.0_dp, 3.0_dp], shape(values))
    weights = [1.0_dp, 2.0_dp, 1.5_dp]
    call initialize_interchange_samples(samples, coordinates, values, weights, &
        "fortfem", "manufactured-unit-square", status)
    call check_condition(status == 0 .and. validate_interchange_samples(samples, status), &
        "facade initializes and validates real physical samples")
    call initialize_interchange_samples(candidate_samples, coordinates, &
        values + reshape([0.0_dp, 0.1_dp, 0.0_dp], shape(values)), weights, &
        "external-fem", "same-physical-samples", status)
    call compare_interchange_samples_facade( &
        samples, candidate_samples, 1.0e-14_dp, absolute_error, relative_error, &
        maximum_error, status)
    call check_condition(status == 0 .and. maximum_error > 0.0_dp .and. &
        absolute_error > 0.0_dp .and. relative_error > 0.0_dp, &
        "facade exposes weighted real comparison")

    complex_values = cmplx(values, reshape([0.2_dp, -0.1_dp, 0.3_dp], &
        shape(values)), dp)
    call initialize_complex_interchange_samples(complex_samples, coordinates, &
        complex_values, weights, "fortfem", "frequency-domain-oracle", status)
    call initialize_complex_interchange_samples(complex_candidate, coordinates, &
        complex_values + cmplx(0.0_dp, reshape([0.0_dp, 0.1_dp, 0.0_dp], &
        shape(values)), dp), weights, "external-bem", "same-physical-samples", status)
    call compare_complex_interchange_samples(complex_samples, complex_candidate, &
        1.0e-14_dp, absolute_error, relative_error, maximum_error, status)
    call check_condition(status == 0 .and. maximum_error > 0.0_dp .and. &
        validate_complex_interchange_samples(complex_samples, status), &
        "facade exposes weighted complex comparison and provenance")

    normalization%normalization_name = "SI"
    normalization%length_unit = "m"
    normalization%time_unit = "s"
    tolerances%coordinate = 1.0e-12_dp
    tolerances%absolute = 1.0e-10_dp
    tolerances%relative = 1.0e-9_dp
    tolerances%residual = 1.0e-8_dp
    timing%total_seconds = 0.25_dp
    timing%peak_memory_bytes = 4096_int64
    call initialize_oracle_manifest(manifest, "external-fem", "1.0", "rev", &
        "BSD-3-Clause", "unit-square", "case-1", "cartesian", "coords", &
        "values", 2, sample_count, normalization, tolerances, timing, "ci", &
        "one-thread", "fortfem-rev", sister_repository_uri="data@rev", &
        success=.true., notes="facade smoke", status=status)
    call check_condition(status == 0 .and. validate_oracle_manifest(manifest, status) .and. &
        trim(manifest%notes) == "facade smoke" .and. &
        trim(manifest%code_license) == "BSD-3-Clause", &
        "facade preserves oracle metadata and validation")

    reference = cmplx(values, 0.0_dp, dp)
    candidates(:, :, 1) = reference
    candidates(:, :, 2) = reference
    candidates(1, 1, 1) = candidates(1, 1, 1) + cmplx(0.005_dp, 0.0_dp, dp)
    candidates(1, 2, 2) = candidates(1, 2, 2) + cmplx(0.1_dp, -0.2_dp, dp)
    reference_dot = cmplx(0.1_dp, -0.1_dp, dp)
    candidates_dot = cmplx(0.0_dp, 0.0_dp, dp)
    candidates_dot(1, :, 2) = cmplx(0.2_dp, 0.1_dp, dp)
    weights_dot = [0.1_dp, -0.2_dp, 0.3_dp]
    backend_kinds = [BOUNDARY_OPERATOR_BACKEND_FEM, BOUNDARY_OPERATOR_BACKEND_BEM]
    do backend = 1, backend_count
        call initialize_boundary_operator_contract(contracts(backend), &
            backend_kinds(backend), &
            "helmholtz", "H1-trace", sample_count, sample_count, .true., .true., &
            .true., .true., .true., .true., "unit", "manufactured", &
            "independent-oracle", "square-boundary", status)
        call check_condition(status == 0 .and. validate_boundary_operator_contract( &
            contracts(backend), status), "facade initializes boundary contract")
        call initialize_boundary_operator_trace_metadata(contracts(backend), &
            BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL, "weighted-tangential", status)
        call check_condition(status == 0, "facade exposes trace metadata contract")
    end do
    call compare_boundary_operator_parity(reference, candidates, weights, contracts, &
        0.01_dp, 0.01_dp, parity, status)
    call check_condition(status == 0 .and. validate_boundary_operator_parity(parity, status) &
        .and. parity%backend_count == backend_count .and. parity%within_tolerance(1) &
        .and. .not. parity%within_tolerance(2), &
        "facade exposes weighted boundary parity with provenance")
    call compare_boundary_operator_parity_jvp(reference, candidates, weights, contracts, &
        0.01_dp, 0.01_dp, reference_dot, candidates_dot, weights_dot, &
        reference_norm_dot, absolute_error_dot, relative_error_dot, status)
    call check_condition(status == 0 .and. abs(reference_norm_dot) < 1.0_dp .and. &
        absolute_error_dot(2) /= 0.0_dp, &
        "facade exposes boundary parity JVP")
    call check_summary("interoperability facade")
end program test_interop_facade
