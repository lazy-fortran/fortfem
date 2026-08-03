program test_surface_trace_parity_ledger
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_boundary_operator_contract, only: &
        BOUNDARY_OPERATOR_BACKEND_BEM, BOUNDARY_OPERATOR_BACKEND_DTN, &
        boundary_operator_contract_t, initialize_boundary_operator_contract
    use fortfem_surface_trace_parity_ledger, only: &
        evaluate_surface_trace_parity_ledger, &
        evaluate_surface_trace_parity_ledger_jvp, &
        evaluate_surface_trace_parity_ledger_vjp
    implicit none

    integer, parameter :: dp = real64, component_count = 2, sample_count = 3
    real(dp), parameter :: epsilon = 1.0e-7_dp
    type(boundary_operator_contract_t) :: contracts(2)
    complex(dp) :: reference_trace(component_count, sample_count)
    complex(dp) :: reference_response(component_count, sample_count)
    complex(dp) :: candidate_trace(component_count, sample_count)
    complex(dp) :: candidate_response(component_count, sample_count)
    complex(dp) :: reference_trace_dot(component_count, sample_count)
    complex(dp) :: reference_response_dot(component_count, sample_count)
    complex(dp) :: candidate_trace_dot(component_count, sample_count)
    complex(dp) :: candidate_response_dot(component_count, sample_count)
    complex(dp) :: reference_trace_plus(component_count, sample_count)
    complex(dp) :: reference_trace_minus(component_count, sample_count)
    complex(dp) :: reference_response_plus(component_count, sample_count)
    complex(dp) :: reference_response_minus(component_count, sample_count)
    complex(dp) :: candidate_trace_plus(component_count, sample_count)
    complex(dp) :: candidate_trace_minus(component_count, sample_count)
    complex(dp) :: candidate_response_plus(component_count, sample_count)
    complex(dp) :: candidate_response_minus(component_count, sample_count)
    complex(dp) :: reference_trace_bar(component_count, sample_count)
    complex(dp) :: reference_response_bar(component_count, sample_count)
    complex(dp) :: candidate_trace_bar(component_count, sample_count)
    complex(dp) :: candidate_response_bar(component_count, sample_count)
    real(dp) :: weights(sample_count), weights_dot(sample_count)
    real(dp) :: weights_plus(sample_count), weights_minus(sample_count)
    real(dp) :: weights_bar(sample_count)
    real(dp) :: reference_norm, absolute_error, relative_error, reciprocity_defect
    real(dp) :: reference_norm_dot, absolute_error_dot, relative_error_dot
    real(dp) :: reciprocity_defect_dot
    real(dp) :: reference_norm_plus, reference_norm_minus
    real(dp) :: absolute_error_plus, absolute_error_minus
    real(dp) :: relative_error_plus, relative_error_minus
    real(dp) :: reciprocity_defect_plus, reciprocity_defect_minus
    real(dp) :: reference_norm_bar, absolute_error_bar, relative_error_bar
    real(dp) :: reciprocity_defect_bar, lhs, rhs
    real(dp) :: reference_squared, error_squared, pairing_one, pairing_two
    real(dp) :: expected_reciprocity, expected_scale
    integer :: status
    logical :: within_tolerance

    reference_trace = reshape([ &
        cmplx(1.0_dp, 0.2_dp, dp), cmplx(-0.4_dp, 0.3_dp, dp), &
        cmplx(0.7_dp, -0.1_dp, dp), cmplx(0.2_dp, 0.8_dp, dp), &
        cmplx(-0.3_dp, 0.5_dp, dp), cmplx(0.6_dp, -0.2_dp, dp)], &
        shape(reference_trace))
    reference_response = reshape([ &
        cmplx(0.4_dp, -0.1_dp, dp), cmplx(0.1_dp, 0.7_dp, dp), &
        cmplx(-0.6_dp, 0.2_dp, dp), cmplx(0.3_dp, -0.5_dp, dp), &
        cmplx(0.8_dp, 0.1_dp, dp), cmplx(-0.2_dp, 0.4_dp, dp)], &
        shape(reference_response))
    candidate_trace = reference_trace + reshape([ &
        cmplx(0.1_dp, -0.03_dp, dp), cmplx(-0.04_dp, 0.08_dp, dp), &
        cmplx(0.06_dp, 0.02_dp, dp), cmplx(-0.02_dp, 0.05_dp, dp), &
        cmplx(0.03_dp, -0.07_dp, dp), cmplx(-0.05_dp, 0.01_dp, dp)], &
        shape(candidate_trace))
    candidate_response = reference_response + reshape([ &
        cmplx(-0.05_dp, 0.06_dp, dp), cmplx(0.02_dp, -0.03_dp, dp), &
        cmplx(0.08_dp, 0.01_dp, dp), cmplx(-0.04_dp, 0.02_dp, dp), &
        cmplx(0.03_dp, 0.05_dp, dp), cmplx(-0.06_dp, -0.02_dp, dp)], &
        shape(candidate_response))
    reference_trace_dot = reshape([ &
        cmplx(0.02_dp, -0.03_dp, dp), cmplx(-0.01_dp, 0.04_dp, dp), &
        cmplx(0.03_dp, 0.02_dp, dp), cmplx(-0.04_dp, 0.01_dp, dp), &
        cmplx(0.01_dp, 0.03_dp, dp), cmplx(-0.02_dp, -0.01_dp, dp)], &
        shape(reference_trace_dot))
    reference_response_dot = reshape([ &
        cmplx(-0.03_dp, 0.02_dp, dp), cmplx(0.04_dp, -0.01_dp, dp), &
        cmplx(0.01_dp, 0.03_dp, dp), cmplx(-0.02_dp, 0.02_dp, dp), &
        cmplx(0.03_dp, -0.04_dp, dp), cmplx(0.02_dp, 0.01_dp, dp)], &
        shape(reference_response_dot))
    candidate_trace_dot = reshape([ &
        cmplx(-0.01_dp, 0.04_dp, dp), cmplx(0.02_dp, 0.01_dp, dp), &
        cmplx(-0.03_dp, 0.02_dp, dp), cmplx(0.01_dp, -0.05_dp, dp), &
        cmplx(0.04_dp, 0.01_dp, dp), cmplx(-0.02_dp, 0.03_dp, dp)], &
        shape(candidate_trace_dot))
    candidate_response_dot = reshape([ &
        cmplx(0.02_dp, 0.01_dp, dp), cmplx(-0.03_dp, 0.02_dp, dp), &
        cmplx(0.01_dp, -0.04_dp, dp), cmplx(0.03_dp, 0.01_dp, dp), &
        cmplx(-0.02_dp, 0.05_dp, dp), cmplx(0.04_dp, -0.01_dp, dp)], &
        shape(candidate_response_dot))
    weights = [1.0_dp, 1.5_dp, 0.75_dp]
    weights_dot = [0.04_dp, -0.03_dp, 0.02_dp]

    call initialize_boundary_operator_contract(contracts(1), &
        BOUNDARY_OPERATOR_BACKEND_BEM, "maxwell", "trace", sample_count, &
        sample_count, .true., .true., .true., .true., .true., .true., "unit", &
        "surface-work", "reference manufactured trace", "surface-fixed-1", status)
    call check_condition(status == 0, "reference contract initializes")
    call initialize_boundary_operator_contract(contracts(2), &
        BOUNDARY_OPERATOR_BACKEND_DTN, "maxwell", "trace", sample_count, &
        sample_count, .true., .true., .true., .true., .true., .true., "unit", &
        "surface-work", "candidate manufactured trace", "surface-fixed-1", status)
    call check_condition(status == 0, "candidate contract initializes")

    call evaluate_surface_trace_parity_ledger( &
        reference_trace, reference_response, candidate_trace, candidate_response, &
        weights, contracts, 0.05_dp, 0.10_dp, reference_norm, absolute_error, &
        relative_error, reciprocity_defect, within_tolerance, status)
    reference_squared = sum(weights*sum(abs(reference_trace)**2, dim=1))
    error_squared = sum(weights*sum(abs(candidate_trace-reference_trace)**2, dim=1))
    pairing_one = real(sum(spread(weights, 1, component_count)* &
        reference_trace*candidate_response), dp)
    pairing_two = real(sum(spread(weights, 1, component_count)* &
        candidate_trace*reference_response), dp)
    expected_scale = max(1.0_dp, abs(pairing_one), abs(pairing_two))
    expected_reciprocity = abs(pairing_one-pairing_two)/expected_scale
    call check_condition(status == 0 .and. abs(reference_norm-sqrt(reference_squared)) < 1.0e-12_dp .and. &
        abs(absolute_error-sqrt(error_squared)) < 1.0e-12_dp .and. &
        abs(relative_error-sqrt(error_squared)/sqrt(reference_squared)) < 1.0e-12_dp .and. &
        abs(reciprocity_defect-expected_reciprocity) < 1.0e-12_dp .and. within_tolerance, &
        "surface parity ledger matches independent weighted work/error oracle")

    call evaluate_surface_trace_parity_ledger_jvp( &
        reference_trace, reference_response, candidate_trace, candidate_response, &
        weights, contracts, 0.05_dp, 0.10_dp, reference_trace_dot, &
        reference_response_dot, candidate_trace_dot, candidate_response_dot, weights_dot, &
        reference_norm_dot, absolute_error_dot, relative_error_dot, &
        reciprocity_defect_dot, status)
    reference_trace_plus = reference_trace + epsilon*reference_trace_dot
    reference_trace_minus = reference_trace - epsilon*reference_trace_dot
    reference_response_plus = reference_response + epsilon*reference_response_dot
    reference_response_minus = reference_response - epsilon*reference_response_dot
    candidate_trace_plus = candidate_trace + epsilon*candidate_trace_dot
    candidate_trace_minus = candidate_trace - epsilon*candidate_trace_dot
    candidate_response_plus = candidate_response + epsilon*candidate_response_dot
    candidate_response_minus = candidate_response - epsilon*candidate_response_dot
    weights_plus = weights + epsilon*weights_dot
    weights_minus = weights - epsilon*weights_dot
    call evaluate_surface_trace_parity_ledger( &
        reference_trace_plus, reference_response_plus, candidate_trace_plus, &
        candidate_response_plus, weights_plus, contracts, 0.05_dp, 0.10_dp, &
        reference_norm_plus, absolute_error_plus, relative_error_plus, &
        reciprocity_defect_plus, within_tolerance, status)
    call evaluate_surface_trace_parity_ledger( &
        reference_trace_minus, reference_response_minus, candidate_trace_minus, &
        candidate_response_minus, weights_minus, contracts, 0.05_dp, 0.10_dp, &
        reference_norm_minus, absolute_error_minus, relative_error_minus, &
        reciprocity_defect_minus, within_tolerance, status)
    call check_condition(status == 0 .and. &
        abs(reference_norm_dot-(reference_norm_plus-reference_norm_minus)/(2.0_dp*epsilon)) < 1.0e-7_dp .and. &
        abs(absolute_error_dot-(absolute_error_plus-absolute_error_minus)/(2.0_dp*epsilon)) < 1.0e-7_dp .and. &
        abs(relative_error_dot-(relative_error_plus-relative_error_minus)/(2.0_dp*epsilon)) < 1.0e-7_dp .and. &
        abs(reciprocity_defect_dot-(reciprocity_defect_plus-reciprocity_defect_minus)/(2.0_dp*epsilon)) < 1.0e-7_dp, &
        "surface parity ledger JVP matches central re-evaluation")
    write(*,*) "DEBUG JVP", reference_norm_dot, (reference_norm_plus-reference_norm_minus)/(2.0_dp*epsilon), &
        absolute_error_dot, (absolute_error_plus-absolute_error_minus)/(2.0_dp*epsilon), &
        relative_error_dot, (relative_error_plus-relative_error_minus)/(2.0_dp*epsilon), &
        reciprocity_defect_dot, (reciprocity_defect_plus-reciprocity_defect_minus)/(2.0_dp*epsilon)

    reference_norm_bar = -0.3_dp
    absolute_error_bar = 0.4_dp
    relative_error_bar = -0.2_dp
    reciprocity_defect_bar = 0.5_dp
    call evaluate_surface_trace_parity_ledger_vjp( &
        reference_trace, reference_response, candidate_trace, candidate_response, &
        weights, contracts, 0.05_dp, 0.10_dp, reference_norm_bar, absolute_error_bar, &
        relative_error_bar, reciprocity_defect_bar, reference_trace_bar, &
        reference_response_bar, candidate_trace_bar, candidate_response_bar, &
        weights_bar, status)
    lhs = reference_norm_bar*reference_norm_dot + absolute_error_bar*absolute_error_dot + &
        relative_error_bar*relative_error_dot + reciprocity_defect_bar*reciprocity_defect_dot
    rhs = real(sum(conjg(reference_trace_bar)*reference_trace_dot) + &
        sum(conjg(reference_response_bar)*reference_response_dot) + &
        sum(conjg(candidate_trace_bar)*candidate_trace_dot) + &
        sum(conjg(candidate_response_bar)*candidate_response_dot), dp) + &
        sum(weights_bar*weights_dot)
    call check_condition(status == 0 .and. abs(lhs-rhs) < 1.0e-7_dp, &
        "surface parity ledger VJP satisfies the real complex adjoint oracle")
    write(*,*) "DEBUG VJP", lhs, rhs, lhs-rhs

    contracts(2)%topology_id = "surface-changed"
    call evaluate_surface_trace_parity_ledger( &
        reference_trace, reference_response, candidate_trace, candidate_response, &
        weights, contracts, 0.05_dp, 0.10_dp, reference_norm, absolute_error, &
        relative_error, reciprocity_defect, within_tolerance, status)
    call check_condition(status /= 0, "surface parity ledger rejects changed topology")
    contracts(2)%topology_id = "surface-fixed-1"
    weights(2) = 0.0_dp
    call evaluate_surface_trace_parity_ledger( &
        reference_trace, reference_response, candidate_trace, candidate_response, &
        weights, contracts, 0.05_dp, 0.10_dp, reference_norm, absolute_error, &
        relative_error, reciprocity_defect, within_tolerance, status)
    call check_condition(status /= 0, "surface parity ledger rejects nonpositive work weights")

    call check_summary("surface trace parity ledger")

end program test_surface_trace_parity_ledger
