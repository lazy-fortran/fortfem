program test_boundary_operator_parity
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_interop, only: &
        BOUNDARY_OPERATOR_BACKEND_BEM, &
        BOUNDARY_OPERATOR_BACKEND_DTN, &
        BOUNDARY_OPERATOR_BACKEND_FEM, &
        BOUNDARY_OPERATOR_BACKEND_PML, &
        boundary_operator_contract_t, &
        boundary_operator_parity_t, &
        compare_boundary_operator_parity, &
        compare_boundary_operator_parity_jvp, &
        compare_boundary_operator_parity_vjp, &
        validate_boundary_operator_parity, &
        initialize_boundary_operator_contract
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: backend_count = 4
    integer, parameter :: sample_count = 3
    type(boundary_operator_contract_t) :: contracts(backend_count)
    type(boundary_operator_parity_t) :: report, report_plus, report_minus, invalid
    complex(dp) :: reference(1, sample_count), candidates(1, sample_count, backend_count)
    complex(dp) :: reference_dot(1, sample_count)
    complex(dp) :: candidates_dot(1, sample_count, backend_count)
    complex(dp) :: reference_plus(1, sample_count), reference_minus(1, sample_count)
    complex(dp) :: candidates_plus(1, sample_count, backend_count)
    complex(dp) :: candidates_minus(1, sample_count, backend_count)
    complex(dp) :: reference_bar(1, sample_count)
    complex(dp) :: candidates_bar(1, sample_count, backend_count)
    real(dp) :: weights(sample_count)
    real(dp) :: weights_dot(sample_count), weights_plus(sample_count)
    real(dp) :: weights_minus(sample_count), reference_norm_dot
    real(dp) :: absolute_error_dot(backend_count), relative_error_dot(backend_count)
    real(dp) :: absolute_error_bar(backend_count), relative_error_bar(backend_count)
    real(dp) :: reference_norm_bar, weights_bar(sample_count)
    real(dp) :: reference_norm_plus, reference_norm_minus
    real(dp) :: absolute_error_plus(backend_count), absolute_error_minus(backend_count)
    real(dp) :: relative_error_plus(backend_count), relative_error_minus(backend_count)
    real(dp) :: lhs, rhs
    real(dp), parameter :: epsilon_fd = 1.0e-7_dp
    integer :: backend_kind(backend_count), status, backend
    logical :: all_passed

    all_passed = .true.
    backend_kind = [BOUNDARY_OPERATOR_BACKEND_FEM, BOUNDARY_OPERATOR_BACKEND_BEM, &
        BOUNDARY_OPERATOR_BACKEND_DTN, BOUNDARY_OPERATOR_BACKEND_PML]
    reference = reshape([cmplx(1.0_dp, 0.0_dp, dp), cmplx(2.0_dp, 0.0_dp, dp), &
        cmplx(-1.0_dp, 0.0_dp, dp)], shape(reference))
    candidates(:, :, 1) = reference
    candidates(1, 3, 1) = reference(1, 3) + cmplx(0.02_dp, -0.01_dp, dp)
    candidates(:, :, 2) = reference
    candidates(1, 1, 2) = reference(1, 1) + cmplx(0.1_dp, 0.0_dp, dp)
    candidates(:, :, 3) = reference
    candidates(1, 2, 3) = reference(1, 2) + cmplx(0.2_dp, 0.0_dp, dp)
    candidates(:, :, 4) = reference
    candidates(1, 2, 4) = reference(1, 2) + cmplx(0.3_dp, 0.0_dp, dp)
    weights = [1.0_dp, 2.0_dp, 1.0_dp]
    reference_dot = reshape([cmplx(0.1_dp, -0.2_dp, dp), &
        cmplx(-0.3_dp, 0.4_dp, dp), cmplx(0.2_dp, 0.1_dp, dp)], &
        shape(reference_dot))
    candidates_dot = 0.0_dp
    candidates_dot(1, :, 1) = [cmplx(0.2_dp, 0.1_dp, dp), &
        cmplx(-0.1_dp, 0.3_dp, dp), cmplx(0.3_dp, -0.2_dp, dp)]
    candidates_dot(1, :, 2) = [cmplx(-0.2_dp, 0.4_dp, dp), &
        cmplx(0.1_dp, -0.3_dp, dp), cmplx(0.2_dp, 0.2_dp, dp)]
    candidates_dot(1, :, 3) = [cmplx(0.3_dp, -0.1_dp, dp), &
        cmplx(-0.2_dp, 0.2_dp, dp), cmplx(0.1_dp, 0.4_dp, dp)]
    candidates_dot(1, :, 4) = [cmplx(-0.1_dp, 0.2_dp, dp), &
        cmplx(0.2_dp, -0.4_dp, dp), cmplx(-0.3_dp, 0.1_dp, dp)]
    weights_dot = [0.1_dp, -0.2_dp, 0.3_dp]

    do backend = 1, backend_count
        call initialize_boundary_operator_contract(contracts(backend), backend_kind(backend), &
            "helmholtz", "H1-trace", sample_count, sample_count, .true., .true., .true., &
            .true., .true., .true., "unit", "manufactured", "parity oracle", &
            "circle-fixed-1", status)
        call record_condition(status == 0, "parity metadata initializes")
    end do

    call compare_boundary_operator_parity(reference, candidates, weights, contracts, &
        0.05_dp, 0.10_dp, report, status)
    call record_condition(status == 0, "boundary operator parity evaluates")
    call record_condition(validate_boundary_operator_parity(report, status), &
        "boundary operator parity validates")
    call record_condition(report%schema_version == "fortfem-boundary-parity-1" .and. &
        trim(report%topology_id) == "circle-fixed-1" .and. &
        trim(report%provenance) == "parity oracle", &
        "canonical compare name preserves parity metadata")
    call record_condition(report%backend_count == backend_count .and. &
        report%sample_count == sample_count .and. report%component_count == 1 .and. &
        abs(report%reference_norm - sqrt(10.0_dp)) < 1.0e-12_dp, &
        "parity report stores the common physical norm")
    call record_condition(report%within_tolerance(1) .and. report%within_tolerance(2) .and. &
        report%within_tolerance(3) .and. .not. report%within_tolerance(4), &
        "parity report distinguishes methods by weighted error")
    call record_condition(abs(report%absolute_error(2) - 0.1_dp) < 1.0e-12_dp .and. &
        abs(report%absolute_error(3) - sqrt(0.08_dp)) < 1.0e-12_dp .and. &
        abs(report%absolute_error(4) - sqrt(0.18_dp)) < 1.0e-12_dp, &
        "parity errors match an independent weighted oracle")

    call compare_boundary_operator_parity_jvp( &
        reference, candidates, weights, contracts, 0.05_dp, 0.10_dp, &
        reference_dot, candidates_dot, weights_dot, reference_norm_dot, &
        absolute_error_dot, relative_error_dot, status)
    reference_norm_plus = 0.0_dp
    reference_norm_minus = 0.0_dp
    reference_plus = reference + epsilon_fd*reference_dot
    reference_minus = reference - epsilon_fd*reference_dot
    candidates_plus = candidates + epsilon_fd*candidates_dot
    candidates_minus = candidates - epsilon_fd*candidates_dot
    weights_plus = weights + epsilon_fd*weights_dot
    weights_minus = weights - epsilon_fd*weights_dot
    call compare_boundary_operator_parity( &
        reference_plus, candidates_plus, weights_plus, contracts, 0.05_dp, 0.10_dp, &
        report_plus, status)
    call compare_boundary_operator_parity( &
        reference_minus, candidates_minus, weights_minus, contracts, 0.05_dp, 0.10_dp, &
        report_minus, status)
    reference_norm_plus = report_plus%reference_norm
    reference_norm_minus = report_minus%reference_norm
    absolute_error_plus = report_plus%absolute_error
    absolute_error_minus = report_minus%absolute_error
    relative_error_plus = report_plus%relative_error
    relative_error_minus = report_minus%relative_error
    call record_condition(status == 0 .and. &
        abs(reference_norm_dot - (reference_norm_plus - reference_norm_minus)/ &
        (2.0_dp*epsilon_fd)) < 1.0e-7_dp .and. &
        maxval(abs(absolute_error_dot - (absolute_error_plus - absolute_error_minus)/ &
        (2.0_dp*epsilon_fd))) < 1.0e-7_dp .and. &
        maxval(abs(relative_error_dot - (relative_error_plus - relative_error_minus)/ &
        (2.0_dp*epsilon_fd))) < 1.0e-7_dp, &
        "canonical compare JVP matches central re-evaluation for all backends")

    reference_norm_bar = -0.4_dp
    absolute_error_bar = [0.2_dp, -0.3_dp, 0.5_dp, -0.1_dp]
    relative_error_bar = [-0.6_dp, 0.7_dp, -0.2_dp, 0.4_dp]
    call compare_boundary_operator_parity_vjp( &
        reference, candidates, weights, contracts, 0.05_dp, 0.10_dp, &
        reference_norm_bar, absolute_error_bar, relative_error_bar, &
        reference_bar, candidates_bar, weights_bar, status)
    lhs = reference_norm_bar*reference_norm_dot + &
        dot_product(absolute_error_bar, absolute_error_dot) + &
        dot_product(relative_error_bar, relative_error_dot)
    rhs = real(sum(conjg(reference_bar)*reference_dot), dp) + &
        real(sum(conjg(candidates_bar)*candidates_dot), dp) + &
        dot_product(weights_bar, weights_dot)
    call record_condition(status == 0 .and. abs(lhs - rhs) < 1.0e-7_dp, &
        "canonical compare VJP satisfies the complex real-part adjoint oracle")

    contracts(2)%topology_id = "different-topology"
    call compare_boundary_operator_parity(reference, candidates, weights, contracts, &
        0.05_dp, 0.10_dp, invalid, status)
    call record_condition(status /= 0, "parity rejects mixed physical topologies")

    call check_summary("boundary operator parity")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_boundary_operator_parity
