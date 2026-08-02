program test_boundary_operator_parity
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        BOUNDARY_OPERATOR_BACKEND_BEM, &
        BOUNDARY_OPERATOR_BACKEND_DTN, &
        BOUNDARY_OPERATOR_BACKEND_FEM, &
        BOUNDARY_OPERATOR_BACKEND_PML, &
        boundary_operator_contract_t, &
        boundary_operator_parity_t, &
        evaluate_boundary_operator_parity, &
        validate_boundary_operator_parity, &
        initialize_boundary_operator_contract
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: backend_count = 4
    integer, parameter :: sample_count = 3
    type(boundary_operator_contract_t) :: contracts(backend_count)
    type(boundary_operator_parity_t) :: report, invalid
    complex(dp) :: reference(1, sample_count), candidates(1, sample_count, backend_count)
    real(dp) :: weights(sample_count)
    integer :: backend_kind(backend_count), status, backend
    logical :: all_passed

    all_passed = .true.
    backend_kind = [BOUNDARY_OPERATOR_BACKEND_FEM, BOUNDARY_OPERATOR_BACKEND_BEM, &
                    BOUNDARY_OPERATOR_BACKEND_DTN, BOUNDARY_OPERATOR_BACKEND_PML]
    reference = reshape([cmplx(1.0_dp, 0.0_dp, dp), cmplx(2.0_dp, 0.0_dp, dp), &
                         cmplx(-1.0_dp, 0.0_dp, dp)], shape(reference))
    candidates(:, :, 1) = reference
    candidates(:, :, 2) = reference
    candidates(1, 1, 2) = reference(1, 1) + cmplx(0.1_dp, 0.0_dp, dp)
    candidates(:, :, 3) = reference
    candidates(1, 2, 3) = reference(1, 2) + cmplx(0.2_dp, 0.0_dp, dp)
    candidates(:, :, 4) = reference
    candidates(1, 2, 4) = reference(1, 2) + cmplx(0.3_dp, 0.0_dp, dp)
    weights = [1.0_dp, 2.0_dp, 1.0_dp]

    do backend = 1, backend_count
        call initialize_boundary_operator_contract(contracts(backend), backend_kind(backend), &
            "helmholtz", "H1-trace", sample_count, sample_count, .true., .true., .true., &
            .true., .true., .true., "unit", "manufactured", "parity oracle", &
            "circle-fixed-1", status)
        call record_condition(status == 0, "parity metadata initializes")
    end do

    call evaluate_boundary_operator_parity(reference, candidates, weights, contracts, &
        0.05_dp, 0.10_dp, report, status)
    call record_condition(status == 0, "boundary operator parity evaluates")
    call record_condition(validate_boundary_operator_parity(report, status), &
        "boundary operator parity validates")
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

    contracts(2)%topology_id = "different-topology"
    call evaluate_boundary_operator_parity(reference, candidates, weights, contracts, &
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
