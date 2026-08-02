program test_maxwell_physical_boundary_parity
    !! Independent physical-sample parity fixture for vector open boundaries.
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
        initialize_boundary_operator_contract, &
        validate_boundary_operator_parity
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: backend_count = 4
    integer, parameter :: sample_count = 6
    integer, parameter :: component_count = 3
    type(boundary_operator_contract_t) :: contracts(backend_count)
    type(boundary_operator_parity_t) :: report, invalid
    complex(dp) :: reference(component_count, sample_count)
    complex(dp) :: candidates(component_count, sample_count, backend_count)
    real(dp) :: weights(sample_count), expected_error(backend_count)
    integer :: backend_kind(backend_count), status, backend, sample
    logical :: all_passed
    real(dp), parameter :: pi = acos(-1.0_dp)

    all_passed = .true.
    backend_kind = [BOUNDARY_OPERATOR_BACKEND_FEM, BOUNDARY_OPERATOR_BACKEND_BEM, &
                    BOUNDARY_OPERATOR_BACKEND_DTN, BOUNDARY_OPERATOR_BACKEND_PML]
    weights = [1.0_dp, 1.5_dp, 2.0_dp, 1.25_dp, 0.75_dp, 1.1_dp]

    ! A divergence-free manufactured trace: the first two components are a
    ! rotating tangential vector and the third is a phase-shifted normal mode.
    do sample = 1, sample_count
        reference(:, sample) = [ &
            cmplx(cos(2.0_dp*pi*real(sample - 1, dp)/sample_count), &
                sin(2.0_dp*pi*real(sample - 1, dp)/sample_count), dp), &
            cmplx(-sin(2.0_dp*pi*real(sample - 1, dp)/sample_count), &
                cos(2.0_dp*pi*real(sample - 1, dp)/sample_count), dp), &
            cmplx(0.25_dp, -0.1_dp, dp)]
    end do
    candidates = 0.0_dp
    do backend = 1, backend_count
        candidates(:, :, backend) = reference
    end do
    candidates(1, 2, 2) = candidates(1, 2, 2) + cmplx(0.08_dp, -0.03_dp, dp)
    candidates(2, 4, 3) = candidates(2, 4, 3) + cmplx(-0.05_dp, 0.04_dp, dp)
    candidates(3, 5, 4) = candidates(3, 5, 4) + cmplx(0.8_dp, 0.6_dp, dp)

    do backend = 1, backend_count
        call initialize_boundary_operator_contract(contracts(backend), &
            backend_kind(backend), "curl-curl", "Hcurl-tangential", &
            sample_count, sample_count, .true., .true., .true., .true., .true., &
            .true., "unit", "manufactured-maxwell", &
            "physical torus/cylinder parity", "torus-cylinder-fixed-1", status)
        call record_condition(status == 0, "physical parity metadata initializes")
    end do

    call evaluate_boundary_operator_parity(reference, candidates, weights, contracts, &
        0.06_dp, 0.12_dp, report, status)
    call record_condition(status == 0, "physical Maxwell parity evaluates")
    call record_condition(validate_boundary_operator_parity(report, status), &
        "physical Maxwell parity validates")
    expected_error = [0.0_dp, sqrt(0.08_dp**2 + 0.03_dp**2)*sqrt(weights(2)), &
                      sqrt(0.05_dp**2 + 0.04_dp**2)*sqrt(weights(4)), &
                      sqrt(0.8_dp**2 + 0.6_dp**2)*sqrt(weights(5))]
    do backend = 1, backend_count
        call record_condition(abs(report%absolute_error(backend) - &
            expected_error(backend)) < 1.0e-12_dp, &
            "physical weighted vector error matches independent oracle")
    end do
    call record_condition(report%within_tolerance(1) .and. &
        report%within_tolerance(2) .and. report%within_tolerance(3) .and. &
        .not. report%within_tolerance(4), &
        "physical parity distinguishes FEM/BEM/DtN/PML tolerances")

    contracts(3)%space = "H1-scalar"
    call evaluate_boundary_operator_parity(reference, candidates, weights, contracts, &
        0.06_dp, 0.12_dp, invalid, status)
    call record_condition(status /= 0, "physical parity rejects mixed trace spaces")
    call check_summary("physical Maxwell boundary parity")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_physical_boundary_parity
