program test_larger_domain_parity
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        BOUNDARY_OPERATOR_BACKEND_PML, &
        boundary_operator_contract_t, &
        initialize_boundary_operator_contract, &
        larger_domain_parity_t, &
        evaluate_larger_domain_parity, &
        validate_larger_domain_parity
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: component_count = 2, sample_count = 3
    complex(dp) :: inner_field(component_count, sample_count)
    complex(dp) :: outer_field(component_count, sample_count)
    real(dp) :: weights(sample_count)
    type(boundary_operator_contract_t) :: inner_contract, outer_contract
    type(larger_domain_parity_t) :: report, invalid
    real(dp) :: expected_absolute, expected_norm, expected_relative
    integer :: status
    logical :: all_passed

    all_passed = .true.
    inner_field(1, :) = [cmplx(1.0_dp, 0.0_dp, dp), &
                          cmplx(2.0_dp, 0.0_dp, dp), &
                          cmplx(-1.0_dp, 0.0_dp, dp)]
    inner_field(2, :) = [cmplx(0.0_dp, 1.0_dp, dp), &
                          cmplx(1.0_dp, 0.0_dp, dp), &
                          cmplx(2.0_dp, 0.0_dp, dp)]
    outer_field = inner_field
    outer_field(1, 1) = cmplx(1.1_dp, 0.2_dp, dp)
    outer_field(2, 2) = cmplx(1.1_dp, 0.0_dp, dp)
    outer_field(2, 3) = cmplx(2.0_dp, -0.1_dp, dp)
    weights = [1.0_dp, 2.0_dp, 1.0_dp]

    call initialize_boundary_operator_contract( &
        inner_contract, BOUNDARY_OPERATOR_BACKEND_PML, "helmholtz", "Hcurl", &
        sample_count, sample_count, .true., .true., .true., .true., .true., .true., &
        "V/m", "unit", "larger-domain independent oracle", "common-interior-1", status)
    call record_condition(status == 0, "inner boundary metadata initializes")
    call initialize_boundary_operator_contract( &
        outer_contract, BOUNDARY_OPERATOR_BACKEND_PML, "helmholtz", "Hcurl", &
        sample_count, sample_count, .true., .true., .true., .true., .true., .true., &
        "V/m", "unit", "larger-domain independent oracle", "common-interior-1", status)
    call record_condition(status == 0, "outer boundary metadata initializes")

    call evaluate_larger_domain_parity( &
        inner_field, outer_field, weights, inner_contract, outer_contract, &
        0.5_dp, 1.5_dp, 0.30_dp, 0.08_dp, report, status)
    call record_condition(status == 0, "larger-domain parity evaluates")
    call record_condition(validate_larger_domain_parity(report, status), &
        "larger-domain parity validates")

    expected_absolute = sqrt(0.08_dp)
    expected_norm = sqrt(17.68_dp)
    expected_relative = expected_absolute/expected_norm
    call record_condition(abs(report%absolute_difference - expected_absolute) < 1.0e-12_dp .and. &
        abs(report%comparison_norm - expected_norm) < 1.0e-12_dp .and. &
        abs(report%relative_difference - expected_relative) < 1.0e-12_dp, &
        "difference and symmetric relative oracle are correct")
    call record_condition(report%farther_boundary .and. &
        abs(report%distance_increase - 1.0_dp) < 1.0e-12_dp .and. &
        abs(report%distance_ratio - 3.0_dp) < 1.0e-12_dp .and. &
        abs(report%relative_difference_per_distance - expected_relative) < 1.0e-12_dp .and. &
        report%within_tolerance, "distance convergence metadata is recorded")

    outer_contract%equation = "laplace"
    call evaluate_larger_domain_parity( &
        inner_field, outer_field, weights, inner_contract, outer_contract, &
        0.5_dp, 1.5_dp, 0.30_dp, 0.08_dp, invalid, status)
    call record_condition(status /= 0, "mixed equation spaces are rejected")

    call check_summary("larger-domain parity")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_larger_domain_parity
