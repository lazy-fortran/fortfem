program test_larger_domain_parity
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        BOUNDARY_OPERATOR_BACKEND_PML, &
        boundary_operator_contract_t, &
        initialize_boundary_operator_contract, &
        larger_domain_parity_t, &
        evaluate_larger_domain_parity, &
        evaluate_larger_domain_parity_jvp, &
        validate_larger_domain_parity
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: component_count = 2, sample_count = 3
    complex(dp) :: inner_field(component_count, sample_count)
    complex(dp) :: outer_field(component_count, sample_count)
    complex(dp) :: inner_field_dot(component_count, sample_count)
    complex(dp) :: outer_field_dot(component_count, sample_count)
    complex(dp) :: inner_field_plus(component_count, sample_count)
    complex(dp) :: outer_field_plus(component_count, sample_count)
    complex(dp) :: inner_field_minus(component_count, sample_count)
    complex(dp) :: outer_field_minus(component_count, sample_count)
    real(dp) :: weights(sample_count)
    real(dp) :: weights_dot(sample_count), weights_plus(sample_count)
    real(dp) :: weights_minus(sample_count)
    type(boundary_operator_contract_t) :: inner_contract, outer_contract
    type(larger_domain_parity_t) :: report, report_plus, report_minus, invalid
    real(dp) :: expected_absolute, expected_norm, expected_relative
    real(dp) :: comparison_norm_dot, absolute_difference_dot
    real(dp) :: relative_difference_dot, relative_per_distance_dot
    real(dp) :: distance_increase_dot, distance_ratio_dot
    real(dp) :: inner_distance_dot, outer_distance_dot
    real(dp) :: expected_comparison_norm_dot, expected_absolute_dot
    real(dp) :: expected_relative_dot, expected_relative_per_distance_dot
    real(dp) :: expected_distance_increase_dot, expected_distance_ratio_dot
    real(dp) :: inner_squared, outer_squared, difference_squared
    real(dp) :: inner_squared_dot, outer_squared_dot, difference_squared_dot
    complex(dp) :: field_difference(component_count, sample_count)
    complex(dp) :: field_difference_dot(component_count, sample_count)
    real(dp), parameter :: epsilon_fd = 1.0e-7_dp
    integer :: status
    integer :: sample, component
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
    inner_field_dot = reshape([ &
        cmplx(0.1_dp, -0.2_dp, dp), cmplx(-0.2_dp, 0.3_dp, dp), &
        cmplx(0.3_dp, 0.1_dp, dp), cmplx(-0.4_dp, 0.2_dp, dp), &
        cmplx(0.5_dp, -0.1_dp, dp), cmplx(-0.1_dp, -0.3_dp, dp)], &
        shape(inner_field_dot))
    outer_field_dot = reshape([ &
        cmplx(-0.2_dp, 0.4_dp, dp), cmplx(0.1_dp, -0.1_dp, dp), &
        cmplx(0.2_dp, 0.3_dp, dp), cmplx(0.3_dp, -0.4_dp, dp), &
        cmplx(-0.2_dp, 0.2_dp, dp), cmplx(0.4_dp, 0.1_dp, dp)], &
        shape(outer_field_dot))
    weights_dot = [0.1_dp, -0.2_dp, 0.3_dp]
    inner_distance_dot = 0.03_dp
    outer_distance_dot = -0.04_dp

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

    field_difference = outer_field - inner_field
    field_difference_dot = outer_field_dot - inner_field_dot
    inner_squared = 0.0_dp
    outer_squared = 0.0_dp
    difference_squared = 0.0_dp
    inner_squared_dot = 0.0_dp
    outer_squared_dot = 0.0_dp
    difference_squared_dot = 0.0_dp
    do sample = 1, sample_count
        do component = 1, component_count
            inner_squared = inner_squared + weights(sample)* &
                abs(inner_field(component, sample))**2
            outer_squared = outer_squared + weights(sample)* &
                abs(outer_field(component, sample))**2
            difference_squared = difference_squared + weights(sample)* &
                abs(field_difference(component, sample))**2
            inner_squared_dot = inner_squared_dot + weights_dot(sample)* &
                abs(inner_field(component, sample))**2 + &
                2.0_dp*weights(sample)*real(conjg(inner_field(component, sample))* &
                inner_field_dot(component, sample), dp)
            outer_squared_dot = outer_squared_dot + weights_dot(sample)* &
                abs(outer_field(component, sample))**2 + &
                2.0_dp*weights(sample)*real(conjg(outer_field(component, sample))* &
                outer_field_dot(component, sample), dp)
            difference_squared_dot = difference_squared_dot + weights_dot(sample)* &
                abs(field_difference(component, sample))**2 + &
                2.0_dp*weights(sample)*real(conjg(field_difference(component, sample))* &
                field_difference_dot(component, sample), dp)
        end do
    end do
    expected_comparison_norm_dot = outer_squared_dot/(2.0_dp*sqrt(outer_squared))
    expected_absolute_dot = difference_squared_dot/(2.0_dp*sqrt(difference_squared))
    expected_relative_dot = (expected_absolute_dot*report%comparison_norm - &
        report%absolute_difference*expected_comparison_norm_dot)/ &
        report%comparison_norm**2
    expected_distance_increase_dot = outer_distance_dot - inner_distance_dot
    expected_distance_ratio_dot = outer_distance_dot/report%inner_boundary_distance - &
        report%outer_boundary_distance*inner_distance_dot/ &
        report%inner_boundary_distance**2
    expected_relative_per_distance_dot = &
        expected_relative_dot/report%distance_increase - &
        report%relative_difference*expected_distance_increase_dot/ &
        report%distance_increase**2
    call evaluate_larger_domain_parity_jvp( &
        inner_field, outer_field, weights, inner_contract, outer_contract, &
        0.5_dp, 1.5_dp, 0.30_dp, 0.08_dp, inner_field_dot, outer_field_dot, &
        weights_dot, inner_distance_dot, outer_distance_dot, comparison_norm_dot, &
        absolute_difference_dot, relative_difference_dot, relative_per_distance_dot, &
        distance_increase_dot, distance_ratio_dot, status)
    call record_condition(status == 0 .and. &
        abs(comparison_norm_dot - expected_comparison_norm_dot) < 1.0e-12_dp .and. &
        abs(absolute_difference_dot - expected_absolute_dot) < 1.0e-12_dp .and. &
        abs(relative_difference_dot - expected_relative_dot) < 1.0e-12_dp .and. &
        abs(relative_per_distance_dot - expected_relative_per_distance_dot) < 1.0e-12_dp .and. &
        abs(distance_increase_dot - expected_distance_increase_dot) < 1.0e-12_dp .and. &
        abs(distance_ratio_dot - expected_distance_ratio_dot) < 1.0e-12_dp, &
        "larger-domain parity JVP matches independent weighted metric oracle")

    inner_field_plus = inner_field + epsilon_fd*inner_field_dot
    outer_field_plus = outer_field + epsilon_fd*outer_field_dot
    weights_plus = weights + epsilon_fd*weights_dot
    inner_field_minus = inner_field - epsilon_fd*inner_field_dot
    outer_field_minus = outer_field - epsilon_fd*outer_field_dot
    weights_minus = weights - epsilon_fd*weights_dot
    call evaluate_larger_domain_parity( &
        inner_field_plus, outer_field_plus, weights_plus, inner_contract, outer_contract, &
        0.5_dp + epsilon_fd*inner_distance_dot, &
        1.5_dp + epsilon_fd*outer_distance_dot, 0.30_dp, 0.08_dp, report_plus, status)
    call evaluate_larger_domain_parity( &
        inner_field_minus, outer_field_minus, weights_minus, inner_contract, outer_contract, &
        0.5_dp - epsilon_fd*inner_distance_dot, &
        1.5_dp - epsilon_fd*outer_distance_dot, 0.30_dp, 0.08_dp, report_minus, status)
    call record_condition( &
        abs(comparison_norm_dot - (report_plus%comparison_norm - &
        report_minus%comparison_norm)/(2.0_dp*epsilon_fd)) < 1.0e-7_dp .and. &
        abs(absolute_difference_dot - (report_plus%absolute_difference - &
        report_minus%absolute_difference)/(2.0_dp*epsilon_fd)) < 1.0e-7_dp .and. &
        abs(relative_per_distance_dot - (report_plus%relative_difference_per_distance - &
        report_minus%relative_difference_per_distance)/(2.0_dp*epsilon_fd)) < 1.0e-7_dp, &
        "larger-domain parity JVP matches complete central re-evaluation")

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
