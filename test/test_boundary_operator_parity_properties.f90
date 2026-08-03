program test_boundary_operator_parity_properties
    use, intrinsic :: iso_fortran_env, only: int32
    use check, only: check_condition, check_property, check_summary, &
        property_random_unit, property_rng_initialize, property_rng_t
    use fortfem_interop, only: &
        BOUNDARY_OPERATOR_BACKEND_BEM, &
        BOUNDARY_OPERATOR_BACKEND_DTN, &
        BOUNDARY_OPERATOR_BACKEND_FEM, &
        boundary_operator_contract_t, boundary_operator_parity_t, &
        compare_boundary_operator_parity, &
        compare_boundary_operator_parity_jvp, &
        compare_boundary_operator_parity_vjp, &
        initialize_boundary_operator_contract
    use fortfem_kinds, only: dp
    implicit none

    logical :: all_passed
    integer(int32) :: first_failed_seed, shrunk_seed

    call check_property( &
        "seeded complex boundary parity value and derivative identities", &
        1357911_int32, 12, parity_case, all_passed, first_failed_seed, &
        shrunk_seed)
    call check_condition(all_passed .and. first_failed_seed == 0_int32 .and. &
        shrunk_seed == 0_int32, &
        "random boundary parity property reports no failure seed")
    call check_summary("random boundary operator parity properties")
    if (.not. all_passed) error stop 1

contains

    logical function parity_case(case_seed)
        integer(int32), intent(in) :: case_seed
        integer, parameter :: component_count = 2
        integer, parameter :: sample_count = 4
        integer, parameter :: backend_count = 3
        real(dp), parameter :: absolute_tolerance = 0.05_dp
        real(dp), parameter :: relative_tolerance = 0.10_dp
        real(dp), parameter :: epsilon_fd = 1.0e-6_dp
        type(property_rng_t) :: rng
        type(boundary_operator_contract_t) :: contracts(backend_count)
        type(boundary_operator_parity_t) :: report, report_plus, report_minus
        complex(dp) :: reference(component_count, sample_count)
        complex(dp) :: candidates(component_count, sample_count, backend_count)
        complex(dp) :: reference_dot(component_count, sample_count)
        complex(dp) :: candidates_dot(component_count, sample_count, backend_count)
        complex(dp) :: reference_plus(component_count, sample_count)
        complex(dp) :: reference_minus(component_count, sample_count)
        complex(dp) :: candidates_plus(component_count, sample_count, backend_count)
        complex(dp) :: candidates_minus(component_count, sample_count, backend_count)
        complex(dp) :: zero_reference(component_count, sample_count)
        complex(dp) :: reference_bar(component_count, sample_count)
        complex(dp) :: candidates_bar(component_count, sample_count, backend_count)
        real(dp) :: weights(sample_count), weights_dot(sample_count)
        real(dp) :: weights_plus(sample_count), weights_minus(sample_count)
        real(dp) :: reference_norm_dot
        real(dp) :: absolute_error_dot(backend_count), relative_error_dot(backend_count)
        real(dp) :: absolute_error_bar(backend_count), relative_error_bar(backend_count)
        real(dp) :: weights_bar(sample_count)
        real(dp) :: expected_reference_norm, expected_absolute(backend_count)
        real(dp) :: expected_relative(backend_count)
        real(dp) :: reference_squared, error_squared
        real(dp) :: reference_norm_plus, reference_norm_minus
        real(dp) :: absolute_error_plus(backend_count)
        real(dp) :: absolute_error_minus(backend_count)
        real(dp) :: relative_error_plus(backend_count)
        real(dp) :: relative_error_minus(backend_count)
        real(dp) :: reference_norm_bar, lhs, rhs
        integer :: backend, component, sample, status
        integer :: backend_kind(backend_count)

        parity_case = .false.
        call property_rng_initialize(rng, case_seed)
        backend_kind = [BOUNDARY_OPERATOR_BACKEND_FEM, &
            BOUNDARY_OPERATOR_BACKEND_BEM, BOUNDARY_OPERATOR_BACKEND_DTN]

        do sample = 1, sample_count
            weights(sample) = 0.5_dp + 1.5_dp*property_random_unit(rng)
            weights_dot(sample) = 0.4_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
            do component = 1, component_count
                reference(component, sample) = random_complex(rng)
                reference_dot(component, sample) = 0.5_dp*random_complex(rng)
            end do
        end do
        reference(1, 1) = reference(1, 1) + cmplx(1.0_dp, 0.5_dp, dp)
        do backend = 1, backend_count
            do sample = 1, sample_count
                do component = 1, component_count
                    candidates(component, sample, backend) = &
                        reference(component, sample) + &
                        0.25_dp*random_complex(rng)
                    candidates_dot(component, sample, backend) = &
                        0.5_dp*random_complex(rng)
                end do
            end do
            candidates(1, 1, backend) = reference(1, 1) + &
                cmplx(0.15_dp*real(backend, dp), -0.12_dp*real(backend, dp), dp)
        end do

        do backend = 1, backend_count
            call initialize_boundary_operator_contract( &
                contracts(backend), backend_kind(backend), "helmholtz", &
                "H1-trace", component_count, sample_count, .true., .true., &
                .true., .true., .true., .true., "unit", "manufactured", &
                "seeded parity property", "property-circle", status)
            if (status /= 0) return
        end do

        call compare_boundary_operator_parity( &
            reference, candidates, weights, contracts, absolute_tolerance, &
            relative_tolerance, report, status)
        if (status /= 0) return
        reference_squared = weighted_squared_norm(reference, weights)
        expected_reference_norm = sqrt(reference_squared)
        if (abs(report%reference_norm - expected_reference_norm) > 2.0e-13_dp) return
        if (expected_reference_norm <= tiny(1.0_dp)) return
        do backend = 1, backend_count
            error_squared = weighted_difference_squared_norm( &
                candidates(:, :, backend), reference, weights)
            expected_absolute(backend) = sqrt(error_squared)
            expected_relative(backend) = expected_absolute(backend)/ &
                expected_reference_norm
        end do
        if (maxval(abs(report%absolute_error - expected_absolute)) > 2.0e-13_dp .or. &
            maxval(abs(report%relative_error - expected_relative)) > 2.0e-13_dp) return

        call compare_boundary_operator_parity_jvp( &
            reference, candidates, weights, contracts, absolute_tolerance, &
            relative_tolerance, reference_dot, candidates_dot, weights_dot, &
            reference_norm_dot, absolute_error_dot, relative_error_dot, status)
        if (status /= 0) return
        reference_plus = reference + epsilon_fd*reference_dot
        reference_minus = reference - epsilon_fd*reference_dot
        candidates_plus = candidates + epsilon_fd*candidates_dot
        candidates_minus = candidates - epsilon_fd*candidates_dot
        weights_plus = weights + epsilon_fd*weights_dot
        weights_minus = weights - epsilon_fd*weights_dot
        call compare_boundary_operator_parity( &
            reference_plus, candidates_plus, weights_plus, contracts, &
            absolute_tolerance, relative_tolerance, report_plus, status)
        if (status /= 0) return
        call compare_boundary_operator_parity( &
            reference_minus, candidates_minus, weights_minus, contracts, &
            absolute_tolerance, relative_tolerance, report_minus, status)
        if (status /= 0) return
        reference_norm_plus = report_plus%reference_norm
        reference_norm_minus = report_minus%reference_norm
        absolute_error_plus = report_plus%absolute_error
        absolute_error_minus = report_minus%absolute_error
        relative_error_plus = report_plus%relative_error
        relative_error_minus = report_minus%relative_error
        if (abs(reference_norm_dot - (reference_norm_plus - reference_norm_minus)/ &
            (2.0_dp*epsilon_fd)) > 2.0e-7_dp .or. &
            maxval(abs(absolute_error_dot - (absolute_error_plus - &
            absolute_error_minus)/(2.0_dp*epsilon_fd))) > 2.0e-7_dp .or. &
            maxval(abs(relative_error_dot - (relative_error_plus - &
            relative_error_minus)/(2.0_dp*epsilon_fd))) > 2.0e-7_dp) return

        reference_norm_bar = 2.0_dp*property_random_unit(rng) - 1.0_dp
        do backend = 1, backend_count
            absolute_error_bar(backend) = 2.0_dp*property_random_unit(rng) - 1.0_dp
            relative_error_bar(backend) = 2.0_dp*property_random_unit(rng) - 1.0_dp
        end do
        call compare_boundary_operator_parity_vjp( &
            reference, candidates, weights, contracts, absolute_tolerance, &
            relative_tolerance, reference_norm_bar, absolute_error_bar, &
            relative_error_bar, reference_bar, candidates_bar, weights_bar, status)
        if (status /= 0) return
        lhs = reference_norm_bar*reference_norm_dot + &
            dot_product(absolute_error_bar, absolute_error_dot) + &
            dot_product(relative_error_bar, relative_error_dot)
        rhs = real(sum(conjg(reference_bar)*reference_dot), dp) + &
            real(sum(conjg(candidates_bar)*candidates_dot), dp) + &
            dot_product(weights_bar, weights_dot)
        if (abs(lhs - rhs) > 2.0e-8_dp) return

        zero_reference = cmplx(0.0_dp, 0.0_dp, dp)
        call compare_boundary_operator_parity_jvp( &
            zero_reference, candidates, weights, contracts, absolute_tolerance, &
            relative_tolerance, reference_dot, candidates_dot, weights_dot, &
            reference_norm_dot, absolute_error_dot, relative_error_dot, status)
        if (status == 0) return
        call compare_boundary_operator_parity_vjp( &
            zero_reference, candidates, weights, contracts, absolute_tolerance, &
            relative_tolerance, reference_norm_bar, absolute_error_bar, &
            relative_error_bar, reference_bar, candidates_bar, weights_bar, status)
        if (status == 0) return

        parity_case = .true.
    end function parity_case

    function random_complex(rng) result(value)
        type(property_rng_t), intent(inout) :: rng
        complex(dp) :: value

        value = cmplx(2.0_dp*property_random_unit(rng) - 1.0_dp, &
            2.0_dp*property_random_unit(rng) - 1.0_dp, dp)
    end function random_complex

    pure real(dp) function weighted_squared_norm(values, weights) result(norm_squared)
        complex(dp), intent(in) :: values(:, :)
        real(dp), intent(in) :: weights(:)
        integer :: component, sample

        norm_squared = 0.0_dp
        do sample = 1, size(values, 2)
            do component = 1, size(values, 1)
                norm_squared = norm_squared + weights(sample)*abs( &
                    values(component, sample))**2
            end do
        end do
    end function weighted_squared_norm

    pure real(dp) function weighted_difference_squared_norm( &
            values, reference, weights) result(norm_squared)
        complex(dp), intent(in) :: values(:, :), reference(:, :)
        real(dp), intent(in) :: weights(:)
        integer :: component, sample

        norm_squared = 0.0_dp
        do sample = 1, size(values, 2)
            do component = 1, size(values, 1)
                norm_squared = norm_squared + weights(sample)*abs( &
                    values(component, sample) - reference(component, sample))**2
            end do
        end do
    end function weighted_difference_squared_norm

end program test_boundary_operator_parity_properties
