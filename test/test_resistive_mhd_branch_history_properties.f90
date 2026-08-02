program test_resistive_mhd_branch_history_properties
    !! Seeded property campaign for neutral branch-history diagnostics.
    !!
    !! All samples and perturbations are caller-owned.  Independent loops
    !! reproduce the path metric, branch-label union, and hysteresis values;
    !! no implementation helper is used as an oracle.
    use, intrinsic :: iso_fortran_env, only: int32
    use check, only: check_condition, check_property, check_summary, &
        property_random_integer, property_random_unit, property_rng_initialize, &
        property_rng_t
    use fortfem_interop, only: &
        compare_resistive_mhd_branch_histories, &
        evaluate_resistive_mhd_branch_diagnostics, &
        evaluate_resistive_mhd_branch_path_metric, &
        evaluate_resistive_mhd_branch_path_metric_jvp, &
        evaluate_resistive_mhd_branch_path_metric_vjp, &
        initialize_resistive_mhd_branch_history, &
        resistive_mhd_branch_diagnostics_t, &
        resistive_mhd_branch_history_t
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    logical :: all_passed
    integer(int32) :: first_failed_seed, shrunk_seed

    call check_property( &
        "random resistive-MHD branch history diagnostics and metric identities", &
        20260802_int32, 20, branch_case, all_passed, first_failed_seed, shrunk_seed)
    call check_condition(all_passed .and. first_failed_seed == 0_int32 .and. &
        shrunk_seed == 0_int32, &
        "random branch-history property reports no failure seed")
    call check_summary("random resistive-MHD branch history properties")
    if (.not. all_passed) error stop 1

contains

    logical function branch_case(case_seed)
        integer(int32), intent(in) :: case_seed
        type(property_rng_t) :: rng
        type(fortsparse_status_t) :: status
        type(resistive_mhd_branch_history_t) :: history, history_plus, history_minus
        type(resistive_mhd_branch_history_t) :: history_second, history_bad
        type(resistive_mhd_branch_diagnostics_t) :: diagnostics, diagnostics_repeat
        real(dp), allocatable :: parameter(:), parameter_dot(:), parameter_bad(:)
        real(dp), allocatable :: state(:, :), state_dot(:, :)
        real(dp), allocatable :: residual(:, :), residual_dot(:, :)
        real(dp), allocatable :: state_second(:, :), residual_second(:, :)
        real(dp), allocatable :: energy(:), energy_dot(:), energy_second(:)
        real(dp), allocatable :: parameter_bar(:), state_bar(:, :)
        real(dp), allocatable :: residual_bar(:, :), energy_bar(:)
        integer, allocatable :: branch_id(:), branch_id_second(:)
        integer :: state_size, sample_count, sample, component, direction, changed
        integer :: expected_multiplicity, expected_hysteresis
        real(dp) :: metric, metric_expected, metric_dot, metric_dot_expected
        real(dp) :: metric_plus, metric_minus, metric_bar, left, right
        real(dp), parameter :: epsilon = 1.0e-6_dp
        real(dp), parameter :: tolerance = 1.0e-12_dp

        branch_case = .false.
        call property_rng_initialize(rng, case_seed)
        state_size = property_random_integer(rng, 1, 3)
        sample_count = property_random_integer(rng, 2, 7)
        direction = merge(1, -1, property_random_integer(rng, 0, 1) == 0)
        allocate(parameter(sample_count), parameter_dot(sample_count), &
            parameter_bad(sample_count), state(state_size, sample_count), &
            state_dot(state_size, sample_count), residual(state_size, sample_count), &
            residual_dot(state_size, sample_count), energy(sample_count), &
            energy_dot(sample_count), branch_id(sample_count))
        parameter(1) = direction*(0.5_dp + property_random_unit(rng))
        do sample = 2, sample_count
            parameter(sample) = parameter(sample - 1) + direction*( &
                0.05_dp + 0.35_dp*property_random_unit(rng))
        end do
        do sample = 1, sample_count
            parameter_dot(sample) = 0.5_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
            energy(sample) = 0.2_dp + 1.6_dp*property_random_unit(rng)
            energy_dot(sample) = 0.6_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
            branch_id(sample) = property_random_integer(rng, 1, 4)
            do component = 1, state_size
                state(component, sample) = &
                    1.4_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
                residual(component, sample) = &
                    1.1_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
                state_dot(component, sample) = &
                    0.5_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
                residual_dot(component, sample) = &
                    0.45_dp*(2.0_dp*property_random_unit(rng) - 1.0_dp)
            end do
        end do

        call initialize_resistive_mhd_branch_history( &
            parameter, state, residual, energy, branch_id, history, status, direction)
        if (status%code /= FORTSPARSE_OK) return
        call evaluate_resistive_mhd_branch_diagnostics(history, diagnostics, status)
        if (status%code /= FORTSPARSE_OK) return
        call evaluate_resistive_mhd_branch_diagnostics( &
            history, diagnostics_repeat, status)
        if (status%code /= FORTSPARSE_OK) return
        expected_multiplicity = count_unique_oracle(branch_id)
        if (diagnostics%branch_multiplicity /= expected_multiplicity .or. &
            diagnostics_repeat%branch_multiplicity /= expected_multiplicity .or. &
            diagnostics%hysteresis_detected .or. &
            diagnostics%hysteresis_sample_count /= 0) return

        call evaluate_resistive_mhd_branch_path_metric(history, metric, status)
        if (status%code /= FORTSPARSE_OK) return
        metric_expected = metric_oracle(parameter, state, residual, energy)
        if (abs(metric - metric_expected) > 3.0e-13_dp) return

        call evaluate_resistive_mhd_branch_path_metric_jvp( &
            history, parameter_dot, state_dot, residual_dot, energy_dot, metric_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        metric_dot_expected = metric_dot_oracle( &
            parameter, parameter_dot, state, state_dot, residual, residual_dot, &
            energy, energy_dot)
        if (abs(metric_dot - metric_dot_expected) > 4.0e-13_dp) return

        call initialize_resistive_mhd_branch_history( &
            parameter + epsilon*parameter_dot, state + epsilon*state_dot, &
            residual + epsilon*residual_dot, energy + epsilon*energy_dot, branch_id, &
            history_plus, status, direction)
        if (status%code /= FORTSPARSE_OK) return
        call evaluate_resistive_mhd_branch_path_metric(history_plus, metric_plus, status)
        if (status%code /= FORTSPARSE_OK) return
        call initialize_resistive_mhd_branch_history( &
            parameter - epsilon*parameter_dot, state - epsilon*state_dot, &
            residual - epsilon*residual_dot, energy - epsilon*energy_dot, branch_id, &
            history_minus, status, direction)
        if (status%code /= FORTSPARSE_OK) return
        call evaluate_resistive_mhd_branch_path_metric(history_minus, metric_minus, status)
        if (status%code /= FORTSPARSE_OK .or. &
            abs(metric_dot - (metric_plus - metric_minus)/(2.0_dp*epsilon)) > 4.0e-8_dp) return

        allocate(parameter_bar(sample_count), state_bar(state_size, sample_count), &
            residual_bar(state_size, sample_count), energy_bar(sample_count))
        metric_bar = 0.4_dp + 1.4_dp*property_random_unit(rng)
        call evaluate_resistive_mhd_branch_path_metric_vjp( &
            history, metric_bar, parameter_bar, state_bar, residual_bar, energy_bar, status)
        left = dot_product(parameter_bar, parameter_dot) + sum(state_bar*state_dot) + &
            sum(residual_bar*residual_dot) + dot_product(energy_bar, energy_dot)
        right = metric_bar*metric_dot
        if (status%code /= FORTSPARSE_OK .or. abs(left - right) > &
            3.0e-12_dp*max(1.0_dp, abs(left), abs(right))) return

        allocate(state_second(state_size, sample_count), &
            residual_second(state_size, sample_count), energy_second(sample_count), &
            branch_id_second(sample_count))
        state_second = state
        residual_second = residual
        energy_second = energy
        changed = property_random_integer(rng, 1, sample_count)
        state_second(:, changed) = state_second(:, changed) + 0.2_dp
        residual_second(:, changed) = residual_second(:, changed) - 0.13_dp
        energy_second(changed) = energy_second(changed) + 0.3_dp
        branch_id_second = branch_id
        branch_id_second(changed) = maxval(branch_id) + 1
        call initialize_resistive_mhd_branch_history( &
            parameter, state_second, residual_second, energy_second, &
            branch_id_second, history_second, status, direction)
        if (status%code /= FORTSPARSE_OK) return
        call compare_resistive_mhd_branch_histories( &
            history, history_second, tolerance, diagnostics, status)
        expected_multiplicity = count_unique_oracle([branch_id, branch_id_second])
        expected_hysteresis = 1
        if (status%code /= FORTSPARSE_OK .or. &
            diagnostics%branch_multiplicity /= expected_multiplicity .or. &
            .not. diagnostics%hysteresis_detected .or. &
            diagnostics%hysteresis_sample_count /= expected_hysteresis .or. &
            abs(diagnostics%max_state_hysteresis - 0.2_dp) > tolerance .or. &
            abs(diagnostics%max_residual_hysteresis - 0.13_dp) > tolerance .or. &
            abs(diagnostics%max_energy_hysteresis - 0.3_dp) > tolerance) return

        parameter_bad = parameter
        parameter_bad(2) = parameter_bad(1)
        call initialize_resistive_mhd_branch_history( &
            parameter_bad, state, residual, energy, branch_id, history_bad, status, direction)
        if (status%code == FORTSPARSE_OK) return
        branch_case = .true.
    end function branch_case

    real(dp) function sample_integrand_oracle(state, residual, energy, sample) result(value)
        real(dp), intent(in) :: state(:, :), residual(:, :), energy(:)
        integer, intent(in) :: sample

        value = energy(sample) + 0.5_dp*( &
            dot_product(state(:, sample), state(:, sample)) + &
            dot_product(residual(:, sample), residual(:, sample)))
    end function sample_integrand_oracle

    real(dp) function metric_oracle(parameter, state, residual, energy) result(metric)
        real(dp), intent(in) :: parameter(:), state(:, :), residual(:, :), energy(:)
        integer :: sample

        metric = 0.0_dp
        do sample = 1, size(parameter) - 1
            metric = metric + 0.5_dp*( &
                sample_integrand_oracle(state, residual, energy, sample) + &
                sample_integrand_oracle(state, residual, energy, sample + 1))* &
                (parameter(sample + 1) - parameter(sample))
        end do
    end function metric_oracle

    real(dp) function metric_dot_oracle( &
            parameter, parameter_dot, state, state_dot, residual, residual_dot, &
            energy, energy_dot) result(metric_dot)
        real(dp), intent(in) :: parameter(:), parameter_dot(:), state(:, :)
        real(dp), intent(in) :: state_dot(:, :), residual(:, :), residual_dot(:, :)
        real(dp), intent(in) :: energy(:), energy_dot(:)
        integer :: sample
        real(dp) :: left, right, left_dot, right_dot
        real(dp) :: parameter_step, parameter_step_dot

        metric_dot = 0.0_dp
        do sample = 1, size(parameter) - 1
            left = sample_integrand_oracle(state, residual, energy, sample)
            right = sample_integrand_oracle(state, residual, energy, sample + 1)
            left_dot = energy_dot(sample) + dot_product( &
                state(:, sample), state_dot(:, sample)) + dot_product( &
                residual(:, sample), residual_dot(:, sample))
            right_dot = energy_dot(sample + 1) + dot_product( &
                state(:, sample + 1), state_dot(:, sample + 1)) + dot_product( &
                residual(:, sample + 1), residual_dot(:, sample + 1))
            parameter_step = parameter(sample + 1) - parameter(sample)
            parameter_step_dot = parameter_dot(sample + 1) - parameter_dot(sample)
            metric_dot = metric_dot + 0.5_dp*(left_dot + right_dot)*parameter_step + &
                0.5_dp*(left + right)*parameter_step_dot
        end do
    end function metric_dot_oracle

    integer function count_unique_oracle(labels) result(count)
        integer, intent(in) :: labels(:)
        integer :: index, previous
        logical :: seen

        count = 0
        do index = 1, size(labels)
            seen = .false.
            do previous = 1, index - 1
                if (labels(previous) == labels(index)) then
                    seen = .true.
                    exit
                end if
            end do
            if (.not. seen) count = count + 1
        end do
    end function count_unique_oracle

end program test_resistive_mhd_branch_history_properties
