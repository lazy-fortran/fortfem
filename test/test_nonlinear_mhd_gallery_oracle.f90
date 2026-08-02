program test_nonlinear_mhd_gallery_oracle
    !! Independent oracle for the neutral island/wall continuation fixture.
    use check, only: check_condition, check_summary
    use fortfem_interop, only: &
        compare_resistive_mhd_branch_histories, &
        evaluate_resistive_mhd_branch_diagnostics, &
        evaluate_resistive_mhd_branch_path_metric, &
        initialize_resistive_mhd_branch_history, &
        resistive_mhd_branch_diagnostics_t, &
        resistive_mhd_branch_history_t
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    integer, parameter :: state_size = 3, sample_count = 5
    real(dp), parameter :: parameter_values(sample_count) = &
        [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
    real(dp), parameter :: tolerance = 1.0e-12_dp
    real(dp) :: forward_state(state_size, sample_count)
    real(dp) :: reverse_state(state_size, sample_count)
    real(dp) :: residual(state_size, sample_count), energy(sample_count)
    real(dp) :: expected_metric, metric
    integer, parameter :: forward_branch(sample_count) = [1, 1, 1, 2, 2]
    integer, parameter :: reverse_branch(sample_count) = [1, 1, 2, 2, 2]
    type(resistive_mhd_branch_history_t) :: forward, reverse
    type(resistive_mhd_branch_diagnostics_t) :: diagnostics
    type(fortsparse_status_t) :: status
    integer :: sample

    do sample = 1, sample_count
        forward_state(:, sample) = manufactured_state(parameter_values(sample), &
            0.0_dp)
        reverse_state(:, sample) = manufactured_state(parameter_values(sample), &
            hysteresis_offset(parameter_values(sample)))
        residual(:, sample) = 0.0_dp
        energy(sample) = 0.5_dp*dot_product(forward_state(:, sample), &
            forward_state(:, sample))
    end do

    call initialize_resistive_mhd_branch_history( &
        parameter_values, forward_state, residual, energy, forward_branch, &
        forward, status)
    call check_condition(status%code == FORTSPARSE_OK, &
        "manufactured forward island/wall path is accepted")
    call evaluate_resistive_mhd_branch_diagnostics(forward, diagnostics, status)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        diagnostics%branch_multiplicity == 2 .and. &
        .not. diagnostics%hysteresis_detected, &
        "forward path has the expected two branch labels")

    call initialize_resistive_mhd_branch_history( &
        parameter_values, reverse_state, residual, energy, reverse_branch, &
        reverse, status)
    call compare_resistive_mhd_branch_histories( &
        forward, reverse, tolerance, diagnostics, status)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        diagnostics%branch_multiplicity == 2 .and. &
        diagnostics%hysteresis_detected .and. &
        diagnostics%hysteresis_sample_count == 3 .and. &
        abs(diagnostics%max_state_hysteresis - &
            hysteresis_offset(0.5_dp)) < tolerance, &
        "forward/reverse paths expose the manufactured island hysteresis")

    call evaluate_resistive_mhd_branch_path_metric(forward, metric, status)
    expected_metric = independent_path_metric(parameter_values, forward_state, energy)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        abs(metric - expected_metric) < tolerance, &
        "branch path metric agrees with the independent trapezoidal oracle")

    call check_summary("nonlinear island/wall gallery oracle")

contains

    pure function manufactured_state(parameter_value, offset) result(state)
        real(dp), intent(in) :: parameter_value, offset
        real(dp) :: state(state_size)

        state(1) = 0.18_dp + 0.24_dp*parameter_value + &
            0.04_dp*parameter_value**2 + offset
        state(2) = 0.65_dp*state(1) + 0.08_dp*parameter_value
        state(3) = 0.05_dp*parameter_value
    end function manufactured_state

    pure real(dp) function hysteresis_offset(parameter_value)
        real(dp), intent(in) :: parameter_value

        hysteresis_offset = 0.12_dp*parameter_value*(1.0_dp - parameter_value)
    end function hysteresis_offset

    pure real(dp) function independent_path_metric(parameter_value, state, energy)
        real(dp), intent(in) :: parameter_value(:), state(:, :), energy(:)
        integer :: index
        real(dp) :: integrand(size(parameter_value))

        integrand = energy + 0.5_dp*sum(state**2, dim=1)
        independent_path_metric = 0.0_dp
        do index = 1, size(parameter_value) - 1
            independent_path_metric = independent_path_metric + 0.5_dp* &
                (parameter_value(index + 1) - parameter_value(index))* &
                (integrand(index + 1) + integrand(index))
        end do
    end function independent_path_metric

end program test_nonlinear_mhd_gallery_oracle
