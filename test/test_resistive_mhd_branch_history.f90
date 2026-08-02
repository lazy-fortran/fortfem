program test_resistive_mhd_branch_history
    use check, only: check_condition, check_summary
    use fortfem_interop, only: &
        compare_resistive_mhd_branch_histories, &
        evaluate_resistive_mhd_branch_diagnostics, &
        evaluate_resistive_mhd_branch_path_metric, &
        evaluate_resistive_mhd_branch_path_metric_jvp, &
        evaluate_resistive_mhd_branch_path_metric_vjp, &
        initialize_resistive_mhd_branch_history, &
        resistive_mhd_branch_diagnostics_t, &
        resistive_mhd_branch_history_t, &
        validate_resistive_mhd_branch_history
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    integer, parameter :: state_size = 2, sample_count = 3
    real(dp), parameter :: epsilon = 1.0e-7_dp
    real(dp) :: parameter(sample_count), state(state_size, sample_count)
    real(dp) :: residual(state_size, sample_count), energy(sample_count)
    real(dp) :: state_dot(state_size, sample_count)
    real(dp) :: residual_dot(state_size, sample_count), energy_dot(sample_count)
    real(dp) :: parameter_dot(sample_count), metric, metric_dot
    real(dp) :: metric_plus, metric_minus, metric_bar
    real(dp) :: parameter_bar(sample_count), state_bar(state_size, sample_count)
    real(dp) :: residual_bar(state_size, sample_count), energy_bar(sample_count)
    integer :: branch_id(sample_count), reverse_branch_id(sample_count)
    type(resistive_mhd_branch_history_t) :: history, reverse_history
    type(resistive_mhd_branch_history_t) :: history_plus_obj, history_minus_obj
    type(resistive_mhd_branch_history_t) :: invalid_history
    type(resistive_mhd_branch_diagnostics_t) :: diagnostics
    type(fortsparse_status_t) :: status
    integer :: index

    parameter = [0.0_dp, 0.5_dp, 1.0_dp]
    state(:, 1) = [1.0_dp, -0.5_dp]
    state(:, 2) = [1.5_dp, -0.25_dp]
    state(:, 3) = [2.0_dp, 0.25_dp]
    residual(:, 1) = [0.2_dp, -0.1_dp]
    residual(:, 2) = [0.1_dp, 0.3_dp]
    residual(:, 3) = [-0.2_dp, 0.4_dp]
    energy = [0.5_dp, 0.75_dp, 1.25_dp]
    branch_id = [1, 1, 2]
    call initialize_resistive_mhd_branch_history( &
        parameter, state, residual, energy, branch_id, history, status)
    call check_condition(status%code == FORTSPARSE_OK, &
        "strictly increasing branch history is accepted")

    call evaluate_resistive_mhd_branch_diagnostics(history, diagnostics, status)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        diagnostics%branch_multiplicity == 2 .and. &
        .not. diagnostics%hysteresis_detected, &
        "branch multiplicity is deterministic for a single history")

    call evaluate_resistive_mhd_branch_path_metric(history, metric, status)
    call check_condition(status%code == FORTSPARSE_OK .and. metric > 0.0_dp, &
        "branch path metric integrates caller-owned samples")

    parameter_dot = [0.1_dp, -0.2_dp, 0.3_dp]
    state_dot(:, 1) = [0.2_dp, -0.1_dp]
    state_dot(:, 2) = [-0.3_dp, 0.4_dp]
    state_dot(:, 3) = [0.1_dp, 0.2_dp]
    residual_dot(:, 1) = [0.1_dp, 0.2_dp]
    residual_dot(:, 2) = [-0.2_dp, 0.1_dp]
    residual_dot(:, 3) = [0.3_dp, -0.1_dp]
    energy_dot = [0.2_dp, -0.1_dp, 0.4_dp]
    call evaluate_resistive_mhd_branch_path_metric_jvp( &
        history, parameter_dot, state_dot, residual_dot, energy_dot, metric_dot, status)
    call initialize_resistive_mhd_branch_history( &
        parameter + epsilon*parameter_dot, state + epsilon*state_dot, &
        residual + epsilon*residual_dot, energy + epsilon*energy_dot, branch_id, &
        history_plus_obj, status)
    call evaluate_resistive_mhd_branch_path_metric(history_plus_obj, metric_plus, status)
    call initialize_resistive_mhd_branch_history( &
        parameter - epsilon*parameter_dot, state - epsilon*state_dot, &
        residual - epsilon*residual_dot, energy - epsilon*energy_dot, branch_id, &
        history_minus_obj, status)
    call evaluate_resistive_mhd_branch_path_metric(history_minus_obj, metric_minus, status)
    call check_condition(abs(metric_dot - (metric_plus - metric_minus) / &
        (2.0_dp*epsilon)) < 2.0e-6_dp, &
        "branch path metric JVP matches a central-difference oracle")

    metric_bar = 1.7_dp
    call evaluate_resistive_mhd_branch_path_metric_vjp( &
        history, metric_bar, parameter_bar, state_bar, residual_bar, energy_bar, status)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        abs(dot_product(parameter_bar, parameter_dot) + &
        sum(state_bar*state_dot) + sum(residual_bar*residual_dot) + &
        dot_product(energy_bar, energy_dot) - metric_bar*metric_dot) < 1.0e-12_dp, &
        "branch path metric VJP satisfies the full dot-product identity")

    reverse_branch_id = [1, 1, 3]
    call initialize_resistive_mhd_branch_history( &
        parameter, state + reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.2_dp, -0.1_dp], [state_size, sample_count]), residual, &
        energy + [0.0_dp, 0.0_dp, 0.3_dp], reverse_branch_id, reverse_history, status)
    call compare_resistive_mhd_branch_histories( &
        history, reverse_history, 1.0e-12_dp, diagnostics, status)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        diagnostics%branch_multiplicity == 3 .and. diagnostics%hysteresis_detected .and. &
        diagnostics%hysteresis_sample_count == 1 .and. &
        diagnostics%max_energy_hysteresis > 0.0_dp, &
        "forward/reverse histories report deterministic multiplicity and hysteresis")

    parameter(2) = parameter(1)
    call initialize_resistive_mhd_branch_history( &
        parameter, state, residual, energy, branch_id, history, status)
    call check_condition(status%code /= FORTSPARSE_OK, &
        "non-monotone continuation parameters are rejected")

    call validate_resistive_mhd_branch_history(invalid_history, status)
    call check_condition(status%code /= FORTSPARSE_OK, &
        "unallocated branch histories are rejected")

    call check_summary("resistive-MHD branch history")

end program test_resistive_mhd_branch_history
