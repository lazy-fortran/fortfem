program test_resistive_wall_midpoint_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        advance_resistive_wall_midpoint, advance_resistive_wall_midpoint_jvp, &
        advance_resistive_wall_midpoint_vjp, evaluate_resistive_wall_energy_balance
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 0.15_dp
    real(dp), parameter :: derivative_step = 1.0e-6_dp
    real(dp) :: inductance(2, 2), resistance(2, 2)
    real(dp) :: inductance_dot(2, 2), resistance_dot(2, 2)
    real(dp) :: inductance_bar(2, 2), resistance_bar(2, 2)
    real(dp) :: current(2), current_dot(2), current_bar(2)
    real(dp) :: voltage_n(2), voltage_next(2), voltage_n_dot(2)
    real(dp) :: voltage_next_dot(2), voltage_n_bar(2), voltage_next_bar(2)
    real(dp) :: current_next(2), current_next_dot(2), current_next_bar(2)
    real(dp) :: current_plus(2), current_minus(2), step_dot, step_bar
    real(dp) :: energy_n, energy_next, input_work, dissipation, balance
    real(dp) :: derivative_error, adjoint_error, lhs, rhs
    integer :: status
    logical :: all_passed

    all_passed = .true.
    inductance = reshape([2.0_dp, 0.2_dp, 0.2_dp, 1.0_dp], [2, 2])
    resistance = reshape([0.4_dp, 0.1_dp, 0.1_dp, 0.3_dp], [2, 2])
    inductance_dot = reshape([0.03_dp, -0.02_dp, -0.02_dp, 0.04_dp], [2, 2])
    resistance_dot = reshape([0.01_dp, 0.02_dp, 0.02_dp, -0.015_dp], [2, 2])
    current = [0.7_dp, -0.4_dp]
    current_dot = [0.03_dp, -0.02_dp]
    voltage_n = [0.8_dp, -0.3_dp]
    voltage_next = [0.5_dp, 0.2_dp]
    voltage_n_dot = [-0.01_dp, 0.04_dp]
    voltage_next_dot = [0.02_dp, -0.03_dp]
    step_dot = 0.011_dp

    call advance_resistive_wall_midpoint( &
        inductance, resistance, step, current, voltage_n, voltage_next, &
        current_next, status)
    call record_condition(status == 0, &
        "resistive-wall midpoint accepts a positive RL step")
    call evaluate_resistive_wall_energy_balance( &
        inductance, resistance, step, current, current_next, voltage_n, &
        voltage_next, energy_n, energy_next, input_work, dissipation, balance, status)
    call record_condition(status == 0 .and. dissipation >= 0.0_dp .and. &
        abs(balance) < 2.0e-13_dp, &
        "midpoint wall update satisfies the discrete energy/passivity ledger")

    call advance_resistive_wall_midpoint_jvp( &
        inductance, resistance, step, current, voltage_n, voltage_next, &
        inductance_dot, resistance_dot, step_dot, current_dot, voltage_n_dot, &
        voltage_next_dot, current_next_dot, status)
    call advance_resistive_wall_midpoint( &
        inductance + derivative_step*inductance_dot, &
        resistance + derivative_step*resistance_dot, step + derivative_step*step_dot, &
        current + derivative_step*current_dot, &
        voltage_n + derivative_step*voltage_n_dot, &
        voltage_next + derivative_step*voltage_next_dot, current_plus, status)
    call advance_resistive_wall_midpoint( &
        inductance - derivative_step*inductance_dot, &
        resistance - derivative_step*resistance_dot, step - derivative_step*step_dot, &
        current - derivative_step*current_dot, &
        voltage_n - derivative_step*voltage_n_dot, &
        voltage_next - derivative_step*voltage_next_dot, current_minus, status)
    derivative_error = maxval(abs(current_next_dot - &
        (current_plus - current_minus)/(2.0_dp*derivative_step)))
    call record_condition(status == 0 .and. derivative_error < 2.0e-8_dp, &
        "resistive-wall midpoint JVP matches central reassembly")

    current_next_bar = [0.7_dp, -0.2_dp]
    call advance_resistive_wall_midpoint_vjp( &
        inductance, resistance, step, current, voltage_n, voltage_next, &
        current_next_bar, inductance_bar, resistance_bar, step_bar, current_bar, &
        voltage_n_bar, voltage_next_bar, status)
    lhs = dot_product(current_next_bar, current_next_dot)
    rhs = sum(inductance_bar*inductance_dot) + sum(resistance_bar*resistance_dot) + &
        step_bar*step_dot + dot_product(current_bar, current_dot) + &
        dot_product(voltage_n_bar, voltage_n_dot) + &
        dot_product(voltage_next_bar, voltage_next_dot)
    adjoint_error = abs(lhs - rhs)
    call record_condition(status == 0 .and. adjoint_error < &
        5.0e-8_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "resistive-wall midpoint VJP satisfies the real adjoint")

    call check_summary("resistive-wall midpoint")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_resistive_wall_midpoint_ad
