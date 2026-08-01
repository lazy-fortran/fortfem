program test_mixed_wave_strang
    use check, only: check_condition, check_summary
    use fortfem_api, only: advance_mixed_wave_strang, &
        advance_mixed_wave_strang_jvp, advance_mixed_wave_strang_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: mass_q(2, 2) = reshape([ &
        2.0_dp, 0.3_dp, 0.3_dp, 1.5_dp], [2, 2])
    real(dp), parameter :: mass_v(2, 2) = reshape([ &
        1.4_dp, -0.2_dp, -0.2_dp, 1.8_dp], [2, 2])
    real(dp), parameter :: coupling_a(2, 2) = reshape([ &
        0.7_dp, -0.3_dp, 0.4_dp, 1.1_dp], [2, 2])
    real(dp), parameter :: coupling_b(2, 2) = reshape([ &
        -0.2_dp, 0.8_dp, -0.6_dp, 0.5_dp], [2, 2])
    real(dp), parameter :: initial_q(2) = [0.7_dp, -0.4_dp]
    real(dp), parameter :: initial_v(2) = [-0.2_dp, 0.9_dp]
    real(dp), parameter :: time_step = 0.17_dp
    real(dp), parameter :: time_step_dot = -0.13_dp
    real(dp), parameter :: finite_difference_step = 1.0e-7_dp
    real(dp), parameter :: q_dot(2) = [0.2_dp, -0.5_dp]
    real(dp), parameter :: v_dot(2) = [0.4_dp, 0.1_dp]
    real(dp), parameter :: q_next_bar(2) = [0.6_dp, -0.1_dp]
    real(dp), parameter :: v_next_bar(2) = [-0.3_dp, 0.8_dp]
    real(dp) :: q(2), v(2), q_reverse(2), v_reverse(2)
    real(dp) :: q_next_dot(2), v_next_dot(2)
    real(dp) :: q_plus(2), v_plus(2), q_minus(2), v_minus(2)
    real(dp) :: q_bar(2), v_bar(2), time_step_bar
    real(dp) :: energy_initial, energy_final, vjp_left, vjp_right
    logical :: all_passed
    type(fortsparse_status_t) :: status

    all_passed = .true.
    q = initial_q
    v = initial_v
    energy_initial = wave_energy(q, v)
    call advance_mixed_wave_strang( &
        mass_q, mass_v, coupling_a, coupling_b, time_step, q, v, status)
    energy_final = wave_energy(q, v)
    call record_condition(status%code == 0, &
        "mixed Strang split accepts compatible wave blocks")
    call record_condition(abs(energy_final - energy_initial) < 2.0e-12_dp, &
        "mixed Strang split preserves the quadratic wave energy")

    q_reverse = q
    v_reverse = v
    call advance_mixed_wave_strang( &
        mass_q, mass_v, coupling_a, coupling_b, -time_step, &
        q_reverse, v_reverse, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(q_reverse - initial_q)) < 2.0e-12_dp .and. &
        maxval(abs(v_reverse - initial_v)) < 2.0e-12_dp, &
        "mixed Strang split is reversible under signed-step reversal")

    call advance_mixed_wave_strang_jvp( &
        mass_q, mass_v, coupling_a, coupling_b, time_step, initial_q, initial_v, &
        time_step_dot, q_dot, v_dot, q_next_dot, v_next_dot, status)
    call record_condition(status%code == 0, &
        "mixed Strang split JVP accepts compatible tangent inputs")

    q_plus = initial_q + finite_difference_step*q_dot
    v_plus = initial_v + finite_difference_step*v_dot
    call advance_mixed_wave_strang( &
        mass_q, mass_v, coupling_a, coupling_b, &
        time_step + finite_difference_step*time_step_dot, q_plus, v_plus, status)
    q_minus = initial_q - finite_difference_step*q_dot
    v_minus = initial_v - finite_difference_step*v_dot
    call advance_mixed_wave_strang( &
        mass_q, mass_v, coupling_a, coupling_b, &
        time_step - finite_difference_step*time_step_dot, q_minus, v_minus, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs((q_plus - q_minus)/(2.0_dp*finite_difference_step) - &
        q_next_dot)) < 3.0e-8_dp .and. &
        maxval(abs((v_plus - v_minus)/(2.0_dp*finite_difference_step) - &
        v_next_dot)) < 3.0e-8_dp, &
        "mixed Strang split JVP matches an independent central difference")

    call advance_mixed_wave_strang_vjp( &
        mass_q, mass_v, coupling_a, coupling_b, time_step, initial_q, initial_v, &
        q_next_bar, v_next_bar, q_bar, v_bar, time_step_bar, status)
    vjp_left = dot_product(q_next_bar, q_next_dot) + &
        dot_product(v_next_bar, v_next_dot)
    vjp_right = dot_product(q_bar, q_dot) + dot_product(v_bar, v_dot) + &
        time_step_bar*time_step_dot
    call record_condition(status%code == 0 .and. &
        abs(vjp_left - vjp_right) < 2.0e-12_dp, &
        "mixed Strang split VJP satisfies the real adjoint identity")

    if (.not. all_passed) error stop 1
    call check_summary("Structure-preserving mixed Strang split")

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

    pure function wave_energy(q_state, v_state) result(energy)
        real(dp), intent(in) :: q_state(:), v_state(:)
        real(dp) :: energy

        energy = 0.5_dp*(dot_product(q_state, matmul(mass_q, q_state)) + &
            dot_product(v_state, matmul(mass_v, v_state)))
    end function wave_energy

end program test_mixed_wave_strang
