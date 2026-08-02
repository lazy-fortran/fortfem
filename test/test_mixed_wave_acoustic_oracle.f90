program test_mixed_wave_acoustic_oracle
    use check, only: check_condition, check_summary
    use fortfem_time, only: advance_mixed_wave_midpoint
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: mode_count = 2, step_count = 400
    real(dp), parameter :: time_step = 0.01_dp
    real(dp), parameter :: frequencies(mode_count) = [1.0_dp, 1.7_dp]
    real(dp), parameter :: initial_q(mode_count) = [1.0_dp, -0.4_dp]
    real(dp), parameter :: initial_v(mode_count) = [0.2_dp, 0.7_dp]
    real(dp) :: mass_q(mode_count, mode_count), mass_v(mode_count, mode_count)
    real(dp) :: coupling(mode_count, mode_count)
    real(dp) :: q(mode_count), v(mode_count), q_exact(mode_count)
    real(dp) :: v_exact(mode_count), energy_initial, energy
    real(dp) :: maximum_error, maximum_energy_drift, time
    real(dp) :: q_reversed(mode_count), v_reversed(mode_count)
    integer :: step, mode
    type(fortsparse_status_t) :: status

    mass_q = 0.0_dp
    mass_v = 0.0_dp
    coupling = 0.0_dp
    do mode = 1, mode_count
        mass_q(mode, mode) = 1.0_dp
        mass_v(mode, mode) = 1.0_dp
        coupling(mode, mode) = frequencies(mode)
    end do

    q = initial_q
    v = initial_v
    energy_initial = 0.5_dp*(dot_product(q, q) + dot_product(v, v))
    maximum_error = 0.0_dp
    maximum_energy_drift = 0.0_dp
    do step = 1, step_count
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling, time_step, q, v, status)
        call check_condition(status%code == 0, &
            "mixed acoustic midpoint accepts compatible modal blocks")
        time = real(step, dp)*time_step
        do mode = 1, mode_count
            q_exact(mode) = initial_q(mode)*cos(frequencies(mode)*time) - &
                initial_v(mode)*sin(frequencies(mode)*time)
            v_exact(mode) = initial_v(mode)*cos(frequencies(mode)*time) + &
                initial_q(mode)*sin(frequencies(mode)*time)
        end do
        maximum_error = max(maximum_error, maxval(abs(q - q_exact)))
        maximum_error = max(maximum_error, maxval(abs(v - v_exact)))
        energy = 0.5_dp*(dot_product(q, q) + dot_product(v, v))
        maximum_energy_drift = max(maximum_energy_drift, &
            abs(energy - energy_initial))
    end do

    call check_condition(maximum_error < 5.0e-4_dp, &
        "mixed midpoint follows the independent acoustic modal solution")
    call check_condition(maximum_energy_drift < 2.0e-12_dp, &
        "mixed midpoint preserves acoustic energy over many steps")

    q_reversed = q
    v_reversed = v
    do step = 1, step_count
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling, -time_step, q_reversed, v_reversed, &
            status)
    end do
    call check_condition(maxval(abs(q_reversed - initial_q)) < 2.0e-12_dp .and. &
        maxval(abs(v_reversed - initial_v)) < 2.0e-12_dp, &
        "mixed midpoint reverses a multi-mode acoustic trajectory")
    call check_summary("Mixed acoustic midpoint analytical oracle")
end program test_mixed_wave_acoustic_oracle
