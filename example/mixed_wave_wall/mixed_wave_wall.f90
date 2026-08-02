program mixed_wave_wall
    !! Manufactured mixed wave / resistive-wall port example.
    use fortfem_time, only: &
        advance_mixed_wave_wall_midpoint, evaluate_mixed_wave_wall_energy_balance
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: steps = 160
    real(dp), parameter :: h = 0.025_dp
    character(*), parameter :: output_directory = "output/example/mixed_wave_wall"
    real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    real(dp), parameter :: green(3) = [0.0_dp, 158.0_dp, 115.0_dp]/255.0_dp
    real(dp) :: mass_q(2, 2), mass_v(2, 2), inductance(1, 1), resistance(1, 1)
    real(dp) :: coupling(2, 2), port(1, 2)
    real(dp) :: q(2), v(2), current(1), q_initial(2), v_initial(2), current_initial(1)
    real(dp) :: time(steps + 1), q_history(2, steps + 1), v_history(2, steps + 1)
    real(dp) :: current_history(steps + 1), energy(steps + 1), dissipation(steps + 1)
    real(dp) :: energy_n, energy_next, step_dissipation, balance
    integer :: step
    type(fortsparse_status_t) :: status

    call execute_command_line("mkdir -p "//output_directory)
    mass_q = reshape([2.0_dp, 0.1_dp, 0.1_dp, 1.5_dp], shape(mass_q))
    mass_v = reshape([1.4_dp, 0.0_dp, 0.0_dp, 1.1_dp], shape(mass_v))
    inductance(1, 1) = 0.9_dp
    resistance(1, 1) = 0.35_dp
    coupling = reshape([1.0_dp, -0.2_dp, 0.3_dp, 0.8_dp], shape(coupling))
    port(1, :) = [0.6_dp, -0.4_dp]
    q_initial = [0.7_dp, -0.3_dp]
    v_initial = [-0.2_dp, 0.5_dp]
    current_initial = [0.4_dp]
    q = q_initial
    v = v_initial
    current = current_initial
    time(1) = 0.0_dp
    q_history(:, 1) = q
    v_history(:, 1) = v
    current_history(1) = current(1)
    energy(1) = 0.5_dp*(dot_product(q, matmul(mass_q, q)) + &
        dot_product(v, matmul(mass_v, v)) + dot_product(current, matmul(inductance, current)))
    dissipation(1) = 0.0_dp
    do step = 1, steps
        call advance_mixed_wave_wall_midpoint( &
            mass_q, mass_v, inductance, resistance, coupling, port, h, q, v, current, status)
        if (status%code /= 0) error stop "mixed wave-wall midpoint failed"
        time(step + 1) = real(step, dp)*h
        q_history(:, step + 1) = q
        v_history(:, step + 1) = v
        current_history(step + 1) = current(1)
        energy(step + 1) = 0.5_dp*(dot_product(q, matmul(mass_q, q)) + &
            dot_product(v, matmul(mass_v, v)) + dot_product(current, matmul(inductance, current)))
        call evaluate_mixed_wave_wall_energy_balance( &
            mass_q, mass_v, inductance, resistance, coupling, port, h, &
            q_history(:, step), v_history(:, step), [current_history(step)], &
            q, v, current, energy_n, energy_next, step_dissipation, balance, status)
        if (status%code /= 0 .or. abs(balance) > 2.0e-11_dp) &
            error stop "mixed wave-wall energy ledger failed"
        dissipation(step + 1) = dissipation(step) + step_dissipation
    end do
    call render_plots()

contains

    subroutine render_plots()
        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(time, q_history(1, :), label="wave coordinate q_1", color=blue)
        call plot(time, q_history(2, :), label="wave coordinate q_2", color=orange)
        call plot(time, v_history(1, :), label="wave velocity v_1", color=green, linestyle="--")
        call plot(time, current_history, label="wall current i", color=[0.8_dp, 0.2_dp, 0.2_dp], linestyle=":")
        call xlabel("time")
        call ylabel("state amplitude")
        call title("Coupled mixed wave--resistive wall solution")
        call legend()
        call savefig(output_directory//"/mixed_wave_wall_solution_1d.png")

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(time, energy, label="coupled energy", color=blue)
        call plot(time, dissipation, label="accumulated resistive dissipation", color=orange)
        call xlabel("time")
        call ylabel("energy / dissipated work")
        call title("Structure-preserving wave--wall energy ledger")
        call legend()
        call savefig(output_directory//"/mixed_wave_wall_energy_1d.png")
    end subroutine render_plots

end program mixed_wave_wall
