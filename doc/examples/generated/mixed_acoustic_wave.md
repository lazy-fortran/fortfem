---
title: mixed_acoustic_wave Example
---

# mixed_acoustic_wave Example

# Mixed acoustic midpoint wave

This fixture is a small physical wave problem for the common mixed
first-order time-step API. Two independent modal pairs satisfy

\[
\dot q=-C^T v,\qquad \dot v=Cq,
\]

with unit mass matrices and modal frequencies one and 1.7 radians per
second. The implicit-midpoint step is the Cayley map of this
port-Hamiltonian system. It preserves the quadratic energy and is exactly
reversible under a signed time-step reversal.

The first plot shows the propagated modal displacement-like and
velocity-like variables. The phase-space and energy plots are secondary
structure diagnostics. The program compares the numerical trajectory with
the independent closed-form oscillator solution and records the wall time in
`benchmark.txt`.

This is a generic acoustic/wave foundation fixture. Pressure, displacement,
electromagnetic, and elasticity clients can supply their own compatible mass
and interconnection blocks without changing the time integrator.

## Usage

```bash
fpm run --example mixed_acoustic_wave
```

## Source Code

```fortran
program mixed_acoustic_wave
    !! Structure-preserving mixed acoustic modal wave.
    !!
    !! The two modal pairs satisfy q_dot = -C^T v and v_dot = C q.
    !! Implicit midpoint is used as the energy-preserving Cayley map.
    use fortfem_api, only: advance_mixed_wave_midpoint
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, set_yscale, title, &
        xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: mode_count = 2, step_count = 400
    real(dp), parameter :: time_step = 0.01_dp
    real(dp), parameter :: frequencies(mode_count) = [1.0_dp, 1.7_dp]
    real(dp), parameter :: initial_q(mode_count) = [1.0_dp, -0.4_dp]
    real(dp), parameter :: initial_v(mode_count) = [0.2_dp, 0.7_dp]
    real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    real(dp), parameter :: green(3) = [0.0_dp, 158.0_dp, 115.0_dp]/255.0_dp
    real(dp), parameter :: purple(3) = [204.0_dp, 121.0_dp, 167.0_dp]/255.0_dp
    character(*), parameter :: output_directory = &
        "output/example/mixed_acoustic_wave"
    real(dp) :: mass_q(mode_count, mode_count), mass_v(mode_count, mode_count)
    real(dp) :: coupling(mode_count, mode_count)
    real(dp) :: q(mode_count), v(mode_count), q_reverse(mode_count)
    real(dp) :: v_reverse(mode_count), energy_initial, energy
    real(dp) :: time(step_count + 1), q_history(mode_count, step_count + 1)
    real(dp) :: v_history(mode_count, step_count + 1)
    real(dp) :: energy_drift(step_count + 1), exact_error(step_count + 1)
    real(dp) :: q_exact(mode_count), v_exact(mode_count)
    real(dp) :: maximum_error, maximum_energy_drift, start_time, elapsed
    integer :: mode, step, unit
    type(fortsparse_status_t) :: status

    call execute_command_line("mkdir -p "//output_directory)
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
    time(1) = 0.0_dp
    q_history(:, 1) = q
    v_history(:, 1) = v
    energy_drift(1) = 0.0_dp
    exact_error(1) = 0.0_dp

    call cpu_time(start_time)
    do step = 1, step_count
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling, time_step, q, v, status)
        if (status%code /= 0) error stop "mixed acoustic midpoint failed"
        time(step + 1) = real(step, dp)*time_step
        q_history(:, step + 1) = q
        v_history(:, step + 1) = v
        energy = 0.5_dp*(dot_product(q, q) + dot_product(v, v))
        energy_drift(step + 1) = energy/energy_initial - 1.0_dp
        call exact_modal_state(time(step + 1), q_exact, v_exact)
        exact_error(step + 1) = max(maxval(abs(q - q_exact)), &
            maxval(abs(v - v_exact)))
    end do
    call cpu_time(elapsed)
    elapsed = elapsed - start_time

    maximum_error = maxval(exact_error)
    maximum_energy_drift = maxval(abs(energy_drift))
    q_reverse = q
    v_reverse = v
    do step = 1, step_count
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling, -time_step, q_reverse, v_reverse, &
            status)
        if (status%code /= 0) error stop "mixed acoustic reverse step failed"
    end do
    if (maxval(abs(q_reverse - initial_q)) > 2.0e-12_dp .or. &
        maxval(abs(v_reverse - initial_v)) > 2.0e-12_dp) &
        error stop "mixed acoustic reversibility oracle failed"
    if (maximum_error > 5.0e-4_dp) &
        error stop "mixed acoustic analytical solution oracle failed"
    if (maximum_energy_drift > 2.0e-12_dp) &
        error stop "mixed acoustic energy invariant failed"

    call render_plots()
    open (newunit=unit, file=output_directory//"/benchmark.txt", &
        status="replace", action="write")
    write (unit, "(a)") "quantity,value"
    write (unit, "(a,es24.16)") "time_step,", time_step
    write (unit, "(a,i0)") "step_count,", step_count
    write (unit, "(a,es24.16)") "maximum_exact_state_error,", maximum_error
    write (unit, "(a,es24.16)") "maximum_relative_energy_drift,", &
        maximum_energy_drift
    write (unit, "(a,es24.16)") "midpoint_wall_time_seconds,", elapsed
    close (unit)

contains

    subroutine exact_modal_state(current_time, q_value, v_value)
        real(dp), intent(in) :: current_time
        real(dp), intent(out) :: q_value(:), v_value(:)
        integer :: local_mode

        do local_mode = 1, mode_count
            q_value(local_mode) = initial_q(local_mode)*cos( &
                frequencies(local_mode)*current_time) - &
                initial_v(local_mode)*sin( &
                frequencies(local_mode)*current_time)
            v_value(local_mode) = initial_v(local_mode)*cos( &
                frequencies(local_mode)*current_time) + &
                initial_q(local_mode)*sin( &
                frequencies(local_mode)*current_time)
        end do
    end subroutine exact_modal_state

    subroutine render_plots()
        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(time, q_history(1, :), label="q mode 1", &
            color=blue, linestyle="-")
        call plot(time, v_history(1, :), label="v mode 1", &
            color=blue, linestyle="--")
        call plot(time, q_history(2, :), label="q mode 2", &
            color=orange, linestyle="-")
        call plot(time, v_history(2, :), label="v mode 2", &
            color=orange, linestyle="--")
        call xlabel("time [s]")
        call ylabel("modal amplitude")
        call title("Mixed acoustic solution: midpoint propagation")
        call legend()
        call savefig(output_directory//"/mixed_acoustic_solution_1d.png")

        call figure(figsize=[7.0_dp, 6.0_dp])
        call plot(q_history(1, :), v_history(1, :), label="mode 1", &
            color=blue)
        call plot(q_history(2, :), v_history(2, :), label="mode 2", &
            color=orange)
        call xlabel("q")
        call ylabel("v")
        call title("Mixed acoustic phase-space trajectories")
        call legend()
        call savefig(output_directory//"/mixed_acoustic_phase_2d.png")

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(time, energy_drift, label="relative energy drift", &
            color=green)
        call xlabel("time [s]")
        call ylabel("E(t)/E(0) - 1")
        call title("Midpoint energy invariant")
        call savefig(output_directory//"/mixed_acoustic_energy_1d.png")

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(time, max(exact_error, tiny(1.0_dp)), &
            label="maximum state error", color=purple)
        call set_yscale("log")
        call xlabel("time [s]")
        call ylabel("absolute error")
        call title("Mixed acoustic analytical comparison")
        call savefig(output_directory//"/mixed_acoustic_exact_error_1d.png")
    end subroutine render_plots

end program mixed_acoustic_wave
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/mixed_acoustic_wave/primary.png)

### mixed_acoustic_energy_1d.png

![mixed_acoustic_energy_1d.png](../../media/examples/mixed_acoustic_wave/mixed_acoustic_energy_1d.png)

### mixed_acoustic_exact_error_1d.png

![mixed_acoustic_exact_error_1d.png](../../media/examples/mixed_acoustic_wave/mixed_acoustic_exact_error_1d.png)

### mixed_acoustic_phase_2d.png

![mixed_acoustic_phase_2d.png](../../media/examples/mixed_acoustic_wave/mixed_acoustic_phase_2d.png)

### mixed_acoustic_solution_1d.png

![mixed_acoustic_solution_1d.png](../../media/examples/mixed_acoustic_wave/mixed_acoustic_solution_1d.png)

---

[← Back to all examples](../index.html)
