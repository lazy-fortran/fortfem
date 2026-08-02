---
title: mixed_wave_3d_structure Example
---

# mixed_wave_3d_structure Example

# Three-dimensional mixed-wave structure

This neutral gallery fixture advances three manufactured Cartesian oscillator
components with the common first-order mixed-wave midpoint block.  The first
plot is the physical (q_x,q_y,q_z) trajectory and overlays the independent
closed-form solution.  Component and energy plots follow as diagnostics.

The example is a time-integration foundation only: it contains no plasma,
material, damping, or boundary physics.  A caller can provide compatible
pressure/velocity, displacement/momentum, Maxwell, elasticity, or tensor-
pressure blocks to the same structure-preserving API.

## Usage

```bash
fpm run --example mixed_wave_3d_structure
```

## Source Code

```fortran
program mixed_wave_3d_structure
    !! Neutral three-dimensional mixed-wave trajectory gallery fixture.
    !!
    !! The three components are independent manufactured Cartesian oscillator
    !! modes.  The midpoint update is the ideal, energy-preserving block; no
    !! damping or absorbing term is hidden in this example.
    use fortfem_api, only: advance_mixed_wave_midpoint
    use fortfem_kinds, only: dp
    use fortplot, only: add_3d_plot, figure, legend, plot, savefig, &
        set_yscale, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: component_count = 3, step_count = 240
    real(dp), parameter :: time_step = 0.02_dp
    real(dp), parameter :: frequencies(component_count) = [ &
        0.75_dp, 1.25_dp, 1.80_dp]
    real(dp), parameter :: initial_q(component_count) = [ &
        0.80_dp, -0.40_dp, 0.20_dp]
    real(dp), parameter :: initial_v(component_count) = [ &
        -0.30_dp, 0.50_dp, 0.70_dp]
    real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    real(dp), parameter :: green(3) = [0.0_dp, 158.0_dp, 115.0_dp]/255.0_dp
    character(*), parameter :: output_directory = &
        "output/example/mixed_wave_3d_structure"
    real(dp) :: mass_q(component_count, component_count)
    real(dp) :: mass_v(component_count, component_count)
    real(dp) :: coupling(component_count, component_count)
    real(dp) :: q(component_count), v(component_count)
    real(dp) :: q_exact(component_count), v_exact(component_count)
    real(dp) :: time(step_count + 1), energy_history(step_count + 1)
    real(dp) :: error_history(step_count + 1)
    real(dp) :: q_history(component_count, step_count + 1)
    real(dp) :: q_exact_history(component_count, step_count + 1)
    real(dp) :: energy_initial, energy, maximum_error, maximum_energy_drift
    real(dp) :: start_time, elapsed
    integer :: command_status, component, step, unit
    type(fortsparse_status_t) :: status

    call execute_command_line( &
        "mkdir -p "//output_directory, exitstat=command_status)
    if (command_status /= 0) error stop "cannot create 3D wave output directory"
    call initialize_gallery_sequence()

    mass_q = 0.0_dp
    mass_v = 0.0_dp
    coupling = 0.0_dp
    do component = 1, component_count
        mass_q(component, component) = 1.0_dp
        mass_v(component, component) = 1.0_dp
        coupling(component, component) = frequencies(component)
    end do

    q = initial_q
    v = initial_v
    energy_initial = wave_energy(q, v)
    time(1) = 0.0_dp
    q_history(:, 1) = q
    q_exact_history(:, 1) = q
    energy_history(1) = energy_initial
    error_history(1) = 0.0_dp

    call cpu_time(start_time)
    do step = 1, step_count
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling, time_step, q, v, status)
        if (status%code /= 0) error stop "3D mixed-wave midpoint failed"
        time(step + 1) = real(step, dp)*time_step
        q_history(:, step + 1) = q
        energy = wave_energy(q, v)
        energy_history(step + 1) = energy
        call exact_state(time(step + 1), q_exact, v_exact)
        q_exact_history(:, step + 1) = q_exact
        error_history(step + 1) = maxval(abs(q - q_exact))
    end do
    call cpu_time(elapsed)
    elapsed = elapsed - start_time

    maximum_error = maxval(error_history)
    maximum_energy_drift = maxval(abs(energy_history/energy_initial - 1.0_dp))
    if (maximum_error > 2.0e-3_dp) &
        error stop "3D mixed-wave analytical solution oracle failed"
    if (maximum_energy_drift > 3.0e-12_dp) &
        error stop "3D mixed-wave energy invariant failed"

    call render_trajectory()
    call record_gallery_stage("physical_solution")
    call render_diagnostics()
    call record_gallery_stage("diagnostics")

    open (newunit=unit, file=output_directory//"/benchmark.csv", &
        status="replace", action="write")
    write (unit, "(a)") "quantity,value"
    write (unit, "(a,es24.16)") "time_step,", time_step
    write (unit, "(a,i0)") "step_count,", step_count
    write (unit, "(a,es24.16)") "maximum_exact_state_error,", maximum_error
    write (unit, "(a,es24.16)") "maximum_relative_energy_drift,", &
        maximum_energy_drift
    write (unit, "(a,es24.16)") "midpoint_wall_time_seconds,", elapsed
    write (unit, "(a)") "primary_plot,mixed_wave_3d_trajectory_3d.png"
    close (unit)

contains

    subroutine initialize_gallery_sequence()
        integer :: sequence_unit

        open (newunit=sequence_unit, &
            file=output_directory//"/gallery_sequence.txt", &
            status="replace", action="write")
        close (sequence_unit)
    end subroutine initialize_gallery_sequence

    subroutine record_gallery_stage(stage)
        character(*), intent(in) :: stage
        integer :: sequence_unit

        open (newunit=sequence_unit, &
            file=output_directory//"/gallery_sequence.txt", &
            status="old", position="append", action="write")
        write (sequence_unit, "(a)") stage
        close (sequence_unit)
    end subroutine record_gallery_stage

    subroutine render_trajectory()
        call figure(figsize=[8.0_dp, 6.5_dp])
        call add_3d_plot( &
            q_history(1, :), q_history(2, :), q_history(3, :), &
            label="midpoint Cartesian trajectory", color=blue, linewidth=2.0_dp)
        call add_3d_plot( &
            q_exact_history(1, :), q_exact_history(2, :), &
            q_exact_history(3, :), label="analytical trajectory", &
            color=orange, linestyle="--", linewidth=1.4_dp)
        call title("3D mixed-wave physical solution trajectory")
        call legend()
        call savefig(output_directory//"/mixed_wave_3d_trajectory_3d.png")
    end subroutine render_trajectory

    subroutine render_diagnostics()
        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(time, q_history(1, :), label="q_x", color=blue)
        call plot(time, q_history(2, :), label="q_y", color=orange)
        call plot(time, q_history(3, :), label="q_z", color=green)
        call xlabel("time [s]")
        call ylabel("Cartesian displacement-like state")
        call title("3D mixed-wave component solution")
        call legend()
        call savefig(output_directory//"/mixed_wave_3d_components_1d.png")

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(time, max(abs(energy_history/energy_initial - 1.0_dp), &
            1.0e-16_dp), label="absolute relative energy drift", color=green)
        call xlabel("time [s]")
        call ylabel("abs(E(t)/E(0) - 1)")
        call title("3D mixed-wave energy invariant")
        call set_yscale("log")
        call legend()
        call savefig(output_directory//"/mixed_wave_3d_energy_1d.png")
    end subroutine render_diagnostics

    pure subroutine exact_state(current_time, q_value, v_value)
        real(dp), intent(in) :: current_time
        real(dp), intent(out) :: q_value(:), v_value(:)
        integer :: local_component

        do local_component = 1, component_count
            q_value(local_component) = initial_q(local_component)*cos( &
                frequencies(local_component)*current_time) - &
                initial_v(local_component)*sin( &
                frequencies(local_component)*current_time)
            v_value(local_component) = initial_v(local_component)*cos( &
                frequencies(local_component)*current_time) + &
                initial_q(local_component)*sin( &
                frequencies(local_component)*current_time)
        end do
    end subroutine exact_state

    pure function wave_energy(q_value, v_value) result(energy)
        real(dp), intent(in) :: q_value(:), v_value(:)
        real(dp) :: energy

        energy = 0.5_dp*(dot_product(q_value, q_value) + &
            dot_product(v_value, v_value))
    end function wave_energy

end program mixed_wave_3d_structure
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/mixed_wave_3d_structure/primary.png)

### mixed_wave_3d_components_1d.png

![mixed_wave_3d_components_1d.png](../../media/examples/mixed_wave_3d_structure/mixed_wave_3d_components_1d.png)

### mixed_wave_3d_energy_1d.png

![mixed_wave_3d_energy_1d.png](../../media/examples/mixed_wave_3d_structure/mixed_wave_3d_energy_1d.png)

### mixed_wave_3d_trajectory_3d.png

![mixed_wave_3d_trajectory_3d.png](../../media/examples/mixed_wave_3d_structure/mixed_wave_3d_trajectory_3d.png)

---

[← Back to all examples](../index.html)
