---
title: nonlinear_resistive_mhd_gallery Example
---

# nonlinear_resistive_mhd_gallery Example

# Neutral nonlinear resistive-MHD island/wall gallery

This bounded manufactured fixture is the physical-first gallery seed for
MHD-14.  It uses the public eight-block nonlinear composition with a small
caller-owned state:

\[
  u=(a_{\rm island},i_{\rm wall},p') .
\]

The first plot shows a smooth manufactured island-flux profile and its
localized wall-current response.  The later plots show forward/reverse
continuation branches and the residual/input/dissipation ledger.  The
machine-readable CSV and JSON files retain the branch multiplicity,
hysteresis, path metric, residual norm, input power, dissipation, and timing.

This is deliberately not an equilibrium or resistive-MHD solver.  The
constitutive callbacks, state selection, geometry, closure, and continuation
policy remain external.  The fixture only verifies that a neutral client can
compose those callbacks and expose a physical state before diagnostics.

Outputs (generated and ignored by git):

- `island_wall_solution_1d.png`
- `island_wall_branches_1d.png`
- `island_wall_ledger_1d.png`
- `island_wall_solution.csv`
- `island_wall_continuation.csv`
- `benchmark.json`

## Usage

```bash
fpm run --example nonlinear_resistive_mhd_gallery
```

## Source Code

```fortran
program nonlinear_resistive_mhd_gallery
    !! Neutral manufactured island/wall continuation gallery fixture.
    !!
    !! The state is deliberately small and closure-neutral: island amplitude,
    !! wall current, and a pressure-like amplitude.  The callbacks exercise
    !! the public eight-block nonlinear composition without selecting a plasma
    !! model or reading an application-specific equilibrium format.
    use fortfem_api, only: &
        RESISTIVE_MHD_AMPERE, &
        RESISTIVE_MHD_FARADAY, &
        RESISTIVE_MHD_PRESSURE, &
        RESISTIVE_MHD_WALL, &
        assemble_nonlinear_resistive_mhd_residual, &
        compare_resistive_mhd_branch_histories, &
        evaluate_resistive_mhd_branch_diagnostics, &
        evaluate_resistive_mhd_branch_path_metric, &
        initialize_resistive_mhd_branch_history, &
        nonlinear_resistive_mhd_energy_ledger_t, &
        resistive_mhd_branch_diagnostics_t, &
        resistive_mhd_branch_history_t
    use fortfem_api, only: RESISTIVE_MHD_BLOCK_COUNT
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    integer, parameter :: state_size = 3, sample_count = 21
    integer, parameter :: profile_count = 161
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    real(dp), parameter :: green(3) = [0.0_dp, 158.0_dp, 115.0_dp]/255.0_dp
    character(*), parameter :: output_directory = &
        "output/example/nonlinear_resistive_mhd_gallery"

    real(dp) :: parameter_values(sample_count)
    real(dp) :: forward_state(state_size, sample_count)
    real(dp) :: reverse_state(state_size, sample_count)
    real(dp) :: forward_residual(state_size, sample_count)
    real(dp) :: reverse_residual(state_size, sample_count)
    real(dp) :: forward_energy(sample_count), reverse_energy(sample_count)
    real(dp) :: residual_norm(sample_count), input_power(sample_count)
    real(dp) :: dissipation(sample_count), balance(sample_count)
    integer :: forward_branch(sample_count), reverse_branch(sample_count)
    real(dp) :: x(profile_count), psi_profile(profile_count)
    real(dp) :: wall_profile(profile_count), state(state_size)
    real(dp) :: residual(state_size), path_metric, start_time, elapsed
    real(dp) :: current_parameter, target_state(state_size)
    real(dp) :: hysteresis, lambda
    integer :: sample, point, unit, command_status
    type(nonlinear_resistive_mhd_energy_ledger_t) :: ledger
    type(resistive_mhd_branch_history_t) :: forward_history, reverse_history
    type(resistive_mhd_branch_diagnostics_t) :: branch_diagnostics
    type(fortsparse_status_t) :: status

    command_status = 0
    call execute_command_line("mkdir -p "//output_directory, exitstat=command_status)
    if (command_status /= 0) error stop "cannot create island/wall output directory"
    call initialize_gallery_sequence()

    do sample = 1, sample_count
        parameter_values(sample) = real(sample - 1, dp)/real(sample_count - 1, dp)
        lambda = parameter_values(sample)
        hysteresis = 0.12_dp*lambda*(1.0_dp - lambda)
        forward_state(:, sample) = manufactured_state(lambda, 0.0_dp)
        reverse_state(:, sample) = manufactured_state(lambda, hysteresis)
        forward_branch(sample) = merge(2, 1, lambda >= 0.5_dp)
        reverse_branch(sample) = merge(2, 1, lambda >= 0.5_dp)

        current_parameter = lambda
        target_state = forward_state(:, sample)
        state = target_state
        call assemble_nonlinear_resistive_mhd_residual( &
            state, nonlinear_value, residual, ledger, status)
        if (status%code /= FORTSPARSE_OK) error stop "island/wall composition failed"
        forward_residual(:, sample) = residual
        forward_energy(sample) = ledger%stored_energy
        residual_norm(sample) = sqrt(dot_product(residual, residual))
        input_power(sample) = ledger%input_power
        dissipation(sample) = ledger%dissipation
        balance(sample) = ledger%balance

        current_parameter = lambda
        target_state = reverse_state(:, sample)
        state = target_state
        call assemble_nonlinear_resistive_mhd_residual( &
            state, nonlinear_value, residual, ledger, status)
        if (status%code /= FORTSPARSE_OK) &
            error stop "reverse island/wall composition failed"
        reverse_residual(:, sample) = residual
        reverse_energy(sample) = ledger%stored_energy
    end do

    call initialize_resistive_mhd_branch_history( &
        parameter_values, forward_state, forward_residual, forward_energy, &
        forward_branch, forward_history, status)
    if (status%code /= FORTSPARSE_OK) error stop "forward branch history failed"
    call initialize_resistive_mhd_branch_history( &
        parameter_values, reverse_state, reverse_residual, reverse_energy, &
        reverse_branch, reverse_history, status)
    if (status%code /= FORTSPARSE_OK) error stop "reverse branch history failed"
    call evaluate_resistive_mhd_branch_diagnostics( &
        forward_history, branch_diagnostics, status)
    if (status%code /= FORTSPARSE_OK) error stop "branch diagnostics failed"
    call compare_resistive_mhd_branch_histories( &
        forward_history, reverse_history, 1.0e-12_dp, branch_diagnostics, status)
    if (status%code /= FORTSPARSE_OK) error stop "branch comparison failed"
    call evaluate_resistive_mhd_branch_path_metric( &
        forward_history, path_metric, status)
    if (status%code /= FORTSPARSE_OK) error stop "branch path metric failed"

    do point = 1, profile_count
        x(point) = -1.0_dp + 2.0_dp*real(point - 1, dp)/real(profile_count - 1, dp)
        psi_profile(point) = forward_state(1, sample_count)*(1.0_dp - x(point)**2)
        wall_profile(point) = forward_state(2, sample_count)* &
            (exp(-((x(point) - 0.82_dp)/0.045_dp)**2) + &
            exp(-((x(point) + 0.82_dp)/0.045_dp)**2))
    end do

    call cpu_time(start_time)
    do sample = 1, 100
        current_parameter = parameter_values(sample_count)
        target_state = forward_state(:, sample_count)
        state = target_state
        call assemble_nonlinear_resistive_mhd_residual( &
            state, nonlinear_value, residual, ledger, status)
    end do
    call cpu_time(elapsed)
    elapsed = elapsed - start_time

    ! The physical solution is intentionally rendered before any diagnostics.
    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(x, psi_profile, label="manufactured island flux psi(x)", &
        color=blue, linewidth=2.4_dp)
    call plot(x, wall_profile, label="wall-current response profile", &
        color=orange, linewidth=2.0_dp)
    call xlabel("normalized radial coordinate x")
    call ylabel("manufactured field amplitude")
    call title("Neutral nonlinear island--wall continuation: physical state")
    call legend()
    call savefig(output_directory//"/island_wall_solution_1d.png")
    call record_gallery_stage("physical_solution")

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(parameter_values, forward_state(1, :), label="forward branch amplitude", &
        color=blue, linewidth=2.0_dp)
    call plot(parameter_values, reverse_state(1, :), label="reverse branch amplitude", &
        color=orange, linestyle="--", linewidth=2.0_dp)
    call xlabel("continuation parameter lambda")
    call ylabel("island amplitude")
    call title("Manufactured continuation branches and hysteresis")
    call legend()
    call savefig(output_directory//"/island_wall_branches_1d.png")
    call record_gallery_stage("branch_diagnostics")

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(parameter_values, input_power, label="signed input power", &
        color=blue, linewidth=2.0_dp)
    call plot(parameter_values, dissipation, label="nonnegative wall dissipation", &
        color=orange, linewidth=2.0_dp)
    call plot(parameter_values, residual_norm, label="residual norm", &
        color=green, linestyle="--", linewidth=2.0_dp)
    call xlabel("continuation parameter lambda")
    call ylabel("power / residual norm")
    call title("Nonlinear composition ledger")
    call legend()
    call savefig(output_directory//"/island_wall_ledger_1d.png")
    call record_gallery_stage("ledger_diagnostics")

    open (newunit=unit, file=output_directory//"/island_wall_solution.csv", &
        status="replace", action="write")
    write (unit, "(a)") "x,psi,wall_current_response"
    do point = 1, profile_count
        write (unit, "(*(es24.16,:,','))") x(point), psi_profile(point), &
            wall_profile(point)
    end do
    close (unit)

    open (newunit=unit, file=output_directory//"/island_wall_continuation.csv", &
        status="replace", action="write")
    write (unit, "(a)") "lambda,forward_amplitude,reverse_amplitude,"// &
        "residual_norm,input_power,dissipation,balance"
    do sample = 1, sample_count
        write (unit, "(*(es24.16,:,','))") parameter_values(sample), &
            forward_state(1, sample), reverse_state(1, sample), residual_norm(sample), &
            input_power(sample), dissipation(sample), balance(sample)
    end do
    close (unit)

    open (newunit=unit, file=output_directory//"/benchmark.json", &
        status="replace", action="write")
    write (unit, "(a)") "{"
    write (unit, "(a)") '  "schema": "fortfem-neutral-island-wall-v1",'
    write (unit, "(a,i0,a)") '  "branch_multiplicity": ', &
        branch_diagnostics%branch_multiplicity, ","
    if (branch_diagnostics%hysteresis_detected) then
        write (unit, "(a)") '  "hysteresis_detected": true,'
    else
        write (unit, "(a)") '  "hysteresis_detected": false,'
    end if
    write (unit, "(a,es24.16,a)") '  "max_state_hysteresis": ', &
        branch_diagnostics%max_state_hysteresis, ","
    write (unit, "(a,es24.16,a)") '  "path_metric": ', path_metric, ","
    write (unit, "(a,es24.16,a)") '  "final_residual_norm": ', &
        residual_norm(sample_count), ","
    write (unit, "(a,es24.16,a)") '  "final_input_power": ', &
        input_power(sample_count), ","
    write (unit, "(a,es24.16,a)") '  "final_dissipation": ', &
        dissipation(sample_count), ","
    write (unit, "(a,es24.16)") '  "composition_seconds_per_call": ', elapsed/100.0_dp
    write (unit, "(a)") "}"
    close (unit)

contains

    subroutine initialize_gallery_sequence()
        integer :: sequence_unit, local_status

        open (newunit=sequence_unit, file=output_directory//"/gallery_sequence.txt", &
            status="replace", action="write", iostat=local_status)
        if (local_status /= 0) error stop "cannot initialize gallery sequence"
        close (sequence_unit)
    end subroutine initialize_gallery_sequence

    subroutine record_gallery_stage(stage)
        character(*), intent(in) :: stage
        integer :: sequence_unit, local_status

        open (newunit=sequence_unit, file=output_directory//"/gallery_sequence.txt", &
            status="old", position="append", action="write", iostat=local_status)
        if (local_status /= 0) error stop "cannot open gallery sequence"
        write (sequence_unit, "(a)", iostat=local_status) stage
        close (sequence_unit)
        if (local_status /= 0) error stop "cannot record gallery sequence"
    end subroutine record_gallery_stage

    pure function manufactured_state(parameter_value, offset) result(state_out)
        real(dp), intent(in) :: parameter_value, offset
        real(dp) :: state_out(state_size)

        state_out(1) = 0.18_dp + 0.24_dp*parameter_value + &
            0.04_dp*parameter_value**2 + offset
        state_out(2) = 0.65_dp*state_out(1) + 0.08_dp*parameter_value
        state_out(3) = 0.05_dp*parameter_value
    end function manufactured_state

    subroutine nonlinear_value(state_in, block_id, block_residual, stored_energy, &
            input_power_out, dissipation_out, local_status)
        real(dp), intent(in) :: state_in(:)
        integer, intent(in) :: block_id
        real(dp), intent(out) :: block_residual(:)
        real(dp), intent(out) :: stored_energy, input_power_out, dissipation_out
        integer, intent(out) :: local_status

        block_residual = 0.0_dp
        stored_energy = 0.0_dp
        input_power_out = 0.0_dp
        dissipation_out = 0.0_dp
        local_status = 1
        if (size(state_in) /= state_size .or. size(block_residual) /= state_size .or. &
            block_id < 1 .or. block_id > RESISTIVE_MHD_BLOCK_COUNT) return
        select case (block_id)
        case (RESISTIVE_MHD_FARADAY)
            block_residual(1) = state_in(1) - target_state(1)
            stored_energy = 0.5_dp*state_in(1)**2
        case (RESISTIVE_MHD_AMPERE)
            block_residual(2) = state_in(2) - target_state(2)
            stored_energy = 0.5_dp*state_in(2)**2
        case (RESISTIVE_MHD_PRESSURE)
            block_residual(3) = state_in(3) - target_state(3)
            stored_energy = 0.5_dp*state_in(3)**2
        case (RESISTIVE_MHD_WALL)
            input_power_out = 0.35_dp*current_parameter*state_in(2)
            dissipation_out = 0.18_dp*state_in(2)**2
        end select
        local_status = 0
    end subroutine nonlinear_value

end program nonlinear_resistive_mhd_gallery
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/nonlinear_resistive_mhd_gallery/primary.png)

### island_wall_branches_1d.png

![island_wall_branches_1d.png](../../media/examples/nonlinear_resistive_mhd_gallery/island_wall_branches_1d.png)

### island_wall_ledger_1d.png

![island_wall_ledger_1d.png](../../media/examples/nonlinear_resistive_mhd_gallery/island_wall_ledger_1d.png)

### island_wall_solution_1d.png

![island_wall_solution_1d.png](../../media/examples/nonlinear_resistive_mhd_gallery/island_wall_solution_1d.png)

---

[← Back to all examples](../index.html)
