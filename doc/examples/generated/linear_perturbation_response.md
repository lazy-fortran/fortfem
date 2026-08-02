---
title: linear_perturbation_response Example
---

# linear_perturbation_response Example

# Manufactured linear perturbation response

This closure-neutral fixture composes the seven public linear-perturbation
blocks

\[
A(\omega)=L+P+V+W+S-\omega^2M+i\omega R
\]

and solves `A state = source` for five coupled Fourier coefficients under an
asymmetric complex drive. The first plot reconstructs the complex state over
physical poloidal angle. It is a
manufactured displacement-like response, not a plasma equilibrium, stability
model, singular-layer closure, or plasma-code file reader.

The frequency sweep forms the dense response matrix `A^-1` and passes it to
FortFEM's public reciprocity/passivity diagnostic. It also compares the exact
frequency JVP of the seven-block composition with an independent centered
difference and checks the forced residual after each solve.

Outputs:

- `linear_response_state_1d.png`: physical complex state first, with magnitude
  and real/imaginary components over poloidal angle;
- `linear_response_phase_1d.png`: unwrapped-free principal phase of that state;
- `linear_response_diagnostics_1d.png`: reciprocity, passivity, derivative,
  and forced-residual diagnostics across frequency;
- `linear_response_state.csv` and `linear_response_diagnostics.csv`:
  reproducible plot sources;
- `benchmark.txt`: small runner-local composition, solve, and diagnostic timing
  record.

The timings are regression context for this tiny manufactured algebraic case,
not production performance claims.

## Usage

```bash
fpm run --example linear_perturbation_response
```

## Source Code

```fortran
program linear_perturbation_response
    !! Physical-first manufactured gallery for seven-block linear response.
    use fortfem_interop, only: &
        assemble_linear_perturbation_operator, &
        assemble_linear_perturbation_operator_jvp, &
        assemble_linear_response_residual, &
        evaluate_linear_response_diagnostics
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    use fortplot, only: figure, legend, plot, savefig, set_yscale, &
        title, xlabel, ylabel
    implicit none

    integer, parameter :: mode_count = 5
    integer, parameter :: sample_count = 241
    integer, parameter :: frequency_count = 17
    integer, parameter :: timing_repetitions = 500
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: difference_step = 1.0e-6_dp
    real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    real(dp), parameter :: black(3) = [0.12_dp, 0.12_dp, 0.12_dp]
    real(dp), parameter :: green(3) = [0.0_dp, 158.0_dp, 115.0_dp]/255.0_dp
    character(*), parameter :: output_directory = &
        "output/example/linear_perturbation_response"
    integer, parameter :: modes(mode_count) = [-2, -1, 0, 1, 2]
    complex(dp), parameter :: imaginary_unit = cmplx(0.0_dp, 1.0_dp, dp)
    complex(dp) :: lorentz(mode_count, mode_count)
    complex(dp) :: pressure(mode_count, mode_count)
    complex(dp) :: inertia(mode_count, mode_count)
    complex(dp) :: vacuum(mode_count, mode_count)
    complex(dp) :: wall(mode_count, mode_count)
    complex(dp) :: resistive(mode_count, mode_count)
    complex(dp) :: singular(mode_count, mode_count)
    complex(dp) :: operator(mode_count, mode_count)
    complex(dp) :: operator_dot(mode_count, mode_count)
    complex(dp) :: operator_plus(mode_count, mode_count)
    complex(dp) :: operator_minus(mode_count, mode_count)
    complex(dp) :: response(mode_count, mode_count)
    complex(dp) :: identity(mode_count, mode_count)
    complex(dp) :: source(mode_count), state(mode_count)
    complex(dp) :: residual(mode_count), field
    complex(dp) :: frequency, frequency_dot
    real(dp) :: theta(sample_count), field_real(sample_count)
    real(dp) :: field_imaginary(sample_count), field_magnitude(sample_count)
    real(dp) :: field_phase(sample_count), frequencies(frequency_count)
    real(dp) :: reciprocity_error(frequency_count)
    real(dp) :: passivity_margin(frequency_count)
    real(dp) :: derivative_error(frequency_count)
    real(dp) :: residual_error(frequency_count)
    real(dp) :: compose_seconds, solve_seconds, diagnostic_seconds
    real(dp) :: start_time, normalization, timing_passivity, timing_reciprocity
    integer :: frequency_index, info, mode, repetition, sample, status, unit

    call execute_command_line("mkdir -p "//output_directory)
    call build_blocks()
    call build_source()
    frequency = cmplx(0.72_dp, 0.0_dp, dp)
    call assemble_linear_perturbation_operator( &
        lorentz, pressure, inertia, vacuum, wall, resistive, singular, &
        frequency, operator, status)
    if (status /= 0) error stop "seven-block operator composition failed"
    call dense_solve(operator, source, state, info)
    if (info /= 0) error stop "manufactured linear response solve failed"
    call assemble_linear_response_residual( &
        operator, state, source, residual, status)
    if (status /= 0 .or. maxval(abs(residual)) > 2.0e-13_dp) &
        error stop "manufactured forced residual failed"

    do sample = 1, sample_count
        theta(sample) = 2.0_dp*pi*real(sample - 1, dp)/ &
            real(sample_count - 1, dp)
        field = cmplx(0.0_dp, 0.0_dp, dp)
        do mode = 1, mode_count
            field = field + state(mode)*exp( &
                imaginary_unit*real(modes(mode), dp)*theta(sample))
        end do
        field_real(sample) = real(field, dp)
        field_imaginary(sample) = aimag(field)
        field_magnitude(sample) = abs(field)
        field_phase(sample) = atan2(aimag(field), real(field, dp))
    end do
    if (maxval(field_magnitude) <= 1.0e-2_dp) &
        error stop "manufactured physical response is trivial"

    identity = cmplx(0.0_dp, 0.0_dp, dp)
    do mode = 1, mode_count
        identity(mode, mode) = cmplx(1.0_dp, 0.0_dp, dp)
    end do
    frequency_dot = cmplx(1.0_dp, 0.0_dp, dp)
    do frequency_index = 1, frequency_count
        frequencies(frequency_index) = 0.25_dp + &
            0.75_dp*real(frequency_index - 1, dp)/ &
            real(frequency_count - 1, dp)
        frequency = cmplx(frequencies(frequency_index), 0.0_dp, dp)
        call assemble_linear_perturbation_operator( &
            lorentz, pressure, inertia, vacuum, wall, resistive, singular, &
            frequency, operator, status)
        if (status /= 0) error stop "frequency-sweep composition failed"
        call dense_solve(operator, identity, response, info)
        if (info /= 0) error stop "frequency-sweep response solve failed"
        call evaluate_linear_response_diagnostics( &
            response, reciprocity_error(frequency_index), &
            passivity_margin(frequency_index), status)
        if (status /= 0) error stop "response diagnostic failed"
        call assemble_linear_perturbation_operator_jvp( &
            lorentz, pressure, inertia, vacuum, wall, resistive, singular, &
            frequency, 0.0_dp*lorentz, 0.0_dp*pressure, 0.0_dp*inertia, &
            0.0_dp*vacuum, 0.0_dp*wall, 0.0_dp*resistive, &
            0.0_dp*singular, frequency_dot, operator_dot, status)
        if (status /= 0) error stop "frequency JVP failed"
        call assemble_linear_perturbation_operator( &
            lorentz, pressure, inertia, vacuum, wall, resistive, singular, &
            frequency + difference_step, operator_plus, status)
        call assemble_linear_perturbation_operator( &
            lorentz, pressure, inertia, vacuum, wall, resistive, singular, &
            frequency - difference_step, operator_minus, status)
        normalization = max(1.0_dp, maxval(abs(operator_dot)))
        derivative_error(frequency_index) = maxval(abs( &
            operator_dot - (operator_plus - operator_minus)/ &
            (2.0_dp*difference_step)))/normalization
        call dense_solve(operator, source, state, info)
        if (info /= 0) error stop "frequency-sweep state solve failed"
        call assemble_linear_response_residual( &
            operator, state, source, residual, status)
        residual_error(frequency_index) = maxval(abs(residual))/ &
            max(1.0_dp, maxval(abs(source)))
    end do
    if (maxval(reciprocity_error) > 5.0e-14_dp) &
        error stop "transpose reciprocity diagnostic regressed"
    if (minval(passivity_margin) <= 0.0_dp) &
        error stop "passivity certificate is nonpositive"
    if (maxval(derivative_error) > 2.0e-9_dp) &
        error stop "frequency JVP central-difference diagnostic failed"

    call measure_timings()
    call write_outputs()
    call make_plots()

contains

    subroutine build_blocks()
        integer :: column, distance, row

        lorentz = cmplx(0.0_dp, 0.0_dp, dp)
        pressure = cmplx(0.0_dp, 0.0_dp, dp)
        inertia = cmplx(0.0_dp, 0.0_dp, dp)
        vacuum = cmplx(0.0_dp, 0.0_dp, dp)
        wall = cmplx(0.0_dp, 0.0_dp, dp)
        resistive = cmplx(0.0_dp, 0.0_dp, dp)
        singular = cmplx(0.0_dp, 0.0_dp, dp)
        do mode = 1, mode_count
            lorentz(mode, mode) = 3.5_dp + &
                0.22_dp*real(modes(mode)*modes(mode), dp)
            pressure(mode, mode) = 0.45_dp + &
                0.04_dp*real(abs(modes(mode)), dp)
            inertia(mode, mode) = 0.85_dp + &
                0.08_dp*real(abs(modes(mode)), dp)
            vacuum(mode, mode) = 0.24_dp
            wall(mode, mode) = 0.11_dp
            resistive(mode, mode) = 0.18_dp + &
                0.025_dp*real(abs(modes(mode)), dp)
            singular(mode, mode) = 0.035_dp + &
                0.045_dp*merge(1.0_dp, 0.0_dp, modes(mode) == 0)
        end do
        do column = 1, mode_count
            do row = column + 1, mode_count
                distance = abs(row - column)
                lorentz(row, column) = -0.16_dp/real(distance, dp)
                pressure(row, column) = 0.035_dp/real(distance, dp)
                vacuum(row, column) = 0.025_dp/real(distance, dp)
                wall(row, column) = 0.010_dp/real(distance, dp)
                singular(row, column) = 0.008_dp/real(distance, dp)
                lorentz(column, row) = lorentz(row, column)
                pressure(column, row) = pressure(row, column)
                vacuum(column, row) = vacuum(row, column)
                wall(column, row) = wall(row, column)
                singular(column, row) = singular(row, column)
            end do
        end do
    end subroutine build_blocks

    subroutine build_source()
        source = [ &
            cmplx(0.08_dp, -0.04_dp, dp), &
            cmplx(0.25_dp, -0.18_dp, dp), &
            cmplx(1.00_dp, 0.15_dp, dp), &
            cmplx(0.12_dp, 0.36_dp, dp), &
            cmplx(-0.05_dp, 0.12_dp, dp)]
    end subroutine build_source

    subroutine measure_timings()
        call cpu_time(start_time)
        do repetition = 1, timing_repetitions
            call assemble_linear_perturbation_operator( &
                lorentz, pressure, inertia, vacuum, wall, resistive, singular, &
                frequency, operator, status)
        end do
        if (status /= 0) error stop "timed operator composition failed"
        call cpu_time(compose_seconds)
        compose_seconds = (compose_seconds - start_time)/ &
            real(timing_repetitions, dp)

        call cpu_time(start_time)
        do repetition = 1, timing_repetitions
            call dense_solve(operator, source, state, info)
        end do
        if (info /= 0) error stop "timed response solve failed"
        call cpu_time(solve_seconds)
        solve_seconds = (solve_seconds - start_time)/ &
            real(timing_repetitions, dp)

        call dense_solve(operator, identity, response, info)
        call cpu_time(start_time)
        do repetition = 1, timing_repetitions
            call evaluate_linear_response_diagnostics( &
                response, timing_reciprocity, timing_passivity, status)
        end do
        if (status /= 0) error stop "timed response diagnostic failed"
        call cpu_time(diagnostic_seconds)
        diagnostic_seconds = (diagnostic_seconds - start_time)/ &
            real(timing_repetitions, dp)
    end subroutine measure_timings

    subroutine write_outputs()
        open (newunit=unit, file=output_directory// &
            "/linear_response_state.csv", status="replace", action="write")
        write (unit, "(a)") "theta_rad,real_state,imag_state,magnitude,phase_rad"
        do sample = 1, sample_count
            write (unit, "(5(es24.16,1x))") theta(sample), &
                field_real(sample), field_imaginary(sample), &
                field_magnitude(sample), field_phase(sample)
        end do
        close (unit)

        open (newunit=unit, file=output_directory// &
            "/linear_response_diagnostics.csv", &
            status="replace", action="write")
        write (unit, "(a)") &
            "frequency,reciprocity_error,passivity_margin,jvp_error,residual_error"
        do frequency_index = 1, frequency_count
            write (unit, "(5(es24.16,1x))") frequencies(frequency_index), &
                reciprocity_error(frequency_index), &
                passivity_margin(frequency_index), &
                derivative_error(frequency_index), &
                residual_error(frequency_index)
        end do
        close (unit)

        open (newunit=unit, file=output_directory//"/benchmark.txt", &
            status="replace", action="write")
        write (unit, "(a,i0)") "state coefficients: ", mode_count
        write (unit, "(a,i0)") "timing repetitions: ", timing_repetitions
        write (unit, "(a,es16.8)") &
            "seven-block composition seconds: ", compose_seconds
        write (unit, "(a,es16.8)") "dense response solve seconds: ", solve_seconds
        write (unit, "(a,es16.8)") &
            "response diagnostic seconds: ", diagnostic_seconds
        write (unit, "(a,es16.8)") "maximum reciprocity error: ", &
            maxval(reciprocity_error)
        write (unit, "(a,es16.8)") "minimum passivity margin: ", &
            minval(passivity_margin)
        write (unit, "(a,es16.8)") "maximum frequency JVP error: ", &
            maxval(derivative_error)
        write (unit, "(a,es16.8)") "maximum forced residual error: ", &
            maxval(residual_error)
        close (unit)
    end subroutine write_outputs

    subroutine make_plots()
        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(theta, field_magnitude, label="magnitude", &
            color=black, linewidth=2.2_dp)
        call plot(theta, field_real, label="real component", &
            color=blue, linestyle="-")
        call plot(theta, field_imaginary, label="imaginary component", &
            color=orange, linestyle="--")
        call xlabel("poloidal angle theta (rad)")
        call ylabel("manufactured response amplitude")
        call title("Physical complex response from seven perturbation blocks")
        call legend()
        call savefig(output_directory//"/linear_response_state_1d.png")

        call figure(figsize=[9.0_dp, 4.8_dp])
        call plot(theta, field_phase, color=blue, linewidth=2.0_dp)
        call xlabel("poloidal angle theta (rad)")
        call ylabel("principal phase (rad)")
        call title("Phase of the manufactured physical response")
        call savefig(output_directory//"/linear_response_phase_1d.png")

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(frequencies, max(reciprocity_error, epsilon(1.0_dp)), &
            label="transpose reciprocity defect", marker="o", color=blue)
        call plot(frequencies, passivity_margin, &
            label="certified passivity margin", marker="s", color=green)
        call plot(frequencies, max(derivative_error, epsilon(1.0_dp)), &
            label="frequency JVP error", marker="^", color=orange, &
            linestyle="--")
        call plot(frequencies, max(residual_error, epsilon(1.0_dp)), &
            label="forced residual", marker="x", color=black, linestyle=":")
        call set_yscale("log")
        call xlabel("dimensionless angular frequency")
        call ylabel("dimensionless diagnostic")
        call title("Reciprocity, passivity, derivative, and residual checks")
        call legend()
        call savefig(output_directory//"/linear_response_diagnostics_1d.png")
    end subroutine make_plots

end program linear_perturbation_response
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/linear_perturbation_response/primary.png)

### linear_response_diagnostics_1d.png

![linear_response_diagnostics_1d.png](../../media/examples/linear_perturbation_response/linear_response_diagnostics_1d.png)

### linear_response_phase_1d.png

![linear_response_phase_1d.png](../../media/examples/linear_perturbation_response/linear_response_phase_1d.png)

### linear_response_state_1d.png

![linear_response_state_1d.png](../../media/examples/linear_perturbation_response/linear_response_state_1d.png)

---

[← Back to all examples](../index.html)
