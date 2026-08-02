---
title: singular_layer_matching_gallery Example
---

# singular_layer_matching_gallery Example

# Neutral singular-layer trace matching

This bounded MHD-08/MHD-09 foundation gallery samples analytical complex inner
and outer traces on a one-dimensional layer coordinate.  The first figure
shows their real profiles and the matched jump prescribed to
`assemble_singular_layer_matching`; complex components are retained in the
CSV source.  A second figure deliberately perturbs that jump and reports the
mismatch residual together with the exact JVP's centered-difference error.

The trace rows, states, jump, and positive weights are caller-owned data.  No
plasma closure, singular asymptotic, equilibrium reader, or production tearing
physics is selected.  Outputs are written to
`output/example/singular_layer_matching_gallery/`:

- `singular_layer_solution_1d.png`: physical outer/inner traces and matched jump;
- `singular_layer_diagnostics_1d.png`: mismatch and derivative diagnostics;
- `singular_layer.csv`: reproducible complex profile and residual samples;
- `provenance.json`: contract, analytical parameters, and diagnostic norms.

## Usage

```bash
fpm run --example singular_layer_matching_gallery
```

## Source Code

```fortran
program singular_layer_matching_gallery
    !! Physical-first neutral inner/outer trace matching gallery.
    !!
    !! The one-dimensional profiles are analytical caller data.  The public
    !! singular-layer contract is applied pointwise to a scalar trace row;
    !! no singular asymptotic, plasma closure, or tearing model is selected.
    use fortfem_api, only: &
        assemble_singular_layer_matching, &
        assemble_singular_layer_matching_jvp
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, set_yscale, title, &
        xlabel, ylabel
    implicit none

    integer, parameter :: sample_count = 161
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: difference_step = 1.0e-6_dp
    real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    real(dp), parameter :: green(3) = [0.0_dp, 158.0_dp, 115.0_dp]/255.0_dp
    character(*), parameter :: output_directory = &
        "output/example/singular_layer_matching_gallery"

    real(dp) :: coordinate(sample_count), outer_operator(sample_count)
    real(dp) :: inner_operator(sample_count), weights(sample_count)
    real(dp) :: outer_operator_dot(sample_count)
    real(dp) :: inner_operator_dot(sample_count), weights_dot(sample_count)
    complex(dp) :: outer_profile(sample_count), inner_profile(sample_count)
    complex(dp) :: matched_jump(sample_count), mismatched_jump(sample_count)
    complex(dp) :: mismatch(sample_count), residual(sample_count)
    complex(dp) :: mismatch_residual(sample_count), residual_jvp(sample_count)
    real(dp) :: jvp_central_error(sample_count)
    complex(dp), parameter :: outer_state = cmplx(0.80_dp, 0.35_dp, dp)
    complex(dp), parameter :: inner_state = cmplx(0.55_dp, -0.20_dp, dp)
    complex(dp), parameter :: outer_state_dot = cmplx(0.03_dp, -0.015_dp, dp)
    complex(dp), parameter :: inner_state_dot = cmplx(-0.02_dp, 0.01_dp, dp)
    integer :: sample, status, status_plus, status_minus, command_status
    real(dp) :: coordinate_value, residual_linf, mismatch_residual_linf
    real(dp) :: jvp_error_linf
    complex(dp) :: matched_jump_dot, mismatch_dot, jump_dot
    complex(dp) :: residual_plus, residual_minus

    command_status = 0
    call execute_command_line("mkdir -p "//output_directory, &
        exitstat=command_status)
    if (command_status /= 0) error stop "cannot create singular-layer output directory"
    call initialize_gallery_sequence()

    do sample = 1, sample_count
        coordinate_value = -1.0_dp + 2.0_dp*real(sample - 1, dp)/ &
            real(sample_count - 1, dp)
        coordinate(sample) = coordinate_value
        outer_operator(sample) = 1.0_dp + 0.15_dp*coordinate_value
        inner_operator(sample) = 0.90_dp - 0.10_dp*coordinate_value
        weights(sample) = 1.0_dp + 0.20_dp*cos(pi*coordinate_value)
        outer_operator_dot(sample) = 0.03_dp*cos(pi*coordinate_value)
        inner_operator_dot(sample) = -0.02_dp*sin(pi*coordinate_value)
        weights_dot(sample) = 0.04_dp*sin(pi*coordinate_value)

        outer_profile(sample) = outer_operator(sample)*outer_state
        inner_profile(sample) = inner_operator(sample)*inner_state
        matched_jump(sample) = outer_profile(sample) - inner_profile(sample)
        mismatch(sample) = 0.025_dp*exp(-4.0_dp*coordinate_value**2)* &
            cmplx(1.0_dp, -0.35_dp, dp)
        mismatched_jump(sample) = matched_jump(sample) + mismatch(sample)

        call evaluate_residual(outer_operator(sample), inner_operator(sample), &
            weights(sample), outer_state, inner_state, matched_jump(sample), &
            residual(sample), status)
        if (status /= 0) error stop "matched singular-layer residual failed"
        call evaluate_residual(outer_operator(sample), inner_operator(sample), &
            weights(sample), outer_state, inner_state, mismatched_jump(sample), &
            mismatch_residual(sample), status)
        if (status /= 0) error stop "mismatched singular-layer residual failed"

        matched_jump_dot = outer_operator_dot(sample)*outer_state + &
            outer_operator(sample)*outer_state_dot - &
            inner_operator_dot(sample)*inner_state - &
            inner_operator(sample)*inner_state_dot
        mismatch_dot = 0.011_dp*sin(pi*coordinate_value)* &
            cmplx(1.0_dp, -0.35_dp, dp)
        jump_dot = matched_jump_dot + mismatch_dot
        call evaluate_jvp(outer_operator(sample), inner_operator(sample), &
            weights(sample), outer_state, inner_state, mismatched_jump(sample), &
            outer_operator_dot(sample), inner_operator_dot(sample), &
            weights_dot(sample), outer_state_dot, inner_state_dot, jump_dot, &
            residual_jvp(sample), status)
        if (status /= 0) error stop "singular-layer JVP failed"

        call evaluate_residual( &
            outer_operator(sample) + difference_step*outer_operator_dot(sample), &
            inner_operator(sample) + difference_step*inner_operator_dot(sample), &
            weights(sample) + difference_step*weights_dot(sample), &
            outer_state + difference_step*outer_state_dot, &
            inner_state + difference_step*inner_state_dot, &
            mismatched_jump(sample) + difference_step*jump_dot, residual_plus, &
            status_plus)
        call evaluate_residual( &
            outer_operator(sample) - difference_step*outer_operator_dot(sample), &
            inner_operator(sample) - difference_step*inner_operator_dot(sample), &
            weights(sample) - difference_step*weights_dot(sample), &
            outer_state - difference_step*outer_state_dot, &
            inner_state - difference_step*inner_state_dot, &
            mismatched_jump(sample) - difference_step*jump_dot, residual_minus, &
            status_minus)
        if (status_plus /= 0 .or. status_minus /= 0) &
            error stop "singular-layer central difference failed"
        jvp_central_error(sample) = abs(residual_jvp(sample) - &
            (residual_plus - residual_minus)/(2.0_dp*difference_step))/ &
            max(1.0_dp, abs(residual_jvp(sample)))
    end do

    residual_linf = maxval(abs(residual))
    mismatch_residual_linf = maxval(abs(mismatch_residual))
    jvp_error_linf = maxval(jvp_central_error)
    if (residual_linf > 2.0e-13_dp) &
        error stop "matched trace residual is not zero"
    if (mismatch_residual_linf < 1.0e-3_dp) &
        error stop "mismatch diagnostic is trivial"
    if (jvp_error_linf > 2.0e-9_dp) &
        error stop "singular-layer JVP central-difference check failed"

    ! Physical state first: caller-owned outer/inner traces and their jump.
    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(coordinate, real(outer_profile, dp), label="outer trace (real)", &
        color=blue, linewidth=2.4_dp)
    call plot(coordinate, real(inner_profile, dp), label="inner trace (real)", &
        color=orange, linestyle="--", linewidth=2.2_dp)
    call plot(coordinate, real(matched_jump, dp), label="matched jump (real)", &
        color=green, linestyle=":", linewidth=2.6_dp)
    call xlabel("normalized layer coordinate s")
    call ylabel("trace amplitude (real part)")
    call title("Neutral inner/outer trace matching: physical profiles")
    call legend()
    call savefig(output_directory//"/singular_layer_solution_1d.png")
    call record_gallery_stage("physical_solution")

    ! Diagnostics second: deliberate jump mismatch and exact JVP parity.
    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(coordinate, abs(mismatch_residual), &
        label="deliberate mismatch residual", color=blue, linewidth=2.4_dp)
    call plot(coordinate, max(1.0e-16_dp, jvp_central_error), &
        label="JVP central-difference error", color=orange, linestyle="--", &
        linewidth=2.2_dp)
    call set_yscale("log")
    call xlabel("normalized layer coordinate s")
    call ylabel("diagnostic magnitude")
    call title("Singular-layer mismatch and derivative diagnostics")
    call legend()
    call savefig(output_directory//"/singular_layer_diagnostics_1d.png")
    call record_gallery_stage("diagnostics")

    call write_field_csv()
    call write_provenance_json()

contains

    subroutine evaluate_residual( &
            outer_operator_value, inner_operator_value, weight_value, &
            outer_state_value, inner_state_value, jump_value, residual_value, &
            status_value)
        real(dp), intent(in) :: outer_operator_value, inner_operator_value
        real(dp), intent(in) :: weight_value
        complex(dp), intent(in) :: outer_state_value, inner_state_value, jump_value
        complex(dp), intent(out) :: residual_value
        integer, intent(out) :: status_value
        complex(dp) :: outer_trace(1, 1), inner_trace(1, 1)
        complex(dp) :: outer_state_vector(1), inner_state_vector(1), jump_vector(1)
        complex(dp) :: residual_vector(1)
        real(dp) :: weight_vector(1)

        outer_trace(1, 1) = cmplx(outer_operator_value, 0.0_dp, dp)
        inner_trace(1, 1) = cmplx(inner_operator_value, 0.0_dp, dp)
        outer_state_vector(1) = outer_state_value
        inner_state_vector(1) = inner_state_value
        jump_vector(1) = jump_value
        weight_vector(1) = weight_value
        call assemble_singular_layer_matching( &
            outer_trace, inner_trace, weight_vector, outer_state_vector, &
            inner_state_vector, jump_vector, residual_vector, status_value)
        residual_value = residual_vector(1)
    end subroutine evaluate_residual

    subroutine evaluate_jvp( &
            outer_operator_value, inner_operator_value, weight_value, &
            outer_state_value, inner_state_value, jump_value, &
            outer_operator_dot_value, inner_operator_dot_value, weight_dot_value, &
            outer_state_dot_value, inner_state_dot_value, jump_dot_value, &
            residual_dot_value, status_value)
        real(dp), intent(in) :: outer_operator_value, inner_operator_value
        real(dp), intent(in) :: weight_value, outer_operator_dot_value
        real(dp), intent(in) :: inner_operator_dot_value, weight_dot_value
        complex(dp), intent(in) :: outer_state_value, inner_state_value, jump_value
        complex(dp), intent(in) :: outer_state_dot_value, inner_state_dot_value
        complex(dp), intent(in) :: jump_dot_value
        complex(dp), intent(out) :: residual_dot_value
        integer, intent(out) :: status_value
        complex(dp) :: outer_trace(1, 1), inner_trace(1, 1)
        complex(dp) :: outer_trace_dot(1, 1), inner_trace_dot(1, 1)
        complex(dp) :: outer_state_vector(1), inner_state_vector(1), jump_vector(1)
        complex(dp) :: outer_state_dot_vector(1), inner_state_dot_vector(1)
        complex(dp) :: jump_dot_vector(1), residual_dot_vector(1)
        real(dp) :: weight_vector(1), weight_dot_vector(1)

        outer_trace(1, 1) = cmplx(outer_operator_value, 0.0_dp, dp)
        inner_trace(1, 1) = cmplx(inner_operator_value, 0.0_dp, dp)
        outer_trace_dot(1, 1) = cmplx(outer_operator_dot_value, 0.0_dp, dp)
        inner_trace_dot(1, 1) = cmplx(inner_operator_dot_value, 0.0_dp, dp)
        outer_state_vector(1) = outer_state_value
        inner_state_vector(1) = inner_state_value
        jump_vector(1) = jump_value
        outer_state_dot_vector(1) = outer_state_dot_value
        inner_state_dot_vector(1) = inner_state_dot_value
        jump_dot_vector(1) = jump_dot_value
        weight_vector(1) = weight_value
        weight_dot_vector(1) = weight_dot_value
        call assemble_singular_layer_matching_jvp( &
            outer_trace, inner_trace, weight_vector, outer_state_vector, &
            inner_state_vector, jump_vector, outer_trace_dot, inner_trace_dot, &
            weight_dot_vector, outer_state_dot_vector, inner_state_dot_vector, &
            jump_dot_vector, residual_dot_vector, status_value)
        residual_dot_value = residual_dot_vector(1)
    end subroutine evaluate_jvp

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

    subroutine write_field_csv()
        integer :: csv_unit, sample_index

        open (newunit=csv_unit, file=output_directory//"/singular_layer.csv", &
            status="replace", action="write")
        write (csv_unit, "(a)") &
            "coordinate,outer_real,outer_imag,inner_real,inner_imag,"// &
            "jump_real,jump_imag,mismatch_real,mismatch_imag,weight,"// &
            "mismatch_residual_real,mismatch_residual_imag,residual_real,"// &
            "residual_imag,residual_jvp_real,residual_jvp_imag,jvp_error"
        do sample_index = 1, sample_count
            write (csv_unit, "(*(es24.16,:,','))") coordinate(sample_index), &
                real(outer_profile(sample_index), dp), &
                aimag(outer_profile(sample_index)), &
                real(inner_profile(sample_index), dp), &
                aimag(inner_profile(sample_index)), &
                real(matched_jump(sample_index), dp), &
                aimag(matched_jump(sample_index)), &
                real(mismatch(sample_index), dp), aimag(mismatch(sample_index)), &
                weights(sample_index), &
                real(mismatch_residual(sample_index), dp), &
                aimag(mismatch_residual(sample_index)), &
                real(residual(sample_index), dp), aimag(residual(sample_index)), &
                real(residual_jvp(sample_index), dp), &
                aimag(residual_jvp(sample_index)), jvp_central_error(sample_index)
        end do
        close (csv_unit)
    end subroutine write_field_csv

    subroutine write_provenance_json()
        integer :: json_unit

        open (newunit=json_unit, file=output_directory//"/provenance.json", &
            status="replace", action="write")
        write (json_unit, "(a)") "{"
        write (json_unit, "(a)") &
            '  "schema": "fortfem-singular-layer-gallery-v1",'
        write (json_unit, "(a)") &
            '  "contract": "singular_layer_matching",'
        write (json_unit, "(a)") &
            '  "closure": "neutral-caller-owned-traces",'
        write (json_unit, "(a)") &
            '  "primary_plot": "singular_layer_solution_1d.png",'
        write (json_unit, "(a,i0,a)") '  "sample_count": ', sample_count, ","
        write (json_unit, "(a,es24.16,a)") '  "matched_residual_linf": ', &
            residual_linf, ","
        write (json_unit, "(a,es24.16,a)") '  "mismatch_residual_linf": ', &
            mismatch_residual_linf, ","
        write (json_unit, "(a,es24.16,a)") '  "jvp_central_error": ', &
            jvp_error_linf, ","
        write (json_unit, "(a,es24.16,a)") '  "outer_operator_slope": ', 0.15_dp, ","
        write (json_unit, "(a,es24.16,a)") '  "inner_operator_slope": ', -0.10_dp, ","
        write (json_unit, "(a,es24.16)") '  "mismatch_amplitude": ', 0.025_dp
        write (json_unit, "(a)") "}"
        close (json_unit)
    end subroutine write_provenance_json

end program singular_layer_matching_gallery
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/singular_layer_matching_gallery/primary.png)

### singular_layer_diagnostics_1d.png

![singular_layer_diagnostics_1d.png](../../media/examples/singular_layer_matching_gallery/singular_layer_diagnostics_1d.png)

### singular_layer_solution_1d.png

![singular_layer_solution_1d.png](../../media/examples/singular_layer_matching_gallery/singular_layer_solution_1d.png)

---

[← Back to all examples](../index.html)
