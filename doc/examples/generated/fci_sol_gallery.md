---
title: fci_sol_gallery Example
---

# fci_sol_gallery Example

---
title: fci_sol_gallery Example
---

This neutral FCI/SOL fixture traces three helical field lines on a toroidal
surface and applies the same fixed-topology maps to FortFEM's conservative
parallel gradient and diffusion support operators. The first image is the
physical three-dimensional field-line solution; the second is a two-dimensional
field-coordinate-plane solution. The diagnostic plot and JSON/CSV ledgers come
after those solution views.

The example deliberately does not choose a plasma species, sheath closure,
transport model, equilibrium reader, or edge-physics code. It provides a small
geometry and operator contract that later SOL applications can consume. The
support pairing follows the PARALLAX field-coordinate-independent construction
described by [Stegmeir et al.](https://doi.org/10.1016/j.cpc.2015.09.016).

Outputs include:

- `fci_sol_field_lines_3d.png`: toroidal surface and traced field lines (the
  physical primary plot);
- `fci_sol_solution_2d.png`: manufactured scalar on the FCI plane grid;
- `fci_sol_gradient_diagnostic_1d.png`: support-gradient diagnostic against the
  analytical line derivative;
- `fci_sol_solution.csv` and `fci_sol_field_lines.csv`: reproducible values;
- `benchmark.json`: timing, conservation, closure, and provenance metadata.

Run it with:

```bash
fpm run --example fci_sol_gallery
```

## Usage

```bash
fpm run --example fci_sol_gallery
```

## Source Code

```fortran
program fci_sol_gallery
    !! Physical-first FCI/SOL gallery fixture.
    !!
    !! This is a deliberately neutral support-operator example: a prescribed
    !! helical field-line family is traced on a torus, while the same geometry
    !! drives a conservative FCI parallel gradient and diffusion action.  No
    !! species, sheath, closure, or plasma-equilibrium model is selected here.
    !! The example is an executable geometry/algebra contract for later SOL
    !! applications, not a production edge solver.
    use fortfem_feec, only: &
        apply_fci_parallel_diffusion, apply_fci_parallel_gradient, &
        trace_fci_field_line_rk4
    use fortfem_kinds, only: dp
    use fortplot, only: &
        add_3d_plot, add_parametric_surface, colorbar, figure, legend, &
        pcolormesh, plot, savefig, title, xlabel, ylabel
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    integer, parameter :: theta_cells = 24
    integer, parameter :: phi_cells = 48
    integer, parameter :: line_count = 3
    integer, parameter :: trace_count = 128
    integer, parameter :: benchmark_repetitions = 64
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: major_radius = 2.4_dp
    real(dp), parameter :: minor_radius = 0.62_dp
    real(dp), parameter :: safety_slope = 0.42_dp
    real(dp), parameter :: parallel_coefficient_value = 1.0_dp
    real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    real(dp), parameter :: green(3) = [0.0_dp, 158.0_dp, 115.0_dp]/255.0_dp
    character(*), parameter :: output_directory = &
        "output/example/fci_sol_gallery"

    real(dp) :: theta_nodes(theta_cells), phi_planes(phi_cells + 1)
    real(dp) :: theta_edges(theta_cells + 1), phi_edges(phi_cells + 1)
    real(dp) :: forward_map(theta_cells, theta_cells, phi_cells)
    real(dp) :: backward_map(theta_cells, theta_cells, phi_cells)
    real(dp) :: line_lengths(theta_cells, phi_cells)
    real(dp) :: parallel_coefficient(theta_cells*phi_cells)
    real(dp) :: canonical_volumes(theta_cells*(phi_cells + 1))
    real(dp) :: staggered_volumes(theta_cells*phi_cells)
    real(dp) :: field(theta_cells*(phi_cells + 1))
    real(dp) :: gradient_field(theta_cells*phi_cells)
    real(dp) :: diffusion_field(theta_cells*(phi_cells + 1))
    real(dp) :: solution_map(phi_cells, theta_cells)
    real(dp) :: gradient_map(phi_cells, theta_cells)
    real(dp) :: line_x(line_count, trace_count + 1)
    real(dp) :: line_y(line_count, trace_count + 1)
    real(dp) :: line_z(line_count, trace_count + 1)
    real(dp) :: surface_x(theta_cells + 1, 37)
    real(dp) :: surface_y(theta_cells + 1, 37)
    real(dp) :: surface_z(theta_cells + 1, 37)
    real(dp) :: endpoint(2), angle, radial, toroidal_angle
    real(dp) :: target_theta, alpha, dtheta, dphi
    real(dp) :: total_mass_rate, dissipation, maximum_gradient
    real(dp) :: action_seconds
    integer :: i, j, k, left, line, sample, unit, command_status
    integer :: repetition, clock_rate, start_count, end_count
    type(fortsparse_status_t) :: status

    call execute_command_line("mkdir -p "//output_directory, &
        exitstat=command_status)
    if (command_status /= 0) error stop "cannot create FCI/SOL output directory"

    dtheta = 2.0_dp*pi/real(theta_cells, dp)
    dphi = 2.0_dp*pi/real(phi_cells, dp)
    do i = 1, theta_cells
        theta_nodes(i) = dtheta*real(i - 1, dp)
        theta_edges(i) = dtheta*real(i - 1, dp)
    end do
    theta_edges(theta_cells + 1) = 2.0_dp*pi
    do j = 1, phi_cells + 1
        phi_planes(j) = dphi*real(j - 1, dp)
        phi_edges(j) = phi_planes(j)
    end do

    call build_support_maps()
    call evaluate_solution()
    call trace_physical_field_lines()
    call render_physical_solution()
    call write_gallery_stage("physical_solution")
    call render_diagnostics()
    call write_gallery_stage("diagnostics")
    call write_outputs()

contains

    subroutine build_support_maps()
        integer :: segment, row, node

        forward_map = 0.0_dp
        backward_map = 0.0_dp
        line_lengths = dphi
        parallel_coefficient = parallel_coefficient_value
        canonical_volumes = 1.0_dp
        staggered_volumes = 1.0_dp
        do segment = 1, phi_cells
            do row = 1, theta_cells
                backward_map(row, row, segment) = 1.0_dp
                target_theta = modulo(theta_nodes(row) + safety_slope*dphi, 2.0_dp*pi)
                left = 1 + int(target_theta/dtheta)
                left = min(theta_cells, max(1, left))
                alpha = (target_theta - theta_nodes(left))/dtheta
                node = 1 + mod(left, theta_cells)
                forward_map(row, left, segment) = 1.0_dp - alpha
                forward_map(row, node, segment) = &
                    forward_map(row, node, segment) + alpha
            end do
        end do
    end subroutine build_support_maps

    subroutine evaluate_solution()
        integer :: segment, row, plane_offset, gradient_offset

        field = 0.0_dp
        do segment = 1, phi_cells + 1
            plane_offset = (segment - 1)*theta_cells
            do row = 1, theta_cells
                field(plane_offset + row) = manufactured_field( &
                    theta_nodes(row), phi_planes(segment))
            end do
        end do

        call apply_fci_parallel_gradient( &
            forward_map, backward_map, line_lengths, field, gradient_field, &
            status)
        if (status%code /= FORTSPARSE_OK) &
            error stop "FCI/SOL gradient action failed"
        call apply_fci_parallel_diffusion( &
            forward_map, backward_map, line_lengths, parallel_coefficient, &
            canonical_volumes, staggered_volumes, field, diffusion_field, &
            status)
        if (status%code /= FORTSPARSE_OK) &
            error stop "FCI/SOL diffusion action failed"

        do segment = 1, phi_cells
            gradient_offset = (segment - 1)*theta_cells
            do row = 1, theta_cells
                solution_map(segment, row) = manufactured_field( &
                    theta_nodes(row) + 0.5_dp*dtheta, &
                    phi_planes(segment) + 0.5_dp*dphi)
                gradient_map(segment, row) = gradient_field(gradient_offset + row)
            end do
        end do
        total_mass_rate = sum(diffusion_field*canonical_volumes)
        dissipation = dot_product(field*canonical_volumes, diffusion_field)
        maximum_gradient = maxval(abs(gradient_field))
    end subroutine evaluate_solution

    subroutine trace_physical_field_lines()
        real(dp) :: initial_point(2)

        do line = 1, line_count
            initial_point = [ &
                2.0_dp*pi*real(line - 1, dp)/real(line_count, dp), &
                minor_radius*(0.45_dp + 0.18_dp*real(line - 1, dp))]
            do sample = 0, trace_count
                toroidal_angle = 2.0_dp*pi*real(sample, dp)/ &
                    real(trace_count, dp)
                call trace_fci_field_line_rk4( &
                    initial_point, toroidal_angle, 4, toroidal_field_line_rhs, &
                    endpoint, status)
                if (status%code /= FORTSPARSE_OK) &
                    error stop "FCI/SOL field-line tracing failed"
                radial = endpoint(2)
                angle = endpoint(1)
                line_x(line, sample + 1) = (major_radius + radial*cos(angle))* &
                    cos(toroidal_angle)
                line_y(line, sample + 1) = (major_radius + radial*cos(angle))* &
                    sin(toroidal_angle)
                line_z(line, sample + 1) = radial*sin(angle)
            end do
        end do

        do j = 1, 37
            toroidal_angle = 2.0_dp*pi*real(j - 1, dp)/36.0_dp
            do i = 1, theta_cells + 1
                angle = 2.0_dp*pi*real(i - 1, dp)/real(theta_cells, dp)
                radial = minor_radius
                surface_x(i, j) = (major_radius + radial*cos(angle))* &
                    cos(toroidal_angle)
                surface_y(i, j) = (major_radius + radial*cos(angle))* &
                    sin(toroidal_angle)
                surface_z(i, j) = radial*sin(angle)
            end do
        end do
    end subroutine trace_physical_field_lines

    subroutine render_physical_solution()
        call figure(figsize=[8.5_dp, 6.5_dp])
        call add_parametric_surface( &
            surface_x, surface_y, surface_z, color="lightsteelblue", &
            alpha=0.30_dp, linewidth=0.2_dp, filled=.true., &
            label="SOL boundary")
        call add_3d_plot(line_x(1, :), line_y(1, :), line_z(1, :), &
            label="field line 1", color=blue, linewidth=2.5_dp)
        call add_3d_plot(line_x(2, :), line_y(2, :), line_z(2, :), &
            label="field line 2", color=orange, linewidth=2.5_dp)
        call add_3d_plot(line_x(3, :), line_y(3, :), line_z(3, :), &
            label="field line 3", color=green, linewidth=2.5_dp)
        call title("FCI/SOL helical field-line solution on a torus")
        call legend()
        call savefig(output_directory//"/fci_sol_field_lines_3d.png")

        call figure(figsize=[8.5_dp, 5.5_dp])
        call pcolormesh(theta_edges, phi_edges, solution_map, cmap="coolwarm")
        call colorbar(label="manufactured SOL potential u")
        call xlabel("poloidal angle theta [rad]")
        call ylabel("toroidal angle phi [rad]")
        call title("FCI/SOL physical solution on field-coordinate planes")
        call savefig(output_directory//"/fci_sol_solution_2d.png")
    end subroutine render_physical_solution

    subroutine render_diagnostics()
        real(dp) :: phi_center(phi_cells), selected_gradient(phi_cells)

        do j = 1, phi_cells
            phi_center(j) = 0.5_dp*(phi_edges(j) + phi_edges(j + 1))
            selected_gradient(j) = gradient_map(j, 1)
        end do
        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(phi_center, selected_gradient, label="FCI parallel gradient", &
            color=blue, linewidth=2.0_dp)
        call plot(phi_center, -0.5_dp*sin(2.0_dp*phi_center), &
            label="analytic line derivative", color=orange, linestyle="--", &
            linewidth=1.6_dp)
        call xlabel("toroidal angle phi [rad]")
        call ylabel("Q u")
        call title("FCI/SOL support-gradient diagnostic")
        call legend()
        call savefig(output_directory//"/fci_sol_gradient_diagnostic_1d.png")
    end subroutine render_diagnostics

    subroutine write_gallery_stage(stage)
        character(*), intent(in) :: stage
        integer :: sequence_unit

        if (stage == "physical_solution") then
            open (newunit=sequence_unit, &
                file=output_directory//"/gallery_sequence.txt", &
                status="replace", action="write")
        else
            open (newunit=sequence_unit, &
                file=output_directory//"/gallery_sequence.txt", &
                status="old", position="append", action="write")
        end if
        write (sequence_unit, "(a)") stage
        close (sequence_unit)
    end subroutine write_gallery_stage

    subroutine write_outputs()
        integer :: segment, row, sample

        open (newunit=unit, file=output_directory//"/fci_sol_solution.csv", &
            status="replace", action="write")
        write (unit, "(a)") "theta,phi,solution,parallel_gradient"
        do segment = 1, phi_cells
            do row = 1, theta_cells
                write (unit, "(4(es24.16,:,','))") &
                    theta_nodes(row) + 0.5_dp*dtheta, &
                    phi_planes(segment) + 0.5_dp*dphi, solution_map(segment, row), &
                    gradient_map(segment, row)
            end do
        end do
        close (unit)

        open (newunit=unit, file=output_directory//"/fci_sol_field_lines.csv", &
            status="replace", action="write")
        write (unit, "(a)") "line,step,x,y,z"
        do line = 1, line_count
            do sample = 1, trace_count + 1
                write (unit, "(2(i0,','),3(es24.16,:,','))") line, sample - 1, &
                    line_x(line, sample), line_y(line, sample), line_z(line, sample)
            end do
        end do
        close (unit)

        call system_clock(count_rate=clock_rate)
        action_seconds = -1.0_dp
        if (clock_rate > 0) then
            call system_clock(start_count)
            do repetition = 1, benchmark_repetitions
                call apply_fci_parallel_gradient( &
                    forward_map, backward_map, line_lengths, field, &
                    gradient_field, status)
                if (status%code /= FORTSPARSE_OK) &
                    error stop "FCI/SOL benchmark gradient failed"
            end do
            call system_clock(end_count)
            action_seconds = real(end_count - start_count, dp)/ &
                real(clock_rate, dp)/real(benchmark_repetitions, dp)
        end if

        open (newunit=unit, file=output_directory//"/benchmark.json", &
            status="replace", action="write")
        write (unit, "(a)") "{"
        write (unit, "(a,i0,a)") '  "theta_cells": ', theta_cells, ","
        write (unit, "(a,i0,a)") '  "phi_cells": ', phi_cells, ","
        write (unit, "(a,i0,a)") '  "field_line_count": ', line_count, ","
        write (unit, "(a,es24.16,a)") '  "mass_rate": ', total_mass_rate, ","
        write (unit, "(a,es24.16,a)") '  "dissipation": ', dissipation, ","
        write (unit, "(a,es24.16,a)") '  "maximum_parallel_gradient": ', &
            maximum_gradient, ","
        write (unit, "(a,es24.16,a)") '  "gradient_action_seconds": ', &
            action_seconds, ","
        write (unit, "(a)") '  "provenance": "analytic-torus-fci-sol-v1",'
        write (unit, "(a)") &
            '  "primary_plot": "fci_sol_field_lines_3d.png",'
        write (unit, "(a)") '  "closure": "neutral-support-operator"'
        write (unit, "(a)") "}"
        close (unit)
    end subroutine write_outputs

    pure function manufactured_field(theta_value, phi_value) result(value)
        real(dp), intent(in) :: theta_value, phi_value
        real(dp) :: value

        value = sin(theta_value - safety_slope*phi_value) + &
            0.35_dp*cos(2.0_dp*phi_value)
    end function manufactured_field

    pure subroutine toroidal_field_line_rhs(phi_value, point_value, derivative)
        real(dp), intent(in) :: phi_value
        real(dp), intent(in) :: point_value(:)
        real(dp), intent(out) :: derivative(:)

        associate (unused_phi => phi_value, unused_point => point_value)
            derivative(1) = safety_slope
            derivative(2) = 0.0_dp
        end associate
    end subroutine toroidal_field_line_rhs

end program fci_sol_gallery
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/fci_sol_gallery/primary.png)

### fci_sol_field_lines_3d.png

![fci_sol_field_lines_3d.png](../../media/examples/fci_sol_gallery/fci_sol_field_lines_3d.png)

### fci_sol_gradient_diagnostic_1d.png

![fci_sol_gradient_diagnostic_1d.png](../../media/examples/fci_sol_gallery/fci_sol_gradient_diagnostic_1d.png)

### fci_sol_solution_2d.png

![fci_sol_solution_2d.png](../../media/examples/fci_sol_gallery/fci_sol_solution_2d.png)

---

[← Back to all examples](../index.html)
