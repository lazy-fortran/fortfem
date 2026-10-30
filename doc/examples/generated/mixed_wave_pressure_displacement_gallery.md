---
title: mixed_wave_pressure_displacement_gallery Example
---

# mixed_wave_pressure_displacement_gallery Example

# Mixed pressure--displacement wave state

This neutral gallery advances a first-order mixed wave state with the existing
structure-preserving midpoint map.  The coordinate block is intentionally
labelled as one pressure amplitude and two displacement amplitudes; the second
block is labelled as the conjugate momentum and two velocities.  It is not a
second-order displacement-only update.

The first three outputs are physical solutions: a three-dimensional
pressure/displacement trajectory, a two-dimensional pressure contour with
displacement arrows, and one-dimensional pressure, displacement, momentum, and
velocity traces.  Energy, exact-state error, reversibility, and symplectic-map
defects are emitted afterward as diagnostics.  Diagonal positive masses make
the manufactured modal solution independently computable without any plasma or
material-specific model.

The corresponding test derives the modal solution and Cayley map independently
of the example implementation.  This fixture is therefore a reusable
structure-preserving time-integration contract for later FEM/IGA wave,
elasticity, acoustics, or anisotropic systems.

## Usage

```bash
fpm run --example mixed_wave_pressure_displacement_gallery
```

## Source Code

```fortran
program mixed_wave_pressure_displacement_gallery
    !! Neutral mixed pressure/displacement wave solution gallery.
    !!
    !! The first-order state is (pressure, displacement) plus
    !! (momentum, velocity).  The midpoint map is energy preserving and
    !! reversible; no second-order-only reduction is used here.
    use fortfem_time, only: advance_mixed_wave_midpoint
    use fortfem_kinds, only: dp
    use fortplot, only: add_3d_plot, colorbar, contourf, figure, legend, plot, &
        quiver, savefig, set_yscale, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: n = 3, step_count = 180
    integer, parameter :: line_count = 121, grid_count = 13
    integer, parameter :: vector_count = grid_count*grid_count
    real(dp), parameter :: time_step = 0.025_dp
    real(dp), parameter :: mass_q(n, n) = reshape([ &
        1.40_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.90_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.70_dp], [n, n])
    real(dp), parameter :: mass_v(n, n) = reshape([ &
        0.80_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.30_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.10_dp], [n, n])
    real(dp), parameter :: coupling(n, n) = reshape([ &
        0.95_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.35_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.75_dp], [n, n])
    real(dp), parameter :: q_initial(n) = [0.65_dp, -0.35_dp, 0.22_dp]
    real(dp), parameter :: v_initial(n) = [-0.18_dp, 0.42_dp, 0.31_dp]
    real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    real(dp), parameter :: green(3) = [0.0_dp, 158.0_dp, 115.0_dp]/255.0_dp
    real(dp), parameter :: red(3) = [213.0_dp, 94.0_dp, 0.0_dp]/255.0_dp
    character(*), parameter :: output_directory = &
        "output/example/mixed_wave_pressure_displacement_gallery"
    real(dp) :: q(n), v(n), q_exact(n), v_exact(n)
    real(dp) :: q_history(n, step_count + 1), v_history(n, step_count + 1)
    real(dp) :: time(step_count + 1), energy_history(step_count + 1)
    real(dp) :: error_history(step_count + 1)
    real(dp) :: q_reverse(n), v_reverse(n)
    real(dp) :: map(2*n, 2*n), symplectic_form(2*n, 2*n)
    real(dp) :: line_coordinate(line_count), pressure_line(line_count)
    real(dp) :: displacement_x_line(line_count), displacement_y_line(line_count)
    real(dp) :: momentum_line(line_count), velocity_x_line(line_count)
    real(dp) :: pressure_grid(grid_count, grid_count)
    real(dp) :: vector_x(vector_count), vector_y(vector_count)
    real(dp) :: displacement_x_grid(vector_count), displacement_y_grid(vector_count)
    real(dp) :: pressure_value, displacement_x, displacement_y
    real(dp) :: energy_initial, maximum_error, maximum_drift
    real(dp) :: reversal_error, symplectic_error, elapsed, start_time
    integer :: command_status, step, point, ix, iy, vector_point, unit
    type(fortsparse_status_t) :: status

    call execute_command_line( &
        "mkdir -p "//output_directory, exitstat=command_status)
    if (command_status /= 0) error stop "cannot create mixed-wave output"
    call initialize_gallery_sequence()

    q = q_initial
    v = v_initial
    energy_initial = wave_energy(q, v)
    time(1) = 0.0_dp
    q_history(:, 1) = q
    v_history(:, 1) = v
    energy_history(1) = energy_initial
    error_history(1) = 0.0_dp
    call cpu_time(start_time)
    do step = 1, step_count
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling, time_step, q, v, status)
        if (status%code /= 0) error stop "mixed-wave midpoint failed"
        time(step + 1) = real(step, dp)*time_step
        q_history(:, step + 1) = q
        v_history(:, step + 1) = v
        energy_history(step + 1) = wave_energy(q, v)
        call exact_state(time(step + 1), q_exact, v_exact)
        error_history(step + 1) = max(maxval(abs(q - q_exact)), &
            maxval(abs(v - v_exact)))
    end do
    call cpu_time(elapsed)
    elapsed = elapsed - start_time

    maximum_error = maxval(error_history)
    maximum_drift = maxval(abs(energy_history/energy_initial - 1.0_dp))
    q_reverse = q
    v_reverse = v
    do step = 1, step_count
        call advance_mixed_wave_midpoint( &
            mass_q, mass_v, coupling, -time_step, q_reverse, v_reverse, status)
        if (status%code /= 0) error stop "mixed-wave reverse step failed"
    end do
    reversal_error = max(maxval(abs(q_reverse - q_initial)), &
        maxval(abs(v_reverse - v_initial)))
    call build_api_map(map, status)
    if (status%code /= 0) error stop "mixed-wave map construction failed"
    call build_symplectic_form(symplectic_form)
    symplectic_error = maxval(abs(matmul(transpose(map), &
        matmul(symplectic_form, map)) - symplectic_form))
    if (maximum_error > 2.0e-3_dp) error stop "mixed-wave exact oracle failed"
    if (maximum_drift > 3.0e-12_dp) error stop "mixed-wave energy failed"
    if (reversal_error > 3.0e-12_dp) error stop "mixed-wave reversibility failed"
    if (symplectic_error > 5.0e-11_dp) error stop "mixed-wave symplectic check failed"

    call reconstruct_fields()
    call render_physical_solution()
    call record_gallery_stage("physical_solution")
    call render_diagnostics()
    call record_gallery_stage("diagnostics")
    call write_outputs(elapsed)

contains

    subroutine exact_state(current_time, q_value, v_value)
        real(dp), intent(in) :: current_time
        real(dp), intent(out) :: q_value(:), v_value(:)
        real(dp) :: frequency, scale_q, scale_v
        integer :: mode

        do mode = 1, n
            frequency = coupling(mode, mode)/sqrt( &
                mass_q(mode, mode)*mass_v(mode, mode))
            scale_q = sqrt(mass_v(mode, mode)/mass_q(mode, mode))
            scale_v = sqrt(mass_q(mode, mode)/mass_v(mode, mode))
            q_value(mode) = q_initial(mode)*cos(frequency*current_time) - &
                scale_q*v_initial(mode)*sin(frequency*current_time)
            v_value(mode) = v_initial(mode)*cos(frequency*current_time) + &
                scale_v*q_initial(mode)*sin(frequency*current_time)
        end do
    end subroutine exact_state

    subroutine reconstruct_fields()
        real(dp) :: x, y

        do point = 1, line_count
            line_coordinate(point) = real(point - 1, dp)/ &
                real(line_count - 1, dp)
            x = line_coordinate(point)
            pressure_line(point) = q(1)*sin(acos(-1.0_dp)*x)
            displacement_x_line(point) = q(2)*sin(2.0_dp*acos(-1.0_dp)*x)
            displacement_y_line(point) = q(3)*cos(acos(-1.0_dp)*x)
            momentum_line(point) = v(1)*sin(acos(-1.0_dp)*x)
            velocity_x_line(point) = v(2)*sin(2.0_dp*acos(-1.0_dp)*x)
        end do
        vector_point = 0
        do iy = 1, grid_count
            y = real(iy - 1, dp)/real(grid_count - 1, dp)
            do ix = 1, grid_count
                x = real(ix - 1, dp)/real(grid_count - 1, dp)
                vector_point = vector_point + 1
                vector_x(vector_point) = x
                vector_y(vector_point) = y
                pressure_value = q(1)*sin(acos(-1.0_dp)*x)* &
                    cos(acos(-1.0_dp)*y)
                displacement_x = q(2)*sin(acos(-1.0_dp)*x)* &
                    sin(acos(-1.0_dp)*y)
                displacement_y = q(3)*cos(acos(-1.0_dp)*x)* &
                    sin(acos(-1.0_dp)*y)
                pressure_grid(ix, iy) = pressure_value
                displacement_x_grid(vector_point) = displacement_x
                displacement_y_grid(vector_point) = displacement_y
            end do
        end do
    end subroutine reconstruct_fields

    subroutine render_physical_solution()
        call figure(figsize=[8.0_dp, 6.5_dp])
        call add_3d_plot( &
            q_history(2, :), q_history(3, :), q_history(1, :), &
            label="midpoint (u_x,u_y,p)", color=blue, linewidth=2.0_dp)
        call exact_trajectory_plot()
        call xlabel("displacement u_x")
        call ylabel("displacement u_y")
        call title("Mixed pressure/displacement physical trajectory")
        call legend()
        call savefig(output_directory//"/solution_3d.png")

        call figure(figsize=[8.0_dp, 6.0_dp])
        call contourf( &
            [(real(ix - 1, dp)/real(grid_count - 1, dp), ix=1,grid_count)], &
            [(real(iy - 1, dp)/real(grid_count - 1, dp), iy=1,grid_count)], &
            pressure_grid, cmap="coolwarm", show_colorbar=.true.)
        call colorbar(label="pressure p(x,y)")
        ! Keep the glyphs within the unit square; CSV output retains the
        ! unscaled physical displacement components.
        call quiver(vector_x, vector_y, 0.06_dp*displacement_x_grid, &
            0.06_dp*displacement_y_grid, color=blue, scale=1.0_dp, width=0.003_dp)
        call xlabel("x")
        call ylabel("y")
        call title("Pressure and displacement vectors")
        call savefig(output_directory//"/solution_2d.png")

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(line_coordinate, pressure_line, label="pressure p", color=blue)
        call plot(line_coordinate, displacement_x_line, &
            label="displacement u_x", color=orange)
        call plot(line_coordinate, displacement_y_line, &
            label="displacement u_y", color=green)
        call plot(line_coordinate, momentum_line, label="momentum pi", &
            color=red, linestyle="--")
        call plot(line_coordinate, velocity_x_line, label="velocity v_x", &
            color=[0.45_dp, 0.25_dp, 0.65_dp], linestyle=":")
        call xlabel("one-dimensional coordinate x")
        call ylabel("field amplitude")
        call title("Mixed pressure, displacement, momentum, and velocity")
        call legend()
        call savefig(output_directory//"/solution_1d.png")
    end subroutine render_physical_solution

    subroutine exact_trajectory_plot()
        real(dp) :: q_value(n), v_value(n)
        real(dp) :: exact_x(step_count + 1), exact_y(step_count + 1)
        real(dp) :: exact_z(step_count + 1)
        integer :: local_step

        do local_step = 1, step_count + 1
            call exact_state(time(local_step), q_value, v_value)
            exact_x(local_step) = q_value(2)
            exact_y(local_step) = q_value(3)
            exact_z(local_step) = q_value(1)
        end do
        call add_3d_plot(exact_x, exact_y, exact_z, &
            label="independent modal solution", color=orange, &
            linestyle="--", linewidth=1.3_dp)
    end subroutine exact_trajectory_plot

    subroutine render_diagnostics()
        real(dp) :: relative_energy(step_count + 1)
        real(dp) :: symplectic_history(step_count + 1)
        real(dp) :: reversal_history(step_count + 1)

        relative_energy = energy_history/energy_initial - 1.0_dp
        symplectic_history = symplectic_error
        reversal_history = reversal_error
        reversal_history(1) = 0.0_dp
        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(time, relative_energy, label="relative energy drift", color=green)
        call plot(time, error_history, label="exact state error", color=orange)
        call xlabel("time")
        call ylabel("ledger value")
        call title("Mixed-wave energy and analytical-state ledgers")
        call legend()
        call savefig(output_directory//"/diagnostics_energy_error.png")

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(time, symplectic_history, label="symplectic map defect", color=blue)
        call plot(time, reversal_history, label="reversibility endpoint", color=red, &
            linestyle="--")
        call xlabel("time")
        call ylabel("absolute defect")
        call title("Mixed-wave reversibility and symplectic ledgers")
        call set_yscale("log")
        call legend()
        call savefig(output_directory//"/diagnostics_structure.png")
    end subroutine render_diagnostics

    subroutine build_api_map(map_value, map_status)
        real(dp), intent(out) :: map_value(:, :)
        type(fortsparse_status_t), intent(out) :: map_status
        real(dp) :: q_work(n), v_work(n)
        integer :: column

        map_value = 0.0_dp
        do column = 1, 2*n
            q_work = 0.0_dp
            v_work = 0.0_dp
            if (column <= n) then
                q_work(column) = 1.0_dp
            else
                v_work(column - n) = 1.0_dp
            end if
            call advance_mixed_wave_midpoint( &
                mass_q, mass_v, coupling, time_step, q_work, v_work, map_status)
            if (map_status%code /= 0) return
            map_value(:, column) = [q_work, v_work]
        end do
    end subroutine build_api_map

    subroutine build_symplectic_form(form)
        real(dp), intent(out) :: form(:, :)
        integer :: mode

        form = 0.0_dp
        do mode = 1, n
            form(mode, n + mode) = mass_q(mode, mode)*mass_v(mode, mode)/ &
                coupling(mode, mode)
            form(n + mode, mode) = -form(mode, n + mode)
        end do
    end subroutine build_symplectic_form

    subroutine write_outputs(wall_time)
        real(dp), intent(in) :: wall_time
        integer :: local_point

        open (newunit=unit, file=output_directory//"/solution_1d.csv", &
            status="replace", action="write")
        write (unit, "(a)") &
            "x,pressure,displacement_x,displacement_y,momentum,velocity_x"
        do local_point = 1, line_count
            write (unit, "(*(es24.16,:,','))") line_coordinate(local_point), &
                pressure_line(local_point), displacement_x_line(local_point), &
                displacement_y_line(local_point), momentum_line(local_point), &
                velocity_x_line(local_point)
        end do
        close (unit)

        open (newunit=unit, file=output_directory//"/solution_2d.csv", &
            status="replace", action="write")
        write (unit, "(a)") "x,y,pressure,displacement_x,displacement_y"
        do local_point = 1, vector_count
            ix = 1 + mod(local_point - 1, grid_count)
            iy = 1 + (local_point - 1)/grid_count
            write (unit, "(*(es24.16,:,','))") vector_x(local_point), &
                vector_y(local_point), pressure_grid(ix, iy), &
                displacement_x_grid(local_point), displacement_y_grid(local_point)
        end do
        close (unit)

        open (newunit=unit, file=output_directory//"/diagnostics.csv", &
            status="replace", action="write")
        write (unit, "(a)") "quantity,value"
        write (unit, "(a,es24.16)") "maximum_exact_state_error,", maximum_error
        write (unit, "(a,es24.16)") "maximum_relative_energy_drift,", maximum_drift
        write (unit, "(a,es24.16)") "reversal_error,", reversal_error
        write (unit, "(a,es24.16)") "symplectic_map_defect,", symplectic_error
        write (unit, "(a,es24.16)") "midpoint_wall_time_seconds,", wall_time
        write (unit, "(a)") "primary_plot,solution_3d.png"
        close (unit)
    end subroutine write_outputs

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

    pure function wave_energy(q_value, v_value) result(value)
        real(dp), intent(in) :: q_value(:), v_value(:)
        real(dp) :: value

        value = 0.5_dp*(dot_product(q_value, matmul(mass_q, q_value)) + &
            dot_product(v_value, matmul(mass_v, v_value)))
    end function wave_energy

end program mixed_wave_pressure_displacement_gallery
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/mixed_wave_pressure_displacement_gallery/primary.png)

### solution_2d.png

![solution_2d.png](../../media/examples/mixed_wave_pressure_displacement_gallery/solution_2d.png)

---

[← Back to all examples](../index.html)
