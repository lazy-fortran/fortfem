---
title: maxwell_mesh_adjoint Example
---

# maxwell_mesh_adjoint Example

# Maxwell PML moving-mesh adjoint

This example benchmarks analytical forward and reverse products for
arbitrary-order-compatible tetrahedral Nedelec PML assembly. The differentiated
inputs include every mesh coordinate, each element's complex Cartesian
stretch, and the wave number.

The JVP propagates one simultaneous mesh/PML perturbation. The VJP consumes a
merged CSC matrix cotangent and returns the full design gradient in one reverse
element sweep, including accumulation at shared vertices. A complete
two-sided reassembly is retained as an independent directional oracle.

Generated gallery artifacts:

- `derivative_timing_1d.png`: measured JVP, VJP, and directional
  finite-difference timings, plus the coordinate-wise finite-difference
  full-gradient cost estimated from measured primal assembly;
- `adjoint_accuracy_1d.png`: relative reverse/central-difference error;
- `mesh_gradient_2d.png`: projected mesh-vertex descent directions;
- `benchmark.txt`: machine-readable sizes, timings, and errors.

No rendered image is committed.

## Method choice

Analytical reverse mode is the production choice for many geometry variables:
its element-sweep count is independent of the number of vertices. FortSym
generates the complex PML coefficient products, while FortNum's generated
determinant/inverse products and analytical Piola products handle geometry.
This path works with GCC. Enzyme remains an optional comparison oracle rather
than a runtime dependency.

## Provenance

- P. Monk, *Finite Element Methods for Maxwell's Equations*, Oxford
  University Press (2003), covariant Piola maps and Nedelec elements.
- J.-M. Jin, *The Finite Element Method in Electromagnetics*, 3rd ed.,
  Wiley (2014), Cartesian perfectly matched layers.
- J. Nocedal and S. Wright, *Numerical Optimization*, 2nd ed., Springer
  (2006), adjoint reduced gradients.

## Usage

```bash
fpm run --example maxwell_mesh_adjoint
```

## Source Code

```fortran
program maxwell_mesh_adjoint
    use fortfem_api, only: assemble_tetra_nedelec_pml_csc, &
        assemble_tetra_nedelec_pml_csc_jvp, &
        assemble_tetra_nedelec_pml_csc_vjp, &
        generate_structured_tetra_box_mesh
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, set_xscale, &
        set_yscale, title, xlabel, ylabel
    use fortsparse, only: csc_z_t, fortsparse_status_t
    implicit none

    integer, parameter :: level_count = 3, repetitions = 3
    real(dp), parameter :: difference_step = 2.0e-7_dp
    real(dp), parameter :: bounds(3, 2) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [3, 2])
    character(*), parameter :: output_directory = &
        "output/example/maxwell_mesh_adjoint"
    real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    real(dp), parameter :: green(3) = [0.0_dp, 158.0_dp, 115.0_dp]/255.0_dp
    real(dp) :: design_variables(level_count), directional_error(level_count)
    real(dp) :: finite_difference_seconds(level_count)
    real(dp) :: full_gradient_fd_seconds(level_count)
    real(dp) :: jvp_seconds(level_count), primal_seconds(level_count)
    real(dp) :: vjp_seconds(level_count)
    real(dp), allocatable :: vertices(:, :), vertices_dot(:, :)
    real(dp), allocatable :: vertices_bar(:, :), display_vertices(:, :)
    real(dp), allocatable :: finest_vertices(:, :), finest_vertices_bar(:, :)
    integer, allocatable :: tetrahedra(:, :)
    complex(dp), allocatable :: stretch(:, :), stretch_dot(:, :)
    complex(dp), allocatable :: stretch_bar(:, :), matrix_values_bar(:)
    type(csc_z_t) :: matrix, matrix_dot, matrix_minus, matrix_plus
    type(fortsparse_status_t) :: sparse_status
    real(dp) :: adjoint_derivative, central_derivative, elapsed, start_time
    real(dp) :: wave_number, wave_number_bar, wave_number_dot
    integer :: entry, level, repetition, unit

    call execute_command_line("mkdir -p "//output_directory)
    wave_number = 1.7_dp
    wave_number_dot = 0.13_dp
    do level = 1, level_count
        call generate_structured_tetra_box_mesh( &
            bounds, [level, level, level], vertices, tetrahedra, entry)
        if (entry /= 0) error stop "Maxwell benchmark mesh generation failed"
        allocate(vertices_dot, mold=vertices)
        vertices_dot = reshape([ &
            (0.02_dp*sin(real(5*entry + 1, dp)), &
            entry=1, size(vertices_dot))], shape(vertices_dot))
        allocate(stretch(3, size(tetrahedra, 2)))
        allocate(stretch_dot(3, size(tetrahedra, 2)))
        do entry = 1, size(tetrahedra, 2)
            stretch(:, entry) = [ &
                cmplx(1.0_dp, 0.15_dp + 0.001_dp*entry, dp), &
                cmplx(1.0_dp, 0.10_dp + 0.0007_dp*entry, dp), &
                cmplx(1.0_dp, 0.20_dp + 0.0005_dp*entry, dp)]
            stretch_dot(:, entry) = [ &
                cmplx(0.01_dp, -0.004_dp, dp), &
                cmplx(-0.006_dp, 0.003_dp, dp), &
                cmplx(0.008_dp, 0.005_dp, dp)]
        end do

        call cpu_time(start_time)
        do repetition = 1, repetitions
            call assemble_tetra_nedelec_pml_csc( &
                vertices, tetrahedra, 1, stretch, wave_number, matrix, &
                sparse_status)
            if (sparse_status%code /= 0) error stop "PML primal assembly failed"
        end do
        call cpu_time(elapsed)
        primal_seconds(level) = &
            (elapsed - start_time)/real(repetitions, dp)
        allocate(matrix_values_bar(matrix%nnz))
        do entry = 1, matrix%nnz
            matrix_values_bar(entry) = cmplx( &
                0.001_dp*cos(real(entry, dp)), &
                0.001_dp*sin(real(2*entry, dp)), dp)
        end do

        call cpu_time(start_time)
        do repetition = 1, repetitions
            call assemble_tetra_nedelec_pml_csc_jvp( &
                vertices, tetrahedra, 1, stretch, wave_number, vertices_dot, &
                stretch_dot, wave_number_dot, matrix_dot, sparse_status)
            if (sparse_status%code /= 0) error stop "PML assembly JVP failed"
        end do
        call cpu_time(elapsed)
        jvp_seconds(level) = (elapsed - start_time)/real(repetitions, dp)

        allocate(vertices_bar, mold=vertices)
        allocate(stretch_bar, mold=stretch)
        call cpu_time(start_time)
        do repetition = 1, repetitions
            call assemble_tetra_nedelec_pml_csc_vjp( &
                vertices, tetrahedra, 1, stretch, wave_number, &
                matrix_values_bar, vertices_bar, stretch_bar, wave_number_bar, &
                sparse_status)
            if (sparse_status%code /= 0) error stop "PML assembly VJP failed"
        end do
        call cpu_time(elapsed)
        vjp_seconds(level) = (elapsed - start_time)/real(repetitions, dp)

        call cpu_time(start_time)
        do repetition = 1, repetitions
            call assemble_tetra_nedelec_pml_csc( &
                vertices + difference_step*vertices_dot, tetrahedra, 1, &
                stretch + difference_step*stretch_dot, &
                wave_number + difference_step*wave_number_dot, matrix_plus, &
                sparse_status)
            call assemble_tetra_nedelec_pml_csc( &
                vertices - difference_step*vertices_dot, tetrahedra, 1, &
                stretch - difference_step*stretch_dot, &
                wave_number - difference_step*wave_number_dot, matrix_minus, &
                sparse_status)
        end do
        call cpu_time(elapsed)
        finite_difference_seconds(level) = &
            (elapsed - start_time)/real(repetitions, dp)
        if (matrix_plus%nnz /= matrix%nnz .or. &
            matrix_minus%nnz /= matrix%nnz) &
            error stop "PML benchmark CSC pattern changed"
        central_derivative = real(sum(conjg(matrix_values_bar)*( &
            matrix_plus%val - matrix_minus%val))/(2.0_dp*difference_step), dp)
        adjoint_derivative = sum(vertices_bar*vertices_dot) + &
            real(sum(conjg(stretch_bar)*stretch_dot), dp) + &
            wave_number_bar*wave_number_dot
        directional_error(level) = &
            abs(adjoint_derivative - central_derivative)/ &
            max(1.0_dp, abs(central_derivative))
        if (directional_error(level) > 2.0e-8_dp) &
            error stop "Maxwell mesh adjoint failed central difference"
        design_variables(level) = real( &
            size(vertices) + 2*size(stretch) + 1, dp)
        full_gradient_fd_seconds(level) = &
            2.0_dp*design_variables(level)*primal_seconds(level)
        if (level == level_count) then
            finest_vertices = vertices
            finest_vertices_bar = vertices_bar
        end if
        deallocate( &
            vertices, vertices_dot, vertices_bar, tetrahedra, stretch, &
            stretch_dot, stretch_bar, matrix_values_bar)
    end do

    call make_plots()
    open (newunit=unit, file=output_directory//"/benchmark.txt", &
        status="replace", action="write")
    write (unit, "(a)") &
        "design_variables primal_s jvp_s vjp_s directional_fd_s full_fd_est_s error"
    do level = 1, level_count
        write (unit, "(f10.0,6es16.7)") design_variables(level), &
            primal_seconds(level), jvp_seconds(level), vjp_seconds(level), &
            finite_difference_seconds(level), &
            full_gradient_fd_seconds(level), directional_error(level)
    end do
    close (unit)

contains

    subroutine make_plots()
        real(dp) :: scale
        real(dp) :: segment_x(2), segment_y(2)
        integer :: vertex

        call figure(figsize=[8.0_dp, 5.5_dp])
        call plot(design_variables, jvp_seconds, label="analytical JVP", &
            marker="o", color=blue)
        call plot(design_variables, vjp_seconds, label="analytical VJP", &
            marker="s", color=orange)
        call plot(design_variables, finite_difference_seconds, &
            label="measured directional FD", marker="^", color=green)
        call plot(design_variables, full_gradient_fd_seconds, &
            label="estimated coordinate FD gradient", linestyle="--")
        call set_xscale("log")
        call set_yscale("log")
        call xlabel("geometry and PML design variables")
        call ylabel("assembly derivative time (s)")
        call title("Maxwell PML mesh derivatives: forward, reverse, and FD")
        call legend()
        call savefig(output_directory//"/derivative_timing_1d.png")

        call figure(figsize=[7.5_dp, 5.0_dp])
        call plot(design_variables, max(directional_error, epsilon(1.0_dp)), &
            marker="o", color=blue)
        call set_xscale("log")
        call set_yscale("log")
        call xlabel("geometry and PML design variables")
        call ylabel("relative adjoint/central-difference error")
        call title("Maxwell mesh-adjoint verification")
        call savefig(output_directory//"/adjoint_accuracy_1d.png")

        scale = 0.12_dp/max(maxval(abs(finest_vertices_bar(1:2, :))), &
            tiny(1.0_dp))
        allocate(display_vertices, mold=finest_vertices)
        display_vertices = finest_vertices - scale*finest_vertices_bar
        call figure(figsize=[7.0_dp, 6.0_dp])
        do vertex = 1, size(finest_vertices, 2)
            segment_x = [ &
                finest_vertices(1, vertex), display_vertices(1, vertex)]
            segment_y = [ &
                finest_vertices(2, vertex), display_vertices(2, vertex)]
            if (vertex == 1) then
                call plot( &
                    segment_x, segment_y, color=orange, linewidth=1.5_dp, &
                    label="adjoint descent displacement")
            else
                call plot( &
                    segment_x, segment_y, color=orange, linewidth=1.5_dp)
            end if
        end do
        call plot( &
            finest_vertices(1, :), finest_vertices(2, :), &
            marker="o", linestyle="none", color=blue, label="mesh vertices")
        call xlabel("x")
        call ylabel("y")
        call title("Projected Maxwell adjoint descent directions")
        call legend()
        call savefig(output_directory//"/mesh_gradient_2d.png")
    end subroutine make_plots

end program maxwell_mesh_adjoint
```

## Generated Plots

### adjoint_accuracy_1d.png

![adjoint_accuracy_1d.png](../../media/examples/maxwell_mesh_adjoint/adjoint_accuracy_1d.png)

### derivative_timing_1d.png

![derivative_timing_1d.png](../../media/examples/maxwell_mesh_adjoint/derivative_timing_1d.png)

### mesh_gradient_2d.png

![mesh_gradient_2d.png](../../media/examples/maxwell_mesh_adjoint/mesh_gradient_2d.png)

---

[← Back to all examples](../index.html)
