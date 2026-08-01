---
title: laplace_exterior_bem_circle Example
---

# laplace_exterior_bem_circle Example

# Exterior Laplace BEM on a circle

This example solves the unbounded exterior Dirichlet problem on the unit
circle with the manufactured decaying solution

\[
u(r,\theta)=\frac{\cos\theta}{r}.
\]

The periodic piecewise-constant single-layer density is assembled and solved
with FortFEM's dense solver, then evaluated directly at exterior points. The
gallery leads with the computed two-dimensional exterior solution, followed by
an analytical radial trace and the recovered boundary density. The interior of
the circle is shown as a neutral zero-valued region in the solution plot; no
artificial far boundary is introduced.

## Provenance

- S. Sauter and C. Schwab, *Boundary Element Methods*, Springer, 2011,
  <https://doi.org/10.1007/978-3-540-68093-2>.
- FortFEM's logarithmic single-layer kernel uses the standard two-dimensional
  Laplace fundamental solution \(G(x,y)=-(2\pi)^{-1}\log|x-y|\).

## Usage

```bash
fpm run --example laplace_exterior_bem_circle
```

## Source Code

```fortran
program laplace_exterior_bem_circle
    !! Exterior Laplace Dirichlet problem on the unit circle.
    !!
    !! The single-layer representation is solved on a periodic piecewise
    !! constant boundary space and evaluated in the unbounded exterior.  The
    !! manufactured solution is u(r,theta) = cos(theta)/r.
    use fortfem_api, only: assemble_laplace_single_layer_constant
    use fortfem_advanced_solvers, only: solver_options, solver_options_t, &
        solver_stats_t, solve
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortfem_kinds, only: dp
    use fortplot, only: add_scatter, colorbar, figure, legend, pcolormesh, &
        plot, savefig, title, xlabel, ylabel
    implicit none

    integer, parameter :: panel_count = 64, quadrature_order = 16
    integer, parameter :: field_nx = 81, field_ny = 81, trace_count = 121
    real(dp), parameter :: exterior_radius = 3.0_dp
    character(*), parameter :: output_directory = &
        "output/example/laplace_exterior_bem_circle"
    real(dp), allocatable :: panel_start(:, :), panel_end(:, :)
    real(dp), allocatable :: single_layer(:, :), density(:), boundary_data(:)
    real(dp), allocatable :: panel_lengths(:), right_hand_side(:)
    real(dp), allocatable :: field_x(:), field_y(:), field(:, :)
    real(dp), allocatable :: exact_field(:, :), radial(:), bem_trace(:)
    real(dp), allocatable :: exact_trace(:), circle_x(:), circle_y(:)
    real(dp) :: matrix_residual, field_error, trace_error, start_time, elapsed
    real(dp) :: boundary_angle, x_value, y_value, radius, pi
    integer :: command_status, i, j, panel, status, unit
    type(solver_options_t) :: options
    type(solver_stats_t) :: statistics

    pi = acos(-1.0_dp)
    call execute_command_line( &
        "mkdir -p "//output_directory, exitstat=command_status)
    if (command_status /= 0) error stop "cannot create BEM output directory"

    allocate( &
        panel_start(2, panel_count), panel_end(2, panel_count), &
        single_layer(panel_count, panel_count), density(panel_count), &
        boundary_data(panel_count))
    do panel = 1, panel_count
        boundary_angle = 2.0_dp*pi*real(panel - 1, dp)/real(panel_count, dp)
        panel_start(:, panel) = [cos(boundary_angle), sin(boundary_angle)]
        boundary_angle = 2.0_dp*pi*real(panel, dp)/real(panel_count, dp)
        panel_end(:, panel) = [cos(boundary_angle), sin(boundary_angle)]
        boundary_angle = 2.0_dp*pi*(real(panel, dp) - 0.5_dp)/ &
            real(panel_count, dp)
        boundary_data(panel) = cos(boundary_angle)
    end do
    allocate(panel_lengths(panel_count), right_hand_side(panel_count))
    do panel = 1, panel_count
        panel_lengths(panel) = norm2(panel_end(:, panel) - panel_start(:, panel))
    end do

    call cpu_time(start_time)
    call assemble_laplace_single_layer_constant( &
        panel_start, panel_end, quadrature_order, single_layer, status)
    if (status /= 0) error stop "exterior circle single-layer assembly failed"
    right_hand_side = panel_lengths*boundary_data
    options = solver_options(method="lapack_lu")
    call solve(single_layer, right_hand_side, density, options, statistics)
    if (.not. statistics%converged) error stop "exterior circle BEM solve failed"
    call cpu_time(elapsed)
    elapsed = elapsed - start_time
    matrix_residual = maxval(abs( &
        matmul(single_layer, density)/panel_lengths - boundary_data))

    allocate( &
        field_x(field_nx), field_y(field_ny), &
        field(field_ny, field_nx), exact_field(field_ny, field_nx))
    do i = 1, field_nx
        field_x(i) = -exterior_radius + 2.0_dp*exterior_radius* &
            real(i - 1, dp)/real(field_nx - 1, dp)
        field_y(i) = -exterior_radius + 2.0_dp*exterior_radius* &
            real(i - 1, dp)/real(field_ny - 1, dp)
    end do
    field = 0.0_dp
    exact_field = 0.0_dp
    field_error = 0.0_dp
    do j = 1, field_ny
        y_value = field_y(j)
        do i = 1, field_nx
            x_value = field_x(i)
            radius = sqrt(x_value*x_value + y_value*y_value)
            if (radius <= 1.02_dp) cycle
            call evaluate_single_layer( &
                [x_value, y_value], density, panel_start, panel_end, &
                field(j, i), status)
            if (status /= 0) error stop "exterior field evaluation failed"
            exact_field(j, i) = x_value/(radius*radius)
            field_error = max(field_error, abs(field(j, i) - exact_field(j, i)))
        end do
    end do

    allocate(radial(trace_count), bem_trace(trace_count), exact_trace(trace_count))
    do i = 1, trace_count
        radial(i) = 1.02_dp + (exterior_radius - 1.02_dp)* &
            real(i - 1, dp)/real(trace_count - 1, dp)
        call evaluate_single_layer( &
            [radial(i), 0.0_dp], density, panel_start, panel_end, &
            bem_trace(i), status)
        if (status /= 0) error stop "exterior trace evaluation failed"
        exact_trace(i) = 1.0_dp/radial(i)
    end do
    trace_error = maxval(abs(bem_trace - exact_trace))

    allocate(circle_x(257), circle_y(257))
    do i = 1, size(circle_x)
        boundary_angle = 2.0_dp*pi*real(i - 1, dp)/real(size(circle_x) - 1, dp)
        circle_x(i) = cos(boundary_angle)
        circle_y(i) = sin(boundary_angle)
    end do
    call figure(figsize=[8.0_dp, 6.5_dp])
    call pcolormesh(field_x, field_y, field, cmap="coolwarm")
    call colorbar(label="computed exterior u_h")
    call plot(circle_x, circle_y, color=[0.0_dp, 0.0_dp, 0.0_dp], &
        linewidth=2.0_dp, &
        label="unit-circle boundary")
    call xlabel("x")
    call ylabel("y")
    call title("Exterior Laplace BEM solution: cos(theta)/r")
    call legend()
    call savefig(output_directory//"/exterior_laplace_solution_2d.png")

    call figure(figsize=[8.0_dp, 5.5_dp])
    call plot(radial, exact_trace, label="analytical r^{-1} cos(theta)", &
        linestyle="-")
    call plot(radial, bem_trace, label="single-layer BEM", linestyle="--")
    call xlabel("radius on the positive x-axis")
    call ylabel("u(r,0)")
    call title("Exterior Laplace BEM radial trace")
    call legend()
    call savefig(output_directory//"/exterior_laplace_trace_1d.png")

    call figure(figsize=[8.0_dp, 6.0_dp])
    call add_scatter( &
        0.5_dp*(panel_start(1, :) + panel_end(1, :)), &
        0.5_dp*(panel_start(2, :) + panel_end(2, :)), c=density, &
        cmap="viridis", marker=".", markersize=7.0_dp, &
        label="single-layer density sigma_h")
    call title("Exterior Laplace BEM boundary density")
    call xlabel("boundary x")
    call ylabel("boundary y")
    call colorbar(label="sigma_h")
    call savefig(output_directory//"/exterior_laplace_density_2d.png")

    open (newunit=unit, file=output_directory//"/benchmark.txt", &
        status="replace", action="write")
    write (unit, "(a,i0)") "circle panels: ", panel_count
    write (unit, "(a,es14.6)") "single-layer solve seconds: ", elapsed
    write (unit, "(a,es14.6)") "boundary matrix residual: ", matrix_residual
    write (unit, "(a,es14.6)") "maximum exterior field error: ", field_error
    write (unit, "(a,es14.6)") "maximum radial trace error: ", trace_error
    close (unit)

contains

    subroutine evaluate_single_layer(point, source_density, starts, ends, &
            value, local_status)
        real(dp), intent(in) :: point(2), source_density(:)
        real(dp), intent(in) :: starts(:, :), ends(:, :)
        real(dp), intent(out) :: value
        integer, intent(out) :: local_status

        real(dp), allocatable :: nodes(:), weights(:)
        real(dp) :: source(2), difference(2), length, distance
        integer :: local_panel, quadrature_point, quadrature_status

        value = 0.0_dp
        local_status = 1
        if (size(source_density) /= size(starts, 2) .or. &
            any(shape(starts) /= shape(ends))) return
        allocate(nodes(quadrature_order), weights(quadrature_order))
        call gauss_legendre_ab( &
            quadrature_order, 0.0_dp, 1.0_dp, nodes, weights)
        do local_panel = 1, size(starts, 2)
            length = norm2(ends(:, local_panel) - starts(:, local_panel))
            do quadrature_point = 1, quadrature_order
                source = starts(:, local_panel) + nodes(quadrature_point)* &
                    (ends(:, local_panel) - starts(:, local_panel))
                difference = point - source
                distance = norm2(difference)
                if (distance <= 1.0e-13_dp) then
                    local_status = 1
                    return
                end if
                value = value - source_density(local_panel)*length* &
                    weights(quadrature_point)*log(distance)/(2.0_dp*acos(-1.0_dp))
            end do
        end do
        local_status = 0
    end subroutine evaluate_single_layer

end program laplace_exterior_bem_circle
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/laplace_exterior_bem_circle/primary.png)

### exterior_laplace_density_2d.png

![exterior_laplace_density_2d.png](../../media/examples/laplace_exterior_bem_circle/exterior_laplace_density_2d.png)

### exterior_laplace_solution_2d.png

![exterior_laplace_solution_2d.png](../../media/examples/laplace_exterior_bem_circle/exterior_laplace_solution_2d.png)

### exterior_laplace_trace_1d.png

![exterior_laplace_trace_1d.png](../../media/examples/laplace_exterior_bem_circle/exterior_laplace_trace_1d.png)

---

[← Back to all examples](../index.html)
