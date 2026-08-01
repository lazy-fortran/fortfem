---
title: cgl_pressure_tensor Example
---

# cgl_pressure_tensor Example

# CGL pressure tensor example

This manufactured profile evaluates the gyrotropic CGL pressure tensor

\[
\mathbf P=p_\perp I+(p_\parallel-p_\perp)\mathbf b\mathbf b^T
\]

for a magnetic direction that rotates across a one-dimensional coordinate.
It also evaluates the generated product-rule force \(\nabla\cdot\mathbf P\),
checks the tensor trace independently, and writes reproducible CSV data.

The generated gallery plots show the diagonal/off-diagonal pressure
components and the two nonzero force-divergence components. This is a
constitutive manufactured-solution example, not a complete equilibrium or
transport solve.

Outputs:

- `cgl_pressure_components_1d.png`
- `cgl_force_divergence_1d.png`
- `cgl_profile.csv`

## Usage

```bash
fpm run --example cgl_pressure_tensor
```

## Source Code

```fortran
program cgl_pressure_tensor
    ! Manufactured CGL pressure and force-divergence profile.
    ! The magnetic direction rotates in the x direction, so both the
    ! gyrotropic pressure components and the product-rule force are visible.
    use fortfem_api, only: &
        evaluate_cgl_pressure_divergence, evaluate_cgl_pressure_tensor
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: point_count = 161
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(*), parameter :: output_directory = &
        "output/example/cgl_pressure_tensor"
    real(dp) :: coordinate(point_count), pressure_11(point_count)
    real(dp) :: pressure_22(point_count), pressure_12(point_count)
    real(dp) :: force_x(point_count), force_y(point_count)
    real(dp) :: direction(3), parallel_gradient(3)
    real(dp) :: perpendicular_gradient(3), direction_gradient(3, 3)
    real(dp) :: pressure_tensor(3, 3), force_divergence(3)
    real(dp) :: x, angle, p_parallel, p_perpendicular, angle_gradient
    real(dp) :: maximum_trace_error
    integer :: point
    type(fortsparse_status_t) :: status
    integer :: unit

    call execute_command_line("mkdir -p "//output_directory)
    maximum_trace_error = 0.0_dp
    do point = 1, point_count
        x = -1.0_dp + 2.0_dp*real(point - 1, dp)/real(point_count - 1, dp)
        coordinate(point) = x
        angle = 0.35_dp*x
        angle_gradient = 0.35_dp
        p_parallel = 3.0_dp + 0.2_dp*x
        p_perpendicular = 1.0_dp + 0.1_dp*cos(pi*x)
        direction = [cos(angle), sin(angle), 0.0_dp]
        parallel_gradient = [0.2_dp, 0.0_dp, 0.0_dp]
        perpendicular_gradient = [-0.1_dp*pi*sin(pi*x), 0.0_dp, 0.0_dp]
        direction_gradient = 0.0_dp
        direction_gradient(1, 1) = -sin(angle)*angle_gradient
        direction_gradient(2, 1) = cos(angle)*angle_gradient

        call evaluate_cgl_pressure_tensor( &
            p_parallel, p_perpendicular, direction, pressure_tensor, status)
        if (status%code /= 0) error stop "CGL tensor evaluation failed"
        call evaluate_cgl_pressure_divergence( &
            p_parallel, p_perpendicular, direction, parallel_gradient, &
            perpendicular_gradient, direction_gradient, force_divergence, status)
        if (status%code /= 0) error stop "CGL force evaluation failed"
        pressure_11(point) = pressure_tensor(1, 1)
        pressure_22(point) = pressure_tensor(2, 2)
        pressure_12(point) = pressure_tensor(1, 2)
        force_x(point) = force_divergence(1)
        force_y(point) = force_divergence(2)
        maximum_trace_error = max(maximum_trace_error, abs( &
            pressure_tensor(1, 1) + pressure_tensor(2, 2) + &
            pressure_tensor(3, 3) - p_parallel - 2.0_dp*p_perpendicular))
    end do
    if (maximum_trace_error > 2.0e-13_dp) &
        error stop "CGL tensor trace oracle failed"

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(coordinate, pressure_11, label="P11")
    call plot(coordinate, pressure_22, label="P22")
    call plot(coordinate, pressure_12, label="P12")
    call xlabel("x")
    call ylabel("pressure-tensor component")
    call title("Manufactured CGL pressure tensor")
    call legend()
    call savefig(output_directory//"/cgl_pressure_components_1d.png")

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(coordinate, force_x, label="div(P)_x")
    call plot(coordinate, force_y, label="div(P)_y")
    call xlabel("x")
    call ylabel("force density")
    call title("CGL product-rule force divergence")
    call legend()
    call savefig(output_directory//"/cgl_force_divergence_1d.png")

    open (newunit=unit, file=output_directory//"/cgl_profile.csv", &
        status="replace", action="write")
    write (unit, "(a)") "x,P11,P22,P12,divP_x,divP_y"
    do point = 1, point_count
        write (unit, "(6(es24.16,1x))") &
            coordinate(point), pressure_11(point), pressure_22(point), &
            pressure_12(point), force_x(point), force_y(point)
    end do
    close (unit)
end program cgl_pressure_tensor
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/cgl_pressure_tensor/primary.png)

### cgl_force_divergence_1d.png

![cgl_force_divergence_1d.png](../../media/examples/cgl_pressure_tensor/cgl_force_divergence_1d.png)

### cgl_pressure_components_1d.png

![cgl_pressure_components_1d.png](../../media/examples/cgl_pressure_tensor/cgl_pressure_components_1d.png)

---

[← Back to all examples](../index.html)
