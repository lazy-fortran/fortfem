---
title: cgl_pressure_tensor Example
---

# cgl_pressure_tensor Example

# CGL pressure tensor example

This manufactured profile evaluates the gyrotropic CGL pressure tensor

\[
\mathbf P=p_\perp I+(p_\parallel-p_\perp)\mathbf b\mathbf b^T
\]

for a magnetic direction that rotates across a two-dimensional sampled field.
It also evaluates the generated product-rule force \(\nabla\cdot\mathbf P\),
checks the tensor against an independent projector oracle, and writes
reproducible CSV data.

The physical-first gallery plot shows the tensor anisotropy and principal
direction field. The original one-dimensional diagonal/off-diagonal pressure
components and the two nonzero force-divergence components remain as
diagnostics. This is a constitutive manufactured-solution example, not a
complete equilibrium or transport solve.

Outputs:

- `cgl_tensor_principal_2d.png`
- `cgl_pressure_components_1d.png`
- `cgl_force_divergence_1d.png`
- `cgl_profile.csv`
- `cgl_tensor_field_2d.csv`

## Usage

```bash
fpm run --example cgl_pressure_tensor
```

## Source Code

```fortran
program cgl_pressure_tensor
    ! Manufactured CGL pressure and force-divergence profile.
    ! The magnetic direction rotates across the sampled x-y field, so both the
    ! gyrotropic pressure components and the product-rule force are visible.
    use fortfem_feec, only: &
        evaluate_cgl_pressure_divergence, evaluate_cgl_pressure_tensor
    use fortfem_kinds, only: dp
    use fortplot, only: colorbar, contourf, figure, legend, plot, quiver, savefig, &
        title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: point_count = 161
    integer, parameter :: tensor_nx = 25, tensor_ny = 21
    integer, parameter :: arrow_nx = 9, arrow_ny = 7
    integer, parameter :: arrow_count = arrow_nx*arrow_ny
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(*), parameter :: output_directory = &
        "output/example/cgl_pressure_tensor"
    real(dp) :: coordinate(point_count), pressure_11(point_count)
    real(dp) :: pressure_22(point_count), pressure_12(point_count)
    real(dp) :: force_x(point_count), force_y(point_count)
    real(dp) :: tensor_x(tensor_nx), tensor_y(tensor_ny)
    real(dp) :: anisotropy_map(tensor_nx, tensor_ny)
    real(dp) :: principal_x(arrow_count), principal_y(arrow_count)
    real(dp) :: principal_u(arrow_count), principal_v(arrow_count)
    real(dp) :: sampled_p11(tensor_nx, tensor_ny)
    real(dp) :: sampled_p22(tensor_nx, tensor_ny)
    real(dp) :: sampled_p12(tensor_nx, tensor_ny)
    real(dp) :: sampled_trace_error(tensor_nx, tensor_ny)
    real(dp) :: direction(3), parallel_gradient(3)
    real(dp) :: perpendicular_gradient(3), direction_gradient(3, 3)
    real(dp) :: pressure_tensor(3, 3), expected_tensor(3, 3)
    real(dp) :: force_divergence(3), x, y, angle
    real(dp) :: p_parallel, p_perpendicular, angle_gradient
    real(dp) :: maximum_trace_error, maximum_tensor_error
    real(dp) :: vector_norm
    integer :: point, ix, iy, arrow
    type(fortsparse_status_t) :: status
    integer :: unit

    call execute_command_line("mkdir -p "//output_directory)
    call initialize_gallery_sequence()
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

    maximum_tensor_error = 0.0_dp
    do iy = 1, tensor_ny
        y = -1.0_dp + 2.0_dp*real(iy - 1, dp)/real(tensor_ny - 1, dp)
        tensor_y(iy) = y
        do ix = 1, tensor_nx
            x = -1.0_dp + 2.0_dp*real(ix - 1, dp)/real(tensor_nx - 1, dp)
            tensor_x(ix) = x
            angle = 0.35_dp*x + 0.22_dp*y
            p_parallel = 3.0_dp + 0.2_dp*x - 0.12_dp*y
            p_perpendicular = 1.0_dp + 0.1_dp*cos(pi*x)*cos(0.5_dp*pi*y)
            direction = [cos(angle), sin(angle), 0.0_dp]

            call evaluate_cgl_pressure_tensor( &
                p_parallel, p_perpendicular, direction, pressure_tensor, status)
            if (status%code /= 0) error stop "CGL 2D tensor evaluation failed"
            expected_tensor = p_perpendicular*identity_tensor()
            expected_tensor = expected_tensor + &
                (p_parallel - p_perpendicular)*outer_product(direction, direction)
            maximum_tensor_error = max(maximum_tensor_error, maxval(abs( &
                pressure_tensor - expected_tensor)))
            sampled_p11(ix, iy) = pressure_tensor(1, 1)
            sampled_p22(ix, iy) = pressure_tensor(2, 2)
            sampled_p12(ix, iy) = pressure_tensor(1, 2)
            anisotropy_map(ix, iy) = p_parallel - p_perpendicular
            sampled_trace_error(ix, iy) = pressure_tensor(1, 1) + &
                pressure_tensor(2, 2) + pressure_tensor(3, 3) - &
                p_parallel - 2.0_dp*p_perpendicular
        end do
    end do
    if (maximum_tensor_error > 2.0e-13_dp) &
        error stop "CGL 2D tensor sampled oracle failed"

    arrow = 0
    do iy = 1, arrow_ny
        y = -1.0_dp + 2.0_dp*real(iy - 1, dp)/real(arrow_ny - 1, dp)
        do ix = 1, arrow_nx
            x = -1.0_dp + 2.0_dp*real(ix - 1, dp)/real(arrow_nx - 1, dp)
            arrow = arrow + 1
            angle = 0.35_dp*x + 0.22_dp*y
            principal_x(arrow) = x
            principal_y(arrow) = y
            principal_u(arrow) = cos(angle)
            principal_v(arrow) = sin(angle)
            vector_norm = sqrt(principal_u(arrow)**2 + principal_v(arrow)**2)
            principal_u(arrow) = 0.08_dp*principal_u(arrow)/vector_norm
            principal_v(arrow) = 0.08_dp*principal_v(arrow)/vector_norm
        end do
    end do

    ! Physical state first: tensor anisotropy with its principal direction.
    call figure(figsize=[8.0_dp, 6.2_dp])
    call contourf(tensor_x, tensor_y, anisotropy_map, cmap="coolwarm", &
        show_colorbar=.true.)
    call colorbar(label="p_parallel - p_perpendicular")
    call quiver(principal_x, principal_y, principal_u, principal_v, scale=1.0_dp, &
        scale_units="xy", angles="xy", color="black", width=0.0025_dp, &
        headwidth=3.0_dp)
    call xlabel("x")
    call ylabel("y")
    call title("CGL tensor anisotropy and principal direction")
    call savefig(output_directory//"/cgl_tensor_principal_2d.png")
    call record_gallery_stage("physical_solution")

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
    call record_gallery_stage("diagnostics")

    open (newunit=unit, file=output_directory//"/cgl_profile.csv", &
        status="replace", action="write")
    write (unit, "(a)") "x,P11,P22,P12,divP_x,divP_y"
    do point = 1, point_count
        write (unit, "(6(es24.16,1x))") &
            coordinate(point), pressure_11(point), pressure_22(point), &
            pressure_12(point), force_x(point), force_y(point)
    end do
    close (unit)

    call write_tensor_field_csv()
contains

    pure function identity_tensor() result(identity)
        real(dp) :: identity(3, 3)

        identity = 0.0_dp
        identity(1, 1) = 1.0_dp
        identity(2, 2) = 1.0_dp
        identity(3, 3) = 1.0_dp
    end function identity_tensor

    pure function outer_product(left, right) result(product)
        real(dp), intent(in) :: left(3), right(3)
        real(dp) :: product(3, 3)

        product = spread(left, dim=2, ncopies=3)* &
            spread(right, dim=1, ncopies=3)
    end function outer_product

    subroutine initialize_gallery_sequence()
        open (newunit=unit, file=output_directory//"/gallery_sequence.txt", &
            status="replace", action="write")
        close (unit)
    end subroutine initialize_gallery_sequence

    subroutine record_gallery_stage(stage)
        character(*), intent(in) :: stage

        open (newunit=unit, file=output_directory//"/gallery_sequence.txt", &
            status="old", position="append", action="write")
        write (unit, "(a)") stage
        close (unit)
    end subroutine record_gallery_stage

    subroutine write_tensor_field_csv()
        integer :: local_unit, local_ix, local_iy
        real(dp) :: local_x, local_y, local_angle
        real(dp) :: local_bx, local_by

        open (newunit=local_unit, &
            file=output_directory//"/cgl_tensor_field_2d.csv", &
            status="replace", action="write")
        write (local_unit, "(a)") &
            "x,y,P11,P22,P12,principal_x,principal_y,trace_error"
        do local_iy = 1, tensor_ny
            local_y = tensor_y(local_iy)
            do local_ix = 1, tensor_nx
                local_x = tensor_x(local_ix)
                local_angle = 0.35_dp*local_x + 0.22_dp*local_y
                local_bx = cos(local_angle)
                local_by = sin(local_angle)
                write (local_unit, "(*(es24.16,:,','))") local_x, local_y, &
                    sampled_p11(local_ix, local_iy), &
                    sampled_p22(local_ix, local_iy), &
                    sampled_p12(local_ix, local_iy), local_bx, local_by, &
                    sampled_trace_error(local_ix, local_iy)
            end do
        end do
        close (local_unit)
    end subroutine write_tensor_field_csv
end program cgl_pressure_tensor
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/cgl_pressure_tensor/primary.png)

### cgl_force_divergence_1d.png

![cgl_force_divergence_1d.png](../../media/examples/cgl_pressure_tensor/cgl_force_divergence_1d.png)

### cgl_pressure_components_1d.png

![cgl_pressure_components_1d.png](../../media/examples/cgl_pressure_tensor/cgl_pressure_components_1d.png)

### cgl_tensor_principal_2d.png

![cgl_tensor_principal_2d.png](../../media/examples/cgl_pressure_tensor/cgl_tensor_principal_2d.png)

---

[← Back to all examples](../index.html)
