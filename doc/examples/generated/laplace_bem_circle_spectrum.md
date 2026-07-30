---
title: laplace_bem_circle_spectrum Example
---

# laplace_bem_circle_spectrum Example

Executable FortFEM laplace_bem_circle_spectrum.f90 example.

## Usage

```bash
fpm run --example laplace_bem_circle_spectrum
```

## Source Code

```fortran
program laplace_bem_circle_spectrum
    use fortfem_api, only: assemble_laplace_hypersingular_linear, &
        assemble_laplace_single_layer_constant
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: mesh_sizes(4) = [12, 24, 48, 96]
    real(dp) :: single_layer_eigenvalue, hypersingular_eigenvalue
    integer :: mesh_id

    print "(a)", "Unit-circle Laplace BEM spectrum for cos(theta)"
    print "(a)", " panels       V_h       error       W_h       error"
    do mesh_id = 1, size(mesh_sizes)
        call circle_mode_eigenvalues( &
            mesh_sizes(mesh_id), single_layer_eigenvalue, &
            hypersingular_eigenvalue)
        print "(i7,4f12.6)", mesh_sizes(mesh_id), &
            single_layer_eigenvalue, abs(single_layer_eigenvalue - 0.5_dp), &
            hypersingular_eigenvalue, &
            abs(hypersingular_eigenvalue - 0.5_dp)
    end do
    print "(a)", "Exact eigenvalues: V_1 = W_1 = 1/2"

contains

    subroutine circle_mode_eigenvalues( &
            panel_count, single_layer_eigenvalue, hypersingular_eigenvalue)
        integer, intent(in) :: panel_count
        real(dp), intent(out) :: single_layer_eigenvalue
        real(dp), intent(out) :: hypersingular_eigenvalue

        integer, allocatable :: panel_nodes(:, :)
        real(dp), allocatable :: hypersingular(:, :), nodal_mode(:)
        real(dp), allocatable :: panel_end(:, :), panel_mode(:)
        real(dp), allocatable :: panel_start(:, :), single_layer(:, :)
        real(dp) :: angle, length, pi
        integer :: panel, status

        allocate(panel_start(2, panel_count), panel_end(2, panel_count))
        allocate(panel_nodes(2, panel_count))
        allocate(single_layer(panel_count, panel_count))
        allocate(hypersingular(panel_count, panel_count))
        allocate(panel_mode(panel_count), nodal_mode(panel_count))

        pi = acos(-1.0_dp)
        do panel = 1, panel_count
            angle = 2.0_dp * pi * real(panel - 1, dp) / real(panel_count, dp)
            panel_start(1, panel) = cos(angle)
            panel_start(2, panel) = sin(angle)
            nodal_mode(panel) = cos(angle)
            angle = 2.0_dp * pi * real(panel, dp) / real(panel_count, dp)
            panel_end(1, panel) = cos(angle)
            panel_end(2, panel) = sin(angle)
            panel_mode(panel) = cos( &
                2.0_dp * pi * (real(panel, dp) - 0.5_dp) / &
                real(panel_count, dp))
            panel_nodes(1, panel) = panel
            panel_nodes(2, panel) = modulo(panel, panel_count) + 1
        end do
        length = norm2(panel_end(:, 1) - panel_start(:, 1))

        call assemble_laplace_single_layer_constant( &
            panel_start, panel_end, 20, single_layer, status)
        if (status /= 0) error stop "Single-layer assembly failed"
        call assemble_laplace_hypersingular_linear( &
            panel_start, panel_end, panel_nodes, panel_count, 20, &
            hypersingular, status)
        if (status /= 0) error stop "Hypersingular assembly failed"

        single_layer_eigenvalue = quadratic_form(single_layer, panel_mode) / &
            (length * dot_product(panel_mode, panel_mode))
        hypersingular_eigenvalue = quadratic_form(hypersingular, nodal_mode) / &
            linear_mass_norm(nodal_mode, length)
    end subroutine circle_mode_eigenvalues

    pure function quadratic_form(matrix, values) result(value)
        real(dp), intent(in) :: matrix(:, :), values(:)
        real(dp) :: value

        integer :: column, row

        value = 0.0_dp
        do column = 1, size(values)
            do row = 1, size(values)
                value = value + &
                    values(row) * matrix(row, column) * values(column)
            end do
        end do
    end function quadratic_form

    pure function linear_mass_norm(values, panel_length) result(norm)
        real(dp), intent(in) :: values(:), panel_length
        real(dp) :: norm

        integer :: panel, successor

        norm = 0.0_dp
        do panel = 1, size(values)
            successor = modulo(panel, size(values)) + 1
            norm = norm + panel_length / 3.0_dp * ( &
                values(panel)**2 + values(panel) * values(successor) + &
                values(successor)**2)
        end do
    end function linear_mass_norm

end program laplace_bem_circle_spectrum
```

## Generated Plots

*No plot artifact is produced by this example.*

---

[← Back to all examples](../index.html)
