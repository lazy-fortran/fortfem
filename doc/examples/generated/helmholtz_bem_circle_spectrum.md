---
title: helmholtz_bem_circle_spectrum Example
---

# helmholtz_bem_circle_spectrum Example

Executable FortFEM helmholtz_bem_circle_spectrum.f90 example.

## Usage

```bash
fpm run --example helmholtz_bem_circle_spectrum
```

## Source Code

```fortran
program helmholtz_bem_circle_spectrum
    use fortfem_boundary, only: assemble_helmholtz_double_layer_constant, &
        assemble_helmholtz_hypersingular_linear, &
        assemble_helmholtz_single_layer_constant
    use fortfem_kinds, only: dp
    use fortplot, only: add_scatter, colorbar, figure, legend, pcolormesh, plot, &
        savefig, set_xscale, set_yscale, title, xlabel, ylabel
    implicit none

    integer, parameter :: mesh_sizes(3) = [12, 24, 48]
    complex(dp), parameter :: exact_single_layer = &
        (0.4497818998782933_dp, 0.4280549908588178_dp)
    complex(dp), parameter :: exact_double_layer = &
        (-0.25522568497710474_dp, 0.2329503859714022_dp)
    complex(dp), parameter :: exact_hypersingular = &
        (0.41099886362254484_dp, -0.1267731564473764_dp)
    complex(dp) :: double_layer, hypersingular, single_layer
    real(dp) :: single_errors(3), double_errors(3), hypersingular_errors(3)
    real(dp) :: mesh_sizes_real(3)
    integer :: command_status, mesh_id

    call execute_command_line( &
        "mkdir -p output/example/helmholtz_bem_circle_spectrum", &
        exitstat=command_status)
    if (command_status /= 0) error stop "cannot create example output directory"

    print "(a)", "Unit-circle outgoing Helmholtz BEM spectrum, k = 1.3, mode 1"
    print "(a)", " panels     rel(V)      rel(K)      rel(W)"
    do mesh_id = 1, size(mesh_sizes)
        call circle_mode_eigenvalues( &
            mesh_sizes(mesh_id), single_layer, double_layer, hypersingular)
        mesh_sizes_real(mesh_id) = real(mesh_sizes(mesh_id), dp)
        single_errors(mesh_id) = relative_error(single_layer, exact_single_layer)
        double_errors(mesh_id) = relative_error(double_layer, exact_double_layer)
        hypersingular_errors(mesh_id) = &
            relative_error(hypersingular, exact_hypersingular)
        print "(i7,3es12.3)", mesh_sizes(mesh_id), &
            single_errors(mesh_id), double_errors(mesh_id), &
            hypersingular_errors(mesh_id)
    end do

    call render_disk_field()
    call render_boundary_mode()

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(mesh_sizes_real, single_errors, label="single layer", marker="o")
    call plot(mesh_sizes_real, double_errors, label="double layer", marker="s")
    call plot(mesh_sizes_real, hypersingular_errors, label="hypersingular", &
        marker="^")
    call set_xscale("log")
    call set_yscale("log")
    call xlabel("circle panels")
    call ylabel("relative spectral error")
    call title("Helmholtz circle BEM spectral convergence")
    call legend()
    call savefig( &
        "output/example/helmholtz_bem_circle_spectrum/convergence_1d.png")

contains

    subroutine render_disk_field()
        integer, parameter :: field_nx = 64, field_ny = 64
        real(dp), parameter :: wave_number = 1.3_dp
        real(dp) :: x_edges(field_nx + 1), y_edges(field_ny + 1)
        real(dp) :: values(field_ny, field_nx)
        real(dp) :: circle_x(field_nx + 1), circle_y(field_nx + 1)
        real(dp) :: x_value, y_value, radius
        integer :: i, j

        do i = 1, field_nx + 1
            x_edges(i) = -1.0_dp + 2.0_dp*real(i - 1, dp)/field_nx
            circle_x(i) = cos(2.0_dp*acos(-1.0_dp)*real(i - 1, dp)/field_nx)
            circle_y(i) = sin(2.0_dp*acos(-1.0_dp)*real(i - 1, dp)/field_nx)
        end do
        y_edges = x_edges
        do j = 1, field_ny
            y_value = 0.5_dp*(y_edges(j) + y_edges(j + 1))
            do i = 1, field_nx
                x_value = 0.5_dp*(x_edges(i) + x_edges(i + 1))
                radius = sqrt(x_value*x_value + y_value*y_value)
                if (radius <= 1.0_dp) then
                    values(j, i) = cos(wave_number*x_value)
                else
                    values(j, i) = 0.0_dp
                end if
            end do
        end do

        call figure(figsize=[8.0_dp, 6.5_dp])
        call pcolormesh(x_edges, y_edges, values, cmap="coolwarm")
        call colorbar(label="Re(exp(i k x))")
        call plot(circle_x, circle_y, color=[0.0_dp, 0.0_dp, 0.0_dp], &
            linewidth=1.5_dp)
        call xlabel("x")
        call ylabel("y")
        call title("Helmholtz BEM incident field inside the circle")
        call savefig( &
            "output/example/helmholtz_bem_circle_spectrum/helmholtz_disk_field_2d.png")
    end subroutine render_disk_field

    subroutine render_boundary_mode()
        integer, parameter :: point_count = 128
        real(dp) :: angle, mode(point_count)
        real(dp) :: x(point_count), y(point_count)
        integer :: point

        do point = 1, point_count
            angle = 2.0_dp*acos(-1.0_dp)*real(point - 1, dp)/ &
                real(point_count, dp)
            x(point) = cos(angle)
            y(point) = sin(angle)
            mode(point) = cos(angle)
        end do
        call figure(figsize=[7.0_dp, 6.5_dp])
        call add_scatter(x, y, c=mode, cmap="coolwarm", marker=".", &
            markersize=5.0_dp, label="cos(theta) spectral mode")
        call title("Helmholtz circle BEM boundary mode")
        call savefig( &
            "output/example/helmholtz_bem_circle_spectrum/boundary_mode_2d.png")
    end subroutine render_boundary_mode

    subroutine circle_mode_eigenvalues( &
            panel_count, single_layer_eigenvalue, double_layer_eigenvalue, &
            hypersingular_eigenvalue)
        integer, intent(in) :: panel_count
        complex(dp), intent(out) :: single_layer_eigenvalue
        complex(dp), intent(out) :: double_layer_eigenvalue
        complex(dp), intent(out) :: hypersingular_eigenvalue

        complex(dp), allocatable :: double_layer_matrix(:, :)
        complex(dp), allocatable :: hypersingular_matrix(:, :)
        complex(dp), allocatable :: single_layer_matrix(:, :)
        integer, allocatable :: panel_nodes(:, :)
        real(dp), allocatable :: nodal_mode(:), panel_end(:, :)
        real(dp), allocatable :: panel_mode(:), panel_start(:, :)
        real(dp) :: angle, length, pi
        integer :: panel, status

        allocate(panel_start(2, panel_count), panel_end(2, panel_count))
        allocate(panel_nodes(2, panel_count))
        allocate(single_layer_matrix(panel_count, panel_count))
        allocate(double_layer_matrix(panel_count, panel_count))
        allocate(hypersingular_matrix(panel_count, panel_count))
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

        call assemble_helmholtz_single_layer_constant( &
            panel_start, panel_end, 1.3_dp, 16, single_layer_matrix, status)
        if (status /= 0) error stop "Single-layer assembly failed"
        call assemble_helmholtz_double_layer_constant( &
            panel_start, panel_end, 1.3_dp, 16, double_layer_matrix, status)
        if (status /= 0) error stop "Double-layer assembly failed"
        call assemble_helmholtz_hypersingular_linear( &
            panel_start, panel_end, panel_nodes, panel_count, 1.3_dp, 16, &
            hypersingular_matrix, status)
        if (status /= 0) error stop "Hypersingular assembly failed"

        single_layer_eigenvalue = quadratic_form( &
            single_layer_matrix, panel_mode) / &
            (length * dot_product(panel_mode, panel_mode))
        double_layer_eigenvalue = quadratic_form( &
            double_layer_matrix, panel_mode) / &
            (length * dot_product(panel_mode, panel_mode))
        hypersingular_eigenvalue = quadratic_form( &
            hypersingular_matrix, nodal_mode) / &
            linear_mass_norm(nodal_mode, length)
    end subroutine circle_mode_eigenvalues

    pure function quadratic_form(matrix, values) result(value)
        complex(dp), intent(in) :: matrix(:, :)
        real(dp), intent(in) :: values(:)
        complex(dp) :: value

        integer :: column, row

        value = (0.0_dp, 0.0_dp)
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

    pure function relative_error(value, reference) result(error)
        complex(dp), intent(in) :: value, reference
        real(dp) :: error

        error = abs(value - reference) / abs(reference)
    end function relative_error

end program helmholtz_bem_circle_spectrum
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/helmholtz_bem_circle_spectrum/primary.png)

### boundary_mode_2d.png

![boundary_mode_2d.png](../../media/examples/helmholtz_bem_circle_spectrum/boundary_mode_2d.png)

### convergence_1d.png

![convergence_1d.png](../../media/examples/helmholtz_bem_circle_spectrum/convergence_1d.png)

### helmholtz_disk_field_2d.png

![helmholtz_disk_field_2d.png](../../media/examples/helmholtz_bem_circle_spectrum/helmholtz_disk_field_2d.png)

---

[← Back to all examples](../index.html)
