---
title: helmholtz_cfie_circle Example
---

# helmholtz_cfie_circle Example

Executable FortFEM helmholtz_cfie_circle.f90 example.

## Usage

```bash
fpm run --example helmholtz_cfie_circle
```

## Source Code

```fortran
program helmholtz_cfie_circle
    use fortfem_api, only: &
        evaluate_helmholtz_combined_potential_constant, &
        solve_helmholtz_cfie_constant
    use fortfem_kinds, only: dp
    use fortplot, only: add_scatter, colorbar, figure, legend, pcolormesh, plot, &
        savefig, set_xscale, set_yscale, title, xlabel, ylabel
    implicit none

    integer, parameter :: mesh_sizes(3) = [16, 32, 64]
    complex(dp), parameter :: exact_scattered = &
        (0.6192993117265686_dp, -0.5636220237912836_dp)
    complex(dp) :: scattered
    real(dp) :: residual, field_errors(3), residuals(3), mesh_sizes_real(3)
    integer :: command_status, mesh_id

    call execute_command_line( &
        "mkdir -p output/example/helmholtz_cfie_circle", &
        exitstat=command_status)
    if (command_status /= 0) error stop "cannot create example output directory"

    print "(a)", "Sound-soft unit-circle scattering, k = eta = 1.3"
    print "(a)", " panels    field error      LU residual"
    do mesh_id = 1, size(mesh_sizes)
        call solve_circle( &
            mesh_sizes(mesh_id), scattered, residual)
        mesh_sizes_real(mesh_id) = real(mesh_sizes(mesh_id), dp)
        field_errors(mesh_id) = abs(scattered - exact_scattered)
        residuals(mesh_id) = residual
        print "(i7,2es16.6)", mesh_sizes(mesh_id), &
            field_errors(mesh_id), residuals(mesh_id)
    end do
    print "(a,2f12.6)", "Mie scattered field at r=2, theta=0.3: ", &
        real(exact_scattered, dp), aimag(exact_scattered)

    call render_field()
    call render_incident_boundary()

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(mesh_sizes_real, field_errors, label="scattered-field error", &
        marker="o")
    call plot(mesh_sizes_real, residuals, label="CFIE residual", marker="s")
    call set_xscale("log")
    call set_yscale("log")
    call xlabel("circle panels")
    call ylabel("error / residual")
    call title("Helmholtz combined-field circle convergence")
    call legend()
    call savefig("output/example/helmholtz_cfie_circle/convergence_1d.png")

contains

    subroutine render_field()
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
        call colorbar(label="Re(incident field)")
        call plot(circle_x, circle_y, color=[0.0_dp, 0.0_dp, 0.0_dp], &
            linewidth=1.5_dp)
        call xlabel("x")
        call ylabel("y")
        call title("Helmholtz CFIE incident field in the disk")
        call savefig( &
            "output/example/helmholtz_cfie_circle/helmholtz_cfie_field_2d.png")
    end subroutine render_field

    subroutine render_incident_boundary()
        integer, parameter :: point_count = 128
        real(dp) :: angle, x(point_count), y(point_count), mode(point_count)
        integer :: point

        do point = 1, point_count
            angle = 2.0_dp*acos(-1.0_dp)*real(point - 1, dp)/ &
                real(point_count, dp)
            x(point) = cos(angle)
            y(point) = sin(angle)
            mode(point) = real( &
                exp(cmplx(0.0_dp, 1.3_dp*x(point), dp)), dp)
        end do
        call figure(figsize=[7.0_dp, 6.5_dp])
        call add_scatter(x, y, c=mode, cmap="coolwarm", marker=".", &
            markersize=5.0_dp, label="Re incident trace")
        call title("Helmholtz CFIE incident field on the circle")
        call savefig( &
            "output/example/helmholtz_cfie_circle/incident_boundary_2d.png")
    end subroutine render_incident_boundary

    subroutine solve_circle(panel_count, scattered, residual)
        integer, intent(in) :: panel_count
        complex(dp), intent(out) :: scattered
        real(dp), intent(out) :: residual

        complex(dp), allocatable :: density(:), incident_trace(:)
        real(dp), allocatable :: panel_end(:, :), panel_start(:, :)
        complex(dp) :: field(1)
        real(dp) :: angle, midpoint(2), pi, point(2, 1)
        integer :: panel, status

        allocate(panel_start(2, panel_count), panel_end(2, panel_count))
        allocate(density(panel_count), incident_trace(panel_count))
        pi = acos(-1.0_dp)
        do panel = 1, panel_count
            angle = 2.0_dp * pi * real(panel - 1, dp) / real(panel_count, dp)
            panel_start(1, panel) = cos(angle)
            panel_start(2, panel) = sin(angle)
            angle = 2.0_dp * pi * real(panel, dp) / real(panel_count, dp)
            panel_end(1, panel) = cos(angle)
            panel_end(2, panel) = sin(angle)
            midpoint = 0.5_dp * ( &
                panel_start(:, panel) + panel_end(:, panel))
            incident_trace(panel) = exp(cmplx(0.0_dp, 1.3_dp * midpoint(1), dp))
        end do

        call solve_helmholtz_cfie_constant( &
            panel_start, panel_end, 1.3_dp, 1.3_dp, incident_trace, 16, &
            density, residual, status)
        if (status /= 0) error stop "Combined-field solve failed"

        point(1, 1) = 2.0_dp * cos(0.3_dp)
        point(2, 1) = 2.0_dp * sin(0.3_dp)
        call evaluate_helmholtz_combined_potential_constant( &
            panel_start, panel_end, 1.3_dp, 1.3_dp, density, point, 16, &
            field, status)
        if (status /= 0) error stop "Combined potential evaluation failed"
        scattered = field(1)
    end subroutine solve_circle

end program helmholtz_cfie_circle
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/helmholtz_cfie_circle/primary.png)

### convergence_1d.png

![convergence_1d.png](../../media/examples/helmholtz_cfie_circle/convergence_1d.png)

### helmholtz_cfie_field_2d.png

![helmholtz_cfie_field_2d.png](../../media/examples/helmholtz_cfie_circle/helmholtz_cfie_field_2d.png)

### incident_boundary_2d.png

![incident_boundary_2d.png](../../media/examples/helmholtz_cfie_circle/incident_boundary_2d.png)

---

[← Back to all examples](../index.html)
