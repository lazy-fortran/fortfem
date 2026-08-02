---
title: nested_surface_solution Example
---

# nested_surface_solution Example

# Nested toroidal surface solution

This physical-first gallery fixture evaluates FortFEM's neutral nested-surface
map for a circular torus and colors the physical surface with a smooth
manufactured scalar solution. Dark parameter lines make both periodic
directions visible. The geometry is produced by the same Fourier/radial map
available to equilibrium-code adapters; the example does not select an
equilibrium model, coordinate convention, or file format.

The first plot is the physical three-dimensional solution. The later plots
show its periodic parameter-space representation and the axis-regular radial
profile used to modulate it. The program checks the evaluated Cartesian
coordinates against the closed-form circular-torus embedding before plotting.

Generated artifacts:

- `nested_surface_solution_3d.png`: scalar solution on the physical torus;
- `nested_surface_solution_parameter_2d.png`: periodic surface field in
  `(theta,zeta)` coordinates;
- `nested_surface_radial_profile_1d.png`: axis-regular radial amplitude;
- `benchmark.txt`: sample counts, oracle error, timings, and provenance.

Generated images and data remain under `output/example/` and are not committed.

## Usage

```bash
fpm run --example nested_surface_solution
```

## Source Code

```fortran
program nested_surface_solution
    !! Physical-first manufactured field on a neutral nested toroidal map.
    use fortfem_core, only: &
        evaluate_axis_regular_radial_basis, evaluate_nested_surface_geometry
    use fortfem_fourier, only: &
        fourier_mode_registry_t, initialize_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortplot, only: &
        add_3d_plot, add_parametric_surface, add_scatter, colorbar, figure, &
        pcolormesh, plot, savefig, title, view_init, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: major_radius = 2.4_dp
    real(dp), parameter :: minor_radius = 0.7_dp
    real(dp), parameter :: geometry_tolerance = 5.0e-12_dp
    real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
    integer, parameter :: theta_count = 49, zeta_count = 73
    integer, parameter :: surface_count = theta_count*zeta_count
    integer, parameter :: radial_count = 81
    integer, parameter :: theta_line_stride = 12, zeta_line_stride = 18
    character(*), parameter :: output_directory = &
        "output/example/nested_surface_solution"
    type(fourier_mode_registry_t) :: registry
    type(fortsparse_status_t) :: status
    complex(dp) :: geometry_coefficients(3, 2)
    complex(dp) :: radial_coefficients(2)
    complex(dp) :: surface_radial(surface_count), radial_values(radial_count)
    real(dp) :: rho(surface_count), theta(surface_count), zeta(surface_count)
    real(dp) :: mapped(3, surface_count), physical(3, surface_count)
    real(dp) :: mapped_jacobian(3, 3, surface_count)
    real(dp) :: physical_jacobian(3, 3, surface_count)
    real(dp) :: metric(3, 3, surface_count), volume(surface_count)
    real(dp) :: surface_x(theta_count, zeta_count)
    real(dp) :: surface_y(theta_count, zeta_count)
    real(dp) :: surface_z(theta_count, zeta_count)
    real(dp) :: solution(theta_count, zeta_count)
    real(dp) :: theta_edges(theta_count), zeta_edges(zeta_count)
    real(dp) :: parameter_solution(zeta_count - 1, theta_count - 1)
    real(dp) :: radial_coordinate(radial_count), radial_profile(radial_count)
    real(dp) :: geometry_start, geometry_seconds, render_start, render_seconds
    real(dp) :: radius, expected(3), maximum_geometry_error
    integer :: i, j, point, unit

    call execute_command_line("mkdir -p "//output_directory)
    call initialize_fourier_mode_registry( &
        registry, [0, 1], [0, 0], 1, 0.0_dp, 0.0_dp, .false., &
        radial_powers=[0, 1], status=status)
    if (status%code /= 0) error stop "nested-surface registry failed"
    geometry_coefficients = cmplx(0.0_dp, 0.0_dp, dp)
    geometry_coefficients(1, 1) = cmplx(major_radius, 0.0_dp, dp)
    geometry_coefficients(1, 2) = cmplx(minor_radius, 0.0_dp, dp)
    geometry_coefficients(2, 2) = cmplx(0.0_dp, -minor_radius, dp)
    radial_coefficients = [ &
        cmplx(0.25_dp, 0.0_dp, dp), cmplx(0.75_dp, 0.0_dp, dp)]

    point = 0
    do j = 1, zeta_count
        zeta_edges(j) = 2.0_dp*pi*real(j - 1, dp)/real(zeta_count - 1, dp)
        do i = 1, theta_count
            if (j == 1) then
                theta_edges(i) = &
                    2.0_dp*pi*real(i - 1, dp)/real(theta_count - 1, dp)
            end if
            point = point + 1
            rho(point) = 1.0_dp
            theta(point) = theta_edges(i)
            zeta(point) = zeta_edges(j)
        end do
    end do

    call cpu_time(geometry_start)
    call evaluate_nested_surface_geometry( &
        registry, geometry_coefficients, rho, theta, zeta, mapped, physical, &
        mapped_jacobian, physical_jacobian, metric, volume, status)
    call cpu_time(geometry_seconds)
    geometry_seconds = geometry_seconds - geometry_start
    if (status%code /= 0) error stop "nested-surface geometry evaluation failed"
    call evaluate_axis_regular_radial_basis( &
        0, [0, 2], radial_coefficients, rho, surface_radial, status)
    if (status%code /= 0) error stop "nested-surface radial basis failed"

    maximum_geometry_error = 0.0_dp
    do point = 1, surface_count
        radius = major_radius + minor_radius*rho(point)*cos(theta(point))
        expected = [ &
            radius*cos(zeta(point)), radius*sin(zeta(point)), &
            minor_radius*rho(point)*sin(theta(point))]
        maximum_geometry_error = max( &
            maximum_geometry_error, maxval(abs(physical(:, point) - expected)))
    end do
    if (maximum_geometry_error > geometry_tolerance) &
        error stop "nested-surface analytical geometry oracle failed"

    point = 0
    do j = 1, zeta_count
        do i = 1, theta_count
            point = point + 1
            surface_x(i, j) = physical(1, point)
            surface_y(i, j) = physical(2, point)
            surface_z(i, j) = physical(3, point)
            solution(i, j) = real(surface_radial(point), dp)*( &
                1.0_dp + 0.30_dp*cos(theta(point) - 2.0_dp*zeta(point)) + &
                0.12_dp*cos(2.0_dp*theta(point) + zeta(point)))
        end do
    end do
    do j = 1, zeta_count - 1
        do i = 1, theta_count - 1
            parameter_solution(j, i) = solution(i, j)
        end do
    end do

    do i = 1, radial_count
        radial_coordinate(i) = real(i - 1, dp)/real(radial_count - 1, dp)
    end do
    call evaluate_axis_regular_radial_basis( &
        0, [0, 2], radial_coefficients, radial_coordinate, radial_values, status)
    if (status%code /= 0) error stop "nested-surface radial profile failed"
    radial_profile = real(radial_values, dp)

    call cpu_time(render_start)
    call render_physical_solution()
    call render_parameter_solution()
    call render_radial_profile()
    call cpu_time(render_seconds)
    render_seconds = render_seconds - render_start
    call write_metadata()

contains

    subroutine render_physical_solution()
        real(dp) :: line_x(max(theta_count, zeta_count))
        real(dp) :: line_y(max(theta_count, zeta_count))
        real(dp) :: line_z(max(theta_count, zeta_count))
        integer :: line

        call figure(figsize=[9.0_dp, 7.0_dp])
        call add_parametric_surface( &
            surface_x, surface_y, surface_z, color="lightgray", &
            alpha=0.22_dp, linewidth=0.0_dp, filled=.true.)
        call add_scatter( &
            reshape(surface_x, [surface_count]), &
            reshape(surface_y, [surface_count]), &
            reshape(surface_z, [surface_count]), &
            c=reshape(solution, [surface_count]), cmap="viridis", &
            marker=".", markersize=3.5_dp)
        do line = 1, theta_count, theta_line_stride
            line_x(:zeta_count) = surface_x(line, :)
            line_y(:zeta_count) = surface_y(line, :)
            line_z(:zeta_count) = surface_z(line, :)
            call add_3d_plot( &
                line_x(:zeta_count), line_y(:zeta_count), line_z(:zeta_count), &
                color="black", linewidth=0.28_dp)
        end do
        do line = 1, zeta_count, zeta_line_stride
            line_x(:theta_count) = surface_x(:, line)
            line_y(:theta_count) = surface_y(:, line)
            line_z(:theta_count) = surface_z(:, line)
            call add_3d_plot( &
                line_x(:theta_count), line_y(:theta_count), line_z(:theta_count), &
                color="black", linewidth=0.28_dp)
        end do
        call colorbar(label="manufactured scalar solution u [1]")
        call xlabel("x [m]")
        call ylabel("y [m]")
        call title("Nested toroidal surface solution (vertical coordinate z [m])")
        call view_init(elev=27.0_dp, azim=-52.0_dp)
        call savefig(output_directory//"/nested_surface_solution_3d.png")
    end subroutine render_physical_solution

    subroutine render_parameter_solution()
        call figure(figsize=[8.5_dp, 5.8_dp])
        call pcolormesh( &
            theta_edges, zeta_edges, parameter_solution, cmap="viridis")
        call colorbar(label="manufactured scalar solution u [1]")
        call xlabel("poloidal angle theta [rad]")
        call ylabel("toroidal angle zeta [rad]")
        call title("Periodic parameter-space view of the same solution")
        call savefig( &
            output_directory//"/nested_surface_solution_parameter_2d.png")
    end subroutine render_parameter_solution

    subroutine render_radial_profile()
        call figure(figsize=[8.0_dp, 5.2_dp])
        call plot( &
            radial_coordinate, radial_profile, color=blue, linewidth=2.0_dp)
        call xlabel("normalized nested-surface radius rho [1]")
        call ylabel("axis-regular amplitude A(rho) [1]")
        call title("Axis-regular radial basis: A = 0.25 + 0.75 rho^2")
        call savefig(output_directory//"/nested_surface_radial_profile_1d.png")
    end subroutine render_radial_profile

    subroutine write_metadata()
        open (newunit=unit, file=output_directory//"/benchmark.txt", &
            status="replace", action="write")
        write (unit, "(a)") "quantity,value"
        write (unit, "(a,i0)") "surface_samples,", surface_count
        write (unit, "(a,es24.16)") &
            "maximum_geometry_error_m,", maximum_geometry_error
        write (unit, "(a,es24.16)") &
            "geometry_evaluation_seconds,", geometry_seconds
        write (unit, "(a,es24.16)") "plot_render_seconds,", render_seconds
        write (unit, "(a)") &
            "provenance,manufactured circular torus R=2.4 m a=0.7 m"
        write (unit, "(a)") &
            "geometry_api,evaluate_nested_surface_geometry"
        write (unit, "(a)") &
            "radial_api,evaluate_axis_regular_radial_basis"
        close (unit)
    end subroutine write_metadata

end program nested_surface_solution
```

## Generated Plots

### primary.png

![primary.png](../../media/examples/nested_surface_solution/primary.png)

### nested_surface_radial_profile_1d.png

![nested_surface_radial_profile_1d.png](../../media/examples/nested_surface_solution/nested_surface_radial_profile_1d.png)

### nested_surface_solution_3d.png

![nested_surface_solution_3d.png](../../media/examples/nested_surface_solution/nested_surface_solution_3d.png)

### nested_surface_solution_parameter_2d.png

![nested_surface_solution_parameter_2d.png](../../media/examples/nested_surface_solution/nested_surface_solution_parameter_2d.png)

---

[← Back to all examples](../index.html)
