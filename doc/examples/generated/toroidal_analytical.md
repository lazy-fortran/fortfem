---
title: toroidal_analytical Example
---

# toroidal_analytical Example

# Toroidal Poisson and Ampère: analytical versus BEM

This example evaluates the real \(n=2,m=1\) toroidal harmonic

\[
\Psi=\sqrt{\cosh\eta-\cos\theta}\,
P_{3/2}^{1}(\cosh\eta)\cos(2\theta)\cos(\phi)
\]

and the current-free magnetostatic field \(\mathbf H=-\nabla\Psi\). The
half-integer Legendre function and derivative come from Fortnum's
Hobson-normalized toroidal special functions.

The numerical path uses exact parametric panels on the constant-\(\eta\)
torus with \(\cosh\eta=2\). It samples only Dirichlet data, solves the dense
Galerkin \(V\phi=(K-\frac12M)u\) system for the unknown P0 Neumann trace, and
then reconstructs the exterior field. The Ampère field is the numerical
gradient of that reconstructed scalar potential. The Helmholtz comparison
uses the corresponding complex curved-torus DtN solve with Laplace
singularity subtraction. Analytical toroidal harmonics and an outgoing point
source are independent oracles.

A second, deliberately small solid-torus run couples a tetrahedral P1 FEM
interior to the same exact-parametric curved BEM surface with the symmetric
Costabel-Han formulation. It reports wall-clock solve time and boundary
error for both Poisson and scalar Helmholtz, making the DtN-only and coupled
paths reproducibly comparable without committing generated media.

The program produces:

- `toroidal_trace_1d.png`: analytical and BEM Poisson/Ampère traces;
- `toroidal_surface_2d.png`: the analytical scalar mode in \((\theta,\phi)\);
- `toroidal_bem_error_2d.png`: BEM relative error over the exterior surface;
- `toroidal_helmholtz_1d.png`: outgoing Helmholtz point source, exact vs BEM;
- `toroidal_fem_bem_1d.png`: analytical vs coupled FEM-BEM boundary traces;
- `toroidal_geometry_3d.png`: the physical toroidal interface;
- `toroidal_trace.csv`: reproducible source data for comparisons.
- `benchmark.txt`: mesh size, separate DtN-solve/evaluation timings, and
  global relative errors.

Rendered outputs are generated in CI and are intentionally not committed.

The separated toroidal harmonic follows
[DLMF §14.19](https://dlmf.nist.gov/14.19). The boundary representation uses
the standard Laplace/Helmholtz Green formulas with outward surface normals.

## Usage

```bash
fpm run --example toroidal_analytical
```

## Source Code

```fortran
program toroidal_analytical
    use fortfem_api, only: &
        cartesian_to_toroidal, &
        evaluate_helmholtz_representation_torus_curved_3d, &
        evaluate_laplace_representation_torus_curved_3d, &
        evaluate_toroidal_harmonic_p, evaluate_toroidal_ampere_field_p, &
        generate_solid_torus_tetra_mesh, generate_torus_surface_mesh, &
        solve_helmholtz_fem_bem_costabel_torus_curved_3d, &
        solve_helmholtz_bem_dtn_torus_curved_3d, &
        solve_laplace_fem_bem_costabel_torus_curved_3d, &
        solve_laplace_bem_dtn_torus_curved_3d, toroidal_point_to_cartesian
    use fortfem_kinds, only: dp
    use fortplot, only: figure, plot, pcolormesh, add_scatter, &
        xlabel, ylabel, title, legend, savefig
    implicit none

    integer, parameter :: trace_points = 61
    integer, parameter :: theta_cells = 24
    integer, parameter :: phi_cells = 16
    integer, parameter :: surface_points = theta_cells*phi_cells
    integer, parameter :: degree_index = 2
    integer, parameter :: order = 1
    real(dp), parameter :: scale = 1.0_dp
    real(dp), parameter :: eta = &
        1.3169578969248167086250463473079684_dp
    real(dp), parameter :: target_eta = 0.35_dp
    real(dp), parameter :: wave_number = 0.2_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(*), parameter :: output_directory = &
        "output/example/toroidal_analytical"
    real(dp) :: theta_trace(trace_points), potential(trace_points)
    real(dp) :: axial_z(trace_points)
    real(dp) :: bem_potential(trace_points)
    real(dp) :: field_norm(trace_points), bem_field_norm(trace_points)
    real(dp) :: theta_edges(theta_cells + 1), phi_edges(phi_cells + 1)
    real(dp) :: surface_value(phi_cells, theta_cells)
    real(dp) :: bem_error(phi_cells, theta_cells)
    real(dp) :: x(surface_points), y(surface_points), z(surface_points)
    real(dp), allocatable :: boundary_trace(:), boundary_flux(:)
    real(dp), allocatable :: parameters(:, :)
    complex(dp), allocatable :: helmholtz_trace(:), helmholtz_flux(:)
    complex(dp), allocatable :: coupled_helmholtz_flux(:)
    complex(dp), allocatable :: coupled_helmholtz_load(:)
    complex(dp), allocatable :: coupled_helmholtz_potential(:)
    complex(dp) :: helmholtz_exact(trace_points)
    complex(dp) :: helmholtz_bem(trace_points)
    integer, allocatable :: volume_tetrahedra(:, :)
    integer, allocatable :: volume_boundary_triangles(:, :)
    real(dp), allocatable :: volume_parameters(:, :), volume_vertices(:, :)
    real(dp), allocatable :: coupled_flux(:), coupled_load(:)
    real(dp), allocatable :: coupled_potential(:)
    real(dp), allocatable :: coupled_theta(:), coupled_exact(:)
    real(dp), allocatable :: coupled_numerical(:)
    complex(dp), allocatable :: coupled_helmholtz_exact(:)
    complex(dp), allocatable :: coupled_helmholtz_numerical(:)
    real(dp), allocatable :: vertices(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp) :: field(3), theta_value, phi_value, denominator, target(3)
    real(dp) :: bem_gradient(3), start_time, end_time
    real(dp) :: bem_seconds, exact_seconds, helmholtz_seconds
    real(dp) :: poisson_solve_seconds, helmholtz_solve_seconds
    real(dp) :: poisson_coupled_seconds, helmholtz_coupled_seconds
    real(dp) :: poisson_coupled_error, helmholtz_coupled_error
    integer :: i, j, point, status, unit

    call execute_command_line("mkdir -p "//output_directory)
    call generate_torus_surface_mesh( &
        2.0_dp/sqrt(3.0_dp), 1.0_dp/sqrt(3.0_dp), 13, 15, &
        vertices, triangles, parameters)
    call build_boundary_data()
    call cpu_time(start_time)
    call solve_laplace_bem_dtn_torus_curved_3d( &
        parameters, triangles, 2.0_dp/sqrt(3.0_dp), 1.0_dp/sqrt(3.0_dp), &
        boundary_trace, 7, boundary_flux, status)
    call cpu_time(end_time)
    if (status /= 0) error stop "Toroidal Laplace DtN solve failed"
    poisson_solve_seconds = end_time - start_time
    call generate_trace()
    call generate_helmholtz_trace()
    call generate_fem_bem_comparison()
    call generate_surface_map()
    call generate_geometry()
    call write_trace_data()

contains

    subroutine generate_trace()
        call cpu_time(start_time)
        do i = 1, trace_points
            theta_trace(i) = 2.0_dp*pi*real(i - 1, dp)/ &
                real(trace_points - 1, dp)
            call evaluate_toroidal_harmonic_p( &
                degree_index, order, target_eta, theta_trace(i), 0.4_dp, &
                potential(i), status)
            if (status /= 0) error stop "Toroidal potential evaluation failed"
            call evaluate_toroidal_ampere_field_p( &
                degree_index, order, scale, target_eta, theta_trace(i), 0.4_dp, &
                field, status)
            if (status /= 0) error stop "Toroidal Ampere evaluation failed"
            field_norm(i) = norm2(field)
        end do
        call cpu_time(end_time)
        exact_seconds = end_time - start_time

        call cpu_time(start_time)
        do i = 1, trace_points
            call toroidal_point_to_cartesian( &
                scale, target_eta, theta_trace(i), 0.4_dp, target)
            call evaluate_laplace_representation_torus_curved_3d( &
                parameters, triangles, 2.0_dp/sqrt(3.0_dp), &
                1.0_dp/sqrt(3.0_dp), boundary_trace, boundary_flux, target, &
                7, bem_potential(i), status)
            if (status /= 0) error stop "Toroidal BEM evaluation failed"
            call numerical_laplace_gradient(target, bem_gradient)
            bem_field_norm(i) = norm2(bem_gradient)
        end do
        call cpu_time(end_time)
        bem_seconds = end_time - start_time

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(theta_trace, potential, label="Poisson analytical", &
            linestyle="-")
        call plot(theta_trace, bem_potential, label="Poisson BEM", &
            linestyle="--")
        call plot(theta_trace, field_norm, label="Ampere analytical |H|", &
            linestyle="-.")
        call plot(theta_trace, bem_field_norm, label="Ampere BEM |H|", &
            linestyle=":")
        call xlabel("poloidal angle theta [rad]")
        call ylabel("normalized mode value")
        call title("Toroidal Poisson and Ampere: analytical vs BEM")
        call legend()
        call savefig(output_directory//"/toroidal_trace_1d.png")
    end subroutine generate_trace

    subroutine generate_surface_map()
        do i = 1, theta_cells + 1
            theta_edges(i) = 2.0_dp*pi*real(i - 1, dp)/ &
                real(theta_cells, dp)
        end do
        do j = 1, phi_cells + 1
            phi_edges(j) = 2.0_dp*pi*real(j - 1, dp)/ &
                real(phi_cells, dp)
        end do
        do i = 1, theta_cells
            theta_value = 0.5_dp*(theta_edges(i) + theta_edges(i + 1))
            do j = 1, phi_cells
                phi_value = 0.5_dp*(phi_edges(j) + phi_edges(j + 1))
                call evaluate_toroidal_harmonic_p( &
                    degree_index, order, target_eta, theta_value, phi_value, &
                    surface_value(j, i), status)
                if (status /= 0) error stop "Toroidal map evaluation failed"
                call toroidal_point_to_cartesian( &
                    scale, target_eta, theta_value, phi_value, target)
                call evaluate_laplace_representation_torus_curved_3d( &
                    parameters, triangles, 2.0_dp/sqrt(3.0_dp), &
                    1.0_dp/sqrt(3.0_dp), boundary_trace, boundary_flux, &
                    target, 7, bem_potential(1), status)
                if (status /= 0) error stop "Toroidal BEM map failed"
                bem_error(j, i) = &
                    abs(bem_potential(1) - surface_value(j, i))
            end do
        end do
        bem_error = bem_error/max(maxval(abs(surface_value)), tiny(1.0_dp))

        call figure(figsize=[8.5_dp, 5.5_dp])
        call pcolormesh(theta_edges, phi_edges, surface_value, cmap="coolwarm")
        call xlabel("poloidal angle theta [rad]")
        call ylabel("toroidal angle phi [rad]")
        call title("Poisson analytical toroidal surface mode")
        call savefig(output_directory//"/toroidal_surface_2d.png")

        call figure(figsize=[8.5_dp, 5.5_dp])
        call pcolormesh(theta_edges, phi_edges, bem_error, cmap="viridis")
        call xlabel("poloidal angle theta [rad]")
        call ylabel("toroidal angle phi [rad]")
        call title("Poisson BEM normalized absolute error")
        call savefig(output_directory//"/toroidal_bem_error_2d.png")
    end subroutine generate_surface_map

    subroutine generate_geometry()
        point = 0
        do j = 1, phi_cells
            phi_value = 2.0_dp*pi*real(j - 1, dp)/real(phi_cells, dp)
            do i = 1, theta_cells
                theta_value = 2.0_dp*pi*real(i - 1, dp)/ &
                    real(theta_cells, dp)
                denominator = cosh(eta) - cos(theta_value)
                point = point + 1
                x(point) = scale*sinh(eta)*cos(phi_value)/denominator
                y(point) = scale*sinh(eta)*sin(phi_value)/denominator
                z(point) = scale*sin(theta_value)/denominator
            end do
        end do

        call figure(figsize=[7.5_dp, 6.5_dp])
        call add_scatter(x, y, z, label="eta=acosh(2) torus", marker=".")
        call title("Toroidal exterior interface")
        call savefig(output_directory//"/toroidal_geometry_3d.png")
    end subroutine generate_geometry

    subroutine generate_helmholtz_trace()
        real(dp) :: radius, source(3)

        source = [2.0_dp/sqrt(3.0_dp), 0.0_dp, 0.0_dp]
        allocate(helmholtz_trace(size(vertices, 2)))
        call build_helmholtz_boundary_data(source)
        call cpu_time(start_time)
        call solve_helmholtz_bem_dtn_torus_curved_3d( &
            parameters, triangles, 2.0_dp/sqrt(3.0_dp), &
            1.0_dp/sqrt(3.0_dp), wave_number, helmholtz_trace, 7, &
            helmholtz_flux, status)
        call cpu_time(end_time)
        if (status /= 0) error stop "Toroidal Helmholtz DtN solve failed"
        helmholtz_solve_seconds = end_time - start_time
        call cpu_time(start_time)
        do i = 1, trace_points
            axial_z(i) = &
                -1.5_dp + 3.0_dp*real(i - 1, dp)/real(trace_points - 1, dp)
            target = [0.0_dp, 0.0_dp, axial_z(i)]
            radius = norm2(target - source)
            helmholtz_exact(i) = &
                exp(cmplx(0.0_dp, wave_number*radius, dp))/radius
            call evaluate_helmholtz_representation_torus_curved_3d( &
                parameters, triangles, 2.0_dp/sqrt(3.0_dp), &
                1.0_dp/sqrt(3.0_dp), helmholtz_trace, helmholtz_flux, &
                target, wave_number, 7, helmholtz_bem(i), status)
            if (status /= 0) error stop "Toroidal Helmholtz BEM failed"
        end do
        call cpu_time(end_time)
        helmholtz_seconds = end_time - start_time

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(axial_z, abs(helmholtz_exact), &
            label="Helmholtz analytical |u|", linestyle="-")
        call plot(axial_z, abs(helmholtz_bem), &
            label="Helmholtz BEM |u|", linestyle="--")
        call xlabel("axial coordinate z")
        call ylabel("field magnitude")
        call title("Outgoing Helmholtz point source through toroidal BEM")
        call legend()
        call savefig(output_directory//"/toroidal_helmholtz_1d.png")
    end subroutine generate_helmholtz_trace

    subroutine generate_fem_bem_comparison()
        integer, parameter :: poloidal_count = 7
        integer, parameter :: toroidal_count = 9
        integer :: first_boundary_vertex, node
        real(dp) :: radius, source(3)

        call generate_solid_torus_tetra_mesh( &
            2.0_dp/sqrt(3.0_dp), 1.0_dp/sqrt(3.0_dp), 2, poloidal_count, &
            toroidal_count, volume_vertices, volume_tetrahedra, &
            volume_boundary_triangles, volume_parameters)
        allocate(coupled_load(size(volume_vertices, 2)))
        allocate(coupled_helmholtz_load(size(volume_vertices, 2)))
        coupled_load = 0.0_dp
        coupled_helmholtz_load = cmplx(0.0_dp, 0.0_dp, dp)
        coupled_load(2) = 1.0_dp
        coupled_helmholtz_load(2) = cmplx(1.0_dp, 0.0_dp, dp)
        source = volume_vertices(:, 2)

        call cpu_time(start_time)
        call solve_laplace_fem_bem_costabel_torus_curved_3d( &
            volume_vertices, volume_tetrahedra, volume_parameters, &
            volume_boundary_triangles, 2.0_dp/sqrt(3.0_dp), &
            1.0_dp/sqrt(3.0_dp), coupled_load, 5, coupled_potential, &
            coupled_flux, status)
        call cpu_time(end_time)
        if (status /= 0) error stop "Toroidal Laplace FEM-BEM solve failed"
        poisson_coupled_seconds = end_time - start_time

        call cpu_time(start_time)
        call solve_helmholtz_fem_bem_costabel_torus_curved_3d( &
            volume_vertices, volume_tetrahedra, volume_parameters, &
            volume_boundary_triangles, 2.0_dp/sqrt(3.0_dp), &
            1.0_dp/sqrt(3.0_dp), wave_number, wave_number, &
            coupled_helmholtz_load, 5, coupled_helmholtz_potential, &
            coupled_helmholtz_flux, status)
        call cpu_time(end_time)
        if (status /= 0) error stop "Toroidal Helmholtz FEM-BEM solve failed"
        helmholtz_coupled_seconds = end_time - start_time

        allocate(coupled_theta(poloidal_count))
        allocate(coupled_exact(poloidal_count))
        allocate(coupled_numerical(poloidal_count))
        allocate(coupled_helmholtz_exact(poloidal_count))
        allocate(coupled_helmholtz_numerical(poloidal_count))
        first_boundary_vertex = 2 + poloidal_count
        do i = 1, poloidal_count
            node = first_boundary_vertex + i - 1
            coupled_theta(i) = volume_parameters(1, node)
            radius = norm2(volume_vertices(:, node) - source)
            coupled_exact(i) = 1.0_dp/(4.0_dp*pi*radius)
            coupled_numerical(i) = coupled_potential(node)
            coupled_helmholtz_exact(i) = &
                exp(cmplx(0.0_dp, wave_number*radius, dp))/(4.0_dp*pi*radius)
            coupled_helmholtz_numerical(i) = coupled_helmholtz_potential(node)
        end do
        poisson_coupled_error = norm2(coupled_numerical - coupled_exact)/ &
            norm2(coupled_exact)
        helmholtz_coupled_error = sqrt(sum(abs( &
            coupled_helmholtz_numerical - coupled_helmholtz_exact)**2))/ &
            sqrt(sum(abs(coupled_helmholtz_exact)**2))

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(coupled_theta, coupled_exact, &
            label="Poisson analytical", linestyle="-")
        call plot(coupled_theta, coupled_numerical, &
            label="Poisson FEM-BEM", linestyle="--")
        call plot(coupled_theta, abs(coupled_helmholtz_exact), &
            label="Helmholtz analytical |u|", linestyle="-.")
        call plot(coupled_theta, abs(coupled_helmholtz_numerical), &
            label="Helmholtz FEM-BEM |u|", linestyle=":")
        call xlabel("poloidal angle theta [rad]")
        call ylabel("boundary potential")
        call title("Solid-torus FEM-BEM coupling against analytical fields")
        call legend()
        call savefig(output_directory//"/toroidal_fem_bem_1d.png")
    end subroutine generate_fem_bem_comparison

    subroutine write_trace_data()
        open (newunit=unit, &
            file=output_directory//"/toroidal_trace.csv", &
            status="replace", action="write")
        write (unit, "(a)") &
            "theta,potential_exact,potential_bem,H_norm_exact,H_norm_bem"
        do i = 1, trace_points
            write (unit, "(*(es24.16,:,','))") &
                theta_trace(i), potential(i), bem_potential(i), &
                field_norm(i), bem_field_norm(i)
        end do
        close (unit)
        open (newunit=unit, file=output_directory//"/benchmark.txt", &
            status="replace", action="write")
        write (unit, "(a,i0)") "boundary triangles: ", size(triangles, 2)
        write (unit, "(a,i0)") "trace targets: ", trace_points
        write (unit, "(a,es14.6)") "analytical seconds: ", exact_seconds
        write (unit, "(a,es14.6)") &
            "Poisson BEM DtN solve seconds: ", poisson_solve_seconds
        write (unit, "(a,es14.6)") "BEM seconds: ", bem_seconds
        write (unit, "(a,es14.6)") &
            "Helmholtz BEM DtN solve seconds: ", helmholtz_solve_seconds
        write (unit, "(a,es14.6)") &
            "Helmholtz BEM seconds: ", helmholtz_seconds
        write (unit, "(a,es14.6)") &
            "Poisson FEM-BEM solve seconds: ", poisson_coupled_seconds
        write (unit, "(a,es14.6)") &
            "Helmholtz FEM-BEM solve seconds: ", helmholtz_coupled_seconds
        write (unit, "(a,es14.6)") "Poisson max relative error: ", maxval( &
            abs(bem_potential - potential))/maxval(abs(potential))
        write (unit, "(a,es14.6)") "Ampere max relative error: ", maxval( &
            abs(bem_field_norm - field_norm))/maxval(field_norm)
        write (unit, "(a,es14.6)") "Helmholtz max relative error: ", maxval( &
            abs(helmholtz_bem - helmholtz_exact))/ &
            maxval(abs(helmholtz_exact))
        write (unit, "(a,es14.6)") &
            "Poisson FEM-BEM boundary relative error: ", poisson_coupled_error
        write (unit, "(a,es14.6)") &
            "Helmholtz FEM-BEM boundary relative error: ", &
            helmholtz_coupled_error
        close (unit)
        if (maxval(abs(bem_potential - potential))/ &
            maxval(abs(potential)) >= 3.0e-1_dp) &
            error stop "Poisson BEM gallery accuracy regression"
        if (maxval(abs(bem_field_norm - field_norm))/ &
            maxval(field_norm) >= 3.0e-1_dp) &
            error stop "Ampere BEM gallery accuracy regression"
        if (maxval(abs(helmholtz_bem - helmholtz_exact))/ &
            maxval(abs(helmholtz_exact)) >= 3.0e-2_dp) &
            error stop "Helmholtz BEM gallery accuracy regression"
        if (poisson_coupled_error >= 3.5e-1_dp) &
            error stop "Poisson FEM-BEM gallery accuracy regression"
        if (helmholtz_coupled_error >= 3.5e-1_dp) &
            error stop "Helmholtz FEM-BEM gallery accuracy regression"
    end subroutine write_trace_data

    subroutine build_boundary_data()
        real(dp) :: eta_at_point

        allocate(boundary_trace(size(vertices, 2)))
        do point = 1, size(vertices, 2)
            call cartesian_to_toroidal( &
                vertices(:, point), scale, eta_at_point, theta_value, phi_value)
            call evaluate_toroidal_harmonic_p( &
                degree_index, order, eta, theta_value, phi_value, &
                boundary_trace(point), status)
            if (status /= 0) error stop "Toroidal boundary trace failed"
        end do
    end subroutine build_boundary_data

    subroutine build_helmholtz_boundary_data(source)
        real(dp), intent(in) :: source(3)
        real(dp) :: displacement(3), radius

        do point = 1, size(vertices, 2)
            displacement = vertices(:, point) - source
            radius = norm2(displacement)
            helmholtz_trace(point) = &
                exp(cmplx(0.0_dp, wave_number*radius, dp))/radius
        end do
    end subroutine build_helmholtz_boundary_data

    subroutine numerical_laplace_gradient(evaluation_point, gradient)
        real(dp), intent(in) :: evaluation_point(3)
        real(dp), intent(out) :: gradient(3)

        real(dp), parameter :: step = 1.0e-4_dp
        real(dp) :: backward_point(3), backward_value
        real(dp) :: forward_point(3), forward_value
        integer :: coordinate

        do coordinate = 1, 3
            forward_point = evaluation_point
            backward_point = evaluation_point
            forward_point(coordinate) = forward_point(coordinate) + step
            backward_point(coordinate) = backward_point(coordinate) - step
            call evaluate_laplace_representation_torus_curved_3d( &
                parameters, triangles, 2.0_dp/sqrt(3.0_dp), &
                1.0_dp/sqrt(3.0_dp), boundary_trace, boundary_flux, &
                forward_point, 7, forward_value, status)
            if (status /= 0) error stop "Forward BEM gradient failed"
            call evaluate_laplace_representation_torus_curved_3d( &
                parameters, triangles, 2.0_dp/sqrt(3.0_dp), &
                1.0_dp/sqrt(3.0_dp), boundary_trace, boundary_flux, &
                backward_point, 7, backward_value, status)
            if (status /= 0) error stop "Backward BEM gradient failed"
            gradient(coordinate) = -(forward_value - backward_value)/(2*step)
        end do
    end subroutine numerical_laplace_gradient

end program toroidal_analytical
```

## Generated Plots

### toroidal_bem_error_2d.png

![toroidal_bem_error_2d.png](../../media/examples/toroidal_analytical/toroidal_bem_error_2d.png)

### toroidal_fem_bem_1d.png

![toroidal_fem_bem_1d.png](../../media/examples/toroidal_analytical/toroidal_fem_bem_1d.png)

### toroidal_geometry_3d.png

![toroidal_geometry_3d.png](../../media/examples/toroidal_analytical/toroidal_geometry_3d.png)

### toroidal_helmholtz_1d.png

![toroidal_helmholtz_1d.png](../../media/examples/toroidal_analytical/toroidal_helmholtz_1d.png)

### toroidal_surface_2d.png

![toroidal_surface_2d.png](../../media/examples/toroidal_analytical/toroidal_surface_2d.png)

### toroidal_trace_1d.png

![toroidal_trace_1d.png](../../media/examples/toroidal_analytical/toroidal_trace_1d.png)

---

[← Back to all examples](../index.html)
