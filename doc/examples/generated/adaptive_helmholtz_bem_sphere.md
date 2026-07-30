---
title: adaptive_helmholtz_bem_sphere Example
---

# adaptive_helmholtz_bem_sphere Example

# Adaptive outgoing Helmholtz BEM on a sphere

This example applies the solve-estimate-mark-refine loop to the complex P0
single-layer equation for a unit sphere. New vertices are projected back to
the analytical surface after conforming red/green refinement.

For constant Dirichlet data and the outgoing convention used by FortFEM, the
radial exterior solution is

\[
u(r)=\frac{\exp(ik(r-1))}{r}.
\]

The field at \(r=2\) is therefore an independent analytical oracle. The
fine-grid Galerkin residual supplies local indicators, and Dörfler marking
selects panels carrying half of their squared norm.

Provenance:

- [Feischl et al., *Convergence of adaptive BEM and adaptive FEM-BEM
  coupling for estimators without h-weighting
  factor*](https://arxiv.org/abs/1405.5306) gives the adaptive loop and
  two-level-estimator setting.
- [NIST DLMF §10.47](https://dlmf.nist.gov/10.47) defines the outgoing
  spherical Hankel functions; its order-zero case reduces to the radial
  expression above after enforcing \(u(1)=1\).

Run with:

```sh
fpm run --example adaptive_helmholtz_bem_sphere
```

## Usage

```bash
fpm run --example adaptive_helmholtz_bem_sphere
```

## Source Code

```fortran
program adaptive_helmholtz_bem_sphere
    use fortfem_api, only: &
        estimate_helmholtz_p0_two_level_residual_3d, &
        evaluate_helmholtz_representation_triangles_3d, &
        generate_sphere_surface_mesh, mark_bem_dorfler, &
        refine_surface_mesh_marked, solve_helmholtz_dirichlet_p0_3d
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, &
        ylabel, yscale
    implicit none

    integer, parameter :: step_count = 5
    real(dp), parameter :: wavenumber = 0.7_dp
    character(*), parameter :: output_directory = &
        "output/example/adaptive_helmholtz_bem_sphere"
    complex(dp), allocatable :: density(:), dirichlet(:)
    integer, allocatable :: parent(:), refined_triangles(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: indicators(:), refined_vertices(:, :)
    real(dp), allocatable :: vertices(:, :)
    logical, allocatable :: marked(:)
    complex(dp) :: exact_field, numerical_field
    real(dp) :: error(step_count), estimator(step_count), panels(step_count)
    integer :: status, step

    call execute_command_line("mkdir -p "//output_directory)
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    exact_field = 0.5_dp*exp(cmplx(0.0_dp, wavenumber, dp))

    do step = 1, step_count
        call solve_helmholtz_dirichlet_p0_3d( &
            vertices, triangles, cmplx(1.0_dp, 0.0_dp, dp), wavenumber, &
            8, density, status)
        if (status /= 0) error stop "adaptive Helmholtz solve failed"
        if (allocated(dirichlet)) deallocate(dirichlet)
        allocate(dirichlet(size(vertices, 2)), source=(0.0_dp, 0.0_dp))
        call evaluate_helmholtz_representation_triangles_3d( &
            vertices, triangles, dirichlet, -density, &
            [0.0_dp, 0.0_dp, 2.0_dp], wavenumber, 8, numerical_field, status)
        if (status /= 0) error stop "adaptive Helmholtz evaluation failed"
        call estimate_helmholtz_p0_two_level_residual_3d( &
            vertices, triangles, density, cmplx(1.0_dp, 0.0_dp, dp), &
            wavenumber, 6, indicators, status)
        if (status /= 0) error stop "adaptive Helmholtz estimate failed"
        panels(step) = real(size(triangles, 2), dp)
        error(step) = abs(numerical_field - exact_field)
        estimator(step) = norm2(indicators)
        print "(i4,2es14.5)", size(triangles, 2), error(step), estimator(step)
        if (step == step_count) exit
        call mark_bem_dorfler(indicators, 0.5_dp, marked, status)
        if (status /= 0) error stop "adaptive Helmholtz marking failed"
        call refine_surface_mesh_marked( &
            vertices, triangles, marked, refined_vertices, &
            refined_triangles, parent, status)
        if (status /= 0) error stop "adaptive Helmholtz refinement failed"
        call project_to_sphere(refined_vertices)
        call move_alloc(refined_vertices, vertices)
        call move_alloc(refined_triangles, triangles)
    end do

    call figure()
    call plot(panels, error, label="field error at r=2")
    call plot(panels, estimator, label="residual estimator")
    call yscale("log")
    call xlabel("surface panels")
    call ylabel("absolute error / estimator")
    call title("Adaptive outgoing Helmholtz BEM on the sphere")
    call legend()
    call savefig(output_directory//"/adaptive_convergence.png")

contains

    subroutine project_to_sphere(points)
        real(dp), intent(inout) :: points(:, :)
        integer :: point

        do point = 1, size(points, 2)
            points(:, point) = points(:, point)/norm2(points(:, point))
        end do
    end subroutine project_to_sphere

end program adaptive_helmholtz_bem_sphere
```

## Generated Plots

*No plot artifact is produced by this example.*

---

[← Back to all examples](../index.html)
