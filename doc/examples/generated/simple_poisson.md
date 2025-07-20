---
title: simple_poisson Example
---

# simple_poisson Example

# Simple Poisson Example

This example demonstrates solving the classic Poisson equation using FortFEM's FEniCS-inspired API.

## Problem Description

We solve the 2D Poisson equation on the unit square:

```
-Δu = 1    in Ω = [0,1] × [0,1]
 u = 0     on ∂Ω
```

Where:
- `Δu = ∂²u/∂x² + ∂²u/∂y²` is the Laplacian operator
- The right-hand side `f = 1` is a constant source term
- Zero Dirichlet boundary conditions are applied on all boundaries

## Features Demonstrated

- **FEniCS-style API**: Natural mathematical notation for expressing weak forms
- **P1 Lagrange elements**: Linear finite elements on triangular mesh
- **Direct solver**: Automatic assembly and solving with LAPACK
- **Mesh visualization**: Automatic mesh plotting showing triangulation
- **Solution visualization**: Scalar field plotting with customizable colormaps

## Mathematical Formulation

The weak formulation seeks `u ∈ H₀¹(Ω)` such that:

```
∫_Ω ∇u · ∇v dx = ∫_Ω f v dx    ∀v ∈ H₀¹(Ω)
```

## Output Files

- `poisson_mesh.png`: Visualization of the finite element mesh (20×20 elements)
- `poisson_solution.png`: Colored contour plot of the solution field

## Expected Results

The solution exhibits the characteristic "bowl" shape with maximum value at the center of the domain, decreasing to zero at the boundaries.
## Usage

```bash
fpm run --example simple_poisson
```

## Source Code

```fortran
program simple_poisson
    use fortfem_kinds
    use fortfem_api
    implicit none

    type(mesh_t) :: mesh
    type(function_space_t) :: Vh
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(function_t) :: f, uh
    type(dirichlet_bc_t) :: bc
    type(form_expr_t) :: a, L

    mesh = unit_square_mesh(20)
    Vh = function_space(mesh, "Lagrange", 1)
    
    u = trial_function(Vh)
    v = test_function(Vh)
    f = constant(1.0_dp)

    a = inner(grad(u), grad(v))*dx
    L = f*v*dx

    bc = dirichlet_bc(Vh, 0.0_dp)
    uh = function(Vh)

    call solve(a == L, uh, bc)
    
    ! Plot mesh
    call plot(mesh, filename="poisson_mesh.png", title="Poisson Mesh (20x20)")
    
    ! Plot solution
    call plot(uh, filename="poisson_solution.png", &
              title="Poisson Solution: -Δu = 1", &
              colormap="viridis")

    write(*,*) "Simple Poisson example completed!"
    write(*,*) "Generated files:"
    write(*,*) "  - Mesh: poisson_mesh.png"
    write(*,*) "  - Solution: poisson_solution.png"

end program simple_poisson```

## Generated Plots

### poisson_solution.png

![poisson_solution.png](../../../artifacts/plots/poisson_solution.png)

### mesh.png

![mesh.png](../../../artifacts/plots/mesh.png)

### mesh_basic.png

![mesh_basic.png](../../../artifacts/plots/mesh_basic.png)

### mesh_coarse.png

![mesh_coarse.png](../../../artifacts/plots/mesh_coarse.png)

### mesh_fine.png

![mesh_fine.png](../../../artifacts/plots/mesh_fine.png)


---

[← Back to Examples](../index.html) | [FortFEM Documentation](../../index.html)
