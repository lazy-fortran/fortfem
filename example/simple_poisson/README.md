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