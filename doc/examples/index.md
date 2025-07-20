---
title: Examples
---

# FortFEM Examples

This page provides an overview of the example programs included with FortFEM. All examples can be found in the [example/](https://github.com/itpplasma/fortfem/tree/main/example) directory.

## Running Examples

To run any example:
```bash
fpm run --example <example_name>
```

To list all available examples:
```bash
fpm run --example --list
```

## Available Examples

### 1. Simple Poisson ([`simple_poisson.f90`](https://github.com/itpplasma/fortfem/blob/main/example/simple_poisson.f90))
Clean, minimal example solving -∆u = 1 on unit square with zero boundary conditions.

**Features**:
- FEniCS-style API with ~10 lines of meaningful code
- P1 Lagrange finite elements
- Direct solver with automatic assembly  
- One-line plotting: `call plot(uh, title="Solution")`

```fortran
mesh = unit_square_mesh(20)
Vh = function_space(mesh, "Lagrange", 1)

u = trial_function(Vh)
v = test_function(Vh)
f = constant(1.0_dp)

a = inner(grad(u), grad(v))*dx
L = f*v*dx

call solve(a == L, uh, bc)
call plot(uh, filename="poisson_solution.png")
```

### 2. Curl-Curl Electromagnetic ([`curl_curl_example.f90`](https://github.com/itpplasma/fortfem/blob/main/example/curl_curl_example.f90))
Vector electromagnetic problem using Nédélec edge elements.

**Features**:
- H(curl) conforming Nédélec edge elements
- GMRES iterative solver for large sparse systems
- Tangential boundary conditions: E×n = 0
- Vector field visualization with streamplots
- Analytical solution comparison: E = [x*y, x²]

```fortran
Vh = vector_function_space(mesh, "Nedelec", 1)
E = vector_trial_function(Vh)
F = vector_test_function(Vh)

a = inner(curl(E), curl(F))*dx + inner(E, F)*dx
L = inner(J, F)*dx

call solve(a == L, Eh, bc)
call plot(Eh, plot_type="streamplot")
```

### 3. Plotting Demonstration ([`plotting_demo.f90`](https://github.com/itpplasma/fortfem/blob/main/example/plotting_demo.f90))
Comprehensive showcase of FortFEM's plotting capabilities.

**Features**:
- Multiple colormap options (viridis, plasma, jet, coolwarm)
- Custom titles and filenames
- Both scalar and vector field plotting
- Fortplotlib integration for publication-quality plots

## Example Structure

Each example follows the FEniCS-style pattern:
1. **Mesh and function space**: `mesh = unit_square_mesh(n)`, `Vh = function_space(mesh, "Lagrange", 1)`
2. **Trial and test functions**: `u = trial_function(Vh)`, `v = test_function(Vh)`
3. **Weak forms**: `a = inner(grad(u), grad(v))*dx`, `L = f*v*dx`
4. **Boundary conditions**: `bc = dirichlet_bc(Vh, 0.0_dp)`
5. **Solve**: `call solve(a == L, uh, bc)`
6. **Visualization**: `call plot(uh, title="Solution")`

## Creating Your Own Examples

To create a new example:

1. Create a new file in `example/` directory
2. Follow the minimal FEniCS-style pattern:

```fortran
program my_example
    use fortfem_api
    implicit none
    
    type(mesh_t) :: mesh
    type(function_space_t) :: Vh
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(function_t) :: uh, f
    type(dirichlet_bc_t) :: bc
    type(form_expr_t) :: a, L
    
    ! Problem setup
    mesh = unit_square_mesh(20)
    Vh = function_space(mesh, "Lagrange", 1)
    
    u = trial_function(Vh)
    v = test_function(Vh)
    f = constant(1.0_dp)
    
    ! Weak form
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx
    
    ! Solve and plot
    uh = function(Vh)
    bc = dirichlet_bc(Vh, 0.0_dp)
    call solve(a == L, uh, bc)
    call plot(uh, title="My Solution")
    
end program
```

## Visualization

FortFEM provides built-in plotting via fortplotlib:

```fortran
! Scalar field plotting
call plot(uh, filename="solution.png", title="Poisson Solution", colormap="viridis")

! Vector field plotting  
call plot(Eh, filename="field.png", plot_type="streamplot", title="E Field")
```

Available colormaps: `viridis`, `plasma`, `jet`, `coolwarm`, `hot`, `gray`

## Benchmarks and Validation

The `benchmark/` directory contains FreeFEM reference solutions for validation:
- Curl-curl convergence studies
- Error norm computations
- Performance comparisons