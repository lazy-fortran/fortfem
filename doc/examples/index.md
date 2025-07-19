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

### 1. Simple Poisson
Minimal example solving -∆u = f on unit square using FEniCS-style API.

```fortran
! Create mesh and function space
mesh = unit_square_mesh(32)
Vh = function_space(mesh, "Lagrange", 1)

! Define trial and test functions
u = trial_function(Vh)
v = test_function(Vh)

! Define weak form
a = inner(grad(u), grad(v))*dx
L = inner(f, v)*dx
```

### 2. Mesh Demo (`mesh_simple.f90`, `mesh_2d_demo.f90`)
Create and analyze triangular meshes.

**Features**:
- Mesh generation
- Quality metrics
- VTK output

### 3. Curl-Curl Problem (`curl_curl_simple.f90`)
Electromagnetic problem using edge elements.

**Features**:
- Nédélec elements
- Vector fields
- Essential BC: E×n = 0

### 4. P2 Elements (`p2_poisson.f90`)
Higher-order finite elements for improved accuracy.

**Features**:
- Quadratic basis functions
- 6 DOFs per triangle
- Better convergence rates

### 5. Mixed Boundary Conditions (`neumann_bc.f90`)
Combining Dirichlet and Neumann conditions.

**Features**:
- Multiple BC types
- Flux conditions
- Natural boundaries

### 6. Convergence Test (`convergence_test.f90`)
Verify theoretical convergence rates.

**Features**:
- Mesh refinement
- Error computation
- Rate calculation

### 7. Basis Functions (`basis_functions.f90`)
Visualize finite element basis functions.

**Features**:
- 1D/2D basis plots
- Shape function evaluation

## Example Structure

Each example follows a similar structure:
1. **Problem setup**: Define geometry and parameters
2. **Mesh creation**: Generate or load mesh
3. **Function space**: Initialize finite element space
4. **Assembly**: Build system matrices and vectors
5. **Boundary conditions**: Apply constraints
6. **Solution**: Solve the linear system
7. **Post-processing**: Visualize or analyze results

## Creating Your Own Examples

To create a new example:

1. Create a new file in `example/` directory
2. Add the program to `fpm.toml` if needed
3. Follow the FEniCS-style pattern:

```fortran
program my_example
    use fortfem_api
    implicit none
    
    ! Declare variables
    type(mesh_t) :: mesh
    type(function_space_t) :: Vh
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(form_expr_t) :: a, L
    
    ! Setup problem
    mesh = unit_square_mesh(32)
    Vh = function_space(mesh, "Lagrange", 1)
    
    ! Define variational form
    u = trial_function(Vh)
    v = test_function(Vh)
    
    a = inner(grad(u), grad(v))*dx
    L = inner(f, v)*dx
    
    ! Solve (to be implemented)
    ! uh = solve(a == L, u, bc)
    
end program
```

## Visualization

Examples use different visualization approaches:
- **fortplotlib**: For 2D plots and basis function visualization
- **VTK output**: For ParaView/VisIt visualization
- **Console output**: For numerical verification

## Performance Examples

Some examples include timing and convergence studies:
- Mesh refinement studies
- Convergence rate verification
- Solver performance comparison