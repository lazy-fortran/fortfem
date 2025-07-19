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

### 1. Simple Poisson (`poisson_simple.f90`)
Minimal example solving -∆u = f on unit square.

```fortran
! Create mesh
call create_unit_square_mesh(mesh, n=20)

! Assemble and solve
call assemble_poisson_2d(mesh, A, f)
call apply_zero_bc(mesh, A, f)
call solve_sparse(A, f, u)
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
3. Follow the pattern:

```fortran
program my_example
    use fortfem
    implicit none
    
    ! Declare variables
    type(mesh_2d_t) :: mesh
    type(p1_space_t) :: V
    ! ... more declarations
    
    ! Setup problem
    call create_mesh(...)
    call V%init(mesh)
    
    ! Assemble system
    call assemble_system(...)
    
    ! Solve
    call solve(...)
    
    ! Post-process
    call output_results(...)
    
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