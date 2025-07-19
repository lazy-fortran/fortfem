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

### 1. Poisson 2D (`poisson_2d.f90`)
Solves the 2D Poisson equation with Dirichlet boundary conditions.

**Problem**: 
- Domain: Unit square [0,1]×[0,1]
- Equation: -∆u = f
- Boundary: u = 0

**Features demonstrated**:
- Mesh generation
- P1 finite element assembly
- Boundary condition application
- Sparse linear system solution

### 2. Basis Function Visualization (`plot_basis.f90`)
Visualizes finite element basis functions using fortplotlib.

**Features demonstrated**:
- P1 and P2 basis function evaluation
- Integration with plotting library
- Reference element transformations

### 3. Mesh 2D Demo (`mesh_2d_demo.f90`)
Demonstrates mesh generation and manipulation capabilities.

**Features demonstrated**:
- Rectangular mesh generation
- Mesh connectivity queries
- Boundary edge identification
- Mesh quality metrics

### 4. Curl-Curl Solver (`curl_curl_example.f90`)
Solves the curl-curl equation using edge elements.

**Problem**:
- Equation: ∇×∇×E + k²E = J
- Boundary: E×n = 0 (perfect electric conductor)

**Features demonstrated**:
- Nédélec edge elements
- H(curl) conforming spaces
- Vector field assembly

### 5. Heat Equation (`heat_equation.f90`)
Solves the time-dependent heat equation.

**Features demonstrated**:
- Time stepping schemes
- Mass matrix assembly
- Mixed boundary conditions

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