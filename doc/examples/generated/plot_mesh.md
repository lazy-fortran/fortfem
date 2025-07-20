---
title: plot_mesh Example
---

# plot_mesh Example

# Mesh Plotting Example

This example demonstrates FortFEM's mesh visualization capabilities, showing how to create and plot finite element meshes at different refinement levels.

## Problem Description

This example creates unit square meshes with varying levels of refinement and generates plots to visualize the triangulation quality and mesh density.

## Features Demonstrated

- **Mesh generation**: Creating unit square meshes with different refinement levels
- **Mesh plotting**: Visualizing triangular finite element meshes
- **Refinement comparison**: Comparing coarse, medium, and fine meshes
- **Plot customization**: Using custom titles and filenames for different mesh plots

## Mesh Configurations

The example generates three different mesh configurations:

1. **Basic mesh** (5 refinements): Standard resolution for typical computations
2. **Fine mesh** (10 refinements): High resolution for accuracy studies  
3. **Coarse mesh** (3 refinements): Low resolution for debugging and prototyping

## Output Files

- `mesh_basic.png`: 5×5 refinement level mesh (25 vertices, 32 triangles)
- `mesh_fine.png`: 10×10 refinement level mesh (100 vertices, 162 triangles)
- `mesh_coarse.png`: 3×3 refinement level mesh (9 vertices, 8 triangles)

## Mesh Statistics

Each plot displays mesh statistics including:
- Number of vertices
- Number of triangular elements
- Number of edges

## Visualization Features

- **Triangle edges**: Blue lines showing element boundaries
- **Vertices**: Red dots marking mesh nodes
- **Axis labels**: Coordinate system reference
- **Automatic scaling**: Proper axis limits with margins
## Usage

```bash
fpm run --example plot_mesh
```

## Source Code

```fortran
program plot_mesh_example
    ! Example demonstrating mesh plotting functionality in FortFEM
    
    use fortfem_api
    implicit none
    
    type(mesh_t) :: mesh
    integer :: n_refinements
    
    ! Create unit square mesh with different refinement levels
    n_refinements = 5
    mesh = unit_square_mesh(n_refinements)
    
    ! Plot basic mesh
    call plot(mesh, filename="mesh_basic.png", title="Unit Square Mesh")
    
    ! Create finer mesh
    n_refinements = 10
    mesh = unit_square_mesh(n_refinements)
    
    ! Plot finer mesh
    call plot(mesh, filename="mesh_fine.png", title="Refined Unit Square Mesh")
    
    ! Create coarse mesh for clarity
    n_refinements = 3
    mesh = unit_square_mesh(n_refinements)
    
    ! Plot with custom title
    call plot(mesh, filename="mesh_coarse.png", title="Coarse Mesh (3x3)")
    
    ! Clean up
    call mesh%destroy()
    
    print *, "Mesh plotting examples completed!"
    print *, "Generated files:"
    print *, "  - mesh_basic.png"
    print *, "  - mesh_fine.png"
    print *, "  - mesh_coarse.png"
    
end program plot_mesh_example```

## Generated Plots

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
