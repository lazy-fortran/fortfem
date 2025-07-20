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