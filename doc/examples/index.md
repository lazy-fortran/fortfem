---
title: Examples
---

# FortFEM Examples

This page provides an overview of the example programs included with FortFEM, complete with source code listings and generated plots.

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

- [curl_curl](generated/curl_curl.html) - Clean FEniCS-style example solving the curl-curl equation
- [fortfem_mesh_benchmark](generated/fortfem_mesh_benchmark.html) - FortFEM mesh generation benchmark to compare with FreeFEM
- [minimal_mesh_example](generated/minimal_mesh_example.html) - Minimal working example of FortFEM mesh generation
- [plot_mesh](generated/plot_mesh.html) - Example demonstrating mesh plotting functionality in FortFEM
- [plotting](generated/plotting.html) - Demonstration of FortFEM plotting capabilities
- [simple_poisson](generated/simple_poisson.html) - 

## Creating Your Own Examples

To create a new example:

1. Create a new file in `example/` directory
2. Follow the minimal FEniCS-style pattern
3. Optionally create a `<example_name>_README.md` file for detailed documentation

## Visualization

FortFEM provides built-in plotting via fortplotlib:

```fortran
! Scalar field plotting
call plot(uh, filename="solution.png", plot_title="Poisson Solution", colormap="viridis")

! Vector field plotting  
call plot(Eh, filename="field.png", plot_type="streamplot", plot_title="E Field")

! Mesh plotting
call plot(mesh, filename="mesh.png", plot_title="FEM Mesh")
```

Available colormaps: `viridis`, `plasma`, `jet`, `coolwarm`, `hot`, `gray`

---

[← Back to Documentation](../index.html)
