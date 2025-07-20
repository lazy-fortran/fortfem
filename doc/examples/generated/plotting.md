---
title: plotting Example
---

# plotting Example

# Plotting Demonstration

This comprehensive example showcases FortFEM's complete plotting and visualization capabilities using multiple colormap options and plot types.

## Problem Description

The example solves multiple Poisson problems with different source terms to demonstrate various plotting features and colormap options available in FortFEM.

## Problems Solved

### 1. Constant Source (f = 1)
Standard Poisson equation with uniform source term:
```
-Δu = 1    in Ω
 u = 0     on ∂Ω
```

### 2. Point Source (f = 10)  
Intensified source to show dynamic range:
```
-Δu = 10   in Ω
 u = 0     on ∂Ω
```

### 3. Same Problem, Different Colormap
Demonstrates colormap flexibility using the same solution data.

## Features Demonstrated

- **Multiple colormaps**: viridis, plasma, coolwarm, jet
- **Custom titles**: Descriptive plot annotations
- **Filename control**: Organized output file naming
- **Solution comparison**: Multiple problems on same mesh
- **Publication quality**: High-resolution PNG output via fortplotlib

## Available Colormaps

FortFEM supports all major scientific colormaps:

- **viridis**: Perceptually uniform, colorblind-friendly (default)
- **plasma**: High contrast, good for presentations
- **coolwarm**: Diverging colormap for signed data
- **jet**: Classic rainbow colormap
- **hot**: Heat map visualization
- **gray**: Grayscale for monochrome publications

## Output Files

- `solution_constant.png`: Constant source with viridis colormap
- `solution_point.png`: Point source with plasma colormap  
- `solution_coolwarm.png`: Same data with cool-warm diverging colors

## API Examples

The example demonstrates the complete plotting API:

```fortran
! Basic plotting
call plot(uh, filename="solution.png")

! With custom title and colormap
call plot(uh, filename="custom.png", plot_title="My Solution", colormap="viridis")

! Different colormaps
call plot(uh, colormap="plasma")    ! High contrast
call plot(uh, colormap="coolwarm")  ! Diverging
call plot(uh, colormap="jet")       ! Rainbow
```

## Integration with fortplotlib

FortFEM's plotting system builds on fortplotlib to provide:
- Publication-quality vector graphics
- Multiple output formats (PNG, PDF, SVG)
- LaTeX-compatible text rendering
- Customizable figure sizes and DPI
## Usage

```bash
fpm run --example plotting
```

## Source Code

```fortran
program plotting_demo
    ! Demonstration of FortFEM plotting capabilities
    ! Shows how to easily visualize FEM solutions with a single plot() command
    
    use fortfem_kinds
    use fortfem_api
    implicit none

    type(mesh_t) :: mesh
    type(function_space_t) :: Vh
    type(function_t) :: uh1, uh2, uh3
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(function_t) :: f
    type(dirichlet_bc_t) :: bc
    type(form_expr_t) :: a, L
    
    write(*,*) "=== FortFEM Plotting Demonstration ==="
    write(*,*) ""
    
    ! Create mesh and function space
    mesh = unit_square_mesh(15)  ! Fine mesh for better plots
    Vh = function_space(mesh, "Lagrange", 1)
    
    u = trial_function(Vh)
    v = test_function(Vh)
    bc = dirichlet_bc(Vh, 0.0_dp)
    
    ! Define weak form for Poisson equation
    a = inner(grad(u), grad(v))*dx
    
    write(*,*) "Creating multiple solutions with different source terms..."
    write(*,*) ""
    
    ! Solution 1: Constant source f = 1
    write(*,*) "1. Constant source f = 1"
    f = constant(1.0_dp)
    L = f*v*dx
    uh1 = function(Vh)
    call solve(a == L, uh1, bc)
    call plot(uh1, filename="solution_constant.png", &
              title="Constant Source: f = 1", &
              colormap="viridis")
    
    ! Solution 2: Point source (approximated by high value in center)
    write(*,*) "2. Point source (approximated)"
    f = constant(10.0_dp)
    L = f*v*dx
    uh2 = function(Vh)
    call solve(a == L, uh2, bc)
    call plot(uh2, filename="solution_point.png", &
              title="Point Source: f = 10", &
              colormap="plasma")
    
    ! Solution 3: Different colormap demonstration
    write(*,*) "3. Same solution with different colormap"
    call plot(uh1, filename="solution_coolwarm.png", &
              title="Constant Source with Cool-Warm Colormap", &
              colormap="coolwarm")
    
    write(*,*) ""
    write(*,*) "Generated plots:"
    write(*,*) "- solution_constant.png   (viridis colormap)"
    write(*,*) "- solution_point.png      (plasma colormap)"
    write(*,*) "- solution_coolwarm.png   (coolwarm colormap)"
    write(*,*) ""
    write(*,*) "All plots created with single plot() command!"
    write(*,*) ""
    write(*,*) "Usage examples:"
    write(*,*) '  call plot(uh)                                   ! Default options'
    write(*,*) '  call plot(uh, "my_solution.png")               ! Custom filename'
    write(*,*) '  call plot(uh, title="My Title")                ! Custom title'
    write(*,*) '  call plot(uh, colormap="plasma")               ! Custom colormap'
    write(*,*) ""
    write(*,*) "Available colormaps: viridis, plasma, coolwarm, jet, etc."

end program plotting_demo```

## Generated Plots

### solution_constant.png

![solution_constant.png](../../../artifacts/plots/solution_constant.png)

### solution_coolwarm.png

![solution_coolwarm.png](../../../artifacts/plots/solution_coolwarm.png)

### solution_point.png

![solution_point.png](../../../artifacts/plots/solution_point.png)

### test_mesh.png

![test_mesh.png](../../../artifacts/plots/test_mesh.png)

### test_mesh_titled.png

![test_mesh_titled.png](../../../artifacts/plots/test_mesh_titled.png)


---

[← Back to Examples](../index.html) | [FortFEM Documentation](../../index.html)
