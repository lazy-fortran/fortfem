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