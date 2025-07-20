# Curl-Curl Electromagnetic Example

This example demonstrates solving electromagnetic problems using Nédélec edge elements for H(curl) conforming finite element spaces.

## Problem Description

We solve the 2D curl-curl equation on the unit square:

```
curl(curl(E)) + E = J    in Ω = [0,1] × [0,1]
E × n = 0                on ∂Ω
```

Where:
- `E` is the electric field vector
- `J` is the current density source term
- `curl(E) = ∂E_y/∂x - ∂E_x/∂y` for 2D vector fields
- Tangential boundary conditions enforce `E × n = 0`

## Analytical Solution

The example uses the analytical solution:
```
E = [x*y, x²]
```

This gives:
- `curl(E) = ∂(x²)/∂x - ∂(x*y)/∂y = 2x - x = x`
- `curl(curl(E)) = ∂x/∂y - ∂(0)/∂x = 0`
- Therefore: `J = curl(curl(E)) + E = [x*y, x²]`

## Features Demonstrated

- **Nédélec edge elements**: H(curl) conforming finite element space
- **Vector problems**: Two-component electromagnetic fields
- **GMRES solver**: Iterative solver for large sparse systems
- **Tangential boundary conditions**: Natural boundary conditions for electromagnetic problems
- **Vector field visualization**: Streamplot representation of electric fields
- **Error analysis**: Comparison with analytical solution

## Mathematical Formulation

The weak formulation seeks `E ∈ H(curl, Ω)` such that:

```
∫_Ω curl(E) · curl(F) dx + ∫_Ω E · F dx = ∫_Ω J · F dx    ∀F ∈ H₀(curl, Ω)
```

## Output Files

- `curlcurl_solution.png`: Vector field plot of numerical solution
- `curlcurl_exact.png`: Vector field plot of analytical solution
- `curl_field.png`: Scalar plot of the curl field

## Convergence Analysis

The example performs convergence studies showing optimal convergence rates for:
- L² error in the electric field
- H(curl) error including curl components

## Applications

This formulation is fundamental for:
- Electromagnetic wave propagation
- Eddy current problems
- Maxwell equations
- Microwave engineering