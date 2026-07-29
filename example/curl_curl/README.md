# Curl-Curl Electromagnetic Example

This is an experimental API and visualization prototype for a future
Nédélec curl-curl solver. It demonstrates the intended form syntax, solver
dispatch, and vector plotting.

It is not currently a numerical validation example. The high-level vector
solver does not yet use the verified global Nédélec assembly or carry the
coefficient and source objects from the symbolic form into assembly.

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

## Intended Analytical Solution

The example uses the analytical solution:
```
E = [x*y, x²]
```

This gives:
- `curl(E) = ∂(x²)/∂x - ∂(x*y)/∂y = 2x - x = x`
- `curl(curl(E)) = ∂x/∂y - ∂(0)/∂x = 0`
- Therefore: `J = curl(curl(E)) + E = [x*y, x²]`

The program plots this field as a reference, but it does not compute an error
against it.

## Features Demonstrated

- Intended FEniCS-style vector form syntax.
- GMRES dispatch through the experimental vector solver.
- Tangential boundary-condition API.
- Vector-field streamplot output.

## Mathematical Formulation

The weak formulation seeks `E ∈ H(curl, Ω)` such that:

```
∫_Ω curl(E) · curl(F) dx + ∫_Ω E · F dx = ∫_Ω J · F dx    ∀F ∈ H₀(curl, Ω)
```

## Output Files

- `curlcurl_solution.png`: Vector field plot of numerical solution
- `curlcurl_exact.png`: Vector field plot of analytical solution
- `curl_field.png`: Scalar plot of the curl field

For verified lowest-order edge moments, covariant Piola mapping, orientations,
and operator assembly, see the Nédélec tests under `test/`.

## Applications

This formulation is fundamental for:
- Electromagnetic wave propagation
- Eddy current problems
- Maxwell equations
- Microwave engineering
