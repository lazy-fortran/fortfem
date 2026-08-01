# Curl-Curl Electromagnetic Example

This example solves a manufactured three-dimensional curl-curl problem with
FortFEM's verified, arbitrary-order tetrahedral Nédélec elements and the
FortSparse direct solver.

## Problem

On the unit cube, solve

```text
curl(curl(E)) + E = J,
n × E = 0 on the boundary,
```

with

```text
E = (0, 0, sin(pi x) sin(pi y)),
J = (0, 0, (2 pi² + 1) sin(pi x) sin(pi y)).
```

Both the field and its analytical curl are sampled after covariant Piola
mapping. Their errors must decrease for Nédélec orders one through five; the
program stops with an error if that independent convergence oracle fails.

## Numerics

- first-kind tetrahedral Nédélec spaces of orders 1–5;
- globally oriented edge, face, and volume degrees of freedom;
- sparse curl-curl-plus-mass assembly and solution through FortSparse;
- homogeneous PEC tangential boundary elimination.

## Output

- `curl_curl_field_slice_2d.png`: the order-five solved Nédélec field
  magnitude with the reconstructed in-plane `curl(E_h)` arrows on a physical
  mid-plane; this is the primary visual result.
- `curl_curl_p_convergence.png`: field and curl errors by polynomial order;
- `convergence.csv`: the numerical values used to generate the plot.

Generated output is written under `output/example/curl_curl/` and is not
checked into the repository.
