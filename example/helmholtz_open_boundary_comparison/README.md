# Helmholtz open-boundary comparison

This example compares three truncations of the same one-dimensional outgoing
Helmholtz wave:

- the exact modal Dirichlet-to-Neumann condition \(u'=iku\);
- a quadratic complex-stretched perfectly matched layer;
- a substantially larger ordinary domain terminated by homogeneous
  Dirichlet data.

The analytical oracle is \(u(x)=e^{ikx}\). The larger-domain result
deliberately demonstrates that distance alone is not a nonreflecting
condition for Helmholtz: an undamped one-dimensional reflection retains unit
magnitude regardless of where the wall is placed.

CI generates `helmholtz_methods_1d.png`, `helmholtz_error_1d.png`, and
`benchmark.txt`. Generated images are not committed.
