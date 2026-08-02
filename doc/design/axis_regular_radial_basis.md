---
title: Axis-regular radial basis
---

# Axis-regular radial basis

`fortfem_axis_regular_radial_basis` evaluates one caller-selected finite
radial polynomial for a scalar Fourier label `m`,

\[
  f(\rho)=\sum_j c_j\rho^{p_j}.
\]

The powers must be strictly increasing and satisfy the same axis contract as
`build_axis_regular_mode_table`:

\[
  p_j\ge |m|, \qquad p_j\equiv |m|\pmod 2.
\]

Coefficients are complex so the result composes directly with Fourier-FEM and
nested-surface representations.  Samples are finite, nonnegative radial
coordinates; their normalization and outer endpoint remain caller-owned.
The evaluator does not select profiles, collocation grids, equilibrium
physics, or a DESC/VMEC/GVEC representation.  Vector and tensor components
may require shifted effective mode labels and remain separate spaces.

## Differentiation

The JVP differentiates both complex coefficients and real sample coordinates,

\[
  \dot f_s=\sum_j \left(
    \dot c_j\rho_s^{p_j}
    +c_jp_j\rho_s^{p_j-1}\dot\rho_s\right),
\]

with the derivative of the constant power defined as zero at the axis.  The
VJP is its exact transpose under the complex real-part pairing

\[
  \langle \bar f,\dot f\rangle
  =\operatorname{Re}\sum_s\overline{\bar f_s}\dot f_s.
\]

The integer mode label and power list are fixed selectors.  Activating a
different mode or power changes the discrete space and is deliberately not a
differentiable operation.

## Verification

`test_axis_regular_radial_basis` compares values with an independently written
polynomial, checks the combined coefficient/coordinate JVP by central
differences, closes the complex real-part adjoint identity, verifies vanishing
at the axis for `m>0`, and rejects bad parity, powers below `abs(m)`, duplicate
powers, negative samples, non-finite coefficients, and incompatible shapes.
