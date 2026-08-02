---
title: Quadratic average-vector-field time step
---

# Quadratic average-vector-field time step

`advance_quadratic_avf` is the neutral ideal-time primitive for a
caller-owned quadratic Hamiltonian

\[
  H(x)=\tfrac12 x^T Kx+g^Tx.
\]

Given a constant interconnection matrix `J`, the update is the exact
average-vector-field/discrete-gradient formula

\[
  x_{n+1}-x_n=hJ\left[K\frac{x_n+x_{n+1}}2+g\right].
\]

It is assembled as one dense linear solve.  If `K` is symmetric and `J` is
skew, the quadratic Hamiltonian is preserved to linear-solve accuracy and a
signed step reversal recovers the input state.  The primitive does not accept
or hide damping, resistivity, viscosity, or PML terms; clients compose those
with `advance_dissipative_cayley` or another separately labelled block.

The companion `advance_quadratic_avf_jvp` and
`advance_quadratic_avf_vjp` routines differentiate the complete map with
respect to `K`, `J`, `g`, the step size, and the input state.  The focused
`test_quadratic_average_vector_field` test uses an independent three-by-three
linear-solve oracle, checks energy and reversibility, compares the JVP with a
central difference, and checks the real adjoint identity for the VJP.

This contract is intentionally physics-neutral: `x` can be a mixed acoustic,
elastic, Maxwell, or other port-Hamiltonian state supplied by FEM, IGA, DG,
or Fourier clients.
