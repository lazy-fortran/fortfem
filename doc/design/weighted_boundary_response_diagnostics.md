# Weighted boundary-response diagnostics

`evaluate_weighted_boundary_response_diagnostics` is a small, neutral
certificate for any square complex boundary map.  A caller supplies a matrix
`A` and positive diagonal work weights `W`; FortFEM inspects `W A` without
constructing a FEM, BEM, DtN, PML, wall, or free-boundary operator.

The normalized reciprocity diagnostic is

\[
  \frac{\lVert W A-(W A)^T\rVert_\max}
       {\max(1,\lVert W A\rVert_\max)}.
\]

The passivity diagnostic is the Gershgorin lower bound of the Hermitian part

\[
  H = \tfrac12\left(W A+(W A)^H\right), \qquad
  \min_i\left(H_{ii}-\sum_{j\ne i}|H_{ij}|\right).
\]

Thus a nonnegative returned lower bound is a conservative certificate for the
weighted matrix oracle.  It is not a claim about a physical model when the
caller has not supplied a physical work pairing, normalization, or radiation
condition.  A negative bound is reported as inconclusive, consistently with
the existing linear-response diagnostics.

The weights are deliberately caller-owned.  They may contain quadrature,
surface metric, impedance, or trace-duality factors from a physical FEM/BEM,
DtN, PML, NESTOR-like, BIEST-like, or virtual-casing adapter.  The diagnostic
therefore gives those backends one common reciprocity/passivity check while
keeping application formats and plasma closures outside FortFEM.

`test_weighted_boundary_response_diagnostics` computes both certificates from
the same matrix formulas independently and checks rejection of nonpositive
weights and nonsquare maps.

## Provenance

The work pairing and Hermitian passivity convention are the same real-part
complex convention used by `linear_response_interchange` and the boundary
trace contracts.  Reciprocity and positive-real/passive response matrices are
standard requirements for reciprocal electromagnetic boundary operators; the
Gershgorin bound is intentionally only a small, deterministic oracle for
this repository.  Production spectral, energy, or far-field certificates
remain owned by the external boundary implementation.
