---
title: FCI boundary-patch mortar
---

# FCI boundary-patch mortar

FortFEM provides a neutral transfer contract for coupling an FCI background
plane to a fitted or cut boundary patch. The patch can be represented by its
own trace basis. The bulk FCI mesh does not have to conform to the vessel or
target surface.

The API consumes two independent trace tables evaluated at the same overlap
quadrature rows:

```text
background_trace(q, i) = B_i(x_q)
patch_trace(q, j)      = P_j(x_q)
overlap_weights(q)    = w_q > 0
ownership_count(q)    = 1
```

The final input is an ownership multiplicity. A value other than one rejects
the overlap before assembly. This prevents a cut or overlapping patch from
contributing the same physical surface twice. Geometry services are
responsible for constructing the rows, weights, and ownership records.

## Transfer algebra

The cross-mass is

\[
C_{ij}=\sum_q w_q B_i(x_q)P_j(x_q).
\]

FortFEM uses its row and column sums as diagonal trace measures,

\[
m^b_i=\sum_j C_{ij},\qquad m^p_j=\sum_i C_{ij}.
\]

The two transfer operators are

\[
T_{b\leftarrow p}=\operatorname{diag}(m^b)^{-1}C,
\qquad
T_{p\leftarrow b}=\operatorname{diag}(m^p)^{-1}C^T.
\]

When both trace tables reproduce constants, these definitions give

\[
T_{b\leftarrow p}{\bf 1}_p={\bf 1}_b,
\qquad
T_{p\leftarrow b}{\bf 1}_b={\bf 1}_p.
\]

They also satisfy the weighted adjoint identity

\[
x_b^T\operatorname{diag}(m^b)T_{b\leftarrow p}x_p
=
(T_{p\leftarrow b}x_b)^T\operatorname{diag}(m^p)x_p.
\]

The implementation rejects nonpositive overlap weights, non-partition-of-unity
trace tables, uncovered trace degrees of freedom, and numerically
rank-deficient cross-mass blocks. `overlap_measure` reports
\(\sum_q w_q\), which gives an independent cut or fitted surface-measure
check.

The JVP and VJP keep ownership and trace topology fixed. The JVP differentiates
the cross-mass, row and column measures, normalized transfers, and overlap
measure. The VJP is checked against the real dot-product identity for all of
those outputs. A change of owner, a new cut-cell row, or a rank change is a
discrete rebuild event and has no derivative under this contract.

## Public entry point

```fortran
call assemble_fci_boundary_patch_mortar( &
    background_trace, patch_trace, overlap_weights, ownership_count, &
    cross_mass, background_mass, patch_mass, background_from_patch, &
    patch_from_background, overlap_measure, status)
call assemble_fci_boundary_patch_mortar_jvp( &
    background_trace, patch_trace, overlap_weights, ownership_count, &
    background_trace_dot, patch_trace_dot, overlap_weights_dot, &
    cross_mass_dot, background_mass_dot, patch_mass_dot, &
    background_from_patch_dot, patch_from_background_dot, &
    overlap_measure_dot, status)
call assemble_fci_boundary_patch_mortar_vjp( &
    background_trace, patch_trace, overlap_weights, ownership_count, &
    cross_mass_bar, background_mass_bar, patch_mass_bar, &
    background_from_patch_bar, patch_from_background_bar, &
    overlap_measure_bar, background_trace_bar, patch_trace_bar, &
    overlap_weights_bar, status)
```

`background_from_patch` maps patch coefficients to the FCI background trace.
`patch_from_background` is its weighted adjoint. The interface is deliberately
agnostic about scalar versus vector trace laws. Tangential projections,
surface currents, sheath models, and material flux laws belong in the caller's
residual.

## Verification

`test_fci_boundary_patch_mortar` uses independent matrices for matching and
nonmatching traces. It checks the cross-mass oracle, row and column measures,
constant reproduction, the weighted adjoint identity, analytical overlap
measure, reversed quadrature orientation, zero-measure rejection, duplicate
ownership rejection, rank-deficient coupling rejection, the fixed-topology JVP
against central differences, and the full real VJP dot-product identity. These
checks use no vessel mesh or plasma-specific boundary model.

The contract is the fitted-patch part of the FCI path described in the
[parallel support-operator design](fci_parallel_operator.html). It implements
the numerical part of issue [#61](https://github.com/lazy-fortran/fortfem/issues/61);
moving Chimera topology, CAD meshing, and application boundary laws remain
outside FortFEM.
